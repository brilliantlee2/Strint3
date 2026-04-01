# This code makes significant use of the UMI-tools package (MIT license).
#
# https://github.com/CGATOxford/UMI-tools
# https://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract
#
# Modified to support a whole BAM as input rather than chromosome-split BAMs.

import argparse
import collections
import itertools
import logging
import multiprocessing
import pathlib
import tempfile

import numpy as np
import pandas as pd
import pysam
from editdistance import eval as edit_distance
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "bam",
        help="BAM file of alignments with tags for gene (GN), corrected barcode (CB) and uncorrected UMI (UR)",
        type=str,
    )

    parser.add_argument(
        "--output",
        help="Output BAM file with corrected UMI tag (UB) [tagged.sorted.bam]",
        type=str,
        default="tagged.sorted.bam",
    )

    parser.add_argument(
        "-i",
        "--ref_interval",
        help="Size of genomic window (bp) to assign as gene name if no gene assigned by feature annotation [1000]",
        type=int,
        default=1000,
    )

    parser.add_argument(
        "--cell_gene_max_reads",
        help="Maximum number of reads to consider for a particular gene + cell barcode combination [20000]",
        type=int,
        default=20000,
    )

    parser.add_argument(
        "-t",
        "--threads",
        help="Threads to use [4]",
        type=int,
        default=4,
    )

    parser.add_argument(
        "--verbosity",
        help="logging level: <=2 logs info, <=3 logs warnings",
        type=int,
        default=2,
    )

    args = parser.parse_args()

    p = pathlib.Path(args.output)
    tempdir = tempfile.TemporaryDirectory(prefix="tmp.", dir=str(p.parents[0]))
    args.tempdir = tempdir.name

    return args


def init_logger(args):
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


def breadth_first_search(node, adj_list):
    searched = set()
    queue = set()
    queue.update((node,))
    searched.update((node,))

    while len(queue) > 0:
        node = queue.pop()
        for next_node in adj_list[node]:
            if next_node not in searched:
                queue.update((next_node,))
                searched.update((next_node,))

    return searched


def get_adj_list_directional(umis, counts, threshold=2):
    """
    Identify all UMIs within the Levenshtein distance threshold
    and where the counts of the first UMI is > (2 * second UMI counts)-1.
    """
    adj_list = {umi: [] for umi in umis}
    iter_umi_pairs = itertools.combinations(umis, 2)
    for umi1, umi2 in iter_umi_pairs:
        if edit_distance(umi1, umi2) <= threshold:
            if counts[umi1] >= (counts[umi2] * 2) - 1:
                adj_list[umi1].append(umi2)
            if counts[umi2] >= (counts[umi1] * 2) - 1:
                adj_list[umi2].append(umi1)

    return adj_list


def get_connected_components_adjacency(umis, graph, counts):
    found = set()
    components = []

    for node in sorted(graph, key=lambda x: counts[x], reverse=True):
        if node not in found:
            component = breadth_first_search(node, graph)
            found.update(component)
            components.append(component)
    return components


def group_directional(clusters, adj_list, counts):
    observed = set()
    groups = []
    for cluster in clusters:
        if len(cluster) == 1:
            groups.append(list(cluster))
            observed.update(cluster)
        else:
            cluster = sorted(cluster, key=lambda x: counts[x], reverse=True)
            temp_cluster = []
            for node in cluster:
                if node not in observed:
                    temp_cluster.append(node)
                    observed.add(node)
            groups.append(temp_cluster)

    return groups


def cluster(counts_dict, threshold=3):
    adj_list = get_adj_list_directional(counts_dict.keys(), counts_dict, threshold)
    clusters = get_connected_components_adjacency(
        counts_dict.keys(), adj_list, counts_dict
    )
    final_umis = [list(x) for x in group_directional(clusters, adj_list, counts_dict)]
    return final_umis


def create_map_to_correct_umi(cluster_list):
    return {y: x[0] for x in cluster_list for y in x}


def correct_umis(umis):
    counts_dict = dict(umis.value_counts())
    umi_map = create_map_to_correct_umi(cluster(counts_dict))
    return umis.replace(umi_map)


def create_region_name(align, args):
    chrom = align.reference_name
    ref_positions = align.get_reference_positions()

    if not ref_positions:
        return "NA"

    start_pos = ref_positions[0]
    end_pos = ref_positions[-1]

    midpoint = int((start_pos + end_pos) / 2)
    interval_start = np.floor(midpoint / args.ref_interval) * args.ref_interval
    interval_end = np.ceil(midpoint / args.ref_interval) * args.ref_interval

    gene = f"{chrom}_{int(interval_start)}_{int(interval_end)}"
    return gene


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def launch_pool(func, func_args, procs=1):
    p = multiprocessing.Pool(processes=procs)
    try:
        results = list(tqdm(p.imap(func, func_args), total=len(func_args)))
        p.close()
        p.join()
    except KeyboardInterrupt:
        p.terminate()
        raise
    return results


def run_groupby(df):
    df["umi_corr"] = df.groupby(["gene_cell"])["umi_uncorr"].transform(correct_umis)
    return df


def collect_records_from_bam(input_bam, args):
    """
    Read all alignments from the BAM, extract read_id / gene / barcode / UR,
    and return a dataframe for UMI clustering.
    """
    bam = pysam.AlignmentFile(input_bam, "rb")

    records = []
    cell_gene_counter = collections.Counter()

    logger.info(f"Collecting read/gene/barcode/UMI records from {input_bam}")
    for align in tqdm(bam.fetch(until_eof=True)):
        if align.is_unmapped:
            continue

        read_id = align.query_name

        for tag in ["UR", "CB", "GN"]:
            assert align.has_tag(tag), f"{tag} tag not found in {read_id}"

        bc_corr = align.get_tag("CB")
        umi_uncorr = align.get_tag("UR")
        gene = align.get_tag("GN")

        if gene == "NA":
            gene = create_region_name(align, args)

        cell_gene_counter[(bc_corr, gene)] += 1
        if cell_gene_counter[(bc_corr, gene)] <= args.cell_gene_max_reads:
            records.append((read_id, gene, bc_corr, umi_uncorr))

    bam.close()

    df = pd.DataFrame.from_records(
        records, columns=["read_id", "gene", "bc", "umi_uncorr"]
    )
    return df


def compute_corrected_umis(df, args):
    """
    Run per-(gene, cell) UMI correction on the dataframe.
    Returns:
      umis: dict read_id -> corrected UMI
      genes: dict read_id -> gene label
    """
    if df.shape[0] == 0:
        return {}, {}

    df["gene_cell"] = df["gene"] + ":" + df["bc"]
    df = df.set_index("gene_cell")

    gene_cell_unique = list(set(df.index))
    gene_cell_per_chunk = 50
    gene_cell_chunks = chunks(gene_cell_unique, gene_cell_per_chunk)

    func_args = []
    for gene_cell_chunk in gene_cell_chunks:
        df_ = df.loc[gene_cell_chunk]
        func_args.append(df_)

    logger.info("Clustering and correcting UMIs")
    results = launch_pool(run_groupby, func_args, args.threads)

    if len(results) > 0:
        df = pd.concat(results, axis=0)
    else:
        df = pd.DataFrame(columns=["read_id", "gene", "bc", "umi_uncorr", "umi_corr"])

    df = df.drop(["bc", "umi_uncorr"], axis=1).set_index("read_id")

    umis = df.to_dict()["umi_corr"]
    genes = df.to_dict()["gene"]
    return umis, genes


def add_tags_all(umis, genes, args):
    """
    Write UB/GN tags back to the whole BAM.
    Existing alignments are preserved; UB/GN are added where a read_id match exists.
    """
    bam = pysam.AlignmentFile(args.bam, "rb")
    bam_out = pysam.AlignmentFile(args.output, "wb", template=bam)

    logger.info(f"Writing corrected tags to {args.output}")
    for align in tqdm(bam.fetch(until_eof=True)):
        read_id = align.query_name

        if (umis.get(read_id) is not None) and (genes.get(read_id) is not None):
            align.set_tag("UB", umis[read_id], value_type="Z")
            align.set_tag("GN", genes[read_id], value_type="Z")

        bam_out.write(align)

    bam.close()
    bam_out.close()

    logger.info(f"Indexing {args.output}")
    pysam.index(args.output)


def main(args):
    init_logger(args)

    df = collect_records_from_bam(args.bam, args)
    logger.info(f"Collected {df.shape[0]} records for UMI correction")

    umis, genes = compute_corrected_umis(df, args)
    logger.info(f"Computed corrected UMI tags for {len(umis)} reads")

    add_tags_all(umis, genes, args)


if __name__ == "__main__":
    args = parse_args()
    main(args)
