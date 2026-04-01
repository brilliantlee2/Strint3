import argparse
import collections
import logging

import pandas as pd
import pysam
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", help="Tagged BAM with CB and UB tags", type=str)
    parser.add_argument("transcripts", help="Read transcript assignments TSV", type=str)
    parser.add_argument(
        "--output",
        help="Output isoform expression matrix TSV",
        type=str,
        default="isoform_expression.tsv",
    )
    parser.add_argument(
        "--verbosity",
        help="logging level: <=2 logs info, <=3 logs warnings",
        type=int,
        default=2,
    )
    return parser.parse_args()


def init_logger(args):
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging.root.setLevel(args.verbosity * 10)


def load_transcript_assignments(path):
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["read_id", "status", "score", "gene", "transcript_id"],
        dtype=str,
    )
    df = df[df["status"] == "Assigned"].copy()
    df = df[df["transcript_id"] != "NA"].copy()
    df = df.drop_duplicates(subset=["read_id"], keep="first")
    return dict(zip(df["read_id"], df["transcript_id"]))


def get_bam_info(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    stats = bam.get_index_statistics()
    n_reads = int(sum([contig.mapped for contig in stats]))
    bam.close()
    return n_reads


def build_isoform_sets(bam_path, read_to_tx, n_reads):
    bam = pysam.AlignmentFile(bam_path, "rb")
    txs = set()
    cells = set()
    umi_sets = collections.defaultdict(set)

    for align in tqdm(bam.fetch(), total=n_reads):
        read_id = align.query_name

        if read_id not in read_to_tx:
            continue
        if not align.has_tag("CB") or not align.has_tag("UB"):
            continue

        tx = read_to_tx[read_id]
        cell = align.get_tag("CB")
        umi = align.get_tag("UB")

        txs.add(tx)
        cells.add(cell)
        umi_sets[(tx, cell)].add(umi)

    bam.close()
    return txs, cells, umi_sets


def populate_matrix(txs, cells, umi_sets):
    rows = sorted(list(txs))
    cols = sorted(list(cells))

    df = pd.DataFrame(0, index=rows, columns=cols)

    for (tx, cell), umis in tqdm(umi_sets.items()):
        df.loc[tx, cell] = len(umis)

    return df


def main(args):
    init_logger(args)

    logger.info(f"Loading transcript assignments from {args.transcripts}")
    read_to_tx = load_transcript_assignments(args.transcripts)

    logger.info(f"Reading BAM tags from {args.bam}")
    n_reads = get_bam_info(args.bam)
    txs, cells, umi_sets = build_isoform_sets(args.bam, read_to_tx, n_reads)

    logger.info("Building isoform expression matrix")
    df = populate_matrix(txs, cells, umi_sets)
    df.to_csv(args.output, sep="\t", index_label="transcript_id")


if __name__ == "__main__":
    args = parse_args()
    main(args)
