import argparse
import logging
import os
import re
import subprocess
import tempfile

import bioframe as bf
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", help="Input BAM file", type=str)
    parser.add_argument("gtf", help="Input GTF file", type=str)
    parser.add_argument(
        "--output",
        help="Output transcript assignment TSV",
        type=str,
        default="read_transcript_assigns.tsv",
    )
    parser.add_argument(
        "-q",
        "--mapq",
        help="Minimum MAPQ to use for transcript assignment [60]",
        type=int,
        default=60,
    )
    parser.add_argument(
        "-c",
        "--chunk_size",
        help="BED chunk size [200000]",
        type=int,
        default=200000,
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


def bam_to_bed(bam_path):
    tmp = tempfile.NamedTemporaryFile(suffix=".bed", delete=False)
    tmp.close()
    cmd = ["bedtools", "bamtobed", "-i", bam_path]
    with open(tmp.name, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)
    return tmp.name


def extract_attr(attr, key):
    m = re.search(rf'{key} "([^"]+)"', str(attr))
    if m:
        return m.group(1)
    return np.nan


def load_gtf(gtf_path):
    cols = [
        "chrom",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]
    df = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=cols,
        dtype=str,
        low_memory=False,
    )
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["chrom", "start", "end", "feature", "attribute"]).copy()

    df = df[df["feature"] == "transcript"].copy()
    if df.shape[0] == 0:
        raise ValueError("No transcript features found in GTF")

    df["transcript_id"] = df["attribute"].apply(lambda x: extract_attr(x, "transcript_id"))
    df["gene_name"] = df["attribute"].apply(lambda x: extract_attr(x, "gene_name"))
    df["gene_id"] = df["attribute"].apply(lambda x: extract_attr(x, "gene_id"))
    df["gene_label"] = df["gene_name"].fillna(df["gene_id"])

    df = df.dropna(subset=["transcript_id"]).copy()

    assert bf.is_bedframe(df), "GTF transcript table is not a valid bedframe"
    return df


def load_bed(bed_path):
    cols = ["chrom", "start", "end", "name", "score", "strand"]
    df = pd.read_csv(bed_path, sep="\t", header=None, names=cols)
    if df.shape[0] > 0:
        assert bf.is_bedframe(df), "BED file not loading as a valid bedframe"
    df["aln_len"] = df["end"] - df["start"]
    return df


def get_overlaps(bed, tx):
    df = bf.overlap(
        bed,
        tx,
        how="left",
        suffixes=("_bed", "_gtf"),
        return_overlap=True,
        return_index=True,
    )
    df = df[
        [
            "index_bed",
            "name_bed",
            "score_bed",
            "transcript_id_gtf",
            "gene_label_gtf",
            "overlap_start",
            "overlap_end",
        ]
    ].fillna(0)

    df = df.rename(
        columns={
            "name_bed": "read_id",
            "score_bed": "score",
            "transcript_id_gtf": "transcript_id",
            "gene_label_gtf": "gene",
        }
    )
    df["score"] = df["score"].astype(int)
    df["overlap_bp"] = df["overlap_end"] - df["overlap_start"]
    return df


def assign_status_low_mapq(df, mapq):
    df.loc[df["score"] < mapq, "status"] = "Unassigned_mapq"
    df.loc[df["score"] < mapq, "transcript_id"] = "NA"
    df.loc[df["score"] < mapq, "gene"] = "NA"
    return df


def assign_status_no_features(df):
    df.loc[df["transcript_id"] == 0, "status"] = "Unassigned_no_features"
    df.loc[df["transcript_id"] == 0, "transcript_id"] = "NA"
    df.loc[df["transcript_id"] == 0, "gene"] = "NA"
    return df


def assign_status_ambiguous_overlap(df):
    is_ambiguous = df[df["status"] == "Unknown"].duplicated(
        subset=["index_bed", "overlap_bp"], keep=False
    )
    ambiguous_idx = is_ambiguous.index[is_ambiguous]
    df.loc[ambiguous_idx, "status"] = "Unassigned_ambiguous"
    df.loc[ambiguous_idx, "transcript_id"] = "NA"
    return df.drop_duplicates(subset=["index_bed", "overlap_bp", "status"])


def find_largest_overlap(df):
    max_idx = df.groupby(["index_bed"])["overlap_bp"].idxmax().sort_values().values
    df = df.loc[max_idx, :]
    unknown_idx = [i for i in max_idx if df.loc[i, "status"] == "Unknown"]
    df.loc[unknown_idx, "status"] = "Assigned"
    return df


def process_bed_chunk(bed_chunk, tx, args):
    df = get_overlaps(bed_chunk, tx)
    df["status"] = "Unknown"
    df = assign_status_low_mapq(df, args.mapq)
    df = assign_status_no_features(df)
    df = assign_status_ambiguous_overlap(df)
    df = find_largest_overlap(df)
    df = df[["read_id", "status", "score", "gene", "transcript_id", "index_bed"]]
    df = df.reset_index(drop=True).drop(columns=["index_bed"])
    return df


def main(args):
    init_logger(args)

    tx = load_gtf(args.gtf)
    bed_path = bam_to_bed(args.bam)
    bed = load_bed(bed_path)

    ext = args.output.split(".")[-1]
    if (bed.shape[0] > 0) and (tx.shape[0] > 0):
        n = int(np.ceil(bed.shape[0] / args.chunk_size))
        chunk_fns = []

        for i, bed_chunk in enumerate(np.array_split(bed, n)):
            df_chunk = process_bed_chunk(bed_chunk, tx, args)
            fn = args.output.replace(ext, f"{i}.{ext}")
            chunk_fns.append(fn)
            df_chunk.to_csv(fn, sep="\t", index=False, header=False)

        with open(args.output, "w") as f_out:
            for fn in chunk_fns:
                with open(fn) as f_in:
                    for line in f_in:
                        f_out.write(line)

        for fn in chunk_fns:
            os.remove(fn)
    else:
        open(args.output, "w").close()

    os.remove(bed_path)


if __name__ == "__main__":
    args = parse_args()
    main(args)
