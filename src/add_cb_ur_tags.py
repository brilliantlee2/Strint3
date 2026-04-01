# 2) add_cb_ur_tags.py
# 用法示例：
# python add_cb_ur_tags.py \
#   --bam filtered.sorted.bam \
#   --tags read_tags.tsv \
#   --output filtered.cb_ur.sorted.bam

import argparse
import pysam
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, help="输入 BAM")
    parser.add_argument("--tags", required=True, help="read_tags.tsv")
    parser.add_argument("--output", required=True, help="输出 BAM")
    parser.add_argument("--read-id-col", default="read_id")
    parser.add_argument("--cb-col", default="cell_id")
    parser.add_argument("--ur-col", default="umi_for_clustering")
    return parser.parse_args()


def load_tag_map(path, read_id_col, cb_col, ur_col):
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = [read_id_col, cb_col, ur_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"缺少列: {missing}")

    df = df[df[read_id_col].str.strip() != ""].copy()
    df = df.drop_duplicates(subset=[read_id_col], keep="first")

    tag_map = {}
    for _, row in df.iterrows():
        rid = row[read_id_col].strip()
        cb = row[cb_col].strip()
        ur = row[ur_col].strip()
        if cb and ur:
            tag_map[rid] = (cb, ur)
    return tag_map


def main():
    args = parse_args()
    tag_map = load_tag_map(args.tags, args.read_id_col, args.cb_col, args.ur_col)

    bam = pysam.AlignmentFile(args.bam, "rb")
    out = pysam.AlignmentFile(args.output, "wb", template=bam)

    total = 0
    tagged = 0

    for aln in bam.fetch(until_eof=True):
        total += 1
        rid = aln.query_name
        if rid in tag_map:
            cb, ur = tag_map[rid]
            aln.set_tag("CB", cb, value_type="Z")
            aln.set_tag("UR", ur, value_type="Z")
            tagged += 1
        out.write(aln)

    bam.close()
    out.close()
    pysam.index(args.output)

    print(f"total_alignments\t{total}")
    print(f"tagged_alignments\t{tagged}")
    print(f"unmatched_alignments\t{total - tagged}")


if __name__ == "__main__":
    main()
