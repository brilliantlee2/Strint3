# 1) prepare_read_tags.py
# 用法示例：
# python prepare_read_tags.py \
#   --input filtered_reads.tsv \
#   --output read_tags.tsv \
#   --read-id-col read_id \
#   --cell-col cell_id \
#   --umi-primary-col putative_umi_5p \
#   --umi-backup-col putative_umi

import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="输入表格，支持 .tsv/.csv")
    parser.add_argument("--output", required=True, help="输出 TSV")
    parser.add_argument("--read-id-col", default="read_id")
    parser.add_argument("--cell-col", default="cell_id")
    parser.add_argument("--umi-primary-col", default="putative_umi_5p")
    parser.add_argument("--umi-backup-col", default="putative_umi")
    return parser.parse_args()


def read_table(path):
    if path.endswith(".csv"):
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def main():
    args = parse_args()
    df = read_table(args.input).copy()

    required = [
        args.read_id_col,
        args.cell_col,
        args.umi_primary_col,
        args.umi_backup_col,
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"缺少列: {missing}")

    df[args.umi_primary_col] = df[args.umi_primary_col].fillna("").astype(str).str.strip()
    df[args.umi_backup_col] = df[args.umi_backup_col].fillna("").astype(str).str.strip()
    df[args.read_id_col] = df[args.read_id_col].astype(str).str.strip()
    df[args.cell_col] = df[args.cell_col].astype(str).str.strip()

    # 主端优先，缺失时回退到备用端
    df["umi_for_clustering"] = df[args.umi_primary_col]
    mask_empty = df["umi_for_clustering"].eq("")
    df.loc[mask_empty, "umi_for_clustering"] = df.loc[mask_empty, args.umi_backup_col]

    out = df[
        [
            args.read_id_col,
            args.cell_col,
            args.umi_primary_col,
            args.umi_backup_col,
            "umi_for_clustering",
        ]
    ].copy()
    out.columns = ["read_id", "cell_id", "umi_primary", "umi_backup", "umi_for_clustering"]

    out = out[out["read_id"] != ""]
    out = out[out["cell_id"] != ""]
    out = out[out["umi_for_clustering"] != ""]

    # 如果同一个 read_id 有重复记录，保留第一条
    out = out.drop_duplicates(subset=["read_id"], keep="first")

    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
