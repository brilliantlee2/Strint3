from collections import Counter
from pathlib import Path

import pandas as pd

from args_parser import set_parser
from utils import (
    assign_read,
    assign_reads_to_cells,
    build_core_cells,
    compute_top1_dominance,
    extract_reads_and_filter_df_by_raw,
    filter_pairs_three_stage,
    get_3p_features_new,
    get_5p_features_new,
    get_bc_whitelist,
    green_msg,
    is_missing,
    norm_bc,
    read_batch_generator,
    revcomp,
    strip_fixed_3p,
    strip_fixed_5p,
)


def write_barcode_lines(path, values):
    with open(path, "w", encoding="utf-8") as handle:
        for value in values:
            handle.write(f"{value}\n")


def count_high_quality_barcodes(df, bc_col, q_col, min_q):
    mask = df[q_col].fillna(-1) >= min_q
    series = df.loc[mask, bc_col].fillna("").astype(str).str.strip()
    series = series[series != ""]
    return Counter(series.value_counts().to_dict())


def build_putative_tables(args):
    read_ids = []
    putative_bcs = []
    putative_bc_min_qs = []
    bc_fixed_locs = []
    umis = []
    umi_fixed_locs = []
    post_umi_flankings = []
    polyA_starts = []
    read_types = []

    read_ids_5p = []
    putative_bcs_5p = []
    putative_bc_min_qs_5p = []
    bc_fixed_locs_5p = []
    umis_5p = []
    umi_fixed_locs_5p = []

    read_batchs = read_batch_generator(args.fastq_fns, args.batch_size)
    for batch in read_batchs:
        for read_info in batch:
            get_3p_features_new(
                read_info,
                read_ids,
                putative_bcs,
                bc_fixed_locs,
                putative_bc_min_qs,
                umis,
                umi_fixed_locs,
                post_umi_flankings,
                polyA_starts,
                read_types,
                FIX6_A=args.FIX6_A_3p,
                FIX5_B=args.FIX5_B_3p,
                FIX6_UMI=args.FIX6_UMI_3p,
            )
            get_5p_features_new(
                read_info,
                read_ids_5p,
                putative_bcs_5p,
                bc_fixed_locs_5p,
                putative_bc_min_qs_5p,
                umis_5p,
                umi_fixed_locs_5p,
                FIX6_A_5p=args.FIX6_A_5p,
                FIX5_B_5p=args.FIX5_B_5p,
                FIX6_UMI_5p=args.FIX6_UMI_5p,
            )

    rst_df_3p = pd.DataFrame(
        {
            "read_id": read_ids,
            "putative_bc": putative_bcs,
            "bc_fixed_locs": bc_fixed_locs,
            "putative_bc_min_qs": putative_bc_min_qs,
            "putative_umi": umis,
            "umi_fixed_locs": umi_fixed_locs,
            "post_umi_flankings": post_umi_flankings,
            "polyA_starts": polyA_starts,
            "read_types": read_types,
        }
    )
    rst_df_5p = pd.DataFrame(
        {
            "read_id": read_ids_5p,
            "putative_bc_5p": putative_bcs_5p,
            "bc_fixed_locs_50": bc_fixed_locs_5p,
            "putative_bc_min_qs_5p": putative_bc_min_qs_5p,
            "putative_umi_5p": umis_5p,
            "umi_fixed_locs_5p": umi_fixed_locs_5p,
        }
    )
    df_merge = rst_df_3p.merge(rst_df_5p, on="read_id", how="inner")

    if args.save_intermediate:
        rst_df_3p.to_csv(args.putative_bc_3p_out, index=False)
        rst_df_5p.to_csv(args.putative_bc_5p_out, index=False)
    df_merge.to_csv(args.putative_bc_out, index=False)

    return rst_df_3p, rst_df_5p, df_merge


def filter_corrected_reads(big_df):
    df = big_df.copy()
    df["BC3c"] = df["BC3_corrected"].map(norm_bc)
    df["BC5c"] = df["BC5_corrected"].map(norm_bc)

    paired = df[(df["BC3c"] != "") & (df["BC5c"] != "")].copy()
    pair_counts = (
        paired.groupby(["BC3c", "BC5c"])
        .size()
        .reset_index(name="pair_n_reads")
    )
    paired2 = paired.merge(pair_counts, on=["BC3c", "BC5c"], how="left")
    bad_read_ids = set(paired2.loc[paired2["pair_n_reads"] < 10, "read_id"])

    df_keep = df[~df["read_id"].isin(bad_read_ids)].copy()
    df_keep = df_keep[(df_keep["BC3c"] != "") | (df_keep["BC5c"] != "")].copy()

    if ("putative_umi" in df_keep.columns) and ("putative_umi_5p" in df_keep.columns):
        miss3 = df_keep["putative_umi"].map(is_missing)
        miss5 = df_keep["putative_umi_5p"].map(is_missing)
        n_removed_both_missing = int((miss3 & miss5).sum())
        df_keep = df_keep[~(miss3 & miss5)].copy()
    else:
        raise KeyError("df中缺少 putative_umi 或 putative_umi_5p，无法按新规则过滤。")

    df_keep = df_keep.drop(columns=["BC3c", "BC5c"], errors="ignore")
    drop_cols = [
        "read_types",
        "putative_bc",
        "bc_fixed_locs",
        "putative_bc_min_qs",
        "umi_fixed_locs",
        "polyA_starts",
        "post_umi_flankings",
        "umi_fixed_locs_5p",
        "bc_fixed_locs_50",
        "putative_bc_min_qs_5p",
        "putative_bc_5p",
        "strand",
    ]
    big_df_filtered = df_keep.drop(columns=drop_cols, errors="ignore")

    stats = {
        "original_rows": int(len(big_df)),
        "bad_paired_removed": int(len(bad_read_ids)),
        "both_empty_barcode_rows": int(((df["BC3c"] == "") & (df["BC5c"] == "")).sum()),
        "both_umi_missing_removed": int(n_removed_both_missing),
        "filtered_rows": int(len(big_df_filtered)),
    }
    return big_df_filtered, stats


def prepare_final_read_table(df):
    df = df.copy()
    df["BC5_30bp"] = df["BC5_corrected"].map(strip_fixed_5p)
    df["BC3_30bp_rc"] = df["BC3_corrected"].map(lambda x: revcomp(strip_fixed_3p(x)))
    df = df.drop(columns=["BC3_corrected", "BC5_corrected", "has_3p", "has_5p"], errors="ignore")
    return df


def build_clean_exports(df_final, pc_final):
    df_out = df_final.copy()
    df_out["BC5n"] = df_out["BC5_30bp"].map(norm_bc) if "BC5_30bp" in df_out.columns else ""
    df_out["BC3n"] = df_out["BC3_30bp_rc"].map(norm_bc) if "BC3_30bp_rc" in df_out.columns else ""

    cols_read = [
        "read_id",
        "putative_umi",
        "putative_umi_5p",
        "BC5n",
        "BC3n",
    ]
    cols_read = [col for col in cols_read if col in df_out.columns]
    df_out = df_out[cols_read].copy()

    pair_counts_kept = pc_final.copy()
    pair_counts_kept = pair_counts_kept[["BC5n", "BC3n", "support_reads"]].copy()
    return df_out, pair_counts_kept


def build_assigned_reads(df_out, pair_counts_kept, dominance_min):
    _, _, core_cells_df, barcode2cell = build_core_cells(pair_counts_kept)
    dom_table = compute_top1_dominance(pair_counts_kept)
    df_assigned, assign_stats = assign_reads_to_cells(
        df_out,
        barcode2cell,
        dom_table,
        dominance_min=dominance_min,
    )

    keep_cols = ["read_id", "putative_umi", "putative_umi_5p", "BC5_30bp", "BC3_30bp_rc", "cell_id"]
    keep_exist = [col for col in keep_cols if col in df_assigned.columns]
    df_assigned_slim = df_assigned[keep_exist].copy()
    df_assigned_slim = df_assigned_slim[
        df_assigned_slim["cell_id"].notna() & (df_assigned_slim["cell_id"] != "")
    ].copy()

    core_type_counts = core_cells_df["type"].value_counts()
    cell_read_stats = (
        df_assigned.dropna(subset=["cell_id"])
        .groupby("cell_id")["read_id"].size()
        .sort_values(ascending=False)
        .reset_index(name="n_reads")
    )

    return (
        core_cells_df,
        df_assigned,
        df_assigned_slim,
        assign_stats,
        core_type_counts,
        cell_read_stats,
        len(barcode2cell),
    )


def run_pipeline(args):
    green_msg("Step 1/7: splitting 3p and 5p putative barcode tables", printit=True)
    rst_df_3p, rst_df_5p, df_merge = build_putative_tables(args)

    green_msg("Step 2/7: counting high-quality putative barcodes and building whitelists", printit=True)
    raw_bc_count_3p = count_high_quality_barcodes(rst_df_3p, "putative_bc", "putative_bc_min_qs", args.minQ)
    raw_bc_count_5p = count_high_quality_barcodes(rst_df_5p, "putative_bc_5p", "putative_bc_min_qs_5p", args.minQ)

    bc_whitelist_3p, ept_bc_3p = get_bc_whitelist(
        raw_bc_count=raw_bc_count_3p,
        full_bc_whitelist=args.full_bc_whitelist_3p,
        exp_cells=args.exp_cells,
        out_plot_fn=args.knee_plot_3p_out,
        DEFAULT_EMPTY_DROP_MIN_ED=args.DEFAULT_EMPTY_DROP_MIN_ED,
        DEFAULT_EMPTY_DROP_NUM=args.DEFAULT_EMPTY_DROP_NUM,
    )
    bc_whitelist_5p, ept_bc_5p = get_bc_whitelist(
        raw_bc_count=raw_bc_count_5p,
        full_bc_whitelist=args.full_bc_whitelist_5p,
        exp_cells=args.exp_cells,
        out_plot_fn=args.knee_plot_5p_out,
        DEFAULT_EMPTY_DROP_MIN_ED=args.DEFAULT_EMPTY_DROP_MIN_ED,
        DEFAULT_EMPTY_DROP_NUM=args.DEFAULT_EMPTY_DROP_NUM,
    )
    write_barcode_lines(args.whitelist_3p_out, bc_whitelist_3p.keys())
    write_barcode_lines(args.whitelist_5p_out, bc_whitelist_5p.keys())
    write_barcode_lines(args.emptydrop_3p_out, ept_bc_3p)
    write_barcode_lines(args.emptydrop_5p_out, ept_bc_5p)

    green_msg("Step 3/7: correcting reads and writing matched FASTQ", printit=True)
    demul_count_tot, count_tot, big_df = assign_read(
        fastq_fns=args.fastq_fns,
        fastq_out=args.fastq_out,
        putative_bc_csv=args.putative_bc_out,
        whitelsit_3p=args.whitelist_3p_out,
        whitelsit_5p=args.whitelist_5p_out,
        max_ed=args.max_ed,
        n_process=args.threads,
        batchsize=args.assign_batchsize,
        minQ=args.minQ,
    )
    if args.save_intermediate:
        big_df.to_csv(args.corrected_bc_out, index=False)

    green_msg("Step 4/7: filtering corrected reads", printit=True)
    big_df_filtered, filter_stats = filter_corrected_reads(big_df)

    green_msg("Step 5/7: stripping fixed sequences and filtering barcode pairs", printit=True)
    df_for_pairs = prepare_final_read_table(big_df_filtered)
    (
        df_final,
        _df_single,
        _df_paired_final,
        _pc_all,
        _pc_min,
        pc_final,
        _pc_dropped,
        pair_stats,
    ) = filter_pairs_three_stage(
        df_for_pairs,
        bc5_col="BC5_30bp",
        bc3_col="BC3_30bp_rc",
        umi3_col="putative_umi",
        umi5_col="putative_umi_5p",
        PAIR_MIN=args.PAIR_MIN,
        TOP1_ALPHA=args.TOP1_ALPHA,
        require_pass_both_ends=args.require_pass_both_ends,
        drop_umiA_ratio_gt=args.drop_umiA_ratio_gt,
    )

    green_msg("Step 6/7: exporting clean read table and pair table", printit=True)
    df_out, pair_counts_kept = build_clean_exports(df_final, pc_final)
    df_out.to_csv(args.clean_reads_out, index=False)
    pair_counts_kept.to_csv(args.pair_counts_out, index=False)

    green_msg("Step 7/7: assigning reads to cells and extracting cell FASTQ", printit=True)
    df_for_assignment = df_final.copy()
    df_for_assignment["BC5n"] = df_for_assignment["BC5_30bp"].map(norm_bc)
    df_for_assignment["BC3n"] = df_for_assignment["BC3_30bp_rc"].map(norm_bc)
    (
        core_cells_df,
        _df_assigned,
        df_assigned_slim,
        assign_stats,
        core_type_counts,
        cell_read_stats,
        core_barcode_count,
    ) = build_assigned_reads(
        df_for_assignment,
        pair_counts_kept,
        dominance_min=args.dominance_min,
    )
    cell_read_stats.to_csv(args.cell_read_stats_out, index=False)

    df_kept, extract_stats = extract_reads_and_filter_df_by_raw(
        df_assigned_slim,
        raw_fastq_gz=args.fastq_out,
        out_fastq_gz=args.cell_fastq_out,
        read_id_col="read_id",
        remove_found=True,
    )
    df_kept.to_csv(args.read_assigned_out, index=False, encoding="utf-8")

    summary = {
        "fastq_files": len(args.fastq_fns),
        "reads_total": int(count_tot),
        "reads_demultiplexed": int(demul_count_tot),
        "putative_rows_3p": int(len(rst_df_3p)),
        "putative_rows_5p": int(len(rst_df_5p)),
        "merged_rows": int(len(df_merge)),
        "filtered_rows": int(len(big_df_filtered)),
        "clean_reads_rows": int(len(df_out)),
        "pair_counts_rows": int(len(pair_counts_kept)),
        "assigned_rows": int(len(df_kept)),
        "core_cells": int(len(core_cells_df)),
        "core_barcodes": int(core_barcode_count),
        "core_type_counts": core_type_counts.to_dict(),
        "cell_read_stats_head": cell_read_stats.head().to_dict(orient="records"),
        "filter_stats": filter_stats,
        "pair_stats": pair_stats,
        "assign_stats": assign_stats,
        "extract_stats": extract_stats,
    }
    return summary


def print_summary(summary, out_dir):
    print("\n=== Strint Summary ===")
    print(f"Output directory: {out_dir}")
    print(f"FASTQ files: {summary['fastq_files']}")
    print(f"Reads total: {summary['reads_total']}")
    print(f"Reads demultiplexed: {summary['reads_demultiplexed']}")
    print(f"Merged putative rows: {summary['merged_rows']}")
    print(f"Filtered rows: {summary['filtered_rows']}")
    print(f"Clean read rows: {summary['clean_reads_rows']}")
    print(f"Pair rows kept: {summary['pair_counts_rows']}")
    print("\nCore cell type counts:")
    for cell_type, count in summary["core_type_counts"].items():
        print(f"{cell_type:16} {count}")
    print(f"core cells: {summary['core_cells']} core barcodes: {summary['core_barcodes']}")
    print(f"Assigned read rows: {summary['assigned_rows']}")
    print("\nassign_stats:")
    print(summary["assign_stats"])
    print("\nTop cell_read_stats:")
    for row in summary["cell_read_stats_head"]:
        print(row)


def main():
    args = set_parser()
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    summary = run_pipeline(args)
    print_summary(summary, args.out_dir)


if __name__ == "__main__":
    main()
