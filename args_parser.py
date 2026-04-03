import argparse
import os
import sys
from pathlib import Path

from utils import err_msg


def check_files_exist(file_list):
    if isinstance(file_list, str):
        file_list = [file_list]
    exit_code = 0
    for fn in file_list:
        if not os.path.exists(fn):
            exit_code = 1
            err_msg(f"Error: can not find {fn}", printit=True)
    if exit_code == 1:
        sys.exit(1)
    return True


def get_files_by_suffix(search_dir, suffix, recursive=True):
    files = []
    if isinstance(suffix, str):
        suffix = [suffix]
    if recursive:
        for item in suffix:
            files.extend(Path(search_dir).rglob(item))
    else:
        for item in suffix:
            files.extend(Path(search_dir).glob(item))
    return sorted(files)


def get_files_from_dir(fastq_dir):
    check_files_exist(fastq_dir)
    if os.path.isdir(fastq_dir):
        fastq_fns = get_files_by_suffix(
            search_dir=fastq_dir,
            suffix=["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"],
            recursive=True,
        )
        if not fastq_fns:
            err_msg(f"No FASTQ files found in directory: {fastq_dir}", printit=True)
            sys.exit(1)
    elif os.path.isfile(fastq_dir):
        fastq_fns = [fastq_dir]
    else:
        err_msg(f"File type of input {fastq_dir} is not supported.", printit=True)
        sys.exit(1)
    return list(map(str, fastq_fns))


def _to_outdir(out_dir, p):
    return p if os.path.isabs(p) else os.path.join(out_dir, os.path.basename(p))


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Strint: barcode splitting, whitelist generation, read correction, "
            "pair filtering, and cell-level assignment for full-length reads."
        )
    )
    parser.add_argument(
        "fastq_fns",
        metavar="<input fastq filename/directory>",
        help=(
            "Full-length sequencing fastq file (.fq/.fastq/.gz) or a directory "
            "containing them."
        ),
        type=get_files_from_dir,
    )

    parser.add_argument(
        "--full-bc-whitelist-3p",
        dest="full_bc_whitelist_3p",
        type=str,
        required=True,
        help="3p full whitelist file path.",
    )
    parser.add_argument(
        "--full-bc-whitelist-5p",
        dest="full_bc_whitelist_5p",
        type=str,
        required=True,
        help="5p full whitelist file path.",
    )

    parser.add_argument(
        "--out_dir",
        type=str,
        default="Strint",
        help="Directory to save all output files (default: Strint).",
    )
    parser.add_argument(
        "--save-intermediate",
        dest="save_intermediate",
        action="store_true",
        help="Write large intermediate tables such as putative_bc_p3/p5 and BC_corrected.",
    )

    parser.add_argument("--batch_size", type=int, default=100000, help="FASTQ batch size.")
    parser.add_argument("--assign_batchsize", type=int, default=10000, help="Batch size used in assign_read.")
    parser.add_argument("--threads", type=int, default=256, help="Worker count for assign_read.")
    parser.add_argument("--minQ", type=int, default=7, help="Minimum barcode quality score for counting.")
    parser.add_argument("--exp_cells", type=int, default=5000, help="Expected number of cells.")
    parser.add_argument("--max_ed", type=int, default=7, help="Maximum edit distance for correction.")
    parser.add_argument(
        "--DEFAULT_EMPTY_DROP_MIN_ED",
        type=int,
        default=5,
        help="Minimum edit distance from empty-drop BC to selected BC.",
    )
    parser.add_argument(
        "--DEFAULT_EMPTY_DROP_NUM",
        type=int,
        default=2000,
        help="Maximum number of empty-drop barcodes to retain.",
    )
    parser.add_argument("--PAIR_MIN", type=int, default=10, help="Minimum read support per barcode pair.")
    parser.add_argument("--TOP1_ALPHA", type=float, default=0.7, help="Top1-anchor alpha threshold.")
    parser.add_argument(
        "--dominance_min",
        type=float,
        default=0.80,
        help="Dominance threshold for absorbing single-end reads into core cells.",
    )
    parser.add_argument(
        "--drop_umiA_ratio_gt",
        type=float,
        default=0.5,
        help="Drop reads when 3p putative UMI A-ratio is above this threshold.",
    )
    parser.add_argument(
        "--require_pass_both_ends",
        action="store_true",
        help="Require both ends to pass the top1-anchor rule in pair filtering.",
    )

    parser.add_argument("--FIX6_A_3p", type=str, default="GCTACC")
    parser.add_argument("--FIX5_B_3p", type=str, default="AGATC")
    parser.add_argument("--FIX6_UMI_3p", type=str, default="TAGGCT")
    parser.add_argument("--FIX6_A_5p", type=str, default="CCTTCC")
    parser.add_argument("--FIX5_B_5p", type=str, default="TGCTG")
    parser.add_argument("--FIX6_UMI_5p", type=str, default="CCACTG")

    parser.add_argument("--putative_bc_out", type=str, default="putative_bc.csv")
    parser.add_argument("--putative_bc_3p_out", type=str, default="putative_bc_p3.csv")
    parser.add_argument("--putative_bc_5p_out", type=str, default="putative_bc_p5.csv")
    parser.add_argument("--whitelist_3p_out", type=str, default="whitelist_3p.csv")
    parser.add_argument("--whitelist_5p_out", type=str, default="whitelist_5p.csv")
    parser.add_argument("--emptydrop_3p_out", type=str, default="empty_bc_list_3p.csv")
    parser.add_argument("--emptydrop_5p_out", type=str, default="empty_bc_list_5p.csv")
    parser.add_argument("--knee_plot_3p_out", type=str, default="knee_plot_3p.png")
    parser.add_argument("--knee_plot_5p_out", type=str, default="knee_plot_5p.png")
    parser.add_argument("--fastq_out", type=str, default="matched_reads.fastq.gz")
    parser.add_argument("--corrected_bc_out", type=str, default="BC_corrected.csv")
    parser.add_argument("--clean_reads_out", type=str, default="ReadIDs_UMI_BC_clean.csv")
    parser.add_argument("--pair_counts_out", type=str, default="pair_counts_kept.csv")
    parser.add_argument("--read_assigned_out", type=str, default="read_assigned_cell.csv")
    parser.add_argument("--cell_fastq_out", type=str, default="cell_reads.fastq.gz")
    parser.add_argument("--cell_read_stats_out", type=str, default="cell_read_stats.csv")

    return parser


def set_parser(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    check_files_exist([args.full_bc_whitelist_3p, args.full_bc_whitelist_5p])

    os.makedirs(args.out_dir, exist_ok=True)

    args.putative_bc_out = _to_outdir(args.out_dir, args.putative_bc_out)
    args.putative_bc_3p_out = _to_outdir(args.out_dir, args.putative_bc_3p_out)
    args.putative_bc_5p_out = _to_outdir(args.out_dir, args.putative_bc_5p_out)
    args.whitelist_3p_out = _to_outdir(args.out_dir, args.whitelist_3p_out)
    args.whitelist_5p_out = _to_outdir(args.out_dir, args.whitelist_5p_out)
    args.emptydrop_3p_out = _to_outdir(args.out_dir, args.emptydrop_3p_out)
    args.emptydrop_5p_out = _to_outdir(args.out_dir, args.emptydrop_5p_out)
    args.knee_plot_3p_out = _to_outdir(args.out_dir, args.knee_plot_3p_out)
    args.knee_plot_5p_out = _to_outdir(args.out_dir, args.knee_plot_5p_out)
    args.fastq_out = _to_outdir(args.out_dir, args.fastq_out)
    args.corrected_bc_out = _to_outdir(args.out_dir, args.corrected_bc_out)
    args.clean_reads_out = _to_outdir(args.out_dir, args.clean_reads_out)
    args.pair_counts_out = _to_outdir(args.out_dir, args.pair_counts_out)
    args.read_assigned_out = _to_outdir(args.out_dir, args.read_assigned_out)
    args.cell_fastq_out = _to_outdir(args.out_dir, args.cell_fastq_out)
    args.cell_read_stats_out = _to_outdir(args.out_dir, args.cell_read_stats_out)

    print(f"[Strint] Collected {len(args.fastq_fns)} FASTQ files.")
    print(f"[Strint] Output directory: {args.out_dir}")

    return args


if __name__ == "__main__":
    parsed_args = set_parser()
    print("\n=== Parsed arguments ===")
    for key, value in vars(parsed_args).items():
        print(f"{key:24}: {value}")
