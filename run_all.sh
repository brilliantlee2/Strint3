#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash run_all.sh \
    --fastq /path/to/full-length.fq.gz \
    --whitelist-3p /path/to/whitelist_3p.txt.gz \
    --whitelist-5p /path/to/whitelist_5p.txt.gz \
    --genome-fa /path/to/genome.fa \
    --junction-bed /path/to/genes.bed \
    --chrom-sizes /path/to/chrom_sizes.tsv \
    --gene-gtf /path/to/genes.gtf \
    --isoform-gtf /path/to/genes.gtf \
    --out-dir /path/to/output_dir \
    [--sample-id sample] \
    [--threads 32] \
    [--cluster-threads 8] \
    [--exp-cells 5000] \
    [--min-q 7] \
    [--max-ed 7] \
    [--pair-min 10] \
    [--top1-alpha 0.7] \
    [--dominance-min 0.8] \
    [--drop-umiA-ratio-gt 0.5] \
    [--gene-assign-mapq 60] \
    [--gene-assign-chunk-size 200000] \
    [--transcript-assign-mapq 60] \
    [--transcript-assign-chunk-size 200000] \
    [--ref-interval 1000] \
    [--cell-gene-max-reads 20000] \
    [--save-intermediate] \
    [--require-pass-both-ends]

Notes:
  1. This script uses the currently active shell environment.
  2. Please activate your conda environment before running.
  3. Outputs are written under:
     <out-dir>/upstream
     <out-dir>/alignment
     <out-dir>/matrix
     <out-dir>/qc
EOF
}

require_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "[ERROR] Command not found: $cmd" >&2
    exit 1
  fi
}

require_file() {
  local path="$1"
  local label="$2"
  if [[ ! -f "$path" ]]; then
    echo "[ERROR] Missing ${label}: $path" >&2
    exit 1
  fi
}

log() {
  echo "[run_all] $*"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOWNSTREAM_DIR="${SCRIPT_DIR}/scripts"

FASTQ=""
WHITELIST_3P=""
WHITELIST_5P=""
GENOME_FA=""
JUNCTION_BED=""
CHROM_SIZES=""
GENE_GTF=""
ISOFORM_GTF=""
OUT_DIR=""
SAMPLE_ID="sample"
THREADS=32
CLUSTER_THREADS=8
EXP_CELLS=5000
MIN_Q=7
MAX_ED=7
PAIR_MIN=10
TOP1_ALPHA=0.7
DOMINANCE_MIN=0.8
DROP_UMIA_RATIO_GT=0.5
GENE_ASSIGN_MAPQ=60
GENE_ASSIGN_CHUNK_SIZE=200000
TRANSCRIPT_ASSIGN_MAPQ=60
TRANSCRIPT_ASSIGN_CHUNK_SIZE=200000
REF_INTERVAL=1000
CELL_GENE_MAX_READS=20000
SAVE_INTERMEDIATE=0
REQUIRE_PASS_BOTH_ENDS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq) FASTQ="$2"; shift 2 ;;
    --whitelist-3p) WHITELIST_3P="$2"; shift 2 ;;
    --whitelist-5p) WHITELIST_5P="$2"; shift 2 ;;
    --genome-fa) GENOME_FA="$2"; shift 2 ;;
    --junction-bed) JUNCTION_BED="$2"; shift 2 ;;
    --chrom-sizes) CHROM_SIZES="$2"; shift 2 ;;
    --gene-gtf) GENE_GTF="$2"; shift 2 ;;
    --isoform-gtf) ISOFORM_GTF="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --cluster-threads) CLUSTER_THREADS="$2"; shift 2 ;;
    --exp-cells) EXP_CELLS="$2"; shift 2 ;;
    --min-q) MIN_Q="$2"; shift 2 ;;
    --max-ed) MAX_ED="$2"; shift 2 ;;
    --pair-min) PAIR_MIN="$2"; shift 2 ;;
    --top1-alpha) TOP1_ALPHA="$2"; shift 2 ;;
    --dominance-min) DOMINANCE_MIN="$2"; shift 2 ;;
    --drop-umiA-ratio-gt) DROP_UMIA_RATIO_GT="$2"; shift 2 ;;
    --gene-assign-mapq) GENE_ASSIGN_MAPQ="$2"; shift 2 ;;
    --gene-assign-chunk-size) GENE_ASSIGN_CHUNK_SIZE="$2"; shift 2 ;;
    --transcript-assign-mapq) TRANSCRIPT_ASSIGN_MAPQ="$2"; shift 2 ;;
    --transcript-assign-chunk-size) TRANSCRIPT_ASSIGN_CHUNK_SIZE="$2"; shift 2 ;;
    --ref-interval) REF_INTERVAL="$2"; shift 2 ;;
    --cell-gene-max-reads) CELL_GENE_MAX_READS="$2"; shift 2 ;;
    --save-intermediate) SAVE_INTERMEDIATE=1; shift ;;
    --require-pass-both-ends) REQUIRE_PASS_BOTH_ENDS=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "[ERROR] Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "${FASTQ}" || -z "${WHITELIST_3P}" || -z "${WHITELIST_5P}" || -z "${GENOME_FA}" || -z "${JUNCTION_BED}" || -z "${CHROM_SIZES}" || -z "${GENE_GTF}" || -z "${ISOFORM_GTF}" || -z "${OUT_DIR}" ]]; then
  echo "[ERROR] Missing required arguments." >&2
  usage
  exit 1
fi

require_cmd python3
require_cmd samtools
require_cmd minimap2
require_cmd bedtools

require_file "${FASTQ}" "FASTQ"
require_file "${WHITELIST_3P}" "3p whitelist"
require_file "${WHITELIST_5P}" "5p whitelist"
require_file "${GENOME_FA}" "genome FASTA"
require_file "${JUNCTION_BED}" "junction BED"
require_file "${CHROM_SIZES}" "chrom sizes"
require_file "${GENE_GTF}" "gene GTF"
require_file "${ISOFORM_GTF}" "isoform GTF"
require_file "${SCRIPT_DIR}/main.py" "main.py"
require_file "${DOWNSTREAM_DIR}/prepare_read_tags.py" "prepare_read_tags.py"
require_file "${DOWNSTREAM_DIR}/add_cb_ur_tags.py" "add_cb_ur_tags.py"
require_file "${DOWNSTREAM_DIR}/assign_genes.py" "assign_genes.py"
require_file "${DOWNSTREAM_DIR}/add_gene_tags.py" "add_gene_tags.py"
require_file "${DOWNSTREAM_DIR}/cluster_umis_allbam.py" "cluster_umis_allbam.py"
require_file "${DOWNSTREAM_DIR}/cell_umi_gene_table.py" "cell_umi_gene_table.py"
require_file "${DOWNSTREAM_DIR}/gene_expression.py" "gene_expression.py"
require_file "${DOWNSTREAM_DIR}/assign_transcripts.py" "assign_transcripts.py"
require_file "${DOWNSTREAM_DIR}/isoform_expression.py" "isoform_expression.py"
require_file "${DOWNSTREAM_DIR}/rna_qc_metrics.py" "rna_qc_metrics.py"
require_file "${DOWNSTREAM_DIR}/Saturation.py" "Saturation.py"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "[WARN] No active conda environment detected. Continuing with current shell PATH." >&2
else
  log "Using conda environment: ${CONDA_PREFIX}"
fi

OUT_DIR="$(cd "$(dirname "${OUT_DIR}")" 2>/dev/null && pwd)/$(basename "${OUT_DIR}")"
UPSTREAM_DIR="${OUT_DIR}/upstream"
ALIGN_DIR="${OUT_DIR}/alignment"
MATRIX_DIR="${OUT_DIR}/matrix"
QC_DIR="${OUT_DIR}/qc"
LOG_DIR="${OUT_DIR}/logs"

mkdir -p "${UPSTREAM_DIR}" "${ALIGN_DIR}" "${MATRIX_DIR}" "${QC_DIR}" "${LOG_DIR}"

UPSTREAM_LOG="${LOG_DIR}/01_upstream.log"
ALIGN_LOG="${LOG_DIR}/02_alignment.log"
GENE_LOG="${LOG_DIR}/03_gene.log"
ISOFORM_LOG="${LOG_DIR}/04_isoform.log"
QC_LOG="${LOG_DIR}/05_qc.log"

log "Step 1/5: running upstream barcode split and cell assignment"
MAIN_ARGS=(
  "${FASTQ}"
  --full-bc-whitelist-3p "${WHITELIST_3P}"
  --full-bc-whitelist-5p "${WHITELIST_5P}"
  --out_dir "${UPSTREAM_DIR}"
  --threads "${THREADS}"
  --batch_size 100000
  --assign_batchsize 10000
  --exp_cells "${EXP_CELLS}"
  --minQ "${MIN_Q}"
  --max_ed "${MAX_ED}"
  --PAIR_MIN "${PAIR_MIN}"
  --TOP1_ALPHA "${TOP1_ALPHA}"
  --dominance_min "${DOMINANCE_MIN}"
  --drop_umiA_ratio_gt "${DROP_UMIA_RATIO_GT}"
)
if [[ "${SAVE_INTERMEDIATE}" -eq 1 ]]; then
  MAIN_ARGS+=(--save-intermediate)
fi
if [[ "${REQUIRE_PASS_BOTH_ENDS}" -eq 1 ]]; then
  MAIN_ARGS+=(--require_pass_both_ends)
fi
python3 "${SCRIPT_DIR}/main.py" "${MAIN_ARGS[@]}" 2>&1 | tee "${UPSTREAM_LOG}"

READ_ASSIGNED_CELL="${UPSTREAM_DIR}/read_assigned_cell.csv"
CELL_READS_FASTQ="${UPSTREAM_DIR}/cell_reads.fastq.gz"
MATCHED_READS_FASTQ="${UPSTREAM_DIR}/matched_reads.fastq.gz"
require_file "${READ_ASSIGNED_CELL}" "read_assigned_cell.csv"
require_file "${CELL_READS_FASTQ}" "cell_reads.fastq.gz"
require_file "${MATCHED_READS_FASTQ}" "matched_reads.fastq.gz"

log "Step 2/5: running alignment and read tagging"
pushd "${ALIGN_DIR}" >/dev/null

python3 "${DOWNSTREAM_DIR}/prepare_read_tags.py" \
  --input "${READ_ASSIGNED_CELL}" \
  --output "${SAMPLE_ID}.read_tags.tsv" 2>&1 | tee "${ALIGN_LOG}"

minimap2 -ax splice -uf --MD -t "${THREADS}" \
  --junc-bed "${JUNCTION_BED}" \
  --secondary=no \
  "${GENOME_FA}" "${CELL_READS_FASTQ}" \
| samtools view --no-PG -b -t "${CHROM_SIZES}" - \
| samtools sort --no-PG -@ "${THREADS}" -o "${SAMPLE_ID}.filtered.sorted.bam" -

samtools index "${SAMPLE_ID}.filtered.sorted.bam"

python3 "${DOWNSTREAM_DIR}/add_cb_ur_tags.py" \
  --bam "${SAMPLE_ID}.filtered.sorted.bam" \
  --tags "${SAMPLE_ID}.read_tags.tsv" \
  --output "${SAMPLE_ID}.filtered.cb_ur.sorted.bam" 2>&1 | tee -a "${ALIGN_LOG}"

bedtools bamtobed -i "${SAMPLE_ID}.filtered.cb_ur.sorted.bam" > "${SAMPLE_ID}.filtered.cb_ur.bed"

popd >/dev/null

log "Step 3/5: gene assignment, UMI clustering, and gene matrix generation"
pushd "${MATRIX_DIR}" >/dev/null

python3 "${DOWNSTREAM_DIR}/assign_genes.py" \
  --output "${SAMPLE_ID}.filtered.read_gene_assigns.tsv" \
  --mapq "${GENE_ASSIGN_MAPQ}" \
  --chunk_size "${GENE_ASSIGN_CHUNK_SIZE}" \
  "${ALIGN_DIR}/${SAMPLE_ID}.filtered.cb_ur.bed" \
  "${GENE_GTF}" 2>&1 | tee "${GENE_LOG}"

python3 "${DOWNSTREAM_DIR}/add_gene_tags.py" \
  --output "${SAMPLE_ID}.filtered.cb_ur.gn.sorted.bam" \
  "${ALIGN_DIR}/${SAMPLE_ID}.filtered.cb_ur.sorted.bam" \
  "${SAMPLE_ID}.filtered.read_gene_assigns.tsv" 2>&1 | tee -a "${GENE_LOG}"

samtools index "${SAMPLE_ID}.filtered.cb_ur.gn.sorted.bam"

python3 "${DOWNSTREAM_DIR}/cluster_umis_allbam.py" \
  --output "${SAMPLE_ID}.filtered.tagged.sorted.bam" \
  --ref_interval "${REF_INTERVAL}" \
  --cell_gene_max_reads "${CELL_GENE_MAX_READS}" \
  --threads "${CLUSTER_THREADS}" \
  "${SAMPLE_ID}.filtered.cb_ur.gn.sorted.bam" 2>&1 | tee -a "${GENE_LOG}"

python3 "${DOWNSTREAM_DIR}/cell_umi_gene_table.py" \
  --output "${SAMPLE_ID}.cell_umi_gene.tsv" \
  "${SAMPLE_ID}.filtered.tagged.sorted.bam" 2>&1 | tee -a "${GENE_LOG}"

python3 "${DOWNSTREAM_DIR}/gene_expression.py" \
  --output "${SAMPLE_ID}.gene_expression.tsv" \
  "${SAMPLE_ID}.filtered.tagged.sorted.bam" 2>&1 | tee -a "${GENE_LOG}"

popd >/dev/null

log "Step 4/5: transcript assignment and isoform matrix generation"
pushd "${MATRIX_DIR}" >/dev/null

python3 "${DOWNSTREAM_DIR}/assign_transcripts.py" \
  --output "${SAMPLE_ID}.read_transcript_assigns.tsv" \
  --mapq "${TRANSCRIPT_ASSIGN_MAPQ}" \
  --chunk_size "${TRANSCRIPT_ASSIGN_CHUNK_SIZE}" \
  "${SAMPLE_ID}.filtered.tagged.sorted.bam" \
  "${ISOFORM_GTF}" 2>&1 | tee "${ISOFORM_LOG}"

python3 "${DOWNSTREAM_DIR}/isoform_expression.py" \
  --output "${SAMPLE_ID}.isoform_expression.tsv" \
  "${SAMPLE_ID}.filtered.tagged.sorted.bam" \
  "${SAMPLE_ID}.read_transcript_assigns.tsv" 2>&1 | tee -a "${ISOFORM_LOG}"

popd >/dev/null

log "Step 5/5: RNA QC and saturation analysis"
pushd "${QC_DIR}" >/dev/null

ln -sf "${MATRIX_DIR}/${SAMPLE_ID}.cell_umi_gene.tsv" cell_umi_gene.tsv
ln -sf "${ALIGN_DIR}/${SAMPLE_ID}.filtered.sorted.bam" filtered.sorted.bam
ln -sf "${ALIGN_DIR}/${SAMPLE_ID}.filtered.sorted.bam.bai" filtered.sorted.bam.bai

python3 "${DOWNSTREAM_DIR}/rna_qc_metrics.py" 2>&1 | tee "${QC_LOG}"
mv -f rna_qc_metrics.tsv "${SAMPLE_ID}.rna_qc_metrics.tsv"
mv -f per_cell_qc.tsv "${SAMPLE_ID}.per_cell_qc.tsv"
mv -f rna_violin_plot.png "${SAMPLE_ID}.rna_violin_plot.png"

python3 "${DOWNSTREAM_DIR}/Saturation.py" 2>&1 | tee "${SAMPLE_ID}.saturation.tsv" | tee -a "${QC_LOG}"
mv -f saturation_curves.png "${SAMPLE_ID}.saturation_curves.png"

popd >/dev/null

log "Pipeline completed successfully."
log "Key outputs:"
echo "  upstream read assignments : ${READ_ASSIGNED_CELL}"
echo "  cell reads fastq          : ${CELL_READS_FASTQ}"
echo "  tagged bam               : ${MATRIX_DIR}/${SAMPLE_ID}.filtered.tagged.sorted.bam"
echo "  gene expression          : ${MATRIX_DIR}/${SAMPLE_ID}.gene_expression.tsv"
echo "  isoform expression       : ${MATRIX_DIR}/${SAMPLE_ID}.isoform_expression.tsv"
echo "  RNA QC                   : ${QC_DIR}/${SAMPLE_ID}.rna_qc_metrics.tsv"
echo "  saturation table         : ${QC_DIR}/${SAMPLE_ID}.saturation.tsv"
