# Strint3 Single-Cell Upstream Pipeline

This repository contains a full upstream single-cell long-read processing workflow.

It supports two usage modes:

1. Run the barcode splitting and cell assignment step directly with `main.py`
2. Run the full upstream workflow with `sc_upstream_pipeline.wdl`
3. Run the full upstream workflow directly in the current conda environment with `run_all.sh`

The full workflow starts from full-length FASTQ and produces:

- cell-assigned reads
- tagged BAM files
- gene expression matrix
- isoform expression matrix
- RNA QC metrics
- sequencing saturation results

## Directory Layout

- [main.py](/Users/liyy/Documents/Documents%20-%20liyy的MacBook%20Air/Projects/CodeX_Space/Strint3/main.py): upstream barcode splitting and cell assignment
- [args_parser.py](/Users/liyy/Documents/Documents%20-%20liyy的MacBook%20Air/Projects/CodeX_Space/Strint3/args_parser.py): command-line arguments for `main.py`
- [utils.py](/Users/liyy/Documents/Documents%20-%20liyy的MacBook%20Air/Projects/CodeX_Space/Strint3/utils.py): helper functions used by `main.py`
- [scripts](/Users/liyy/Documents/Documents%20-%20liyy的MacBook%20Air/Projects/CodeX_Space/Strint3/scripts): downstream scripts for BAM tagging, gene assignment, UMI clustering, matrix generation, and QC
- [sc_upstream_pipeline.wdl](/Users/liyy/Documents/Documents%20-%20liyy的MacBook%20Air/Projects/CodeX_Space/Strint3/sc_upstream_pipeline.wdl): one-command workflow definition
- [sc_upstream_pipeline.inputs.json](/Users/liyy/Documents/Documents%20-%20liyy的MacBook%20Air/Projects/CodeX_Space/Strint3/sc_upstream_pipeline.inputs.json): example WDL inputs
- [run_all.sh](/Users/liyy/Documents/Documents%20-%20liyy的MacBook%20Air/Projects/CodeX_Space/Strint3/run_all.sh): no-container full pipeline runner using the current shell/conda environment

## Prerequisites

The WDL is designed to use software already installed on the server.

Required executables:

- `python3`
- `samtools`
- `minimap2`
- `bedtools`
- `docker`

Required Python packages:

- `pandas`
- `numpy`
- `pysam`
- `bioframe`
- `tqdm`
- `matplotlib`
- `seaborn`
- `networkx`
- `editdistance`
- `fast_edit_distance`

You also need a WDL runner, for example:

- `miniwdl`
- `cromwell`

## Build the Workflow Image

The WDL tasks run inside containers. Build the provided image first:

```bash
cd /home/liyy/2.project/C4_V3/Strint3
docker build -t strint3-sc:latest .
```

Quick checks:

```bash
docker run --rm strint3-sc:latest python3 --version
docker run --rm strint3-sc:latest samtools --version
docker run --rm strint3-sc:latest minimap2 --version
docker run --rm strint3-sc:latest bedtools --version
docker run --rm strint3-sc:latest python3 - <<'PY'
import pandas, numpy, pysam, bioframe, tqdm, matplotlib, seaborn, networkx, editdistance
print("python dependencies ok")
PY
```

## Step 1: Run Upstream Split Only

This step starts from full-length FASTQ and produces the key file:

- `read_assigned_cell.csv`

Example:

```bash
python3 main.py \
  /home/liyy/2.project/C4_V3/20260326/250F701400011_3/260F100027011.full-length.fq.gz \
  --full-bc-whitelist-3p /home/liyy/2.project/C4_V3/Strint3/whitelist_3p_3seg_rc.txt.gz \
  --full-bc-whitelist-5p /home/liyy/2.project/C4_V3/Strint3/whitelist_5p_3seg.txt.gz \
  --out_dir test_simple
```

Main outputs under `test_simple/`:

- `read_assigned_cell.csv`
- `cell_reads.fastq.gz`
- `matched_reads.fastq.gz`
- `ReadIDs_UMI_BC_clean.csv`
- `pair_counts_kept.csv`
- `cell_read_stats.csv`

Optional large intermediate outputs can be enabled with:

```bash
--save-intermediate
```

## Step 2: Run Full WDL Workflow

The WDL runs the full pipeline:

1. barcode splitting and cell assignment
2. read tag preparation
3. splice-aware alignment with minimap2
4. CB/UR tag injection
5. gene assignment
6. GN tag injection
7. UMI clustering
8. gene matrix generation
9. transcript assignment
10. isoform matrix generation
11. RNA QC
12. saturation analysis

### Update Inputs

Edit [sc_upstream_pipeline.inputs.json](/Users/liyy/Documents/Documents%20-%20liyy的MacBook%20Air/Projects/CodeX_Space/Strint3/sc_upstream_pipeline.inputs.json) and set the real server paths for:

- `fastq`
- `whitelist_3p`
- `whitelist_5p`
- `genome_fa`
- `junction_bed`
- `chrom_sizes`
- `gene_gtf`
- `isoform_gtf`
- `main_py`
- `prepare_read_tags_py`
- `add_cb_ur_tags_py`
- `assign_genes_py`
- `add_gene_tags_py`
- `cluster_umis_allbam_py`
- `cell_umi_gene_table_py`
- `gene_expression_py`
- `assign_transcripts_py`
- `isoform_expression_py`
- `rna_qc_metrics_py`
- `saturation_py`

Important:

- `junction_bed` should be the splice junction BED used by minimap2
- `chrom_sizes` should match the genome FASTA
- `gene_gtf` is used for gene assignment
- `isoform_gtf` is used for transcript assignment
- each `*_py` input should point to the actual script file on the server

### Run with miniwdl

```bash
cd /home/liyy/2.project/C4_V3/Strint3
miniwdl run --cfg miniwdl.cfg sc_upstream_pipeline.wdl -i sc_upstream_pipeline.inputs.json
```

### Run with Cromwell

```bash
java -jar cromwell.jar run sc_upstream_pipeline.wdl -i sc_upstream_pipeline.inputs.json
```

## Major Workflow Outputs

The workflow returns these main outputs:

- `read_assigned_cell.csv`
- `cell_reads.fastq.gz`
- `sample.filtered.sorted.bam`
- `sample.filtered.cb_ur.sorted.bam`
- `sample.filtered.cb_ur.gn.sorted.bam`
- `sample.filtered.tagged.sorted.bam`
- `sample.cell_umi_gene.tsv`
- `sample.gene_expression.tsv`
- `sample.read_transcript_assigns.tsv`
- `sample.isoform_expression.tsv`
- `rna_qc_metrics.tsv`
- `per_cell_qc.tsv`
- `rna_violin_plot.png`
- `res.tsv`
- `saturation_curves.png`

## Notes on Current Script Behavior

Two downstream scripts currently use fixed input filenames internally:

- [20260402/rna_qc_metrics.py](/Users/liyy/Documents/Documents%20-%20liyy的MacBook%20Air/Projects/CodeX_Space/Strint3/20260402/rna_qc_metrics.py)
- [20260402/Saturation.py](/Users/liyy/Documents/Documents%20-%20liyy的MacBook%20Air/Projects/CodeX_Space/Strint3/20260402/Saturation.py)

The WDL handles this by creating local symbolic links before running them.

## Recommended Server Checks

Before running the WDL, it is a good idea to verify:

```bash
python3 --version
samtools --version
minimap2 --version
bedtools --version
docker --version
```

And check that Python can import the required modules:

```bash
python3 - <<'PY'
import pandas, numpy, pysam, bioframe, tqdm, matplotlib, seaborn, networkx, editdistance
print("python dependencies ok")
PY
```

## Troubleshooting

If the workflow fails early:

- check all input paths in `sc_upstream_pipeline.inputs.json`
- check that each `*_py` script path is correct
- check that `main.py` can run standalone first
- check that `docker build -t strint3-sc:latest .` completed successfully
- check that `miniwdl.cfg` points to the same image tag you built

If alignment fails:

- confirm `genome_fa`, `junction_bed`, and `chrom_sizes` belong to the same reference
- confirm `minimap2` and `samtools` are on `PATH`

If gene or transcript assignment fails:

- confirm `gene_gtf` and `isoform_gtf` are valid GTF files
- confirm `bedtools` is installed
- confirm `bioframe` is available in the Python environment

If QC or saturation fails:

- check that `cell_umi_gene.tsv` was produced correctly
- check that matplotlib/seaborn are available

## Step 3: Run Full Pipeline Without WDL

If Docker or WDL execution is inconvenient on the server, you can run the full pipeline directly with [run_all.sh](/Users/liyy/Documents/Documents%20-%20liyy的MacBook%20Air/Projects/CodeX_Space/Strint3/run_all.sh).

This mode uses the currently active shell environment, so please activate your conda environment first.

Example:

```bash
cd /home/liyy/2.project/C4_V3/Strint3
conda activate Cyclone

bash run_all.sh \
  --fastq /home/liyy/2.project/C4_V3/20260326/250F701400011_3/260F100027011.full-length.fq.gz \
  --whitelist-3p /home/liyy/2.project/C4_V3/Strint3/whitelist_3p_3seg_rc.txt.gz \
  --whitelist-5p /home/liyy/2.project/C4_V3/Strint3/whitelist_5p_3seg.txt.gz \
  --genome-fa /home/liyy/1.data/REF/GRCm39/genome.fa \
  --junction-bed /home/liyy/1.data/REF/GRCm39/genes.bed \
  --chrom-sizes /home/liyy/1.data/REF/GRCm39/chrom_sizes.tsv \
  --gene-gtf /home/liyy/1.data/REF/GRCm39/genes.gtf \
  --isoform-gtf /home/liyy/1.data/REF/GRCm39/genes.gtf \
  --out-dir /home/liyy/2.project/C4_V3/Strint3/test_simple_all \
  --sample-id test_simple \
  --threads 32 \
  --cluster-threads 8
```

Outputs are organized under:

- `test_simple_all/upstream`
- `test_simple_all/alignment`
- `test_simple_all/matrix`
- `test_simple_all/qc`
- `test_simple_all/logs`

Key outputs:

- `upstream/read_assigned_cell.csv`
- `upstream/cell_reads.fastq.gz`
- `matrix/test_simple.gene_expression.tsv`
- `matrix/test_simple.isoform_expression.tsv`
- `qc/test_simple.rna_qc_metrics.tsv`
- `qc/test_simple.saturation.tsv`
