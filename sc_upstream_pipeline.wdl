version 1.1

workflow sc_upstream_pipeline {
  input {
    File fastq
    File whitelist_3p
    File whitelist_5p

    File genome_fa
    File junction_bed
    File chrom_sizes
    File gene_gtf
    File isoform_gtf

    File main_py
    File prepare_read_tags_py
    File add_cb_ur_tags_py
    File assign_genes_py
    File add_gene_tags_py
    File cluster_umis_allbam_py
    File cell_umi_gene_table_py
    File gene_expression_py
    File assign_transcripts_py
    File isoform_expression_py
    File rna_qc_metrics_py
    File saturation_py

    String sample_id = "sample"

    Int strint_threads = 32
    Int strint_batch_size = 100000
    Int strint_assign_batchsize = 10000
    Int minimap2_threads = 32
    Int cluster_umi_threads = 8

    Int exp_cells = 5000
    Int min_q = 7
    Int max_ed = 7
    Int empty_drop_min_ed = 5
    Int empty_drop_num = 2000
    Int pair_min = 10
    Float top1_alpha = 0.7
    Float dominance_min = 0.80
    Float drop_umiA_ratio_gt = 0.5

    Int gene_assign_mapq = 60
    Int gene_assign_chunk_size = 200000
    Int transcript_assign_mapq = 60
    Int transcript_assign_chunk_size = 200000
    Int ref_interval = 1000
    Int cell_gene_max_reads = 20000

    Boolean save_intermediate = false
    Boolean require_pass_both_ends = false
  }

  call RunStrintSplit {
    input:
      main_py = main_py,
      fastq = fastq,
      whitelist_3p = whitelist_3p,
      whitelist_5p = whitelist_5p,
      threads = strint_threads,
      batch_size = strint_batch_size,
      assign_batchsize = strint_assign_batchsize,
      exp_cells = exp_cells,
      min_q = min_q,
      max_ed = max_ed,
      empty_drop_min_ed = empty_drop_min_ed,
      empty_drop_num = empty_drop_num,
      pair_min = pair_min,
      top1_alpha = top1_alpha,
      dominance_min = dominance_min,
      drop_umiA_ratio_gt = drop_umiA_ratio_gt,
      save_intermediate = save_intermediate,
      require_pass_both_ends = require_pass_both_ends
  }

  call PrepareReadTags {
    input:
      prepare_read_tags_py = prepare_read_tags_py,
      read_assigned_cell_csv = RunStrintSplit.read_assigned_cell_csv,
      sample_id = sample_id
  }

  call MapCellReads {
    input:
      fastq = RunStrintSplit.cell_reads_fastq_gz,
      genome_fa = genome_fa,
      junction_bed = junction_bed,
      chrom_sizes = chrom_sizes,
      sample_id = sample_id,
      threads = minimap2_threads
  }

  call AddCBURTags {
    input:
      add_cb_ur_tags_py = add_cb_ur_tags_py,
      bam = MapCellReads.filtered_sorted_bam,
      tags_tsv = PrepareReadTags.read_tags_tsv,
      sample_id = sample_id
  }

  call BamToBed {
    input:
      bam = AddCBURTags.cb_ur_bam,
      sample_id = sample_id
  }

  call AssignGenes {
    input:
      assign_genes_py = assign_genes_py,
      bed = BamToBed.cb_ur_bed,
      gene_gtf = gene_gtf,
      sample_id = sample_id,
      mapq = gene_assign_mapq,
      chunk_size = gene_assign_chunk_size
  }

  call AddGeneTags {
    input:
      add_gene_tags_py = add_gene_tags_py,
      bam = AddCBURTags.cb_ur_bam,
      gene_assigns_tsv = AssignGenes.read_gene_assigns_tsv,
      sample_id = sample_id
  }

  call ClusterUmis {
    input:
      cluster_umis_allbam_py = cluster_umis_allbam_py,
      bam = AddGeneTags.cb_ur_gn_bam,
      sample_id = sample_id,
      ref_interval = ref_interval,
      cell_gene_max_reads = cell_gene_max_reads,
      threads = cluster_umi_threads
  }

  call CellUmiGeneTable {
    input:
      cell_umi_gene_table_py = cell_umi_gene_table_py,
      bam = ClusterUmis.tagged_bam,
      sample_id = sample_id
  }

  call GeneExpression {
    input:
      gene_expression_py = gene_expression_py,
      bam = ClusterUmis.tagged_bam,
      sample_id = sample_id
  }

  call AssignTranscripts {
    input:
      assign_transcripts_py = assign_transcripts_py,
      bam = ClusterUmis.tagged_bam,
      isoform_gtf = isoform_gtf,
      sample_id = sample_id,
      mapq = transcript_assign_mapq,
      chunk_size = transcript_assign_chunk_size
  }

  call IsoformExpression {
    input:
      isoform_expression_py = isoform_expression_py,
      bam = ClusterUmis.tagged_bam,
      transcript_assigns_tsv = AssignTranscripts.read_transcript_assigns_tsv,
      sample_id = sample_id
  }

  call RNAQCMetrics {
    input:
      rna_qc_metrics_py = rna_qc_metrics_py,
      cell_umi_gene_tsv = CellUmiGeneTable.cell_umi_gene_tsv,
      mapped_bam = MapCellReads.filtered_sorted_bam,
      mapped_bam_bai = MapCellReads.filtered_sorted_bam_bai
  }

  call Saturation {
    input:
      saturation_py = saturation_py,
      cell_umi_gene_tsv = CellUmiGeneTable.cell_umi_gene_tsv
  }

  output {
    File read_assigned_cell_csv = RunStrintSplit.read_assigned_cell_csv
    File cell_reads_fastq_gz = RunStrintSplit.cell_reads_fastq_gz
    File read_tags_tsv = PrepareReadTags.read_tags_tsv
    File filtered_sorted_bam = MapCellReads.filtered_sorted_bam
    File filtered_sorted_bam_bai = MapCellReads.filtered_sorted_bam_bai
    File cb_ur_bam = AddCBURTags.cb_ur_bam
    File cb_ur_bam_bai = AddCBURTags.cb_ur_bam_bai
    File cb_ur_bed = BamToBed.cb_ur_bed
    File read_gene_assigns_tsv = AssignGenes.read_gene_assigns_tsv
    File cb_ur_gn_bam = AddGeneTags.cb_ur_gn_bam
    File cb_ur_gn_bam_bai = AddGeneTags.cb_ur_gn_bam_bai
    File tagged_bam = ClusterUmis.tagged_bam
    File tagged_bam_bai = ClusterUmis.tagged_bam_bai
    File cell_umi_gene_tsv = CellUmiGeneTable.cell_umi_gene_tsv
    File gene_expression_tsv = GeneExpression.gene_expression_tsv
    File read_transcript_assigns_tsv = AssignTranscripts.read_transcript_assigns_tsv
    File isoform_expression_tsv = IsoformExpression.isoform_expression_tsv
    File rna_qc_metrics_tsv = RNAQCMetrics.rna_qc_metrics_tsv
    File per_cell_qc_tsv = RNAQCMetrics.per_cell_qc_tsv
    File rna_violin_plot_png = RNAQCMetrics.rna_violin_plot_png
    File saturation_table_tsv = Saturation.saturation_table_tsv
    File saturation_curves_png = Saturation.saturation_curves_png
  }
}

task RunStrintSplit {
  input {
    File main_py
    File fastq
    File whitelist_3p
    File whitelist_5p
    Int threads
    Int batch_size
    Int assign_batchsize
    Int exp_cells
    Int min_q
    Int max_ed
    Int empty_drop_min_ed
    Int empty_drop_num
    Int pair_min
    Float top1_alpha
    Float dominance_min
    Float drop_umiA_ratio_gt
    Boolean save_intermediate
    Boolean require_pass_both_ends
  }

  command <<<
    set -euo pipefail
    mkdir -p strint_out
    python3 "~{main_py}" "~{fastq}" \
      --full-bc-whitelist-3p "~{whitelist_3p}" \
      --full-bc-whitelist-5p "~{whitelist_5p}" \
      --out_dir strint_out \
      --threads ~{threads} \
      --batch_size ~{batch_size} \
      --assign_batchsize ~{assign_batchsize} \
      --exp_cells ~{exp_cells} \
      --minQ ~{min_q} \
      --max_ed ~{max_ed} \
      --DEFAULT_EMPTY_DROP_MIN_ED ~{empty_drop_min_ed} \
      --DEFAULT_EMPTY_DROP_NUM ~{empty_drop_num} \
      --PAIR_MIN ~{pair_min} \
      --TOP1_ALPHA ~{top1_alpha} \
      --dominance_min ~{dominance_min} \
      --drop_umiA_ratio_gt ~{drop_umiA_ratio_gt} \
      ~{if save_intermediate then "--save-intermediate" else ""} \
      ~{if require_pass_both_ends then "--require_pass_both_ends" else ""}
  >>>

  output {
    File read_assigned_cell_csv = "strint_out/read_assigned_cell.csv"
    File cell_reads_fastq_gz = "strint_out/cell_reads.fastq.gz"
    File matched_reads_fastq_gz = "strint_out/matched_reads.fastq.gz"
    File clean_reads_csv = "strint_out/ReadIDs_UMI_BC_clean.csv"
    File pair_counts_csv = "strint_out/pair_counts_kept.csv"
    File cell_read_stats_csv = "strint_out/cell_read_stats.csv"
  }

  runtime {
    cpu: threads
    memory: "32G"
  }
}

task PrepareReadTags {
  input {
    File prepare_read_tags_py
    File read_assigned_cell_csv
    String sample_id
  }

  command <<<
    set -euo pipefail
    python3 "~{prepare_read_tags_py}" \
      --input "~{read_assigned_cell_csv}" \
      --output "~{sample_id}.read_tags.tsv"
  >>>

  output {
    File read_tags_tsv = "~{sample_id}.read_tags.tsv"
  }
}

task MapCellReads {
  input {
    File fastq
    File genome_fa
    File junction_bed
    File chrom_sizes
    String sample_id
    Int threads
  }

  command <<<
    set -euo pipefail
    minimap2 -ax splice -uf --MD -t ~{threads} \
      --junc-bed "~{junction_bed}" \
      --secondary=no \
      "~{genome_fa}" "~{fastq}" \
    | samtools view --no-PG -b -t "~{chrom_sizes}" - \
    | samtools sort --no-PG -@ ~{threads} -o "~{sample_id}.filtered.sorted.bam" -
    samtools index "~{sample_id}.filtered.sorted.bam"
  >>>

  output {
    File filtered_sorted_bam = "~{sample_id}.filtered.sorted.bam"
    File filtered_sorted_bam_bai = "~{sample_id}.filtered.sorted.bam.bai"
  }

  runtime {
    cpu: threads
    memory: "32G"
  }
}

task AddCBURTags {
  input {
    File add_cb_ur_tags_py
    File bam
    File tags_tsv
    String sample_id
  }

  command <<<
    set -euo pipefail
    python3 "~{add_cb_ur_tags_py}" \
      --bam "~{bam}" \
      --tags "~{tags_tsv}" \
      --output "~{sample_id}.filtered.cb_ur.sorted.bam"
  >>>

  output {
    File cb_ur_bam = "~{sample_id}.filtered.cb_ur.sorted.bam"
    File cb_ur_bam_bai = "~{sample_id}.filtered.cb_ur.sorted.bam.bai"
  }
}

task BamToBed {
  input {
    File bam
    String sample_id
  }

  command <<<
    set -euo pipefail
    bedtools bamtobed -i "~{bam}" > "~{sample_id}.filtered.cb_ur.bed"
  >>>

  output {
    File cb_ur_bed = "~{sample_id}.filtered.cb_ur.bed"
  }
}

task AssignGenes {
  input {
    File assign_genes_py
    File bed
    File gene_gtf
    String sample_id
    Int mapq
    Int chunk_size
  }

  command <<<
    set -euo pipefail
    python3 "~{assign_genes_py}" \
      --output "~{sample_id}.filtered.read_gene_assigns.tsv" \
      --mapq ~{mapq} \
      --chunk_size ~{chunk_size} \
      "~{bed}" "~{gene_gtf}"
  >>>

  output {
    File read_gene_assigns_tsv = "~{sample_id}.filtered.read_gene_assigns.tsv"
  }
}

task AddGeneTags {
  input {
    File add_gene_tags_py
    File bam
    File gene_assigns_tsv
    String sample_id
  }

  command <<<
    set -euo pipefail
    python3 "~{add_gene_tags_py}" \
      --output "~{sample_id}.filtered.cb_ur.gn.sorted.bam" \
      "~{bam}" "~{gene_assigns_tsv}"
    samtools index "~{sample_id}.filtered.cb_ur.gn.sorted.bam"
  >>>

  output {
    File cb_ur_gn_bam = "~{sample_id}.filtered.cb_ur.gn.sorted.bam"
    File cb_ur_gn_bam_bai = "~{sample_id}.filtered.cb_ur.gn.sorted.bam.bai"
  }
}

task ClusterUmis {
  input {
    File cluster_umis_allbam_py
    File bam
    String sample_id
    Int ref_interval
    Int cell_gene_max_reads
    Int threads
  }

  command <<<
    set -euo pipefail
    python3 "~{cluster_umis_allbam_py}" \
      --output "~{sample_id}.filtered.tagged.sorted.bam" \
      --ref_interval ~{ref_interval} \
      --cell_gene_max_reads ~{cell_gene_max_reads} \
      --threads ~{threads} \
      "~{bam}"
  >>>

  output {
    File tagged_bam = "~{sample_id}.filtered.tagged.sorted.bam"
    File tagged_bam_bai = "~{sample_id}.filtered.tagged.sorted.bam.bai"
  }

  runtime {
    cpu: threads
    memory: "24G"
  }
}

task CellUmiGeneTable {
  input {
    File cell_umi_gene_table_py
    File bam
    String sample_id
  }

  command <<<
    set -euo pipefail
    python3 "~{cell_umi_gene_table_py}" \
      --output "~{sample_id}.cell_umi_gene.tsv" \
      "~{bam}"
  >>>

  output {
    File cell_umi_gene_tsv = "~{sample_id}.cell_umi_gene.tsv"
  }
}

task GeneExpression {
  input {
    File gene_expression_py
    File bam
    String sample_id
  }

  command <<<
    set -euo pipefail
    python3 "~{gene_expression_py}" \
      --output "~{sample_id}.gene_expression.tsv" \
      "~{bam}"
  >>>

  output {
    File gene_expression_tsv = "~{sample_id}.gene_expression.tsv"
  }
}

task AssignTranscripts {
  input {
    File assign_transcripts_py
    File bam
    File isoform_gtf
    String sample_id
    Int mapq
    Int chunk_size
  }

  command <<<
    set -euo pipefail
    python3 "~{assign_transcripts_py}" \
      --output "~{sample_id}.read_transcript_assigns.tsv" \
      --mapq ~{mapq} \
      --chunk_size ~{chunk_size} \
      "~{bam}" "~{isoform_gtf}"
  >>>

  output {
    File read_transcript_assigns_tsv = "~{sample_id}.read_transcript_assigns.tsv"
  }
}

task IsoformExpression {
  input {
    File isoform_expression_py
    File bam
    File transcript_assigns_tsv
    String sample_id
  }

  command <<<
    set -euo pipefail
    python3 "~{isoform_expression_py}" \
      --output "~{sample_id}.isoform_expression.tsv" \
      "~{bam}" "~{transcript_assigns_tsv}"
  >>>

  output {
    File isoform_expression_tsv = "~{sample_id}.isoform_expression.tsv"
  }
}

task RNAQCMetrics {
  input {
    File rna_qc_metrics_py
    File cell_umi_gene_tsv
    File mapped_bam
    File mapped_bam_bai
  }

  command <<<
    set -euo pipefail
    ln -s "~{cell_umi_gene_tsv}" cell_umi_gene.tsv
    ln -s "~{mapped_bam}" filtered.sorted.bam
    ln -s "~{mapped_bam_bai}" filtered.sorted.bam.bai
    python3 "~{rna_qc_metrics_py}"
  >>>

  output {
    File rna_qc_metrics_tsv = "rna_qc_metrics.tsv"
    File per_cell_qc_tsv = "per_cell_qc.tsv"
    File rna_violin_plot_png = "rna_violin_plot.png"
  }
}

task Saturation {
  input {
    File saturation_py
    File cell_umi_gene_tsv
  }

  command <<<
    set -euo pipefail
    ln -s "~{cell_umi_gene_tsv}" cell_umi_gene.tsv
    python3 "~{saturation_py}" | tee res.tsv
  >>>

  output {
    File saturation_table_tsv = "res.tsv"
    File saturation_curves_png = "saturation_curves.png"
  }
}
