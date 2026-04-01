#!/bin/bash
fq=/home/liyy/2.project/C4_V3/20260204/Debug/20260302/cell_reads.fastq.gz
genome=/home/liyy/1.data/REF/GRCH38/genome.fa
bed=/home/liyy/1.data/REF/GRCH38/genes.bed
chrome_size=/home/liyy/1.data/REF/GRCH38/genes.bed/chrom_sizes.tsv
minimap2 -ax splice -uf --MD -t 32 \
  --junc-bed ${bed} \
  --secondary=no \
  ${genome} ${fq} > tmp.sam

samtools view --no-PG tmp.sam \
  -t ${chrome_size} \
  -o tmp.unsort.bam

samtools sort --no-PG tmp.unsort.bam -o filtered.sorted.bam
samtools index filtered.sorted.bam


#前期reference处理工作
#samtools faidx genome.fa
#paftools.js gff2bed -j genes.gtf > genes.bed
#cut -f1,2 genome.fa.fai | sort -V > chrom_sizes.tsv
