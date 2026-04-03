#!/bin/bash
fq=../cell_reads.fastq.gz
genome=/home/liyy/1.data/REF/GRCm39/genome.fa
bed=/home/liyy/1.data/REF/GRCm39/genes.bed
chrome_size=/home/liyy/1.data/REF/GRCm39/chrom_sizes.tsv
threads=96

minimap2 -ax splice -uf --MD -t ${threads} \
  --junc-bed ${bed} \
  --secondary=no \
  ${genome} ${fq} \
| samtools view --no-PG -b -t ${chrome_size} - \
| samtools sort --no-PG -@ ${threads} -o filtered.sorted.bam -

samtools index filtered.sorted.bam


#前期reference处理工作
#samtools faidx genome.fa
#paftools.js gff2bed -j genes.gtf > genes.bed
#cut -f1,2 genome.fa.fai | sort -V > chrom_sizes.tsv
