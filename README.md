# Strint3
A prerocess tools for demultiplexing barcodes in MGI full-length single cell RNA-seq.



脚本  输入  输出
1. prepare_read_tags.py  read_assigned_cell.csv  read_tags.tsv  python prepare_read_tags.py --input read_assigned_cell.csv  --output read_tags.tsv

2. add_cb_ur_tags.py  filtered.sorted.bam  read_tags.tsv python ../add_cb_ur_tags.py  --bam filtered.sorted.bam --tags read_tags.tsv  --output filtered.cb_ur.sorted.bam

3. 比对参考辅助文件准备 bedtools bamtobed -i filtered.cb_ur.sorted.bam > filtered.cb_ur.bed &

4. scripts/assign_genes.py python  ../assign_genes.py  --output filtered.read_gene_assigns.tsv filtered.cb_ur.bed /home/liyy/1.data/REF/GRCH38/genes.gtf

5. scripts/add_gene_tags.py python ../add_gene_tags.py  --output filtered.cb_ur.gn.sorted.bam  filtered.cb_ur.sorted.bam filtered.read_gene_assigns.tsv 

6. scripts/cluster_umis_allbam.py python ../cluster_umis_allbam.py  --output filtered.tagged.sorted.bam filtered.cb_ur.gn.sorted.bam

7. scripts/cell_umi_gene_table.py   python ../cell_umi_gene_table.py  filtered.tagged.sorted.bam --output cell_umi_gene.tsv

8. scripts/gene_expression.py  python ../gene_expression.py  --output gene_expression.tsv filtered.tagged.sorted.bam


python3 main.py  --full-bc-whitelist-3p /home/liyy/2.project/C4_V3/Strint3/SampleData2/  --full-bc-whitelist-5p /path/to/whitelist_5p.txt.gz



python3 main.py /home/liyy/2.project/C4_V3/20260326/250F701400011_3/260F100027011.full-length.fq.gz \
  --full-bc-whitelist-3p /home/liyy/2.project/C4_V3/Strint3/whitelist_3p_3seg_rc.txt.gz \
  --full-bc-whitelist-5p /home/liyy/2.project/C4_V3/Strint3/whitelist_5p_3seg.txt.gz \
  --out_dir test_simple
