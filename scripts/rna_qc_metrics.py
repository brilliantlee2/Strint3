import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# =========================
# Input files
# =========================
CELL_UMI_GENE_TSV = "cell_umi_gene.tsv"
TOTAL_BAM = "filtered.sorted.bam"

# =========================
# Helpers
# =========================
def count_bam_reads(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    stats = bam.get_index_statistics()
    n_reads = int(sum(contig.mapped for contig in stats))
    bam.close()
    return n_reads


def is_mito_gene(gene_name):
    gene_name = str(gene_name)
    return gene_name.startswith("MT-") or gene_name.startswith("mt-") or gene_name.startswith("Mt-")


# =========================
# Load tables
# =========================
df = pd.read_csv(CELL_UMI_GENE_TSV, sep="\t")

# 这里假设 cell_umi_gene.tsv 这几列存在：
# read_id, gene, barcode, umi
required_cols = {"read_id", "gene", "barcode", "umi"}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns in {CELL_UMI_GENE_TSV}: {missing}")

# =========================
# Basic counts
# =========================
estimated_cells = df["barcode"].nunique()
total_reads = count_bam_reads(TOTAL_BAM)
cell_associated_reads = df["read_id"].nunique()

# =========================
# Per-cell stats
# =========================
reads_per_cell = df.groupby("barcode")["read_id"].nunique()
umis_per_cell = df.groupby("barcode")["umi"].nunique()
genes_per_cell = df.groupby("barcode")["gene"].nunique()

# mito percent by reads
df["is_mito"] = df["gene"].apply(is_mito_gene)
mito_reads_per_cell = df.groupby("barcode")["is_mito"].sum()
mito_percent_per_cell = (mito_reads_per_cell / reads_per_cell) * 100

# =========================
# Summary metrics
# =========================
metrics = {
    "Estimated number of cells": estimated_cells,
    "Mean reads per cell": total_reads / estimated_cells if estimated_cells > 0 else 0,
    "Mean cell-associated reads per cell": reads_per_cell.mean(),
    "Mean UMI counts per cell": umis_per_cell.mean(),
    "Median UMI counts per cell": umis_per_cell.median(),
    "Mean Genes per cell": genes_per_cell.mean(),
    "Median Genes per cell": genes_per_cell.median(),
    "Total genes detected": df["gene"].nunique(),
    "Fraction reads in cells": cell_associated_reads / total_reads if total_reads > 0 else 0,
}

metrics_df = pd.DataFrame(
    {"Metric": list(metrics.keys()), "Value": list(metrics.values())}
)

# 格式化输出
formatted_metrics_df = metrics_df.copy()
formatted_metrics_df["Value"] = formatted_metrics_df.apply(
    lambda row: f"{row['Value']:.2%}" if row["Metric"] == "Fraction reads in cells"
    else f"{row['Value']:,.2f}" if isinstance(row["Value"], float)
    else f"{row['Value']:,}",
    axis=1,
)

formatted_metrics_df.to_csv("rna_qc_metrics.tsv", sep="\t", index=False)

print("\nRNA QC Metrics")
print(formatted_metrics_df.to_string(index=False))

# =========================
# Save per-cell QC table
# =========================
per_cell_qc = pd.DataFrame({
    "barcode": reads_per_cell.index,
    "reads": reads_per_cell.values,
    "umis": umis_per_cell.reindex(reads_per_cell.index).values,
    "genes": genes_per_cell.reindex(reads_per_cell.index).values,
    "mito_percent": mito_percent_per_cell.reindex(reads_per_cell.index).values,
})

per_cell_qc.to_csv("per_cell_qc.tsv", sep="\t", index=False)

# =========================
# Violin plot
# =========================
plot_df = per_cell_qc.melt(
    id_vars="barcode",
    value_vars=["genes", "umis", "mito_percent"],
    var_name="metric",
    value_name="value"
)

label_map = {
    "genes": "Genes",
    "umis": "UMIs",
    "mito_percent": "Mito percent",
}
plot_df["metric"] = plot_df["metric"].map(label_map)

sns.set(style="whitegrid", font_scale=1.2)

fig, axes = plt.subplots(1, 3, figsize=(16, 6))

for ax, metric in zip(axes, ["Genes", "UMIs", "Mito percent"]):
    sub = plot_df[plot_df["metric"] == metric]
    sns.violinplot(
        data=sub,
        x="metric",
        y="value",
        inner="box",
        color="#4C9BD5",
        cut=0,
        ax=ax
    )
    ax.set_xlabel("")
    ax.set_ylabel(metric)
    ax.set_title(metric)

fig.suptitle("RNA Violin Plot", fontsize=18)
plt.tight_layout()
plt.savefig("rna_violin_plot.png", dpi=300, bbox_inches="tight")
plt.show()
