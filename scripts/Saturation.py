import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("cell_umi_gene.tsv", sep="\t")

fractions = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

records = []

for frac in fractions:
    df_sub = df.sample(frac=frac, random_state=42)

    n_reads = df_sub.shape[0]
    df_sub["gene_bc_umi"] = (
        df_sub["gene"].astype(str) + "_" +
        df_sub["barcode"].astype(str) + "_" +
        df_sub["umi"].astype(str)
    )
    n_deduped_reads = df_sub["gene_bc_umi"].nunique()
    saturation = 1 - (n_deduped_reads / n_reads)

    genes_per_cell = df_sub.groupby("barcode")["gene"].nunique().median()
    umis_per_cell = df_sub.groupby("barcode")["umi"].nunique().median()
    reads_per_cell = df_sub.groupby("barcode")["read_id"].nunique().median()

    records.append((frac, n_reads, reads_per_cell, genes_per_cell, umis_per_cell, saturation))

res = pd.DataFrame(
    records,
    columns=["fraction", "reads", "reads_per_cell", "genes_per_cell", "umis_per_cell", "saturation"]
)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

axes[0].plot(res["reads_per_cell"], res["genes_per_cell"], marker="o")
axes[0].set_xlabel("Median reads per cell")
axes[0].set_ylabel("Median genes per cell")
axes[0].set_title("Genes per cell")

axes[1].plot(res["reads_per_cell"], res["umis_per_cell"], marker="o")
axes[1].set_xlabel("Median reads per cell")
axes[1].set_ylabel("Median UMIs per cell")
axes[1].set_title("UMIs per cell")

axes[2].plot(res["reads_per_cell"], res["saturation"], marker="o")
axes[2].set_xlabel("Median reads per cell")
axes[2].set_ylabel("Sequencing saturation")
axes[2].set_title("Saturation")

plt.tight_layout()
plt.savefig("saturation_curves.png", dpi=300)
plt.show()

print(res)
