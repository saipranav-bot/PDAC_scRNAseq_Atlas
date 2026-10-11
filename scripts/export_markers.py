import os
import scanpy as sc

adata = sc.read_h5ad("data/processed/PDAC_scanpy_final.h5ad")

os.makedirs("results/Markers/Cell_Type_Markers", exist_ok=True)

groups = adata.uns["rank_genes_groups"]["names"].dtype.names

print("Groups found:")
print(groups)

for g in groups:
    print(f"Saving {g}")

    df = sc.get.rank_genes_groups_df(
        adata,
        group=g
    )

    df.to_csv(
        f"results/Markers/Cell_Type_Markers/{g}_markers.csv",
        index=False
    )

print("\nFinished exporting all marker genes.")
