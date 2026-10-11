import os
import pandas as pd

input_dir = "results/Markers/Cell_Type_Markers"
output_dir = "results/Markers/Filtered_Cell_Type_Markers"

os.makedirs(output_dir, exist_ok=True)

for file in os.listdir(input_dir):

    if not file.endswith(".csv"):
        continue

    df = pd.read_csv(os.path.join(input_dir, file))

    df = df[
        (df["pvals_adj"] < 0.05) &
        (df["logfoldchanges"] > 1)
    ].sort_values(
        by="logfoldchanges",
        ascending=False
    )

    df.to_csv(
        os.path.join(output_dir, file),
        index=False
    )

print("Filtered marker genes exported.")
