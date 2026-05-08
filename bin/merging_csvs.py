import glob
import pandas as pd

dfs_path = "../output/features/results/*.csv"

print(dfs_path)

df = pd.concat([pd.read_csv(f) for f in glob.glob(dfs_path)], ignore_index=True)

df_selected = pd.concat([
    df.iloc[:, :1],
    df.iloc[:, [9]],
    df.iloc[:, 21:24],
    df.iloc[:, [25]],
    df.iloc[:, 50:]
], axis=1)
df_selected.to_csv("../output/features/merged_results.csv", index=False)
