import pandas as pd

# Load files
stats_df = pd.read_csv("/Users/awiener/projects_yay/Early-diagnosis-final/significant_phosphositess.csv")
metadata_df = pd.read_csv("/Users/awiener/projects_yay/Early-diagnosis-final/merged_result.csv")

# Merge on Phosphosite ID
merged_df = pd.merge(
    stats_df,
    metadata_df,
    on="Phosphosite",
    how="inner"  # Keep only mapped phosphosites
)

print(f"Merged data shape: {merged_df.shape}")
print(merged_df.head())

# Check for unmapped entries
unmapped_stats = stats_df[~stats_df['Phosphosite'].isin(metadata_df['Phosphosite'])]
print(f"Unmapped stats entries: {len(unmapped_stats)}")

unmapped_meta = metadata_df[~metadata_df['Phosphosite'].isin(stats_df['Phosphosite'])]
print(f"Unmapped metadata entries: {len(unmapped_meta)}")

merged_df.to_csv("merged_results_final.csv", index=False)
print("\nMerge complete. Check merged_results.csv")