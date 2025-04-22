import pandas as pd

# Read 1.csv and filter CPT rows --------------------------------------------
df1 = pd.read_csv("/Users/awiener/Downloads/CPTAC3_Pancreatic_Ductal_Adenocarcinoma_Phosphoproteome.phosphosite.tmt11 (1).tsv", index_col=0)

# More flexible regex pattern for CPT rows
cpt_pattern = r'^CPT\d+.*Log Ratio.*$'  # Case-insensitive, allows variations
df1 = df1[df1.index.str.contains(cpt_pattern, regex=True, case=False)]

# Clean CPT names (remove "Log Ratio" and whitespace)
df1.index = df1.index.str.replace(r'\s*Log Ratio\s*', '', regex=True).str.strip()

# Clean phosphosite column names in 1.csv
df1.columns = df1.columns.str.strip()

# Debug: Check phosphosites in 1.csv
print("", df1.columns.tolist()[:5])

# Transpose and reshape data -----------------------------------------------
df1_transposed = df1.T  # Columns become CPTs, rows become phosphosites

# Convert to long format
df1_long = df1_transposed.reset_index().melt(
    id_vars="index",
    var_name="CPT",
    value_name="Log_Ratio"
).rename(columns={"index": "Phosphosite"})

df1_long["Phosphosite"] = df1_long["Phosphosite"].str.strip().str.lower()  # Clean phosphosite names
# Debug: Check sample log ratios
print("\nSample data from CPTAC_raw_data.csv after reshaping:")
print(df1_long.head())

# Read and clean 2.csv -----------------------------------------------------
df2 = pd.read_csv("/Users/awiener/projects_yay/Early-diagnosis-final/merged_result.csv")
df2["Phosphosite"] = df2["Phosphosite"].str.strip().str.lower()

# Debug: Check phosphosites in 2.csv
print("\nPhosphosites in mergedresult.csv:", df2["Phosphosite"].unique()[:5])

# 2. (Optional) Debug: list a few unmatched sites
unmatched = set(df2["Phosphosite"]) - set(df1_long["Phosphosite"])
print("Unmatched phosphosites:", list(unmatched)[:10])

# Merge data ---------------------------------------------------------------
merged = df2.merge(
    df1_long[["Phosphosite", "Log_Ratio"]],
    on="Phosphosite",
    how="left"
).rename(columns={"Log_Ratio": "Log Ratio"})

# Debug: Check merge results
print("\nMerge validation:")
print("Rows in merged data:", len(merged))
print("Non-empty log ratios:", merged["Log Ratio"].notna().sum())

# Save result --------------------------------------------------------------
merged.to_csv("merged_results.csv", index=False)
print("\nMerge complete. Check merged_results.csv")