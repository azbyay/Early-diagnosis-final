import pandas as pd

# Load the CSV files
df1 = pd.read_csv('/Users/awiener/Downloads/PDC_biospecimen_manifest_04132025_193114.csv')
df2 = pd.read_csv('/Users/awiener/Downloads/phospho_with_gene_mapping.csv')

# Ensure no duplicates in df1 (to avoid unintended row multiplication)
df1 = df1.drop_duplicates(subset='Aliquot Submitter ID', keep='first')

# Merge (inner join to enforce matching IDs in both files)
merged_df = df2.merge(
    df1[['Aliquot Submitter ID', 'Tissue Type']],
    on='Aliquot Submitter ID',
    how='inner'  # Use inner join to exclude non-matching IDs
)

# Save the merged result
merged_df.to_csv('merged_result.csv', index=False)

# Verify no null values in 'tissue type'
assert merged_df['Tissue Type'].isna().sum() == 0, "Null values detected in 'tissue type'"
print("Merge successful with no null values!")