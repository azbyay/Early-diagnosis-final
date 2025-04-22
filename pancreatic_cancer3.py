import pandas as pd

# Load data 

import pandas as pd

# Step 1: Filter biospecimen by clinical Case IDs
clinical = pd.read_csv('/Users/awiener/Downloads/PDC_clinical_manifest_04132025_193659.csv', sep=',', header=0)
clinical_case_ids = set(clinical['Case ID'].unique())

biospecimen = pd.read_csv('/Users/awiener/Downloads/PDC_biospecimen_manifest_04132025_193114.csv', sep=',',header=0)
filtered_biospecimen = biospecimen[biospecimen['Case ID'].isin(clinical_case_ids)]
# In case of any duplicates in the Aliquot Submitter ID, remove it.
clean_biospecimen = filtered_biospecimen.drop_duplicates(
    subset=['Aliquot Submitter ID'],
    keep=False
)

if clean_biospecimen.empty:
    raise ValueError("No valid biospecimen data found for the provided clinical Case IDs.")

aliquot_ids = clean_biospecimen['Aliquot Submitter ID'].unique().tolist()

# Step 2: Filter proteome and phosphoproteome data
def filter_columns(df, id_col):
    # Keep identifier column and valid log ratio columns
    valid_columns = [id_col] + [
        col for col in df.columns
        if ' Log Ratio' in col
        and not any(x in col for x in ['unshared', 'withdrawn'])
    ]
    return df[valid_columns]

# Load and filter proteome
proteome = pd.read_csv('/Users/awiener/Downloads/CPTAC3_Pancreatic_Ductal_Adenocarcinoma_Proteome.tmt11.tsv', sep='\t')
proteome = filter_columns(proteome, 'Gene')

# Load and filter phosphoproteome
phosphoproteome = pd.read_csv('/Users/awiener/Downloads/CPTAC3_Pancreatic_Ductal_Adenocarcinoma_Phosphoproteome.phosphosite.tmt11 (1).tsv', sep='\t')
phosphoproteome = filter_columns(phosphoproteome, 'Phosphosite')

# Step 3-5: Find valid aliquots present in both datasets
valid_aliquots = []
for aliquot in aliquot_ids:
    col_name = f"{aliquot} Log Ratio"
    proteome_has = col_name in proteome.columns
    phosphoproteome_has = col_name in phosphoproteome.columns
    if proteome_has and phosphoproteome_has:
        valid_aliquots.append(col_name)
    else: print(f"Excluding {aliquot} - present in proteome: {proteome_has}, {phosphoproteome_has}")
    

# Step 6: Save filtered data
if valid_aliquots:
    # For proteome
    proteome_filtered = proteome[['Gene'] + valid_aliquots]
    proteome_filtered.to_csv('filtered_proteome(1).csv', index=False)
    
    # For phosphoproteome
    phospho_filtered = phosphoproteome[['Phosphosite'] + valid_aliquots]
    phospho_filtered.to_csv('filtered_phosphoproteome(1).csv', index=False)
    print(f"Successfully saved data with {len(valid_aliquots)} valid aliquots")
else:
    print("No valid aliquots found in both datasets")

print("Processing complete!")