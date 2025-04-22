import pandas as pd
from pancreatic_cancer3 import proteome_filtered, phospho_filtered, filtered_biospecimen, aliquot
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
#Load the respective filtered datasets: clinical and biospecimen

clinical = pd.read_csv('/Users/awiener/Downloads/PDC_clinical_manifest_04132025_193659.csv', sep=',', header=0)
biospecimen = pd.read_csv('/Users/awiener/Downloads/PDC_biospecimen_manifest_04132025_193114.csv', sep=',', header=0)
filtered_proteome = pd.read_csv('filtered_proteome(1).csv', sep=',', header=0)
filtered_phosphoproteome = pd.read_csv('filtered_phosphoproteome(1).csv', sep=',', header=0)

# Debugging: Print the columns to confirm they are parsed correctly

print("Biospecimen columns:", biospecimen.columns)
print("Biospecimen columns after loading:", biospecimen.columns.tolist())
print("Filtered biospecimen columns:", filtered_biospecimen.columns)
print(f"aliquot: {aliquot}")

# Debugging: Print column names and sample rows
print("Biospecimen columns:", biospecimen.columns)
print("Filtered biospecimen columns:", filtered_biospecimen.columns)
print("Clinical columns:", clinical.columns)
print("Proteome filtered columns:", proteome_filtered.columns)
print("Phosphoproteome filtered columns:", phospho_filtered.columns)


# Validate required columns in biospecimen
required_columns = ['Aliquot Submitter ID', 'Case ID']
for col in required_columns:
    if col not in biospecimen.columns:
        raise KeyError(f"Missing required column in biospecimen: {col}")

aliquot_column = f"{aliquot} Log Ratio"
if aliquot_column not in proteome_filtered.columns:
    print("Available columns in proteome_filtered:", proteome_filtered.columns)
    raise KeyError(f"Aliquot column '{aliquot_column}' not found in proteome_filtered columns.")

# Debugging: Print the aliquot column being used
print(f"Using aliquot column for merge: {aliquot_column}")


filtered_biospecimen.columns = filtered_biospecimen.columns.str.strip()
biospecimen.columns = biospecimen.columns.str.strip()


#validate the aliquot ids
proteome_aliquots = set([col.replace(' Log Ratio', '') for col in proteome_filtered.columns if 'Log Ratio' in col])
phospho_aliquots = set([col.replace(' Log Ratio', '') for col in phospho_filtered.columns if 'Log Ratio' in col])

assert proteome_aliquots == phospho_aliquots, \
    f"Aliquot mismatch! Proteome: {proteome_aliquots - phospho_aliquots} | Phospho: {phospho_aliquots - proteome_aliquots}"

# Debugging: Print aliquot sets
print("Proteome aliquots:", proteome_aliquots)
print("Phosphoproteome aliquots:", phospho_aliquots)



# Debugging: Print aliquot sets
print("Proteome aliquots:", proteome_aliquots)
print("Phosphoproteome aliquots:", phospho_aliquots)



#after the bunch of debugging steps above, proceed to reshape the proteome data to long format.
proteome_long = proteome_filtered.melt(
    id_vars=['Gene'], 
    var_name='Aliquot_LogRatio_Column',  # Temporary column to hold "CPTXXXX Log Ratio" strings
    value_name='Log Ratio'
)
# Extract Aliquot Submitter ID by stripping " Log Ratio" from column names
proteome_long['Aliquot Submitter ID'] = proteome_long['Aliquot_LogRatio_Column'].str.replace(' Log Ratio', '')
proteome_long = proteome_long.drop(columns=['Aliquot_LogRatio_Column'])  # Cleanup

# Convert to string if needed
proteome_long['Aliquot Submitter ID'] = proteome_long['Aliquot Submitter ID'].astype(str)
filtered_biospecimen['Aliquot Submitter ID'] = filtered_biospecimen['Aliquot Submitter ID'].astype(str)

# Merge proteome + biospecimen (linking both)
# Step 2: Merge with Biospecimen Data (Using Aliquot Submitter ID)
# ----------------------------------------------------------------------------
# Ensure biospecimen has Aliquot Submitter ID (not Aliquot_ID)
merged_proteome = pd.merge(
    proteome_long,
    filtered_biospecimen[['Aliquot Submitter ID', 'Case ID']],  # Use Submitter ID, not Aliquot ID
    on='Aliquot Submitter ID',
    how='inner'  # Keep only matching aliquots
)


# Merge with clinical
final_data = pd.merge(
    merged_proteome,
    clinical[[
        'Case ID', 'AJCC Pathologic Stage', 'Days to Recurrence', 
                       'Days to Last Follow Up', 'Vital Status', 'Progression Free Survival', 'Tumor Grade', 
                       'Tumor Stage', 'Primary Diagnosis', 'Morphology', 'Residual Disease', 
                       'AJCC Pathologic M', 'AJCC Pathologic N', 'AJCC Pathologic T', 'Age at Diagnosis', 'Sex', 'Ethnicity', 'Race'
                       ]],
    on='Case ID'
)

print("Merge successful! Final columns:", final_data.columns.tolist())

# Save final_data to CSV
output_path = "merged_clinical_proteome_data.csv"
#final_data.to_csv(output_path, index=False)  # index=False avoids adding row numbers
print(f"Saved merged data to: {output_path}")


# Check 1: Verify Aliquot Submßitter ID matches between datasets
missing_aliquots = proteome_long[~proteome_long['Aliquot Submitter ID'].isin(filtered_biospecimen['Aliquot Submitter ID'])]
if not missing_aliquots.empty:
    print(f"Warning: {len(missing_aliquots)} proteome aliquots missing in biospecimen data")

# Check 2: Ensure no confusion between Aliquot_ID and Aliquot_Submitter_ID
assert 'Aliquot ID' not in merged_proteome.columns, \
    "Aliquot ID should not exist in merged data - use Aliquot Submitter ID instead"




# #2---------Process Phosphoproteome Dataset------------
# Reshape the phosphoproteome data to long format
phosphoproteome_long = filtered_phosphoproteome.melt(
    id_vars=['Phosphosite'], 
    var_name='Aliquot_LogRatio_Column',  # Temporary column to hold "CPTXXXX Log Ratio" strings
    value_name='Log Ratio'
)

phosphoproteome_long['Aliquot Submitter ID'] = phosphoproteome_long['Aliquot_LogRatio_Column'].str.replace(' Log Ratio', '')
phosphoproteome_long = phosphoproteome_long.drop(columns=['Aliquot_LogRatio_Column']) 

# Step 3: Merge with Biospecimen Data (Using Aliquot Submitter ID)
# ----------------------------------------------------------------------------
# Ensure biospecimen has Aliquot Submitter ID (not Aliquot_ID)
merged_phospho = pd.merge(
    phosphoproteome_long,
    filtered_biospecimen[['Aliquot Submitter ID', 'Case ID']],  # Use Submitter ID, not Aliquot ID
    on='Aliquot Submitter ID',
    how='inner'  # Keep only matching aliquots
)


# Merge with clinical
final_data = pd.merge(
    merged_phospho,
    clinical[[
        'Case ID', 'AJCC Pathologic Stage', 'Days to Recurrence', 
                       'Days to Last Follow Up', 'Vital Status', 'Progression Free Survival', 'Tumor Grade', 
                       'Tumor Stage', 'Primary Diagnosis', 'Morphology', 'Residual Disease', 
                       'AJCC Pathologic M', 'AJCC Pathologic N', 'AJCC Pathologic T', 'Age at Diagnosis', 'Sex', 'Ethnicity', 'Race'
                       ]],
    on='Case ID'
)

print("Merge successful! Final columns:", final_data.columns.tolist())

# Save final_data to CSV
output_path = "merged_clinical_phosphoproteome_data.csv"
final_data.to_csv(output_path, index=False)  # index=False avoids adding row numbers
print(f"Saved merged data to: {output_path}")