import pandas as pd
import re
# Define the phospho DataFrame
phospho = pd.read_csv('merged_clinical_phosphoproteome_data.csv', sep=',', header=0)

# Check both NaN and empty strings
def is_null_or_empty(x):
    return pd.isna(x) or (isinstance(x, str) and x.strip() == '')

null_mask = phospho['Phosphosite'].apply(is_null_or_empty)
print(f"True invalid entries: {null_mask.sum()}")

def process_uniprot_refseq(refseq_str):
    """Process RefSeq Protein ID column entries"""
    # Split into individual accessions
    if pd.isna(refseq_str):
        return []
    accessions = []
    for part in str(refseq_str).split(';'):
        # Remove any brackets and their contents (isoform annotations)
        base = part.split('[')[0].strip()
        # Remove version numbers
        base = re.sub(r'\.\d+$', '', base)
        if base.startswith(('NP_', 'XP_')):  # Keep only RefSeq protein accessions
            accessions.append(base)
    return accessions

# Load and process UniProt data
uniprot = pd.read_csv(
    '/Users/awiener/Downloads/uniprotkb_proteome_UP000005640_2025_04_16.tsv', 
    sep='\t',
    usecols=['RefSeq', 'Gene Names'],
    dtype={'RefSeq': str, 'Gene Names': 'string' },
    keep_default_na=False,
    na_values = ['']
).replace('', pd.NA)

# Explode RefSeq IDs into individual rows
uniprot['Protein_Accession'] = uniprot['RefSeq'].apply(process_uniprot_refseq)
uniprot = uniprot.explode('Protein_Accession').dropna(subset=['Protein_Accession'])

# Process Gene Names (take first gene if multiple)
uniprot['Gene'] = uniprot['Gene Names'].str.split(';').str[0].str.strip()
print("Null values in RefSeq:", uniprot['RefSeq'].isna().sum())
print("Sample processed accessions:", uniprot['Protein_Accession'].head(10).tolist())
# Process phosphosite accessions
phospho['Protein_Accession'] = (
    phospho['Phosphosite']
    .str.split(':').str[0]
    .replace(r'^\s*$', pd.NA, regex=True)              # Extract NP_XXXXX.2 from NP_XXXXX.2:s76
    .str.split('.').str[0]              # Remove version numbers (NP_XXXXX)
)

print(f"Final nulls: {phospho['Protein_Accession'].isna().sum()}")


#Full Null Reporting
def report_nulls(df, name):
    total = len(df)
    nulls = df.isna().sum().sum()
    print(f"\n{name} null report:")
    print(f"- Total entries: {total}")
    print(f"- Null entries: {nulls} ({nulls/total:.1%})")
    print("- Columns with nulls:")
    print(df.isna().sum())

report_nulls(phospho, "Phosphoproteome Data")
report_nulls(uniprot, "UniProt Mappings")



original_refseq = "NP_001123497.1;NP_001123498.2 [A6NFQ2-3];NP_775949.2 [A6NFQ2-2];XP_006715990.1;XP_006715991.1"
processed = process_uniprot_refseq(original_refseq)
print(processed)
# Output: ['NP_001123497', 'NP_001123498', 'NP_775949', 'XP_006715990', 'XP_006715991']

# Create mapping dictionary
refseq_to_gene = uniprot[['Protein_Accession', 'Gene']].drop_duplicates().set_index('Protein_Accession')['Gene'].to_dict()

# Map phosphosite data
phospho['Gene'] = phospho['Protein_Accession'].map(refseq_to_gene)

# Check example phosphosite
test_phospho = phospho[phospho['Phosphosite'].str.contains('NP_997229.2:s76')]
print(f"Test phosphosite mapped to: {test_phospho['Gene'].values[0]}")

# Check coverage
mapped_percent = phospho['Gene'].notna().mean() * 100
print(f"Mapped phosphosites: {mapped_percent:.1f}%")


# Map gene names to phospho DataFrame
phospho['Gene'] = phospho['Protein_Accession'].map(refseq_to_gene)

# Save to CSV
phospho.to_csv('phospho_with_gene_mapping.csv', index=False)

