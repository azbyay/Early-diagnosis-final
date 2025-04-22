import pandas as pd
import re

# Define the phospho DataFrame
phospho = pd.read_csv('/Users/awiener/projects_yay/Early-diagnosis-final/merged_results_final.csv', sep=',', header=0)

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
    '/Users/awiener/Downloads/uniprotkb_proteome_UP000005640_2025_04_16 (1).tsv', 
    sep='\t',
    usecols=['RefSeq', 'Gene Names', 'Pathway'],
    dtype={'RefSeq': str, 'Gene Names': 'string', 'Pathway':'string' },
    keep_default_na=False,
    na_values = ['']
).replace('', pd.NA)

# Explode RefSeq IDs into individual rows
uniprot = uniprot.rename(columns={'Gene Names':'Gene'})
uniprot['Protein_Accession'] = uniprot['RefSeq'].apply(process_uniprot_refseq)
uniprot = uniprot.explode('Protein_Accession').dropna(subset=['Protein_Accession'])


# Process Gene Names (take first gene if multiple)
def _clean_pathway(L):
    # only join if L is a list of strings
    if isinstance(L, list):
        return '; '.join(p.strip() for p in L if p and p.strip())
    else:
        return pd.NA

uniprot['Pathway'] = (
    uniprot['Pathway']
      # remove leading “PATHWAY:” if present
      .str.replace(r'^\s*PATHWAY:\s*', '', regex=True)
      # drop trailing “{…}.” evidence tags
      .str.replace(r'\s*\{.*\}\.?\s*$', '', regex=True)
      # split on semicolon into list
      .str.split(';')
      # strip whitespace and drop empties, then re‑join with “; ”
      .apply(_clean_pathway)
)


refseq_to_gene = (
    uniprot[['Protein_Accession', 'Gene']]
    .drop_duplicates()
    .set_index('Protein_Accession')['Gene']
    .to_dict()
)
refseq_to_pathway = (
    uniprot[['Protein_Accession', 'Pathway']]
    .drop_duplicates()
    .set_index('Protein_Accession')['Pathway']
    .to_dict()
)



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


phospho['Gene'] = phospho['Protein_Accession'].map(refseq_to_gene)
phospho['Pathway'] = phospho['Protein_Accession'].map(refseq_to_pathway)

# …existing code saving…
phospho.to_csv('phospho_with_gene_mapping_and_pathwaynew.csv', index=False)



