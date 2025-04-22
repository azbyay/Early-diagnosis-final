""" ---------- Still Preprocess Data To Match and Merge Proteome and Phosphoproteome columns ---------- """
import pandas as pd
from mygene import MyGeneInfo
import psutil #import library for log memory 
import os

dtype_dict = {
    'Log Ratio': 'float32',
    'AJCC Pathologic Stage': 'category',
    'Gene': 'category',
    'Phosphosite': 'category',
    'Protein_Accession': 'category'
}

def log_memory():
    print(f"Memory used: {psutil.Process().memory_info().rss / 1024**2:.2f} MB")

# Load phosphoproteome data
phospho = pd.read_csv('merged_clinical_phosphoproteome_data.csv', usecols=['Phosphosite', 'Log Ratio', 'AJCC Pathologic Stage'], dtype=dtype_dict)

phospho['Protein_Accession'] = phospho['Phosphosite'].str.split(':').str[0].astype('category')
log_memory()  # ~500 MB
print('Mapping genes to protein accessions...')

# Extract NP_XXXXX from Phosphosite column
phospho['Protein_Accession'] = phospho['Phosphosite'].str.split(':').str[0]


# ---------------------------
# ---------------------------
# 3. Hybrid Mapping Approach
# ---------------------------

cache_file = 'protein_gene_mapping.feather'
uniprot_mapping_file = 'uniprotkb_proteome_UP000005640_2025_04_16.tsv'  # Pre-download from UniProt

# Try to load cached mappings first
if os.path.exists(cache_file):
    mapping = pd.read_feather(cache_file)
else:
    # Initialize empty mapping
    mapping = pd.DataFrame(columns=['Protein_Accession', 'Gene'])
    

    # Add UniProt fallback mapping
    if os.path.exists(uniprot_mapping_file):
        uniprot_map = pd.read_csv(uniprot_mapping_file, usecols=['RefSeq Protein ID', 'Gene Names'])
        uniprot_map = uniprot_map.rename(columns={
            'RefSeq Protein ID': 'Protein_Accession',
            'Gene Names': 'symbol'
        })
        mapping = pd.concat([mapping, uniprot_map])
    
    # Clean and save mappings
    if not mapping.empty:
        mapping = mapping.dropna().drop_duplicates()
        mapping['symbol'] = mapping['symbol'].str.split(';').explode().str.strip()
        mapping.to_feather(cache_file)
    else:
        print("Warning: No mappings found. Continuing with partial data.")
        mapping = pd.DataFrame(columns=['Protein_Accession', 'symbol'])

# ---------------------------
# 4. Process Mappings
# ---------------------------
if not mapping.empty:
    mapping.columns = ['Protein_Accession', 'Gene']
    mapping = mapping.dropna().drop_duplicates()
    
    # Merge with phospho data
    phospho_mapped = phospho.merge(
        mapping,
        on='Protein_Accession',
        how='left'
    )
    
    # Save unmapped entries for inspection
    unmapped = phospho_mapped[phospho_mapped['Gene'].isna()]
    if not unmapped.empty:
        unmapped[['Protein_Accession']].to_csv('unmapped_accessions.csv', index=False)
        print(f"{len(unmapped)} unmapped accessions saved for review")
    
    # Proceed with mapped entries
    phospho_mapped = phospho_mapped.dropna(subset=['Gene'])
else:
    phospho_mapped = phospho.copy()
    phospho_mapped['Gene'] = pd.NA


# Phase 1: Merge phospho with mappings
phospho_mapped = phospho.merge(
    mapping,
    on='Protein_Accession',
    how='inner'
).drop(columns='Protein_Accession')

log_memory()  # ~800 MB

# Phase 2: Load and merge proteome
proteome = pd.read_csv(
    'merged_clinical_proteome_data.csv',
    usecols=['Gene', 'Log Ratio', 'AJCC Pathologic Stage'],
    dtype=dtype_dict
)

# Split into chunks for merging
merged_chunks = []
chunk_size = 10000  # Adjust based on memory

if not phospho_mapped.empty:
    for i in range(0, len(phospho_mapped), chunk_size):
        chunk = phospho_mapped.iloc[i:i+chunk_size]
        merged = chunk.merge(
            proteome,
            on='Gene',
            suffixes=('_phospho', '_proteome'),
            how='inner'
        )
        merged_chunks.append(merged)
        log_memory()
else:
    print("""
    No mapped phosphoproteome data available.
    Possible reasons:
    - All phosphosite mappings failed
    - No matching genes between proteome and phospho data
    - Input files contain invalid accession IDs
    """)
    exit(1)

# ---------------------------
# 5. Final Assembly (with validation)
# ---------------------------
print("Combining results...")

if len(merged_chunks) == 0:
    raise ValueError("""
    No data to merge. Possible causes:
    1. Gene mappings failed for all phosphosites
    2. No overlap between proteome and phospho genes
    3. Input files are empty
    Check 'unmapped_accessions.csv' and mapping files.
    """)

try:
    final = pd.concat(merged_chunks, ignore_index=True)
except ValueError as e:
    print(f"Merge failed: {str(e)}")
    print(f"Number of chunks: {len(merged_chunks)}")
    print(f"Chunk sizes: {[len(c) for c in merged_chunks]}")
    exit(1)

print("Phospho mapped sample:", phospho_mapped.head())
print("Proteome sample:", proteome.head())
common_genes = set(phospho_mapped['Gene']).intersection(set(proteome['Gene']))
print(f"Shared genes: {len(common_genes)}")
if os.path.exists('unmapped_accessions.csv'):
    unmapped = pd.read_csv('unmapped_accessions.csv')
    print("Unmapped accessions:", unmapped['Protein_Accession'].unique())

# ---------------------------
# 6. Save Results
# ---------------------------
print("Saving...")
final.to_csv('integrated_proteome_phospho.csv', index=False)

# Save
merged.to_csv('integrated_proteome_phospho.csv', index=False)

# Split multi-gene entries
mapping['Gene'] = mapping['Gene'].str.split(';')
mapping = mapping.explode('Gene')

"""

Verify gene symbols match between datasets

"""
assert set(proteome['Gene']).issuperset(set(phospho_mapped['Gene']))

# Check merged dimensions
print(f"Proteome rows: {len(proteome)}")
print(f"Phospho rows: {len(phospho_mapped)}")
print(f"Merged rows: {len(merged)}")

