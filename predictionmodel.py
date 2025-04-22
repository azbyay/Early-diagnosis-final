import pandas as pd
import numpy as np
import re
from sklearn.preprocessing import MultiLabelBinarizer
from scipy.sparse import vstack, csr_matrix
from sklearn.decomposition import IncrementalPCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split

# Load data and initial processing
full_data = pd.read_csv('/Users/awiener/Downloads/phospho_with_gene_mapping.csv', low_memory=False)
print("Unique AJCC Stages:", full_data['AJCC Pathologic Stage'].unique())
print("Count per stage:\n", full_data['AJCC Pathologic Stage'].value_counts())

def process_phosphosite(row):
    """Extract phosphorylation sites with protein and gene info."""
    try:
        protein_part, sites_part = row['Phosphosite'].split(':')
    except ValueError:
        return []
    protein_accession = protein_part.split('.')[0]
    sites = re.findall(r'[sty]\d+', sites_part)
    gene = 'Unknown_Gene' if pd.isna(row['Gene']) else str(row['Gene']).split()[0].strip()
    return [f"{protein_accession}_{gene}_{site}" for site in sites]


# First pass: Collect all unique phosphorylation sites
all_labels = set()
chunk_size = 10000  # Adjust based on available RAM
all_features = None
all_targets = []

column_types = {
    'Phosphosite': 'category',
    'Gene': 'category',
    'AJCC Pathologic Stage': 'category',

}
for chunk in pd.read_csv(
    '/Users/awiener/Downloads/phospho_with_gene_mapping.csv',
    chunksize=chunk_size,
    dtype=column_types,
    low_memory=False
):
    chunk['Protein_Gene_Site'] = chunk.apply(process_phosphosite, axis=1)
    for sites in chunk['Protein_Gene_Site']:
        all_labels.update(sites)

all_labels = sorted(all_labels)
print(f"Total unique phosphorylation sites: {len(all_labels)}")

# Initialize MultiLabelBinarizer
mlb = MultiLabelBinarizer(classes=all_labels, sparse_output=True) #already set the classes during initialization
mlb.fit([])

# Second pass: Process data and encode features
#second load:

#all_features = None
all_features = csr_matrix((0, len(all_labels)))
all_targets = []
all_stages = full_data['AJCC Pathologic Stage'].astype('category').cat.categories.tolist()

total_chunks = 0 
for chunk in pd.read_csv('/Users/awiener/Downloads/phospho_with_gene_mapping.csv',
    chunksize=chunk_size,
    dtype=column_types,
    low_memory=False
):
    #{
    #    'Phosphosite': 'category',
    #    'Gene': 'category',
    #    'AJCC Pathologic Stage': pd.CategoricalDtype(categories=all_stages)
    #},
    #low_memory=False
    try:
        total_chunks += 1
        # Encode target
        chunk['AJCC Pathologic Stage'] = chunk['AJCC Pathologic Stage'].astype(
            pd.CategoricalDtype(categories=all_stages))

        # Get targets
        targets = chunk['AJCC Pathologic Stage'].cat.codes.values
        num_samples = len(chunk)

        chunk['Protein_Gene_Site'] = chunk.apply(process_phosphosite, axis=1)
        chunk_features = mlb.transform(chunk['Protein_Gene_Site'])

        if chunk_features.shape[0] != num_samples:
            print(f" Feature/target mismatch in chunk {total_chunks}: {chunk_features.shape[0]} vs {num_samples}")
            continue
        all_features = vstack([all_features, chunk_features])
        all_targets.append(targets)

        #if len(targets) != num_samples:
        #    print(f"Mismatch in target length: {len(targets)} vs {num_samples}")
        #    continue
        #all_targets.append(targets)
        if total_chunks % 50 == 0:
            print(f"Processed {total_chunks} chunks ({all_features.shape[0]} samples)")

        # Process features (ensure list format)
        ##chunk['Protein_Gene_Site'] = chunk.apply(process_phosphosite, axis=1)
        ##chunk_features = mlb.transform(chunk['Protein_Gene_Site'])
        #chunk['Protein_Gene_Site'] = chunk['Protein_Gene_Site'].apply(
        #    lambda x: x if isinstance(x, list) else []
        #)

        
        # Stack features
        ##all_features = vstack([all_features, chunk_features])
    #moved except up 
    except Exception as e:
        print(f"Error processing chunk: {str(e)}")
        if not chunk.empty and hasattr(chunk, 'index'):
            print(f"Chunk starts at index: {chunk.index.min()}")
        continue
    """try:
        # Encode target
        all_targets.append(chunk['AJCC Pathologic Stage'].cat.codes.values)
        
        # Process features
        chunk['Protein_Gene_Site'] = chunk.apply(process_phosphosite, axis=1)
        
        # Transform with MultiLabelBinarizer
        if len(chunk) > 0:
            chunk_features = mlb.transform(chunk['Protein_Gene_Site'])
        else:
            chunk_features = csr_matrix((0, len(all_labels)))

        # Safeguard stacking
        if chunk_features.shape[1] != all_features.shape[1]:
            print(f"Mismatch! Expected {all_features.shape[1]} features, got {chunk_features.shape[1]}")
            raise ValueError("Feature dimension mismatch in chunk")
           
        all_features = vstack([all_features, chunk_features])
    """    

    # Encode target
    #all_targets.append(chunk['AJCC Pathologic Stage'].cat.codes.values)
    
    # Process features
    #chunk['Protein_Gene_Site'] = chunk.apply(process_phosphosite, axis=1)
    #chunk['Protein_Gene_Site'] = chunk['Protein_Gene_Site'].apply(lambda x: x if isinstance(x, list) else [])

     # Transform with MultiLabelBinarizer
    #if len(chunk) > 0:
    #    chunk_features = mlb.transform(chunk['Protein_Gene_Site'])
    #else:
    #    chunk_features = csr_matrix((0, len(all_labels)))

        # Safeguard stacking
    #    if all_features is None:
    #        all_features = chunk_features
    #    else:
            # Assert and check shapes before stacking
            #if chunk_features.shape[1] != all_features.shape[1]:
                #print(f"Mismatch! all_features shape: {all_features.shape}, chunk_features shape: {chunk_features.shape}")
                #raise ValueError("Feature dimension mismatch before stacking.")

            #if len(chunk_features.shape) != 2:
                #print(f"Invalid chunk_features dimension: {chunk_features.shape}")
                #raise ValueError("chunk_features is not 2D")

            #if len(all_features.shape) != 2:
                #print(f"Invalid all_features dimension: {all_features.shape}")
                #raise ValueError("all_features is not 2D")

    #        all_features = vstack([all_features, chunk_features])

# Combine all target arrays
y = np.concatenate(all_targets)
# Final validation
if all_features.shape[0] != len(y):
    raise  ValueError (f"Critical mismatch: {all_features.shape[0]} features vs {len(y)} targets")
#if len(y) == 0:
    #raise ValueError("No targets accumulated - check data processing steps")

print(f"Final feature matrix shape: {all_features.shape}")
print(f"Target vector length: {len(y)}")
print(f"\nFinal dataset:")
print(f"- Samples: {all_features.shape[0]}")
print(f"- Features: {all_features.shape[1]}")
print(f"- Targets: {len(y)}")


# Dimensionality Reduction with IncrementalPCA
n_components = 100
ipca = IncrementalPCA(n_components=n_components, batch_size=chunk_size)

# Fit PCA incrementally on dense chunks
for i in range(0, all_features.shape[0], chunk_size):
    chunk_sparse = all_features[i:i+chunk_size]
    ipca.partial_fit(chunk_sparse.toarray())

# Transform features incrementally
X_pca = []
for i in range(0, all_features.shape[0], chunk_size):
    chunk_sparse = all_features[i:i+chunk_size]
    X_pca_chunk = ipca.transform(chunk_sparse.toarray())
    X_pca.append(X_pca_chunk)

X_pca = np.vstack(X_pca)

# Train-test split and model training
X_train, X_test, y_train, y_test = train_test_split(X_pca, y, test_size=0.2, stratify=y, random_state=42)

model = LogisticRegression(class_weight='balanced', max_iter=1000, random_state=42)
model.fit(X_train, y_train)

# Evaluation
print(f"Model Accuracy: {model.score(X_test, y_test):.2f}")
