# 1. Start with merged data (stats + phospho-gene mapping)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gseapy as gp

# --------------------------
# 1. Load Merged Data
# --------------------------
def load_merged_data():
    """Load your combined dataset"""
    # Replace with your actual file path
    df = pd.read_csv("/Users/awiener/projects_yay/Early-diagnosis-final/phospho_with_gene_mapping_and_pathwaynew.csv")
    
    # Ensure these columns exist in your data:
    # Phosphosite (unique ID), log2FC, q_value, Gene, AJCC_Stage, Phospho_Level, t_stat
    required_columns = [
        'Phosphosite', 'log2FC', 'q_value', 
        'Gene', 'AJCC Pathologic Stage','t_stat'
    ]
    assert all(col in df.columns for col in required_columns), "Missing required columns!"
    
    return df.set_index('Phosphosite')  # Set phosphosite as index

# --------------------------
# 2. Filter Significant Sites
# --------------------------
def filter_significant(df):
    """Filter sites with q < 0.05 and |log2FC| > 1"""
    return df[
        (df['q_value'] < 0.05) & 
        (abs(df['log2FC']) > 1)
    ]

# --------------------------
# 3. Biological Weighting
# --------------------------
def apply_weights(df):
    """Calculate weights from q-values"""
    # Clip q_value to avoid division by zero (minimum 0.00001)
    df['weight'] = 1 / df['q_value'].clip(lower=1e-5)
    return df


# --------------------------
# 4. Multi-Site Analysis
# --------------------------
def calculate_gene_level(df):
    """Calculate weighted average phosphorylation per gene"""
    return df.groupby('Gene').apply(
        lambda group: np.average(group['log2FC'], weights=group['weight'])
    )

# --------------------------
# 5. Enhanced Residuals
# --------------------------
def calculate_residuals(original_df, gene_means):
    """Calculate confidence-weighted residuals"""
    residuals = original_df.copy()
    residuals['residual'] = np.nan  # Initialize column
    
    for gene in gene_means.index:
        # Get all sites for this gene
        gene_mask = residuals['Gene'] == gene
        
        # Calculate residuals weighted by t-statistic
        residuals.loc[gene_mask, 'residual'] = (
            residuals.loc[gene_mask, 'log2FC'] - gene_means[gene]
        ) * np.sqrt(abs(residuals.loc[gene_mask, 't_stat']))
    
    return residuals

# --------------------------
# 6. Pathway Analysis
# --------------------------
def run_pathway_analysis(gene_list):
    """Run pathway enrichment using gseapy"""
    return gp.enrichr(
        gene_list=gene_list,
        gene_sets='KEGG_2021_Human',  
        organism='Human'
    )


if __name__ == "__main__":
    # Step 1: Load data
    merged_data = load_merged_data()
    significant_data = filter_significant(merged_data)
    weighted_data=apply_weights(significant_data)
    # Step 2: Filter significant phosphosites
    significant_data = filter_significant(merged_data)
    print(f"Working with {len(significant_data)} significant phosphosites")
    
    long_data = weighted_data
    
    
    # Step 5: Gene-level analysis
    gene_levels = calculate_gene_level(long_data)
    
    # Step 6: Calculate residuals
    residual_data = calculate_residuals(weighted_data, gene_levels)
    
    # Step 7: Pathway analysis
    pathway_results = run_pathway_analysis(gene_levels.index.tolist())
    
    # Save results
    gene_levels.to_csv("gene_level_predictions.csv")
    residual_data.to_csv("residual_analysis.csv")
    pathway_results.results.to_csv("pathway_enrichment.csv")
    print("Results saved to CSV files!")

    # --------------------------
    # Validation Plots
    # --------------------------
    plt.figure(figsize=(15, 5))
    
    # Plot 1: Weight distribution
    plt.subplot(1, 2, 1)
    plt.hist(weighted_data['weight'], bins=50, color='skyblue')
    plt.xlabel('1/q_value Weight')
    plt.ylabel('Frequency')
    plt.title('Weight Distribution')
    
    # Plot 2: Residual vs t-stat
    plt.subplot(1, 2, 2)
    plt.scatter(residual_data['t_stat'], residual_data['residual'], alpha=0.5, color='salmon')
    plt.xlabel('t-statistic')
    plt.ylabel('Residual')
    plt.title('Residual Dependence on t-stat')
    
    plt.tight_layout()
    plt.savefig("validation_plots.png")
    plt.show()