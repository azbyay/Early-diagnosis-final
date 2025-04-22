import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import t, f
from statsmodels.stats.multitest import fdrcorrection

# 1. Load Phosphoproteomics Data -------------------------------------------------
file_path = "/Users/awiener/projects_yay/Early-diagnosis-final/merged_result.csv"
def load_phospho_data(file_path):
    #Load long-format data
    phospho_long = pd.read_csv(file_path)
    # Replace 'sample' with the actual column name for sample IDs
    sample_column = 'Aliquot Submitter ID' 
    if sample_column not in phospho_long.columns:
        raise KeyError(f"Column '{sample_column}' not found in the dataset. Available columns: {list(phospho_long.columns)}")
    
    #Pivot to wide format (samples as rows, phosphosites as columns)
    valid_sample_pattern = r'^CPT\d+$'
    phospho_long = phospho_long[phospho_long[sample_column].str.match(valid_sample_pattern)]
    
    #Validate log2 transformation
    phospho_wide = phospho_long.pivot(index=sample_column, columns='Phosphosite', values='Log Ratio')
    print(f"Data range: [{phospho_wide.min().min():.2f}, {phospho_wide.max().max():.2f}] log2 units")
    
    #Filter low-observation sites
    min_obs = 5
    valid_sites = phospho_wide.columns[phospho_wide.count() >= min_obs]
    return phospho_wide[valid_sites]

# Example usage:
phospho_log2 = load_phospho_data("/Users/awiener/projects_yay/Early-diagnosis-final/merged_result.csv")

# 2. Create Design Matrix --------------------------------------------------------
# Sample metadata (must match phospho data index)
num_samples = len(phospho_log2.index)  # Get the number of samples
half_samples = num_samples // 2  # Split samples into two groups (adjust as needed)
sample_info = pd.DataFrame({
    'sample_id': phospho_log2.index,
    'group': ['Normal']*half_samples + ['Tumor']*(num_samples - half_samples) 
}).set_index('sample_id')

# Create design matrix (add intercept for limma-style analysis)
design = pd.get_dummies(sample_info['group'], drop_first=True)
design.insert(0, 'Intercept', 1)  # Add intercept column

# 3. LIMMA-Style Analysis Function -----------------------------------------------
def limma_py(phospho_data, design_matrix):
    """LIMMA-style analysis with robust type handling."""
    n_sites = phospho_data.shape[1]
    coefs = np.full(n_sites, np.nan)  # Initialize with NaNs
    stderrs = np.full(n_sites, np.nan)
    dfs = np.full(n_sites, np.nan)
    sigmas = np.full(n_sites, np.nan)

    # Get the coefficient name for the group (e.g., 'group_Tumor')
    group_coef = [col for col in design_matrix.columns if col != 'Intercept'][0]

    for i, site in enumerate(phospho_data.columns):
        # Extract response variable and ensure numeric type
        y = phospho_data[site].dropna().astype(np.float64)
        if y.empty:
            continue  # Skip if no data

        # Align design matrix with available samples
        X = design_matrix.loc[y.index].astype(np.float64)
        if X.empty:
            continue

        # Check for sufficient samples and matrix rank
        if len(y) < X.shape[1] + 1 or np.linalg.matrix_rank(X) < X.shape[1]:
            continue

        # Fit model and extract coefficients
        try:
            model = sm.OLS(y, X).fit()
            coefs[i] = model.params[group_coef]  # Use column name
            stderrs[i] = model.bse[group_coef]
            sigmas[i] = np.sqrt(model.mse_resid)
            dfs[i] = model.df_resid
        except:
            continue  # Skip sites causing errors

    # Empirical Bayes moderation and FDR
    valid = ~np.isnan(coefs)
    if not np.any(valid):
        return pd.DataFrame()  # No valid results

    s2_prior = np.nanmedian(sigmas[valid]**2)
    df_prior = np.sum(valid)

    s2_post = (df_prior * s2_prior + dfs[valid] * sigmas[valid]**2) / (df_prior + dfs[valid])
    mod_t = coefs[valid] / (stderrs[valid] * np.sqrt(s2_post / sigmas[valid]**2))

    pvals = 2 * t.cdf(-np.abs(mod_t), df_prior + dfs[valid])
    qvals = fdrcorrection(pvals)[1]

    return pd.DataFrame({
        'log2FC': coefs[valid],
        't_stat': mod_t,
        'p_value': pvals,
        'q_value': qvals,
        'sigma': sigmas[valid],
        'df_total': df_prior + dfs[valid]
    }, index=phospho_data.columns[valid])

# 4. Run Analysis -----------------------------------------------------------------
results = limma_py(phospho_log2, design)

# 5. Save and Visualize Results ---------------------------------------------------
# Save significant results (q < 0.05, |log2FC| > 1)
significant = results[(results['q_value'] < 0.05) & (results['log2FC'].abs() > 1)]
significant.to_csv("significant_phosphositess.csv")

# Volcano plot
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(10,6))
sns.scatterplot(
    x='log2FC', 
    y=-np.log10(results['q_value']),
    data=results,
    hue=results['q_value'] < 0.05,
    palette={True: 'red', False: 'gray'},
    alpha=0.7,
    s=50
)


print("Data Summary:\n", phospho_log2.describe())

# If all values are ≥0, it’s expected but not an error.
plt.title("Differential Phosphorylation Analysis")
plt.xlabel("log2(Fold Change)")
plt.ylabel("-log10(q-value)")
plt.axhline(-np.log10(0.05), linestyle='--', color='grey')
plt.axvline(-1, color='blue', linestyle=':')
plt.axvline(1, color='blue', linestyle=':')
plt.show()