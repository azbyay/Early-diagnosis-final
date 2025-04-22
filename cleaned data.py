import pandas as pd

# Load data (replace with your actual file paths)
proteome_df = pd.read_csv("/Users/awiener/Downloads/CPTAC3_Pancreatic_Ductal_Adenocarcinoma_Proteome.tmt11.tsv")
phosphoproteome_df = pd.read_csv("/Users/awiener/Downloads/CPTAC3_Pancreatic_Ductal_Adenocarcinoma_Phosphoproteome.phosphosite.tmt11 (1).tsv")

# Remove "Unshared Log Ratio" columns from proteome_df
proteome_df = proteome_df.drop(columns=[col for col in proteome_df.columns if "Unshared Log Ratio" in col])

# Remove "Withdrawn Log Ratio" column from phosphoproteome_df
phosphoproteome_df = phosphoproteome_df.drop(columns=["Withdrawn Log Ratio"], errors='ignore')  # errors='ignore' handles cases where the column might not exist


# Save the cleaned data to new files
proteome_df.to_csv("cleaned_proteome.csv", index=False)
phosphoproteome_df.to_csv("cleaned_phosphoproteome.csv", index=False)
