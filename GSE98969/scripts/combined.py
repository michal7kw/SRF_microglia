# %%
import os
import pandas as pd
import numpy as np
import gzip
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
import gc

# %%
def load_geo_data(file_path):
    with gzip.open(file_path, 'rt') as f:
        data = pd.read_csv(f, sep='\t', index_col=0)
    return data

# %%
# Load DAM data from GSE98969 and aggregate to bulk-like samples
dam_files = [os.path.join("./data", f) for f in os.listdir("./data") if f.endswith('.txt.gz') and f.startswith('GSM')]
dam_data_list = []
for file in dam_files:
    sample_data = load_geo_data(file)
    sample_name = os.path.basename(file).split('_')[0]  # Use GSM ID as sample name
    # Sum expression across all cells to create a bulk-like sample
    bulk_like_sample = sample_data.sum(axis=1)
    dam_data_list.append(pd.DataFrame(bulk_like_sample, columns=[sample_name]))

dam_data = pd.concat(dam_data_list, axis=1)
# Save the concatenated data to a CSV file
output_file = "dam_data_bulk_like.csv"
dam_data.to_csv(output_file)
print(f"Data saved to {output_file}")

# %%
dam_data.head()

# %%
def create_distribution_plots(data, title_prefix):
    # Boxplot
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=data)
    plt.title(f"{title_prefix} - Boxplot of Gene Expression")
    plt.xlabel("Samples")
    plt.ylabel("Log2 Expression")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f"{title_prefix.lower().replace(' ', '_')}_boxplot.png")
    plt.close()

    # Density plot
    plt.figure(figsize=(12, 6))
    for column in data.columns:
        sns.kdeplot(data=data[column], label=column)
    plt.title(f"{title_prefix} - Density Plot of Gene Expression")
    plt.xlabel("Log2 Expression")
    plt.ylabel("Density")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"{title_prefix.lower().replace(' ', '_')}_density.png")
    plt.close()

# %%
create_distribution_plots(np.log2(dam_data + 1), "Unfiltered Bulk-like Data")


# %%
gc.collect()

# %%
print("Summary statistics for unfiltered bulk-like data:")
print(dam_data.describe())

# %%
def preprocess_data(data):
    # Log2 transform
    data_norm = np.log2(data + 1)
    # Filter low-expression genes
    data_filtered = data_norm[data_norm.mean(axis=1) > 1]
    return data_filtered

# aged_data_filtered = preprocess_data(aged_data)
dam_data_filtered = preprocess_data(dam_data)

# %%
create_distribution_plots(dam_data_filtered, "Filtered Bulk-like Data")

# %%
print("\nSummary statistics for filtered bulk-like data:")
print(dam_data_filtered.describe())

# %%
# common_genes = set(aged_data_filtered.index) & set(dam_data_filtered.index)
# aged_common = aged_data_filtered.loc[common_genes]
# dam_common = dam_data_filtered.loc[common_genes]

# %%
def calculate_similarity(data1, data2):
    correlations = []
    for col in data2.columns:
        corr, _ = pearsonr(data1.iloc[:, 0], data2[col])
        correlations.append(corr)
    return pd.Series(correlations, index=data2.columns)

similarity_index = calculate_similarity(aged_common, dam_common)

# %%
plt.figure(figsize=(10, 6))
sns.histplot(similarity_index, kde=True)
plt.title("Distribution of Similarity Index (Pearson Correlation)")
plt.xlabel("Correlation Coefficient")
plt.ylabel("Frequency")
plt.savefig("similarity_distribution.png")
plt.close()

print("\nSummary of Similarity Index:")
print(similarity_index.describe())

# %%
plt.figure(figsize=(12, 8))
combined_data = pd.concat([aged_common.iloc[:, 0], dam_common], axis=1)
sns.clustermap(combined_data, cmap="viridis", z_score=0, standard_scale=1)
plt.title("Gene Expression Heatmap")
plt.savefig("expression_heatmap.png")
plt.close()

similarity_index.to_csv("similarity_index.csv")
print("Analysis complete. Results saved to CSV and plots generated.")



