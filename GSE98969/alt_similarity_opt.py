# %% [markdown]
# # Environment

# %%
import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob
from scipy import stats
import seaborn as sns
import gc
from tqdm import tqdm

# Set random seed for reproducibility
np.random.seed(42)

cluster = True

# %%
import rpy2
from rpy2.robjects import pandas2ri
pandas2ri.activate()

# %%
if cluster:
    wd_dir = '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_microglia/GSE98969'
else:
    wd_dir = '/home/michal/WSL_GitHub/SRF_microglia/GSE98969'
os.chdir(wd_dir)

# %%
if cluster:
    sharon_dir = '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_microglia/Sharon_RNA/samples/all/runs/all/fastq/merge-by-read/trimmed/trimmomatic/mapped/STAR/merged/featureCounts/merged'
else:
    sharon_dir = '/home/michal/WSL_GitHub/SRF_microglia/data/sharon_rna'

# %%
%load_ext rpy2.ipython

# %%
%%R
x <- c(1, 2, 3, 4, 5)
mean(x)

# %% [markdown]
# # Load GSE98969 data

# %%
# Get list of all MARS-seq files
mars_seq_files = glob('./data/extracted/GSM*.txt.gz')
print(mars_seq_files)

# %%
# Read metadata
metadata = pd.read_csv('./data/metadata_csv.csv')
metadata_dict = metadata.set_index('geo_accession').to_dict('index')
metadata_dict[list(metadata_dict.keys())[0]]

# %%
# def read_mars_seq_file_in_chunks(filename, chunksize=10000):
#     total_rows = sum(1 for _ in pd.read_csv(filename, sep='\t', compression='gzip', chunksize=chunksize))
    
#     pbar = tqdm(total=total_rows, desc="Reading file", unit="rows")
#     for chunk in pd.read_csv(filename, sep='\t', compression='gzip', index_col=0, chunksize=chunksize):
#         pbar.update(len(chunk))
#         yield chunk
#     pbar.close()

# def process_mars_seq_files(mars_seq_files, metadata_dict, chunksize=10000):
#     adata = None
#     obs_names_counter = {}
    
#     for file in tqdm(mars_seq_files, desc="Processing files", unit="file"):
#         print(f"\nProcessing file: {file}")
#         gsm_id = os.path.basename(file).split('_')[0]
        
#         for i, chunk in enumerate(read_mars_seq_file_in_chunks(file, chunksize)):
#             # Make observation names unique within the file
#             chunk.index = [f"{gsm_id}_{idx}" for idx in chunk.index]
            
#             # Ensure observation names are unique across all files
#             new_obs_names = []
#             for name in chunk.index:
#                 if name in obs_names_counter:
#                     obs_names_counter[name] += 1
#                     new_name = f"{name}_{obs_names_counter[name]}"
#                 else:
#                     obs_names_counter[name] = 0
#                     new_name = name
#                 new_obs_names.append(new_name)
            
#             chunk.index = new_obs_names
            
#             current_adata = ad.AnnData(chunk.T)
            
#             # Add metadata
#             if gsm_id in metadata_dict:
#                 for key, value in metadata_dict[gsm_id].items():
#                     current_adata.obs[key] = value
#             else:
#                 print(f"Warning: No metadata found for {gsm_id}")
            
#             if adata is None:
#                 adata = current_adata
#             else:
#                 adata = ad.concat([adata, current_adata], join='outer', fill_value=0)
            
#             # Clear temporary variables
#             del chunk, current_adata
#             gc.collect()
        
#         gc.collect()
    
#     return adata

# %%
# adata = process_mars_seq_files(mars_seq_files, metadata_dict)
# print(adata)

# %%
# Function to read a single MARS-seq file
def read_mars_seq_file(filename):
    df = pd.read_csv(filename, sep='\t', compression='gzip', index_col=0)
    return df


# %%
# Read all MARS-seq files and store them in a list
adatas = []
for file in mars_seq_files:
    df = read_mars_seq_file(file)
    adata = ad.AnnData(df.T)
    adata.var_names_make_unique()
    
    # Extract GSM ID from filename
    gsm_id = os.path.basename(file).split('_')[0]
    
    # Add metadata
    if gsm_id in metadata_dict:
        for key, value in metadata_dict[gsm_id].items():
            adata.obs[key] = value
    else:
        print(f"Warning: No metadata found for {gsm_id}")
    
    adatas.append(adata)

# %%
gc.collect()

# %%
# Concatenate all AnnData objects
adata = ad.concat(adatas, join='outer', fill_value=0)
adata

# %%
adata.X

# %%
adata.obs.head()

# %%
# Remove prefixes
adata.obs['Treatment'] = adata.obs['Treatment'].str.replace('treatment: ', '')
adata.obs['age'] = adata.obs['age'].str.replace('mouse age: ', '')
adata.obs['strain'] = adata.obs['strain'].str.replace('strain: ', '')
adata.obs['organ'] = adata.obs['organ'].str.replace('organ: ', '')


# %%
# Basic preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)


# %%
# Calculate quality control metrics
adata.var['mt'] = adata.var_names.str.startswith('mt-')  # Identify mitochondrial genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# %%
# Plot QC metrics
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

# %%
# Filter cells based on QC metrics (adjust these thresholds as needed)
adata = adata[adata.obs.n_genes_by_counts < 3000, :]
adata = adata[adata.obs.total_counts < 10000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]

# %%
# Plot QC metrics
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

# %%
# Save raw counts in a new layer
adata.layers['counts'] = adata.X.copy()

# Normalize data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# %%
# Identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Plot highly variable genes
sc.pl.highly_variable_genes(adata)

# %%
# Set the .raw attribute of AnnData object to the normalized and logarithmized raw gene expression
adata.raw = adata

# Scale data
sc.pp.scale(adata, max_value=10)

# %%
# Perform PCA
sc.tl.pca(adata, svd_solver='arpack')

# %%
sc.pl.pca_variance_ratio(adata, n_pcs=50)

# %%
# Plot PCA results
sc.pl.pca(adata, color=['Treatment', 'age', 'region', 'strain', 'organ'])

# %%
# Compute neighborhood graph

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
# Perform UMAP Embedding
sc.tl.umap(adata)

# %%
# Plot UMAP
sc.pl.umap(adata, color=['Treatment', 'age', 'region', 'strain', 'organ'])

# %%
# Perform clustering
sc.tl.leiden(adata)

# %%
# Plot clustering results
sc.pl.umap(adata, color=['leiden', 'Treatment', 'age', 'region', 'strain', 'organ'])

# %%
%%capture
# Find marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')

# %%
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# %%
# Save the results
adata.write('./output/GSE98969_microglia_results.h5ad')

print("Analysis complete. Results saved to 'GSE98969_microglia_results.h5ad'.")

# %%
adata.X

# %%
# Import necessary libraries
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import matplotlib.pyplot as plt
import seaborn as sns

# Create a count matrix from the AnnData object
# We need to use the raw counts, not the normalized data
if 'counts' in adata.layers:
    count_matrix = pd.DataFrame(adata.layers['counts'].T, index=adata.var_names, columns=adata.obs_names)
else:
    print("Warning: Raw counts not found. Attempting to reconstruct from normalized data.")
    # Reverse log-normalization (assuming log1p was used) and round to integers
    count_matrix = pd.DataFrame(np.round(np.expm1(adata.X)).T, index=adata.var_names, columns=adata.obs_names)

# Ensure all values are non-negative integers
count_matrix = count_matrix.clip(lower=0).astype(int)

# Create a sample information dataframe
sample_info = adata.obs[['leiden']].copy()
sample_info.columns = ['condition']

# Ensure that count_matrix and sample_info have the same number of samples
common_samples = count_matrix.columns.intersection(sample_info.index)
count_matrix = count_matrix[common_samples]
sample_info = sample_info.loc[common_samples]

print(f"Number of samples in count matrix: {count_matrix.shape[1]}")
print(f"Number of samples in sample info: {len(sample_info)}")

# Create DeseqDataSet object
dds = DeseqDataSet(
    counts=count_matrix,
    metadata=sample_info,
    design_factors="condition",
    refit_cooks=True
)

# Run DESeq2 analysis
dds.deseq2()

# Get results
stat_res = DeseqStats(dds, contrast=["condition", "1", "0"])
stat_res.summary()
results = stat_res.results_df

# Filter for significantly differentially expressed genes
significant_genes = results[(results['padj'] < 0.05) & (abs(results['log2FoldChange']) > 1)]

print(f"Number of significantly differentially expressed genes: {len(significant_genes)}")

# Save results to CSV
results.to_csv('deseq2_results.csv')
significant_genes.to_csv('deseq2_significant_genes.csv')

# Create a volcano plot
plt.figure(figsize=(10, 8))
sns.scatterplot(data=results, x='log2FoldChange', y='-log10(pvalue)', 
                hue=(results['padj'] < 0.05) & (abs(results['log2FoldChange']) > 1))
plt.title('Volcano Plot of DESeq2 Results')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 P-value')
plt.savefig('deseq2_volcano_plot.png')
plt.close()

print("DESeq2 analysis complete. Results saved to 'deseq2_results.csv' and 'deseq2_significant_genes.csv'.")
print("Volcano plot saved as 'deseq2_volcano_plot.png'.")


# %%


# %%


# %%


# %%


# %%


# %%


# %%


# %%


# %%
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def load_data(file_path):
    """Load the saved results."""
    adata = sc.read(file_path)
    print(f"Data loaded successfully. Shape: {adata.shape}")
    print(f"Observations: {list(adata.obs_keys())}")
    print(f"Variables: {list(adata.var_keys())}")
    return adata

def prepare_data(adata):
    """Prepare the data for analysis."""
    adata.layers['log_norm'] = adata.X.copy()
    adata.X = adata.layers['counts'].copy()
    return adata

def create_pseudobulk(adata, group_by):
    """Create pseudo-bulk data from single-cell data."""
    adata.obs[group_by] = pd.Categorical(adata.obs[group_by])
    indicator = pd.get_dummies(adata.obs[group_by])
    
    pseudobulk = ad.AnnData(
        X=indicator.values.T @ adata.X,
        obs=pd.DataFrame(index=indicator.columns),
        var=adata.var.copy()
    )
    
    pseudobulk.obs = adata.obs.groupby(group_by).first()
    
    for layer in adata.layers.keys():
        pseudobulk.layers[layer] = indicator.values.T @ adata.layers[layer]
    
    return pseudobulk

def normalize_pseudobulk(pseudobulk):
    """Normalize pseudo-bulk data."""
    pseudobulk_norm = pseudobulk.copy()
    pseudobulk_norm.X = pseudobulk_norm.X / pseudobulk_norm.X.sum(axis=1, keepdims=True) * 1e6
    pseudobulk_norm.layers['counts'] = pseudobulk_norm.layers['counts'] / pseudobulk_norm.layers['counts'].sum(axis=1, keepdims=True) * 1e6
    pseudobulk_norm.layers['log_norm'] = np.log1p(pseudobulk_norm.X)
    return pseudobulk_norm

def log_transform(pseudobulk_norm):
    """Log transform the normalized data."""
    pseudobulk_log = pseudobulk_norm.copy()
    pseudobulk_log.X = np.log2(pseudobulk_log.X + 1)
    return pseudobulk_log

def perform_de(adata, group1, group2, min_samples=3):
    """
    Perform differential expression analysis.
    
    Parameters:
    - adata: AnnData object
    - group1, group2: lists of sample names for each group
    - min_samples: minimum number of samples required in each group to perform t-test
    """
    genes = []
    pvalues = []
    log2fc = []
    
    for gene in adata.var_names:
        group1_data = adata[:, gene].X[adata.obs.index.isin(group1)].flatten()
        group2_data = adata[:, gene].X[adata.obs.index.isin(group2)].flatten()
        
        # Remove zero values
        group1_data = group1_data[group1_data != 0]
        group2_data = group2_data[group2_data != 0]
        
        # Check if we have enough non-zero samples
        if len(group1_data) >= min_samples and len(group2_data) >= min_samples:
            # Perform t-test
            t_stat, p_value = stats.ttest_ind(group1_data, group2_data)
            
            # Calculate log2 fold change
            mean1 = np.mean(group1_data) if len(group1_data) > 0 else 1e-9
            mean2 = np.mean(group2_data) if len(group2_data) > 0 else 1e-9
            log2fc_value = np.log2(mean1 / mean2)
            
            genes.append(gene)
            pvalues.append(p_value)
            log2fc.append(log2fc_value)
        else:
            # If not enough samples, add NaN values
            genes.append(gene)
            pvalues.append(np.nan)
            log2fc.append(np.nan)
    
    return pd.DataFrame({'gene': genes, 'pvalue': pvalues, 'log2fc': log2fc})

def plot_volcano(de_results):
    """Create and save a volcano plot."""
    plt.figure(figsize=(10, 8))
    plt.scatter(de_results['log2fc'], -np.log10(de_results['pvalue']), alpha=0.5)
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 P-value')
    plt.title('Volcano Plot: 5XFAD vs C57BL/6')
    plt.show()

def plot_heatmap(adata, genes):
    """Create and save a heatmap of top differentially expressed genes."""
    plt.figure(figsize=(12, 10))
    sc.pl.heatmap(adata, var_names=genes, groupby='strain', show_gene_labels=True, cmap='viridis', dendrogram=False)
    plt.show()
    

def save_results(pseudobulk, pseudobulk_norm, pseudobulk_log, de_results):
    """Save analysis results to CSV files."""
    pseudobulk.write_csvs('pseudobulk_raw_counts.csv')
    pseudobulk_norm.write_csvs('pseudobulk_normalized.csv')
    pd.DataFrame(pseudobulk_log.X.T, index=pseudobulk_log.var_names, columns=pseudobulk_log.obs_names).to_csv('pseudobulk_log_transformed.csv')
    de_results.to_csv('differential_expression_results.csv')


# %%
# Load and prepare data
adata = load_data('./output/GSE98969_microglia_results.h5ad')

# %%
adata

# %%
adata = prepare_data(adata)

# %%
# Create and process pseudo-bulk data
pseudobulk = create_pseudobulk(adata, 'strain')
pseudobulk_norm = normalize_pseudobulk(pseudobulk)
pseudobulk_log = log_transform(pseudobulk_norm)

# %%
# Perform differential expression analysis
de_results = perform_de(pseudobulk_log, ['5XFAD'], ['C57BL/6'])
de_results = de_results.sort_values('pvalue')

# %%
de_results

# %%
plot_volcano(de_results)

# %%
# plot_heatmap(pseudobulk_log, de_results.head(50)['gene'])

# %%
# Create pseudo-bulk data
def create_pseudobulk(adata, group_by):
    """
    This function creates pseudo-bulk data from single-cell data.
    
    Parameters:
    - adata: AnnData object containing single-cell data
    - group_by: String, the column name in adata.obs to group cells by
    
    The function does the following:
    1. Converts the grouping column to categorical data type.
    2. Creates an indicator matrix for each group.
    3. Sums the counts for each group to create pseudo-bulk data.
    4. Maintains relevant metadata from the original data.
    5. Preserves all layers from the original data in the pseudo-bulk data.
    
    Returns:
    - pseudobulk: AnnData object containing the pseudo-bulk data
    """
    # Convert group_by column to categorical if it's not already
    adata.obs[group_by] = pd.Categorical(adata.obs[group_by])
    
    # Create indicator matrix
    indicator = pd.get_dummies(adata.obs[group_by])
    
    # Sum the counts for each group
    pseudobulk = ad.AnnData(
        X=indicator.values.T @ adata.X,
        obs=pd.DataFrame(index=indicator.columns),
        var=adata.var.copy()
    )
    
    # Maintain relevant metadata
    pseudobulk.obs = adata.obs.groupby(group_by).first()
    
    # Ensure the layers are preserved
    for layer in adata.layers.keys():
        pseudobulk.layers[layer] = indicator.values.T @ adata.layers[layer]
    
    return pseudobulk

# %%
pseudobulk = create_pseudobulk(adata, 'strain')

# %%
pseudobulk

# %%
pseudobulk.X

# %%
# # Normalize pseudo-bulk data
# pseudobulk_norm = pseudobulk.copy()
# pseudobulk_norm.X = pseudobulk_norm.X / pseudobulk_norm.X.sum(axis=1, keepdims=True) * 1e6  # CPM normalization
# pseudobulk_norm.layers['counts'] = pseudobulk_norm.layers['counts'] / pseudobulk_norm.layers['counts'].sum(axis=1, keepdims=True) * 1e6  # CPM normalization for counts layer
# pseudobulk_norm.layers['log_norm'] = np.log1p(pseudobulk_norm.X)  # Log-normalize the main matrix
# pseudobulk_norm.X

# %%
# Log transform
pseudobulk_log = np.log2(pseudobulk_norm + 1)

# %%
pseudobulk_log.head()

# %%
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
data = pd.read_csv('pseudobulk_log.tsv', sep='\t', index_col=0)

# Create AnnData object
adata = sc.AnnData(data.T)

# Add gene names to var
adata.var_names = data.index

# Add condition information to obs
adata.obs['condition'] = adata.obs.index

# Perform differential expression analysis
sc.tl.rank_genes_groups(adata, 'condition', method='wilcoxon')

# Get results
results = adata.uns['rank_genes_groups']
groups = results['names'].dtype.names

# Create a DataFrame with the results
def get_df(key):
    return pd.DataFrame({group + '_' + key: results[key][group] for group in groups})

results_df = pd.concat([get_df(key) for key in ['names', 'scores', 'pvals', 'pvals_adj']], axis=1)

# Save results to CSV
results_df.to_csv('differential_expression_results.csv')

# Visualize top differentially expressed genes
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
plt.savefig('top_differentially_expressed_genes.png')

# Create a heatmap of top differentially expressed genes
top_genes = results_df.iloc[:50, results_df.columns.get_level_values(1)=='names'].values.flatten()
sc.pl.heatmap(adata, top_genes, groupby='condition', show_gene_labels=True, figsize=(12, 8))
plt.savefig('heatmap_top_differentially_expressed_genes.png')

print("Differential expression analysis completed. Results saved to 'differential_expression_results.csv'")
print("Visualizations saved as 'top_differentially_expressed_genes.png' and 'heatmap_top_differentially_expressed_genes.png'")

# %%
# Perform differential expression analysis
def perform_de(data, group1, group2):
    genes = []
    pvalues = []
    log2fc = []
    for gene in data.index:
        t_stat, p_value = stats.ttest_ind(data.loc[gene, group1], data.loc[gene, group2])
        genes.append(gene)
        pvalues.append(p_value)
        log2fc.append(np.log2(data.loc[gene, group1].mean() / data.loc[gene, group2].mean()))
    return pd.DataFrame({'gene': genes, 'pvalue': pvalues, 'log2fc': log2fc})


# %%
# Example: Differential expression between 5XFAD and C57BL/6
de_results = perform_de(pseudobulk_log, ['5XFAD'], ['C57BL/6'])
de_results = de_results.sort_values('pvalue')

# %%
de_results.head()

# %%
# Volcano plot
plt.figure(figsize=(10, 8))
plt.scatter(de_results['log2fc'], -np.log10(de_results['pvalue']), alpha=0.5)
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 P-value')
plt.title('Volcano Plot: 5XFAD vs C57BL/6')
plt.show()

# Heatmap of top differentially expressed genes
top_genes = de_results.head(50)['gene']
plt.figure(figsize=(12, 10))
sns.heatmap(pseudobulk_log.loc[top_genes], cmap='viridis', center=0)
plt.title('Top 50 Differentially Expressed Genes')
plt.show()

# Save pseudo-bulk data
pseudobulk.to_csv('pseudobulk_raw_counts.csv')
pseudobulk_norm.to_csv('pseudobulk_normalized.csv')
pseudobulk_log.to_csv('pseudobulk_log_transformed.csv')
de_results.to_csv('differential_expression_results.csv')

print("Pseudo-bulk analysis complete. Results saved to CSV files and plots.")

# %%



