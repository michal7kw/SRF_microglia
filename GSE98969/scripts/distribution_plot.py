# %%
import os
import pandas as pd
import numpy as np
import gzip
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
import gc

output_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/GSE98969/output"
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
    plt.savefig(os.path.join(output_dir, f"{title_prefix.lower().replace(' ', '_')}_boxplot.png"))
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
    plt.savefig(os.path.join(output_dir, f"{title_prefix.lower().replace(' ', '_')}_density.png"))
    plt.close()

# %%
# Load the CSV file
dam_data = pd.read_csv(os.path.join(output_dir, "dam_data_bulk_like.csv"), index_col=0)

# Create distribution plots for the loaded data
create_distribution_plots(dam_data, "Bulk-like Data")

# %%
create_distribution_plots(np.log2(dam_data + 1), "Unfiltered Bulk-like Data")