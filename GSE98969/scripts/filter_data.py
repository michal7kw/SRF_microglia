# %%
import os
import pandas as pd
import numpy as np
import gzip
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
import gc

output_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/GSE98969"

# Load the CSV file
dam_data = pd.read_csv(os.path.join(output_dir, "dam_data_bulk_like.csv"), index_col=0)

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
create_distribution_plots(os.path.join(output_dir, "output", dam_data_filtered), "Filtered Bulk-like Data")
