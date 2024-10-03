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

# %%
def load_geo_data(file_path):
    with gzip.open(file_path, 'rt') as f:
        data = pd.read_csv(f, sep='\t', index_col=0)
    return data

# %%
# Load DAM data from GSE98969 and aggregate to bulk-like samples
dam_files = [os.path.join(os.path.join(output_dir, "data"), f) for f in os.listdir(os.path.join(output_dir, "data")) if f.endswith('.txt.gz') and f.startswith('GSM')]
dam_data_list = []
for file in dam_files:
    sample_data = load_geo_data(file)
    sample_name = os.path.basename(file).split('_')[0]  # Use GSM ID as sample name
    # Sum expression across all cells to create a bulk-like sample
    bulk_like_sample = sample_data.sum(axis=1)
    dam_data_list.append(pd.DataFrame(bulk_like_sample, columns=[sample_name]))

dam_data = pd.concat(dam_data_list, axis=1)
# Save the concatenated data to a CSV file
output_file = os.path.join(output_dir, "output", "dam_data_bulk_like.csv")
dam_data.to_csv(output_file)
print(f"Data saved to {output_file}")