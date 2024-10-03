#!/bin/bash
#SBATCH --job-name=load_data
#SBATCH --account kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mail-type=END ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=kubacki.michal@hst.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/logs/test_script.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/logs/test_script.out"

# Load the appropriate conda environment (if needed)
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Run your Python script
python /beegfs/scratch/ric.broccoli/kubacki.michal/GSE98969/load_data.py