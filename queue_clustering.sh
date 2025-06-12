#!/bin/bash
#PBS -N classification_batch
#PBS -l nodes=1:ppn=64
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -o logs/obs_error_batch.log
#PBS -V

cd /home/jorge.gacitua/experimentos/Sensitivity_Experiments

# Load conda into the environment
source ~/miniconda3/etc/profile.d/conda.sh  # Adjust if using Anaconda or Mambaforge

# Activate your conda environment
conda activate wrf_python

# Confirm environment is active and packages are available
which python
python -c "import netCDF4, wrf; print('netCDF4 and wrf-python OK')"

# Run the script
python -u classification.py > log.log 2>&1

echo "=== All experiments finished ==="