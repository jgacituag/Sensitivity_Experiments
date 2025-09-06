#!/bin/bash
#PBS -N wrf_exp_1
#PBS -q larga
#PBS -l nodes=1:ppn=48
#PBS -o logs/output_1.log
#PBS -e logs/error_1.log

source /opt/load-libs.sh 1
cd /home/jorge.gacitua/salidas/Sensitivity_Experiments

# Load conda into the environment
source /home/jorge.gacitua/salidas/miniconda3/etc/profile.d/conda.sh  # Adjust if using Anaconda or Mambaforge

# Activate your conda environment
conda activate wrf_python

# Confirm environment is active and packages are available
which python
python -c "import netCDF4, wrf; print('netCDF4 and wrf-python OK')"

export OMP_NUM_THREADS=48

# Also set MPI environment if WRF was compiled with MPI support
export OMP_PROC_BIND=true
export OMP_PLACES=cores

python -u ./run_multiple_conf.py 200 400