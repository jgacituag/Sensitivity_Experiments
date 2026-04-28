#!/bin/bash
#PBS -N wrf_batch
#PBS -q larga
#PBS -l nodes=1:ppn=48
#PBS -o path_to_output_log
#PBS -e path_to_output_log

REPO=path_to_repo_root
LOGDIR="$REPO"/logs

# Capturar todo stdout/stderr al log independientemente de PBS
exec > >(tee -a $LOGDIR/wrf_batch.log) 2>&1
# Load system libraries
source /opt/load-libs.sh 1

# Move to the directory where the job was submitted
cd "$REPO"

# Activate Conda environment

source path_to_conda/etc/profile.d/conda.sh
conda activate Sensitivity_Experiments_env

# Set OpenMP threads
export OMP_NUM_THREADS=48
export OMP_PROC_BIND=true
export OMP_PLACES=cores

python -u src/run_batch.py --start $START --end $END

