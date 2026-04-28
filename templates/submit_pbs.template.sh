#!/bin/bash
#PBS -N wrf_batch
#PBS -q larga
#PBS -l nodes=1:ppn=48
#PBS -o output_${PBS_JOBID}.log
#PBS -e error_${PBS_JOBID}.log

# Load system libraries
source /opt/load-libs.sh 1

# Move to the directory where the job was submitted
cd $PBS_O_WORKDIR

# Activate Conda environment
eval "$(conda shell.bash hook)"
conda activate Sensitivity_Experiments_env

# Set OpenMP threads
export OMP_NUM_THREADS=48
export OMP_PROC_BIND=true
export OMP_PLACES=cores

# The arguments $1 and $2 correspond to the start and end row index of the CSV
# Usage: qsub -v "1=0,2=200" submit_pbs.sh
python -u src/run_batch.py --start $1 --end $2

