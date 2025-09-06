#!/bin/bash  
source /opt/load-libs.sh 3
export PATH=/opt/intel/oneapi/intelpython/latest/bin:/home/jorge.gacitua/salidas/miniconda3/bin:/home/jorge.gacitua/salidas/miniconda3/condabin:/opt/mpich/mpich-3.4.2/intel/2021.4.0/bin:/opt/intel/oneapi/vtune/2021.7.1/bin64:/opt/intel/oneapi/vpl/2021.6.0/bin:/opt/intel/oneapi/mpi/2021.4.0/libfabric/bin:/opt/intel/oneapi/mpi/2021.4.0/bin:/opt/intel/oneapi/mkl/2021.4.0/bin/intel64:/opt/intel/oneapi/itac/2021.4.0/bin:/opt/intel/oneapi/inspector/2021.4.0/bin64:/opt/intel/oneapi/dpcpp-ct/2021.4.0/bin:/opt/intel/oneapi/dev-utilities/2021.4.0/bin:/opt/intel/oneapi/debugger/10.2.4/gdb/intel64/bin:/opt/intel/oneapi/compiler/2021.4.0/linux/lib/oclfpga/llvm/aocl-bin:/opt/intel/oneapi/compiler/2021.4.0/linux/lib/oclfpga/bin:/opt/intel/oneapi/compiler/2021.4.0/linux/bin/intel64:/opt/intel/oneapi/compiler/2021.4.0/linux/bin:/opt/intel/oneapi/clck/2021.4.0/bin/intel64:/opt/intel/oneapi/advisor/2021.4.0/bin64:/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/snap/bin
export LD_LIBRARY_PATH=/opt/intel/oneapi/compiler/2021.4.0/linux/compiler/lib/intel64_lin:$LD_LIBRARY_PATH

#PBS -N classification_batch
#PBS -l nodes=1:ppn=64
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -o logs/obs_error_batch.log
#PBS -V

cd /home/jorge.gacitua/salidas/Sensitivity_Experiments

# Load conda into the environment
source /home/jorge.gacitua/salidas/miniconda3/etc/profile.d/conda.sh  # Adjust if using Anaconda or Mambaforge

# Activate your conda environment
conda activate wrf_python

# Confirm environment is active and packages are available
which python
python -c "import netCDF4, wrf; print('netCDF4 and wrf-python OK')"

# Run the script
conda run -n wrf_python python -u run_multiple_conf.py  > log_run_multiple.log 2>&1

echo "=== All experiments finished==="
