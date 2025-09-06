#!/bin/bash
#PBS -N wrf_exps
#PBS -t 0-4              # 5 jobs only
#PBS -l nodes=1:ppn=48
#PBS -q larga
#PBS -o logs/output_$PBS_ARRAYID.log
#PBS -e logs/error_$PBS_ARRAYID.log

source /opt/load-libs.sh 3
export PATH=/opt/intel/oneapi/intelpython/latest/bin:/home/jorge.gacitua/salidas/miniconda3/bin:/home/jorge.gacitua/salidas/miniconda3/condabin:/opt/mpich/mpich-3.4.2/intel/2021.4.0/bin:/opt/intel/oneapi/vtune/2021.7.1/bin64:/opt/intel/oneapi/vpl/2021.6.0/bin:/opt/intel/oneapi/mpi/2021.4.0/libfabric/bin:/opt/intel/oneapi/mpi/2021.4.0/bin:/opt/intel/oneapi/mkl/2021.4.0/bin/intel64:/opt/intel/oneapi/itac/2021.4.0/bin:/opt/intel/oneapi/inspector/2021.4.0/bin64:/opt/intel/oneapi/dpcpp-ct/2021.4.0/bin:/opt/intel/oneapi/dev-utilities/2021.4.0/bin:/opt/intel/oneapi/debugger/10.2.4/gdb/intel64/bin:/opt/intel/oneapi/compiler/2021.4.0/linux/lib/oclfpga/llvm/aocl-bin:/opt/intel/oneapi/compiler/2021.4.0/linux/lib/oclfpga/bin:/opt/intel/oneapi/compiler/2021.4.0/linux/bin/intel64:/opt/intel/oneapi/compiler/2021.4.0/linux/bin:/opt/intel/oneapi/clck/2021.4.0/bin/intel64:/opt/intel/oneapi/advisor/2021.4.0/bin64:/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/snap/bin
export LD_LIBRARY_PATH=/opt/intel/oneapi/compiler/2021.4.0/linux/compiler/lib/intel64_lin:$LD_LIBRARY_PATH

cd $PBS_O_WORKDIR

module load mpi  # adjust to your HPC modules

# Each job runs 200 experiments
CHUNK=200
START=$(( PBS_ARRAY_INDEX * CHUNK ))
END=$(( START + CHUNK ))

python run_single_conf.py $START $END