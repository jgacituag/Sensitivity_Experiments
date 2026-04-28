#!/bin/bash
#PBS -N wrf_exp_batch
#PBS -q larga
#PBS -l nodes=1:ppn=48
#PBS -o logs/output_${PBS_JOBID}.log
#PBS -e logs/error_${PBS_JOBID}.log

# Cargar librerías del sistema
source /opt/load-libs.sh 1

# Moverse al directorio desde donde se envió el trabajo
cd $PBS_O_WORKDIR

# Activar el entorno conda
eval "$(conda shell.bash hook)"
conda activate wrf_python

# Verificación de dependencias
python -c "import netCDF4, wrf; print('netCDF4 and wrf-python OK')"

# Configuración de paralelización
export OMP_NUM_THREADS=48
export OMP_PROC_BIND=true
export OMP_PLACES=cores

# Ejecución del script principal
# Los argumentos (ej. 0 200) pueden pasarse como variables al hacer qsub, 
# o dejarse fijos si se crearán múltiples scripts queue_X.sh
python -u ./run_multiple_conf.py 0 200