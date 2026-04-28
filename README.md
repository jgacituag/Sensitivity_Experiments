
# WRF Supercell Sensitivity Experiments

This repository contains an automated, modular, and scalable workflow for running idealized supercell simulations using the Weather Research and Forecasting (WRF) model. It is designed to take parameter samples (e.g., generated via SALib), dynamically modify the WRF `input_sounding` and `namelist.input`, and execute batches of simulations in a High-Performance Computing (HPC) environment via PBS.

## Repository Structure

```text
Sensitivity_Experiments/
├── README.md                           # This guide
├── .gitignore                          # Ignored files (logs, NetCDF outputs, local configs)
├── environment.yml                     # Conda environment definition
├── Exploracion_Simulaciones.ipynb      # Jupyter notebook for ensemble sanity check and visualization
├── src/                      
│   ├── generate_sample.py              # Script to generate the SALib experimental design
│   ├── run_batch.py                    # Main orchestrator (loops through CSV and runs WRF)
│   └── wrf_utils.py                    # Physics and utility functions (sounding & namelist editors)
└── templates/                     
    ├── config.template.yaml            # Template for the master configuration file
    └── submit_pbs.template.sh          # Template for the PBS job submission script
``` 

## Setup & Installation

### 1. Clone the repository
```bash
git clone https://github.com/jgacituag/Sensitivity_Experiments.git
cd Sensitivity_Experiments
```

### 2. Environment Configuration
This project uses a Conda environment to manage Python dependencies. You must create and activate this environment before running the scripts.

Create the environment using the provided `environment.yml` file:
```bash
conda env create -f environment.yml
```

Once the installation is complete, activate the environment:
```bash
conda activate Sensitivity_Experiments_env
```

### 3. Local Configuration (Important)
To prevent local paths from being tracked by Git, this repository uses template files. You must create your local copies before running the model.

**Copy the configuration template:**
```bash
cp templates/config.template.yaml config.yaml
```
Open `config.yaml` and update the `paths` section with the absolute paths to your compiled WRF ideal/real executables and your desired output directory.

**Copy the PBS submission template:**
```bash
cp templates/submit_pbs.template.sh submit_pbs.sh
```
Open `submit_pbs.sh` and ensure the PBS directives (queues, nodes) match your HPC cluster's architecture.

## Usage

### Step 1: Generate the Experimental Design
Before running the simulations, generate the SALib parameter space. You can configure the number of parameters and total samples by editing the header of `src/generate_sample.py`.

Run the generator:
```bash
python src/generate_sample.py
```
This will create a CSV file inside the `data/` directory. Ensure the `csv_file` path in your `config.yaml` matches this new file.

### Step 2: Submit the Batch Job
This workflow is designed to run in batches inside a compute node. Submit the job to your cluster's queue system specifying the start and end rows of the CSV file you want to process.

For example, to run simulations from row 0 to 199 in the CSV:
```bash
qsub -v START=0,END=200 submit_pbs.sh
```

### How it works under the hood:
1. The PBS script requests the node and executes `src/run_batch.py`.
2. The script reads `config.yaml` and the specified rows from your CSV file.
3. For each row, it creates an isolated temporary directory.
4. It maps the normalized parameters to physical variables (wind shear, moisture, stability) using `src/wrf_utils.py`.
5. It runs `./ideal.exe` and `./wrf.exe`.
6. The final `wrfout` NetCDF file is moved to your output directory, and the temporary execution folder is deleted to save disk space.