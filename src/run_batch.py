import sys
import os
import csv
import yaml
import shutil
import argparse
import subprocess
import wrf_utils

def main():
    parser = argparse.ArgumentParser(description="Batch runner for idealized WRF simulations.")
    parser.add_argument("--start", type=int, required=True, help="Start index of the CSV file")
    parser.add_argument("--end", type=int, required=True, help="End index of the CSV file (exclusive)")
    args = parser.parse_args()

    # Load configuration
    with open('config.yaml', 'r') as file:
        conf = yaml.safe_load(file)

    csv_path = conf['paths']['csv_file']
    wrf_base = conf['paths']['wrf_base']
    out_dir_base = conf['paths']['output_dir']
    
    if not os.path.exists(csv_path):
        print(f"Error: CSV file not found at {csv_path}")
        sys.exit(1)

    with open(csv_path, newline='') as csvfile:
        csv_reader = csv.reader(csvfile)
        next(csv_reader)  # Skip header
        
        for row_number, row in enumerate(csv_reader):
            if args.start <= row_number < args.end:
                print(f"\n--- Starting simulation {row_number:03d} ---")
                
                # 1. Setup temporary directory for this specific run
                tmp_dir = os.path.join(out_dir_base, 'tmp', f"run_{row_number:03d}")
                os.makedirs(tmp_dir, exist_ok=True)
                
                # Setup final output directory
                final_out_dir = os.path.join(out_dir_base, f"run_{row_number:03d}")
                os.makedirs(final_out_dir, exist_ok=True)

                # 2. Copy WRF executables and base files
                # Added -r just in case there are necessary subdirectories or hidden tables
                os.system(f"cp -r {wrf_base}/* {tmp_dir}/")
                
                # 3. Map SALib row to physical WRF parameters
                current_params = wrf_utils.map_salib_to_params(row, conf['wrf_params'])
                
                # 4. Modify inputs
                wrf_utils.modify_input_sounding(os.path.join(tmp_dir, 'input_sounding'), current_params)
                wrf_utils.edit_namelist(os.path.join(tmp_dir, 'namelist.input'), current_params)
                
                # 5. Execute WRF
                os.chdir(tmp_dir)
                print("Running ideal.exe...")
                subprocess.run(["./ideal.exe"], stdout=open('out_ideal.log', 'w'), stderr=subprocess.STDOUT)
                
                print("Running wrf.exe...")
                # Ensure OpenMP threads are set for this specific subprocess
                my_env = os.environ.copy()
                my_env["OMP_NUM_THREADS"] = str(conf['wrf_params'].get('nthreads', 48))
                subprocess.run(["./wrf.exe"], stdout=open('out_wrf.log', 'w'), stderr=subprocess.STDOUT, env=my_env)
                
                # 6. Move outputs and logs
                out_nc = f"wrfout_{row_number:03d}.nc"
                wrf_output_name = "./wrfout_d01_0001-01-01_00:00:00"

                if os.path.exists(wrf_output_name):
                    os.system(f"mv {wrf_output_name} {os.path.join(final_out_dir, out_nc)}")
                    print(f"Success: WRF output generated for {row_number:03d}.")
                else:
                    os.system(f"cp out_ideal.log out_wrf.log {final_out_dir}/")
                    if os.path.exists("rsl.error.0000"):
                        os.system(f"cp rsl.* {final_out_dir}/")
                    print(f"ERROR: WRF failed. Check logs in {final_out_dir}")
                # Go back to base dir and remove temp folder to save space
                os.chdir(out_dir_base)
                shutil.rmtree(tmp_dir)
                
                print(f"--- Completed simulation {row_number:03d} ---")

if __name__ == "__main__":
    main()