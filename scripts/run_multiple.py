import sys
import csv
import sensiwrf.wrf_module as wrf
from sensiwrf import config as conf
import os

def run_wrf_single_conf(row, row_number):
    param1, param2, param3, param4, param5 = map(float, row)
    conf_dict = conf.conf
    conf_dict['expname']   = 'Multi_5_vars_saltelli'
    conf_dict['run_num']   = row_number
    conf_dict['modelpath'] = '/home/jorge.gacitua/salidas/em_quarter_ss'
    conf_dict['datapath']  = '/home/jorge.gacitua/salidas'

    # Wind profile
    conf_dict['modify_wind_profile'] = True
    conf_dict['remove_mean_wind']    = True
    conf_dict['shear_type']          = 'Curved'
    conf_dict['total_shear_depth']   = 4000+(param1*(8000-4000))
    conf_dict['int_total_shear']     = 0.00001+(param2*(30-0))
    conf_dict['curved_shear_per']    = param3

    conf_dict['shear_depth_u']    = 0.0
    conf_dict['shear_strength_u'] = 5.0e-3
    conf_dict['shear_depth_v']    = 8000.0
    conf_dict['shear_strehgth_v'] = 0.0

    conf_dict['llj_amp']   = 0.0
    conf_dict['llj_h']     = 1500.0
    conf_dict['llj_width'] = 500.0
    conf_dict['llj_dir']   = 360.0

    # Stability
    conf_dict['modify_stability'] = True
    conf_dict['stability_factor'] = -0.5+(param4*(3.0--0.5))
    conf_dict['stability_factor_height'] = 5000 #2500+(param5*(10000-2500))#Height of maximum warming / cooling. (tipical range [2500 , 10000 ] )

    # Moisture
    conf_dict['modify_moisture_profile'] = True
    conf_dict['dry_run'] = False
    conf_dict['low_level_moisture_height'] = 2000.0    #Low level moisture modification will take effect below this level.2000+(param8*(3000-2000))
    conf_dict['low_level_moisture_mult_factor'] = -15+(param5*(30--15))
    conf_dict['mid_level_moisture_height'] = 2000.0    #Mid level moisture modification will take effect above this level. 2000+(param8*(3000-2000))
    conf_dict['mid_level_moisture_mult_factor'] = 0.0  #Moisture modification factor for mid levels (tipical range [-10.0 , 10.0 ] ) -10+(param7*(10--10))

    # Sensitivity function
    conf_dict['sf_type'] = 'mean_precipitation'
    conf_dict['sf_percentile'] = 99.0

    # Model integration
    conf_dict['run_model'] = True
    conf_dict['plot_exp']  = False
    conf_dict['model_xdomain'] = 2000.0 * 42
    conf_dict['model_nx'] = 80
    conf_dict['model_dx'] = 2000.0
    conf_dict['model_dt'] = 12
    conf_dict['model_nz'] = 41
    conf_dict['model_dy'] = conf_dict['model_dx']
    conf_dict['model_ny'] = conf_dict['model_nx']

    # Set the number of threads for WRF to use all available cores
    # This will be used in the wrf_module.py when calling wrf.exe
    conf_dict['nthreads'] = int(os.environ.get('OMP_NUM_THREADS', '48'))
    
    print(f"Starting simulation {row_number} with {conf_dict['nthreads']} threads")

    wrf.run_wrf(conf_dict)
    
    print(f"Completed simulation {row_number}")

if __name__ == "__main__":
    start_idx = int(sys.argv[1])
    end_idx   = int(sys.argv[2])

    #csv_file_path = 'sampling_batch_1000_8_variables.csv'
    csv_file_path = 'saltelli_d5_1st_le1000.csv'
    with open(csv_file_path, newline='') as csvfile:
        csv_reader = csv.reader(csvfile)
        next(csv_reader)  # skip header
        for row_number, row in enumerate(csv_reader):
            if start_idx <= row_number < end_idx:
                run_wrf_single_conf(row, row_number)