"""
Utility functions for modifying WRF idealized supercell initial conditions.
Handles the mapping of SALib parameters, modifications to the input_sounding,
and dynamic editing of the namelist.input.
"""

import os
import numpy as np

# Physical constants
CP = 1004.0       # Specific heat of dry air at constant pressure (J/K)
RD = 287.0        # Gas constant for dry air (J/kg/K)
G = 9.81          # Gravity (m/s^2)
P0 = 100000.0     # Reference pressure (Pa)


def map_salib_to_params(row, base_params):
    """
    Maps normalized SALib output row values [0, 1] to physical WRF parameters.
    Returns an updated configuration dictionary.
    """
    param1, param2, param3, param4, param5 = map(float, row)
    conf = base_params.copy()

    # Convert YAML lists to NumPy arrays to avoid math errors in sounding modifications
    conf['theta_layer_limits'] = np.array(conf.get('theta_layer_limits', [0.0, 3000.0, 15000.0]))
    conf['dthetadz_layer'] = np.array(conf.get('dthetadz_layer', [0.004, 0.004, 0.004]))
    # Wind Profile
    conf['total_shear_depth'] = 4000.0 + (param1 * (8000.0 - 4000.0))
    conf['int_total_shear']   = 0.00001 + (param2 * (30.0 - 0.0))
    conf['curved_shear_per']  = param3
    
    # Stability
    conf['stability_factor'] = -0.5 + (param4 * (3.0 - -0.5))
    
    # Moisture
    conf['low_level_moisture_mult_factor'] = -15.0 + (param5 * (30.0 - -15.0))
    
    return conf


def read_input_sounding(file_path):
    """Reads the original WRF input_sounding file format into a dictionary."""
    sounding = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
    lines = [line.strip() for line in lines]
    line_split = lines[0].split()
    
    sounding['surf_pressure'] = float(line_split[0])
    sounding['surf_theta']    = float(line_split[1])
    sounding['surf_qv']       = float(line_split[2])
    sounding['nlevs']         = len(lines) - 1
    
    sounding['height'] = np.zeros(sounding['nlevs'])
    sounding['theta']  = np.zeros(sounding['nlevs'])
    sounding['qv']     = np.zeros(sounding['nlevs'])
    sounding['u']      = np.zeros(sounding['nlevs'])
    sounding['v']      = np.zeros(sounding['nlevs'])
    sounding['p']      = np.zeros(sounding['nlevs'])
    
    tmp_press = np.zeros(sounding['nlevs'])
    tmp_press[0] = (sounding['surf_pressure'] * 100.0) ** (RD / CP)

    for i in range(sounding['nlevs']):
        line_split = lines[i + 1].split()
        sounding['height'][i] = float(line_split[0])
        sounding['theta'][i]  = float(line_split[1])
        sounding['qv'][i]     = float(line_split[2])
        sounding['u'][i]      = float(line_split[3])
        sounding['v'][i]      = float(line_split[4])

        if i >= 1:
            theta_mean = 0.5 * (sounding['theta'][i] + sounding['theta'][i-1])
            dz = sounding['height'][i] - sounding['height'][i-1]
            tmp_press[i] = tmp_press[i-1] - G * (P0 ** (RD / CP)) * dz / (CP * theta_mean)
            
    sounding['p'] = tmp_press ** (CP / RD) / 100.0
    sounding['t'] = sounding['theta'] * tmp_press / (P0 ** (RD / CP)) 

    return sounding


def write_input_sounding(file_path, sounding):
    """Writes the modified dictionary back to the WRF input_sounding format."""
    with open(file_path, 'w') as f:
        f.write(f"{sounding['surf_pressure']} {sounding['surf_theta']} {sounding['surf_qv']}\n")
        for i in range(sounding['nlevs']):
            f.write(f"{sounding['height'][i]} {sounding['theta'][i]} {sounding['qv'][i]} {sounding['u'][i]} {sounding['v'][i]}\n")


def modify_input_sounding(file_path, conf):
    """Applies wind, stability, and moisture perturbations to the sounding."""
    sounding = read_input_sounding(file_path)

    # --- 1. Modify Wind Profile ---
    if conf.get('modify_wind_profile', False):
        sounding['u'][:] = 0.0
        sounding['v'][:] = 0.0

        shear_type = conf.get('shear_type', 'Linear')

        if shear_type == 'Linear':
            max_u = conf['shear_depth_u'] * conf['shear_strength_u']
            sounding['u'] = sounding['height'] * conf['shear_strength_u']
            sounding['u'][sounding['height'] > conf['shear_depth_u']] = max_u

            max_v = conf['shear_depth_v'] * conf['shear_strength_v']
            sounding['v'] = sounding['height'] * conf['shear_strength_v']
            sounding['v'][sounding['height'] > conf['shear_depth_v']] = max_v

        elif shear_type == 'Curved':
            curved_shear_theta = conf['curved_shear_per'] * np.pi

            if conf['curved_shear_per'] > 0.0:
                curved_int_shear = conf['int_total_shear'] * conf['curved_shear_per']
                curved_dep_shear = conf['total_shear_depth'] * conf['curved_shear_per']
                curved_amp_shear = curved_int_shear / curved_shear_theta

                mask = sounding['height'] <= curved_dep_shear
                frac = (sounding['height'][mask] / curved_dep_shear) * curved_shear_theta
                sounding['v'][mask] = -1.0 * curved_amp_shear * np.sin(frac)
                sounding['u'][mask] = -1.0 * curved_amp_shear * np.cos(frac)
            else:
                curved_amp_shear = 0.0
                curved_dep_shear = 0.0

            u_ini = -1.0 * curved_amp_shear * np.cos(curved_shear_theta)
            v_ini = -1.0 * curved_amp_shear * np.sin(curved_shear_theta)

            linear_int_shear = conf['int_total_shear'] * (1.0 - conf['curved_shear_per'])
            u_fin = u_ini + linear_int_shear * np.sin(curved_shear_theta)
            v_fin = v_ini - linear_int_shear * np.cos(curved_shear_theta)

            mask = (sounding['height'] > curved_dep_shear) & (sounding['height'] <= conf['total_shear_depth'])
            num = sounding['height'][mask] - curved_dep_shear
            den = (conf['total_shear_depth'] - curved_dep_shear) if (conf['total_shear_depth'] > curved_dep_shear) else 1.0
            w = num / den
            sounding['u'][mask] = u_ini + (u_fin - u_ini) * w
            sounding['v'][mask] = v_ini + (v_fin - v_ini) * w

            mask = sounding['height'] >= conf['total_shear_depth']
            sounding['u'][mask] = u_fin
            sounding['v'][mask] = v_fin

        # Low-level jet and surface wind modifications
        llj_gauss = np.exp(-0.5 * (conf['llj_h'] - sounding['height'])**2 / (conf['llj_width']**2))
        sounding['u'] -= conf['llj_amp'] * llj_gauss * np.sin(np.deg2rad(conf['llj_dir']))
        sounding['v'] -= conf['llj_amp'] * llj_gauss * np.cos(np.deg2rad(conf['llj_dir']))
        
        if 'surf_u' in conf: sounding['u'] += conf['surf_u']
        if 'surf_v' in conf: sounding['v'] += conf['surf_v']

        if conf.get('remove_mean_wind', False):
            mask = sounding['height'] < 6000.0
            sounding['u'] -= np.mean(sounding['u'][mask])
            sounding['v'] -= np.mean(sounding['v'][mask])

    # --- 2. Modify Stability (Theta) ---
    if conf.get('modify_stability', False):
        H = conf['stability_factor_height']
        factor = np.sin(sounding['height'] * np.pi / (2 * H))
        factor[factor < 0.0] = 0.0
        factor[sounding['height'] > 2 * H] = 0.0
        factor = 1.0 + 0.01 * conf['stability_factor'] * factor

        sounding['theta'] = factor * sounding['theta']
        pi = ((sounding['p'] * 100.0) / P0) ** (RD / CP)
        sounding['t'] = sounding['theta'] * pi

    # --- 3. Modify Moisture Profile ---
    if conf.get('modify_moisture_profile', False):
        if conf.get('dry_run', False):
            sounding['qv'][:] = 0.0
        else:
            tmp_z = (-sounding['height'] + conf['low_level_moisture_height']) / 250.0
            tmp_factor = 1.0 + 0.01 * conf['low_level_moisture_mult_factor'] / (1.0 + np.exp(-tmp_z))
            sounding['qv'] *= tmp_factor

            tmp_z = (sounding['height'] - conf['mid_level_moisture_height']) / 250.0
            tmp_factor = 1.0 + 0.01 * conf['mid_level_moisture_mult_factor'] / (1.0 + np.exp(-tmp_z))
            sounding['qv'] *= tmp_factor

    # Saturation Correction
    pi = ((sounding['p'] * 100.0) / P0) ** (RD / CP)
    sounding['t'] = sounding['theta'] * pi
    epsilon = 0.622
    es = 6.112 * np.exp(17.67 * (sounding['t'] - 273.16) / (sounding['t'] - 273.16 + 243.5))
    sounding['qvs'] = es * epsilon / (sounding['p'] - (1 - epsilon) * es) * 1000.0
    mask_sat = sounding['qv'] > sounding['qvs']
    sounding['qv'][mask_sat] = sounding['qvs'][mask_sat]

    write_input_sounding(file_path, sounding)


def edit_namelist(file_path, conf):
    """Dynamically edits the parameters of the namelist.input file."""
    temp_file = file_path + '.tmp'
    with open(temp_file, 'w') as new_file:
        with open(file_path, 'r') as old_file:
            for line in old_file:
                if '@@DT@@' in line:
                    new_file.write(line.replace('@@DT@@', str(conf['model_dt'])))
                elif '@@DT_FRACT_NUM@@' in line:
                    new_file.write(line.replace('@@DT_FRACT_NUM@@', str(conf['model_dt_fract_num'])))
                elif '@@DT_FRACT_DEN@@' in line:
                    new_file.write(line.replace('@@DT_FRACT_DEN@@', str(conf['model_dt_fract_den'])))
                elif '@@DX@@' in line:
                    new_file.write(line.replace('@@DX@@', str(conf['model_dx'])))
                elif '@@DY@@' in line:
                    new_file.write(line.replace('@@DY@@', str(conf['model_dy'])))
                elif '@@NX@@' in line:
                    new_file.write(line.replace('@@NX@@', str(conf['model_nx'])))
                elif '@@NY@@' in line:
                    new_file.write(line.replace('@@NY@@', str(conf['model_ny'])))
                elif '@@NZ@@' in line:
                    new_file.write(line.replace('@@NZ@@', str(conf['model_nz'])))
                else:
                    new_file.write(line)

    os.rename(temp_file, file_path)