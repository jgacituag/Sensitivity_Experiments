#This script run the wrf model. 
import pickle
import numpy as np
import os 
import wrf_module as wrf

#############################################################################################
# CONFIGURATION PARAMETERS (in a future version this can be an input to a function)
#############################################################################################

conf = dict()
conf['expname']   = 'TEST'
conf['run_num']   = 2
conf['modelpath'] = '/home/jruiz/share/EXPERIMENTS/EXPERIMENTS_UNWEATHER/comp_model/em_squall_2d'
conf['datapath']  = '/home/jruiz/share/EXPERIMENTS/EXPERIMENTS_UNWEATHER/DATA/'
conf['nthreads']  = 3 #Number of threads that will be used for running WRF

#Parameters controling the shape of the wind profile.
conf['modify_wind_profile'] = True  #Are we going to modify the original wind profile?
conf['remove_mean_wind']    = False #True for convection experiments, false for mountain wave experiments.
conf['surf_u']           = 10.0     #U-component wind at the surface 
conf['surf_v']           = 0.0      #V-component wind at the surface.
conf['shear_depth_u']    = 15000.0  #Heigth of no shear (for the u component)
conf['shear_strength_u'] = 0.0e-3   #Shear in m/s^2 (for the u component)
conf['shear_depth_v']    = 0.0      #Heigth of no shear (for the v component)
conf['shear_strength_v'] = 0.0e-3   #Shear in m/s^2 (for the v component)
conf['llju_amp']         = 0.0      #Low level jet local maximum (for the u component)
conf['llju_h']           = 1500.0   #Low level jet heigth (for the u component)
conf['llju_width']       = 500.0    #Low level jet width  (for the u component)
conf['lljv_amp']         = 0.0      #Low level jet local maximum (for the v component)
conf['lljv_h']           = 0.0      #Low level jet heigth (for the v component)
conf['lljv_width']       = 500.0    #Low level jet width  (for the v component)

#Parameters controling the stability (temperature profile)
conf['modify_stability'] = True          #Are we going to modify the original stability?
conf['stability_factor'] = -1.0           #Factor controling the stability change (tipical range [-1.5 , 1.5] )
conf['stability_factor_height'] = 10000.0 #Height of maximum warming / cooling. (tipical range [2500 , 10000 ] )

#Parameters controling the shape of the potential temperature profile.
conf['modify_theta_profile'] = False  #Are we going to modify the original temperature?
conf['surf_theta']           = 288.0 #Potential temperature at the lowest model level.
#Temperature profile will be computed as costant d(theta)/dz in user defined layers. 
conf['theta_layer_limits']  = np.array([0.0,3000.0,15000.0])  #Layer limits.
conf['dthetadz_layer']      = np.array([4.0,4.0,4.0])/1000.0 #dthetadz corresponding to each layer.

#Parameters controling the shape of the moisture profile.
conf['modify_moisture_profile'] = True        #Are we going to modify the original moisture profile?
conf['dry_run']                 = False       #Assume 0.0 moisture content at all levels?
conf['low_level_moisture_height'] = 2000.0    #Low level moisture modification will take effect below this level.
conf['low_level_moisture_mult_factor'] = 0.0  #Moisture modification for low levels (tipical range [-15.0 , 15.0] )   
conf['mid_level_moisture_height'] = 2000.0    #Mid level moisture modification will take effect above this level.
conf['mid_level_moisture_mult_factor'] = 0.0  #Moisture modification factor for mid levels (tipical range [-10.0 , 10.0 ] )

#Parameters controling the model parameters (parameters controlling the model integration)
conf['model_xdomain'] = 250.0 * 202  #Do not change (total extension of the domain in the x direction, meters)
conf['model_nx'] = 202               #Number of grid points in the X direction (integer) 
conf['model_dx'] = conf['model_xdomain'] / conf['model_nx']  #Model resolution is computed from the number of grid points to keep the model domain unchanged (meters) 
conf['model_dt_fract_num'] = 0       #Numerator for fractional time step (from 1 to 10)
conf['model_dt_fract_den'] = 10      #Denominator for fractional time step (for the moment we kept this equal to 10)
#From the WRF documentation> Example, if you want to use 60.3 sec as your time step,
#                            set model_dt = 60, 
#                            model_dt_fract_num = 3, and 
#                            model_dt_fract_den = 10
conf['model_dt'] = 3                 #Integration time step (seconds, integer)            
conf['model_nz'] = 60                #Number of vertical levels (integer)
conf['model_dy'] = conf['model_dx']  #Model resolution in the y-direction (for the moment we kept this equal to the resolution in the x-direction). 
conf['model_ny'] = 3                 #Number of grid points in the y direction (for the moment we keep this constant). 
#Other possible parameters (difussion, etc)



#Parameters controlling the sensitivity function to be used
conf['sf_type'] = 'mean_precipitation' 
conf['sf_percentile'] = 99.0 

conf['run_model'] = True  #Are we going to run the model?
conf['plot_exp']  = True  #Are we going to do detailed plots of the model solution?


#############################################################################################
# PREPARATORY STEPS
#############################################################################################

#Create output directory
expoutpath = conf['datapath'] + '/' + conf['expname']  
os.makedirs(expoutpath,exist_ok=True)
#Copy configuration file to the output directory
with open( expoutpath + '/conf_' + str( conf['run_num'] ) + '.pkl' , 'wb') as handle:
    pickle.dump( conf , handle, protocol=pickle.HIGHEST_PROTOCOL)

#Create temporary directory 
#We need a temporary directory to run wrf model 
#we create a random one assuming that several experiments will be running at the same time. 

tmp_dir = conf['datapath'] + '/tmp/tmpdir_' + str( np.round( np.random.rand() * 100000 ) )
os.makedirs( tmp_dir , exist_ok = True )

#############################################################################################
# STEP 1 - GENERATE THE STORM "ENVIRONMENT" 
#############################################################################################

if conf['run_model'] :
   #Copy wrf executables and data to the temp directory.
   os.system( 'cp ' + conf['modelpath'] + '/* ' + tmp_dir + '/' )
   input_sounding_file = tmp_dir + '/input_sounding' 
   sounding_ori = wrf.read_input_sounding( input_sounding_file )
   wrf.modify_input_sounding( input_sounding_file , conf )
   sounding_mod = wrf.read_input_sounding( input_sounding_file )

#############################################################################################
# STEP 2 - EDIT THE MODEL NAMELIST (MODEL INTEGRATION PARAMETERS AMONG OTHERS) 
#############################################################################################

#if conf['run_model'] :
#   #Modify the namelist.input
#   wrf.edit_namelist( tmp_dir + '/namelist.input', conf )

#############################################################################################
# STEP 3 - RUN THE MODEL 
#############################################################################################
outdatafile = expoutpath + '/wrfout_' + str( conf['run_num'] ) + '.nc'

#if conf['run_model'] :
#   os.chdir( tmp_dir ) #Change working directory to temp_dir.name
#   #Run the model (to executables ideal.exe prepares the initial conditions for the model
#   #integration, wrf.exe integrates the Navier-Stokes equation from the initial conditions
#   #and for the period of the simulation). 
#   print('Running ideal.exe')
#   os.system('./ideal.exe > out_ideal.log')
#   print('Running wrf.exe')
#   os.system('export OMP_NUM_THREADS=' + str(conf['nthreads'])  + ';./wrf.exe > out_wrf.log')
#   #Move the output (netcdf file containing 4D arrays representing spatio-temporal variation
#   #of different physical variables such as temperatur, pressure, wind components, clouds)
#   os.system('mv ./wrfout_d01_0001-01-01_00:00:00 ' + outdatafile )

#############################################################################################
# STEP 4 - OBTAIN SIMULATION METRICS  
#############################################################################################

#if conf['run_model'] :

#   sensitivity_score =  wrf.sensitivity_function_single_file( outdatafile , sf_type= conf['sf_type'], sf_percentile = conf['sf_percentile'] , inilev = 0 , endlev = 10 ) 

#   print('The sensitivity score is : ', sensitivity_score )

#############################################################################################
# STEP 5 (OPTIONAL) - PLOT THE OUTPUT 
#############################################################################################

#if conf['plot_exp'] :
#   #Import modules for ploting
#   import wrf_plot   as wrfp
#   #Create a directory for the figures.
#   plot_path = expoutpath + '/figs_' + str(conf['run_num']) + '/'
#   os.makedirs(plot_path,exist_ok=True)
#   #Plot storm "environment" 
#   wrf_data = wrf.get_profile_single_file( outdatafile , xp=0 , yp=0 , tp=0 ) #Read the data from netcdf file.
#   wrfp.plot_sounding( wrf_data , plot_path , show=False )                    #Plot the storm "environment" 
#
#   #Plot time evolution of "reflectivity" 
#   wrf_data_2 = wrf.get_data_vslice_single_file( outdatafile , slice_width=1 , slice_index=1 , slice_type='vy', t_start = 0 , t_end = 7 , write_pkl = False )
#   wrfp.plot_ref_vertvel( wrf_data_2 , plot_path , scale_factor = 1.0 , arrow_scale_factor = 1.0 , ybound = [0,15000] ) 






