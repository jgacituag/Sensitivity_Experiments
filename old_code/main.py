import numpy as np
import sensiwrf.wrf_module as wrf
from sensiwrf import config as conf


#Load the configuration from conf.py
conf = conf.conf 


conf['expname']   = 'TEST_SENS'

conf['run_num']   = 2
conf['surf_u']    = 10.0 


#Call run wrf with the given configuration. 
wrf.run_wrf( conf ) 





