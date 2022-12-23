import params
from pykDict import *
from cosmolopy import fidcosmo

def getparams(filename):
    
    dict=pykDict()
    dict.read_from_file(filename)

    params.verbose         = dict['verbose']

    params.format          = dict['format']
    params.folder          = dict['folder']
    params.inputfile       = dict['inputfile']
    params.mapfile         = dict['mapfile']

    params.H0              = dict['H0']
    params.omegab          = dict['omegab']
    params.omegac          = dict['omegac']
    params.omegal          = dict['omegal']
    params.scalar_index    = dict['scalar_index']
    params.sigma8          = dict['sigma8']
    params.min_redshift    = dict['min_redshift']
    params.max_redshift    = dict['max_redshift']

    params.box_size        = dict['box_size']
    params.fov             = dict['fov']
    params.flat            = dict['flat']
    params.numdens         = dict['numdens']

    params.hod             = dict['hod']
    params.LM             = dict['LM']

    params.shang_L0        = dict['shang_L0']
    params.shang_Td        = dict['shang_Td']
    params.shang_beta      = dict['shang_beta']
    params.shang_eta       = dict['shang_eta']
    params.shang_I0        = dict['shang_I0']
    params.shang_Mmin      = dict['shang_Mmin']
    params.shang_Msmin     = dict['shang_Msmin']
    params.shang_Mpeak     = dict['shang_Mpeak']
    params.shang_sigmaM    = dict['shang_sigmaM']
    params.shang_nc2ns     = dict['shang_nc2ns']

    params.nu_obs          = dict['nu_obs'] * 1e9 # GHz to Hz
    params.nside           = dict['nside']

    params.min_mass        = dict['min_mass']
    params.M1              = dict['M1']
    params.M0              = dict['M0']
    params.sigma_logM      = dict['sigma_logM']
    params.logM_min        = dict['logM_min']
    params.alpha           = dict['alpha']
    params.mass_type       = dict['mass_type']

    params.omegam          = params.omegab+params.omegac
    params.h               = float(params.H0)/100
    
    # for cosmolopy 
    params.cosmo                   = fidcosmo
    params.cosmo['h']              = params.h
    params.cosmo['n']              = params.scalar_index
    params.cosmo['omega_M_0']      = params.omegam
    params.cosmo['omega_b_0']      = params.omegab
    params.cosmo['omega_lambda_0'] = params.omegal
    params.cosmo['sigma_8']        = params.sigma8
