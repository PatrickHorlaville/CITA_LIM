#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 14:42:17 2017

@author: rkerr31
"""
import numpy as np
import scipy as scp
###CONSTANTS####################################################################
Omega_M = 0.25
h       = 0.7
rsmn    = 2.3935904705839488
rsmx    = 3.442482686918404
nzbn    = 10000.
nuco    = 115.3
del_mf  = 1.0
alpha   = 1.37
beta    = -1.74
amtsr   = 8.461595e-8     #convert arcmin**2 to steradians
pixbm   = 3.              #in arcminutes
rhomean = 2.775E11*Omega_M*(h**2)
#################################################################################

za=np.linspace(0,4,1000)
def hubble(za):
    return h*100*np.sqrt(Omega_M*(1+za)**3+1-Omega_M)
def drdz(za):
    return 3e5 / hubble(za)
ra  = np.cumsum(drdz(za)*(za[1]-za[0]))
ra -= ra[0]
z_to_r = scp.interpolate.interp1d(za,ra)
r_to_z = scp.interpolate.interp1d(ra,za)

halo = np.array([ -2.27656738e+02,  -1.11031670e+02,   6.82393750e+03,
        -6.19773560e+01,   2.53635387e+01,   5.94999695e+01,
         4.08865547e+00,  -2.26925354e+02,  -1.11384079e+02,
         6.82308496e+03], dtype='float32')
    
x   = halo[0]
y   = halo[1]
z   = halo[2]
r   = np.sqrt(x**2+y**2+z**2)
vx  = halo[3]
vy  = halo[4]
vz  = halo[5]
Rth = halo[6]
rs = r_to_z(r)
M   = 4.0/3*np.pi*(Rth**3)*rhomean

bins = np.linspace(rsmn, rsmx, nzbn+1)
bctr = (bins[:-1] + bins[1:]) / 2
mpcb = z_to_r(bins)
zmpc = (mpcb[:-1] + mpcb[1:]) / 2
###FUNCTION##TO##GRAB##SFR######################################################       
dat_zp1, dat_logm, dat_logsfr, _ = np.loadtxt("sfr_release.dat", unpack=True) # Columns are: z+1, logmass, logsfr, logstellarmass
# Intermediate processing of tabulated data                                                                         
dat_logzp1 = np.log10(dat_zp1)
dat_sfr    = 10.**dat_logsfr

# Reshape arrays                                                                                                    
dat_logzp1  = np.unique(dat_logzp1)  # log(z+1), 1D                                                                 
dat_logm    = np.unique(dat_logm)  # log(Mhalo), 1D                                                                 
dat_sfr     = np.reshape(dat_sfr, (dat_logm.size, dat_logzp1.size))

# Get interpolated SFR value(s)                                                                                     
rbv         = scp.interpolate.RectBivariateSpline(dat_logm, dat_logzp1, dat_sfr, kx=1, ky=1)
Mco = np.logspace(9,13.6,1000)   
###FUNCTION##FOR##COMPLETE##CONVERSION##########################################
def LCO(Mhal, z, chi):
    sfr      = rbv.ev(np.log10(Mhal+1E2), np.log10(z+1))
    zeroes   = sfr < 2E-4
    sfr[zeroes]=0.0
    lir      = sfr * 1e10 / del_mf
    alphainv = 1./alpha
    lcop     = lir**alphainv * 10**(-beta * alphainv)
    Lco      =  4.9e-5 * lcop
    numin    = nuco/(rsmx+1.)
    numax    = nuco/(rsmn+1.)
    dnu      = (numax-numin)/nzbn
    pxsz     = pixbm**2*amtsr
    convfac  = 2.63083e-6
    nuobs    = nuco/(z+1.)
    I_line   = (Lco)/(4*np.pi*chi**2*(1+z)**2*dnu*pxsz)
    T_line   = (1./2)*convfac*I_line/nuobs**2
    return T_line*1E6
################################################################################
T_line = LCO(M, r, rs)

