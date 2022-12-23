#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 11:32:30 2017

@author: rkerr31
"""
import numpy as np
#import sympy as sp
import matplotlib.pyplot as plt
import scipy.interpolate as scp

Omega_M = 0.25
h       = 0.7
rhomean = 2.775E11*Omega_M*(h**2)
pixbm   = 1./2              #in arcminutes
amtdeg  = 1./60           #convert arcmin to degrees
dgtrad  = np.pi/180.      #convert degrees to radians
nzbins  = 10000.          #number of redshift bins 
c       = 2.998E5
del_mf  = 1.0
alpha   = 1.37
beta    = -1.74
nuco    = 115.3
amtsr   = 8.461595e-8     #convert arcmin**2 to steradians
G       = 6.67E-11
mpctm   = 3.086E22
convfac  = 2.63083e-6
##LOAD##DATA####################################################################
data = open('1160Mpc_n4096_nb25_nt17_merge.pksc.13579', 'rb')
nhalo = np.fromfile(data, dtype='int32', count=1)
rthmax = np.fromfile(data, dtype='float32', count=1)
redshift = np.fromfile(data, dtype='float32', count=1)
nhalobytes=nhalo*10
hdata = np.fromfile(data, dtype='float32', count=nhalobytes)
hdata=np.reshape(hdata, (nhalo, 10))
print "start list generation"
X   = hdata[:,0]
print '1'
Y   = hdata[:,1]
print '2'
Z   = hdata[:,2]
print '3'
R   = np.sqrt(X**2+Y**2+Z**2)
print '4'
VX  = hdata[:,3]
print '5'
VY  = hdata[:,4]
print '6'
VZ  = hdata[:,5]
print '7'
Rth = hdata[:,6]
print '8'
M   = 4.0/3*np.pi*(Rth**3)*rhomean
print '9'
S   = (np.sqrt(4.0/9*np.pi*G*rhomean*6.7702543e-38)*Rth*mpctm)/1000.
print 'done'
pntanx  = 0.0
pntany  = 0.0
pbrad   = pixbm*amtdeg*dgtrad  #pixel beam in radians 
slope   = np.tan(pbrad/2)      #slope of beam relative to beam centre
bmctrx  = np.tan(pntanx*dgtrad)
bmctry  = np.tan(pntany*dgtrad)
mxz = np.max(Z)
mnz = np.min(Z)
mxx = np.max(X)
print 'rs start' 
za=np.linspace(0,4,1000)
def hubble(za):
    return h*100*np.sqrt(Omega_M*(1+za)**3+1-Omega_M)
def drdz(za):
    return 3e5 / hubble(za)
ra  = np.cumsum(drdz(za)*(za[1]-za[0]))
ra -= ra[0]
z_to_r = scp.interpolate.interp1d(za,ra)
r_to_z = scp.interpolate.interp1d(ra,za)

rs = r_to_z(R)
print "rs done"
def vcorr(rs, x, y, z, vx, vy, vz):
        xyzv = np.asarray([x, y, z])
        vvec = (np.asarray([vx, vy, vz])).T
        nrmf = np.sqrt(x**2+y**2+z**2) #make array of normalization factors
        uvec = (xyzv/nrmf).T 
        vrad = np.einsum('ij,ij->i', vvec, uvec)
        znew = ((1+rs)*(1+(vrad/c)))-1
        return znew
print "vcorr start"
rsvc = vcorr(rs, X, Y, Z, VX, VY, VZ)
print "corr distance"
rvc  = z_to_r(rsvc)
print "done"
print "distance filter"
M   = M[(rvc < (mxz)) & (rvc > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
X   = X[(rvc < (mxz)) & (rvc > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
Y   = Y[(rvc < (mxz)) & (rvc > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
Z   = Z[(rvc < (mxz)) & (rvc > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
R   = R[(rvc < (mxz)) & (rvc > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
VX  = VX[(rvc < (mxz)) & (rvc > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
VY  = VY[(rvc < (mxz)) & (rvc > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
VZ  = VZ[(rvc < (mxz)) & (rvc > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
S   = S[(rvc < (mxz)) & (rvc > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
RSVC   = rsvc[(rvc < (mxz)) & (rvc > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
print "done"
bins   = np.linspace(np.min(RSVC), np.max(RSVC), nzbins+1)
binctr = (bins[:-1] + bins[1:]) / 2
def pltdata(pntanx=pntanx, pntany=pntany):
    print 'start'
    ML     = M[(X > (bmctrx-slope)*Z) & (X< (bmctrx+slope)*Z) & (Y>(bmctry-slope)*Z) & (Y< (bmctry+slope)*Z)]
    print '1'
    RL     = R[(X > (bmctrx-slope)*Z) & (X< (bmctrx+slope)*Z) & (Y>(bmctry-slope)*Z) & (Y< (bmctry+slope)*Z)]
    print '2'
    RSVCL  = RSVC[(X > (bmctrx-slope)*Z) & (X< (bmctrx+slope)*Z) & (Y>(bmctry-slope)*Z) & (Y< (bmctry+slope)*Z)]
    print '3'
    SL      = S[(X > (bmctrx-slope)*Z) & (X< (bmctrx+slope)*Z) & (Y>(bmctry-slope)*Z) & (Y< (bmctry+slope)*Z)]
    ##BIN###########################################################################
#    bins   = np.linspace(min(RSVC), max(RSVC), nzbins+1)
#    mpcb   = np.linspace(min(RL) , max(RL) , nzbins+1)
#    zmpc   = (mpcb[:-1] + mpcb[1:]) / 2
#    zspect = np.histogram(RSVC, weights=ML, bins=bins)[0]
#    binctr = (bins[:-1] + bins[1:]) / 2
    ###START##CONVERSION##TO##LCO###################################################
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
    ###MAKE##A##FUNCTION##FOR##IT###################################################
    def LCO(Mhal, z, chi):
        sfr      = rbv.ev(np.log10(Mhal+1E2), np.log10(z+1))
        zeroes   = sfr < 2E-4
        sfr[zeroes]=0.0
        lir      = sfr * 1e10 / del_mf
        alphainv = 1./alpha
        lcop     = lir**alphainv * 10**(-beta * alphainv)
        Lco      =  4.9e-5 * lcop
        return Lco
            
    return RSVCL, LCO(ML, RSVCL, RL), SL

RS,LCO,SF = pltdata()

def gaussian(x,s,LCO,mu):
    return LCO/(s/c*np.sqrt(2*np.pi))*np.e**(-(x-mu)**2/(2*(s/c)**2))

spect = np.zeros(binctr.shape)
print 'start'
for i in xrange(len(RS)):
    spect = spect + gaussian(binctr, SF[i], LCO[i], RS[i])
print 'end'

def Tline(Lco,z):
    chi      = z_to_r(z)
    numin    = nuco/(3.4+1.)
    numax    = nuco/(2.4+1.)
    dnu      = (numax-numin)/nzbins
    pxsz     = pixbm**2*amtsr
    nuobs    = nuco/(z+1.)
    I_line   = (Lco)/(4*np.pi*chi**2*(1+z)**2*dnu*pxsz)
    T_line   = (1./2)*convfac*I_line/nuobs**2
    return T_line*1E6

SpctT = Tline(LCO, RS, )

plt.plot(binctr, spect)

