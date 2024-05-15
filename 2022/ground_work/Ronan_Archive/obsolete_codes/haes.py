#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 11:32:30 2017

@author: rkerr31
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as scp

Omega_M = 0.25
h       = 0.7
rhomean = 2.775E11*Omega_M*(h**2)
pixbm   = 1./2            #in arcminutes
amtdeg  = 1./60           #convert arcmin to degrees
dgtrad  = np.pi/180.      #convert degrees to radians
nzbins  = 100000.          #number of redshift bins 
c       = 2.998E5
del_mf  = 1.0
alpha   = 1.37
beta    = -1.74
nuco    = 4.571E14
amtsr   = 8.461595e-8     #convert arcmin**2 to steradians
G       = 6.67E-11
mpctm   = 3.086E22
convfac = 2.63083e-6
mnz     = 2.4
mxz     = 3.4
##LOAD##DATA####################################################################
data = np.load('fullcathae.npy')
print "start list generation"
X   = data[:,0]
print '1'
Y   = data[:,1]
print '2'
Z   = data[:,2]
print '3'
rs   = np.sqrt(X**2+Y**2+Z**2)
print '7'
L   = data[:,3]/3.826E33
VR  = data[:,4]
RH  = data[:,5]
S   = (np.sqrt(4.0/9*np.pi*G*rhomean*6.7702543e-38)*RH*mpctm)/1000.
print 'done'
pntanx  = -4.049+(0.5/60.)
pntany  = 3.236
pbrad   = pixbm*amtdeg*dgtrad  #pixel beam in radians 
slope   = np.tan(pbrad/2)      #slope of beam relative to beam centre
bmctrx  = np.tan(pntanx*dgtrad)
bmctry  = np.tan(pntany*dgtrad)
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

rs = r_to_z(rs)
#mnrs    = z_to_r(2.4)
#mxrs    = z_to_r(3.4)
print "rs done"
def vcorr(rs, vr):
        znew = ((1+rs)*(1+(vr/c)))-1
        return znew
print "vcorr start"
rsvc = vcorr(rs, VR)
print "done"
print "distance filter"
cond1 = (rsvc < mxz) & (rsvc > mnz)
L    = L[cond1]
X    = X[cond1]
Y    = Y[cond1]
Z    = Z[cond1]
VR   = VR[cond1]
S    = S[cond1]
RSVC = rsvc[cond1]
cond2 = (X > (bmctrx-slope)*Z) & (X< (bmctrx+slope)*Z) & (Y>(bmctry-slope)*Z) & (Y< (bmctry+slope)*Z)
print "done"
bins   = np.linspace(np.min(RSVC), np.max(RSVC), nzbins+1)
binctr = (bins[:-1] + bins[1:]) / 2
def pltdata(pntanx=pntanx, pntany=pntany):
    print 'start'
    LL     = L[cond2]
    print '1'
    RSVCL  = RSVC[cond2]
    print '3'
    SL      = S[cond2]
    return RSVCL, LL, SL

RS,LCO,SF = pltdata()

def gaussian(x,s,mu):
    return np.e**(-(x-mu)**2/(2*(s/c)**2))

def nrmfac(lco, spec):
    return lco/sum(spec)

spect = np.zeros(binctr.shape)
print 'start'
for i in xrange(len(RS)):
    gussn = gaussian(binctr, SF[i], RS[i])
    norm  = nrmfac(LCO[i], gussn)
    spect = spect + (gussn*norm)
print 'end'

chi = z_to_r(binctr)
def IHA(Lha,z,chi):
    numin    = nuco/(mxz+1.)
    numax    = nuco/(mnz+1.)
    dnu      = (numax-numin)/nzbins
    pxsz     = pixbm**2*amtsr
    nuobs    = nuco/(z+1.)
    I_line   = (Lha)/(4*np.pi*dnu*chi**2*(1+z)**2*pxsz)
#    T_line   = (1./2)*convfac*I_line/nuobs**2               
    return I_line*1E6

SpctT = IHA(spect, binctr, chi)

plt.semilogy(binctr, SpctT+1e-3)

#bns = np.linspace(0.1, np.log(1500), 30)
#bns = np.exp(bns)
#nvbnsc = (bns[:-1] + bns[1:]) / 2
#nvhl = plt.hist(SpctT, bins=bns)[0]
'''
f, axarr = plt.subplots(3, sharex=True)
axarr[0].semilogy(binctr, Sf1+1E-3)
axarr[2].set_ylabel('$T_{CO} [\mu K]$')
axarr[1].semilogy(binctr, Sf2+1E-3)
axarr[1].set_ylabel('$T_{CO} [\mu K]$')
axarr[2].semilogy(binctr, Sf3+1E-3)
plt.xlabel('$Z$')
axarr[0].set_ylabel('$T_{CO} [\mu K]$')
'''