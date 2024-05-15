#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 13:25:25 2017

@author: rkerr31
"""
import numpy as np
import math
import scipy.interpolate as scp

G       = 6.67E-11
msol    = 1.989E30
mpctm   = 3.086E22
Omega_M = 0.25
h       = 0.7
rhomean = 2.775E11*Omega_M*(h**2)

satx = np.load('./1160Mpc_n4096_nb25_nt17_merge.pksc.13579_censat.npz_FILES/x_sat.npy')
saty = np.load('./1160Mpc_n4096_nb25_nt17_merge.pksc.13579_censat.npz_FILES/y_sat.npy')
satz = np.load('./1160Mpc_n4096_nb25_nt17_merge.pksc.13579_censat.npz_FILES/z_sat.npy')
sams = np.load('./1160Mpc_n4096_nb25_nt17_merge.pksc.13579_censat.npz_FILES/M_sat.npy')
nmss = np.load('./1160Mpc_n4096_nb25_nt17_merge.pksc.13579_censat.npz_FILES/ns_sat.npy')

cenx = np.load('./1160Mpc_n4096_nb25_nt17_merge.pksc.13579_censat.npz_FILES/x_cen.npy')
ceny = np.load('./1160Mpc_n4096_nb25_nt17_merge.pksc.13579_censat.npz_FILES/y_cen.npy')
cenz = np.load('./1160Mpc_n4096_nb25_nt17_merge.pksc.13579_censat.npz_FILES/z_cen.npy')
cems = np.load('./1160Mpc_n4096_nb25_nt17_merge.pksc.13579_censat.npz_FILES/M_cen.npy')
nmsc = np.load('./1160Mpc_n4096_nb25_nt17_merge.pksc.13579_censat.npz_FILES/ns_cen.npy')

Rvir = ((3*sams)/(4*np.pi*rhomean*200))**(1./3)
rsol = np.sqrt(satx**2 + saty**2 + satz**2)
rsolc = np.sqrt(cenx**2 + ceny**2 + cenz**2)

za=np.linspace(0,4,1000)
def hubble(za):
    return h*100*np.sqrt(Omega_M*(1+za)**3+1-Omega_M)
def drdz(za):
    return 3e5 / hubble(za)
ra  = np.cumsum(drdz(za)*(za[1]-za[0]))
ra -= ra[0]
z_to_r = scp.interpolate.interp1d(za,ra)
r_to_z = scp.interpolate.interp1d(ra,za)

rs = r_to_z(rsol)

RS  = Rvir/(7.85*(sams/(2E12/h))**-0.081*1/(1+rs)**0.71)
crr = Rvir/RS

cencox = []
cencoy = []
cencoz = []
for i in xrange(len(cenx)):
    if int(nmsc[i])>0:
        for k in xrange(int(nmsc[i])):
            cencox.append(cenx[i])
            cencoy.append(ceny[i])
            cencoz.append(cenz[i])

cencox = np.asarray(cencox)
cencoy = np.asarray(cencoy)
cencoz = np.asarray(cencoz)     

dxyz = np.r_['0,2',(satx-cencox), (saty-cencoy), (satz-cencoz)]
rsat = np.sqrt(dxyz[0]**2+dxyz[1]**2+dxyz[2]**2)
xrr = rsat/RS
rhos = 3*sams/(4*np.pi*RS**3)

def consts(cr):
    return sams/(np.log(1+cr)-(cr/(1+cr)))

cs = consts(crr)
            
def Min(r, cst):
    return cst*(np.log(1+r)-(r/(1+r)))

Minr = Min(xrr,cs)

vtan = np.sqrt(G*Minr*msol/(rsat*mpctm))/1000.

def perpvect(xc, yc, zc, xs, ys, zs):
    xyzc = (np.asarray([xc, yc, zc])/np.sqrt(xc**2+yc**2+zc**2)).T
    xyzs = (np.asarray([xs, ys, zs])/np.sqrt(xs**2+ys**2+zs**2)).T
    angles = np.asarray(map(math.acos, (np.einsum('ij,ij->i', xyzc, xyzs))))
    velang = (angles-(np.pi/2))
    rands = np.cos(np.random.rand(len(xc))*2*np.pi)
    rv = np.cos(velang)*rands
    return np.asarray(rv)

uvvec = perpvect(cencox, cencoy, cencoz, dxyz[0], dxyz[1], dxyz[2])
vrads  = vtan*uvvec
vradc  = np.zeros(len(cenx))

sm = sams/(nmss+1)
cm = cems/(nmsc+1)

cendatarr = np.r_['0,2', cenx, ceny, cenz, rsolc, cm, vradc].T
satdatarr = np.r_['0,2', satx, saty, satz, rsol,  sm, vrads].T
                 
arrf = np.r_[cendatarr, satdatarr]
np.save('fullcat.npy', arrf)

