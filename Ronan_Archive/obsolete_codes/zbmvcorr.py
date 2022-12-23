import numpy as np
#import sympy as sp
import matplotlib.pyplot as plt
import scipy.interpolate as scp

Omega_M = 0.25
h       = 0.7
rhomean = 2.775E11*Omega_M*(h**2)
pixbm   = 0.5             #in arcminutes
amtdeg  = 1./60           #convert arcmin to degrees
dgtrad  = np.pi/180.      #convert degrees to radians
nzbins  = 10000.          #number of redshift bins 
c       = 2.998E5
del_mf  = 1.0
alpha   = 1.37
beta    = -1.74
nuco    = 115.3
amtsr   = 8.461595e-8     #convert arcmin**2 to steradians
f3      = 0.000391003842023 
convfa1 = 6.215e-9
convfa2 = 2.63083e-6
print "saved constants"
##LOAD##DATA####################################################################
data = open('1160Mpc_n4096_nb25_nt17_merge.pksc.13579', 'rb')
print "got data"
nhalo = np.fromfile(data, dtype='int32', count=1)
rthmax = np.fromfile(data, dtype='float32', count=1)
redshift = np.fromfile(data, dtype='float32', count=1)
nhalobytes=nhalo*10
hdata = np.fromfile(data, dtype='float32', count=nhalobytes)
hdata=np.reshape(hdata, (nhalo, 10))
X   = hdata[:,0]
Y   = hdata[:,1]
Z   = hdata[:,2]
print "done xyz"
R   = np.sqrt(X**2+Y**2+Z**2)
VX  = hdata[:,3]
VY  = hdata[:,4]
VZ  = hdata[:,5]
Rth = hdata[:,6]
M   = 4.0/3*np.pi*(Rth**3)*rhomean
print "done loading"

data = []
###R##TO##Z##CONVERSION#########################################################
#D = 6500.
#z1 = sp.Symbol('z1')
#print sp.solve(((c/(100.*h*D)*((z1+1)**2-1)/((z1+1)**2+1))-(h*100.*D))-1, z1)
za=np.linspace(0,4,1000)
def hubble(za):
    return h*100*np.sqrt(Omega_M*(1+za)**3+1-Omega_M)
def drdz(za):
    return 3e5 / hubble(za)
ra  = np.cumsum(drdz(za)*(za[1]-za[0]))
ra -= ra[0]
z_to_r = scp.interpolate.interp1d(za,ra)
r_to_z = scp.interpolate.interp1d(ra,za)
RS = r_to_z(R)
###GET#MASSES#OF#H#I############################################################
MminHI = 1.e10 * (30./60)**3 * ((1+RS)/4)**(-3./2)
MmaxHI = 1.e10 * (200./60)**3 * ((1+RS)/4)**(-3./2)

MmaxHI=MmaxHI[M>MminHI]
M=M[M>MminHI]
HImass = f3*M/(1+M/MmaxHI)

L_HI = convfa1 * HImass
##SET##MORE##USEFUL##VALUES#####################################################
pntanx  = 0.0
pntany  = 0.0
pbrad   = pixbm*amtdeg*dgtrad  #pixel beam in radians 
slope   = np.tan(pbrad/2)      #slope of beam relative to beam centre
bmctrx  = np.tan(pntanx*dgtrad)
bmctry  = np.tan(pntany*dgtrad)
mxz=np.max(Z)
mnz=np.min(Z)
mxx=np.max(X)
ptdir = [0,0, 6378]
ptr   = np.sqrt(ptdir[0]**2+ptdir[1]**2+ptdir[2]**2)
xyzg  = np.asarray([X, Y, Z])
gnrm  = np.sqrt(xyzg[0]**2+xyzg[1]**2+xyzg[2]**2)
xyzg  = xyzg/gnrm
lparl = np.dot(ptdir, xyzg)
lperp = np.sqrt(abs(ptr**2-lparl**2))

#ML  = M[(X > (bmctrx-slope)*Z) & (X< (bmctrx+slope)*Z) & (Y>(bmctry-slope)*Z) & (Y< (bmctry+slope)*Z) & (R < mxz) & (R > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
#XL  = X[(X > (bmctrx-slope)*Z) & (X< (bmctrx+slope)*Z) & (Y>(bmctry-slope)*Z) & (Y< (bmctry+slope)*Z) & (R < mxz) & (R > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
#YL  = Y[(X > (bmctrx-slope)*Z) & (X< (bmctrx+slope)*Z) & (Y>(bmctry-slope)*Z) & (Y< (bmctry+slope)*Z) & (R < mxz) & (R > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
#ZL  = Z[(X > (bmctrx-slope)*Z) & (X< (bmctrx+slope)*Z) & (Y>(bmctry-slope)*Z) & (Y< (bmctry+slope)*Z) & (R < mxz) & (R > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
#RL  = R[(X > (bmctrx-slope)*Z) & (X< (bmctrx+slope)*Z) & (Y>(bmctry-slope)*Z) & (Y< (bmctry+slope)*Z) & (R < mxz) & (R > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
#VXL  = VX[(X > (bmctrx-slope)*Z) & (X< (bmctrx+slope)*Z) & (Y>(bmctry-slope)*Z) & (Y< (bmctry+slope)*Z) & (R < mxz) & (R > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
#VYL  = VY[(X > (bmctrx-slope)*Z) & (X< (bmctrx+slope)*Z) & (Y>(bmctry-slope)*Z) & (Y< (bmctry+slope)*Z) & (R < mxz) & (R > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
#VZL  = VZ[(X > (bmctrx-slope)*Z) & (X< (bmctrx+slope)*Z) & (Y>(bmctry-slope)*Z) & (Y< (bmctry+slope)*Z) & (R < mxz) & (R > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
#RSL  = RS[(X > (bmctrx-slope)*Z) & (X< (bmctrx+slope)*Z) & (Y>(bmctry-slope)*Z) & (Y< (bmctry+slope)*Z) & (R < mxz) & (R > (mnz*np.sqrt((mxz**2+mxx**2))/mxz))]
#def vcorr(rs, x, y, z, vx, vy, vz):
#    xyzv = np.asarray([XL, YL, ZL])
#    vvec = (np.asarray([VXL, VYL, VZL])).T
#    nrmf = np.apply_along_axis(np.linalg.norm, 0, xyzv) #make array of normalization factors
#    uvec = (xyzv/nrmf).T 
#    vrad = np.einsum('ij,ij->i', vvec, uvec)
#    znew = ((1+rs)*(1+(vrad/c)))-1
#    return znew
#RSVC = vcorr(RSL, XL, YL, ZL, VXL, VYL, VZL)
###BIN###########################################################################
#bins   = np.linspace(min(RSVC), max(RSVC), nzbins+1)
#mpcb   = np.linspace(min(RL) , max(RL) , nzbins+1)
#zmpc   = (mpcb[:-1] + mpcb[1:]) / 2
#zspect = np.histogram(RSVC, weights=ML, bins=bins)[0]
#binctr = (bins[:-1] + bins[1:]) / 2
####START##CONVERSION##TO##LCO###################################################
#dat_zp1, dat_logm, dat_logsfr, _ = np.loadtxt("sfr_release.dat", unpack=True) # Columns are: z+1, logmass, logsfr, logstellarmass
## Intermediate processing of tabulated data                                                                         
#dat_logzp1 = np.log10(dat_zp1)
#dat_sfr    = 10.**dat_logsfr
#
## Reshape arrays                                                                                                    
#dat_logzp1  = np.unique(dat_logzp1)  # log(z+1), 1D                                                                 
#dat_logm    = np.unique(dat_logm)  # log(Mhalo), 1D                                                                 
#dat_sfr     = np.reshape(dat_sfr, (dat_logm.size, dat_logzp1.size))
#
## Get interpolated SFR value(s)                                                                                     
#rbv         = scp.interpolate.RectBivariateSpline(dat_logm, dat_logzp1, dat_sfr, kx=1, ky=1)
#Mco = np.logspace(9,13.6,1000)      
####MAKE##A##FUNCTION##FOR##IT###################################################
#def LCO(Mhal, z, chi):
#    sfr      = rbv.ev(np.log10(Mhal+1E2), np.log10(z+1))
#    zeroes   = sfr < 2E-4
#    sfr[zeroes]=0.0
#    lir      = sfr * 1e10 / del_mf
#    alphainv = 1./alpha
#    lcop     = lir**alphainv * 10**(-beta * alphainv)
#    Lco      =  4.9e-5 * lcop
#    numin    = nuco/(3.4+1.)
#    numax    = nuco/(2.4+1.)
#    dnu      = (numax-numin)/nzbins
#    pxsz     = pixbm**2*amtsr
#    nuobs    = nuco/(z+1.)
#    I_line   = (Lco)/(4*np.pi*chi**2*(1+z)**2*dnu*pxsz)
#    T_line   = (1./2)*convfa2*I_line/nuobs**2
#    return T_line*1E6
#
#        
#T_line = LCO(zspect, binctr, zmpc)
#plt.semilogy(binctr, T_line+1E-3)
#plt.xlabel('$Z$')
#plt.ylabel('$T_{co} [\mu K]$')