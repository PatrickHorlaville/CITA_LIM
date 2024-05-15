import numpy as np
import astropy.cosmology as cosmo

Omega_M = 0.25
h       = 0.7
rhomean = 2.775E11*Omega_M*(h**2)
pixbm   = 4.              #in arcminutes
amtdeg  = 1./60           #convert arcmin to degrees
dgtrad  = np.pi/180.      #convert degrees to radians
nzbins  = 10000.          #number of redshift bins 
c       = 2.998E5
##LOAD##DATA####################################################################
data = open('1160Mpc_n4096_nb25_nt17_merge.pksc.13579', 'rb')
nhalo = np.fromfile(data, dtype='int32', count=1)
rthmax = np.fromfile(data, dtype='float32', count=1)
redshift = np.fromfile(data, dtype='float32', count=1)
nhalobytes=nhalo*10
hdata = np.fromfile(data, dtype='float32', count=nhalobytes)
hdata=np.reshape(hdata, (nhalo, 10))
x   = hdata[:,0]
y   = hdata[:,1]
z   = hdata[:,2]
r   = np.sqrt(x**2+y**2+z**2)
vx  = hdata[:,3]
vy  = hdata[:,4]
vz  = hdata[:,5]
Rth = hdata[:,6]
M   = 4.0/3*np.pi*(Rth**3)*rhomean
################################################################################
pntanx  = 0.5
pntany  = 0.5
pbrad   = pixbm*amtdeg*dgtrad  #pixel beam in radians 
slope   = np.tan(pbrad/2)      #slope of beam relative to beam centre
bmctrx  = np.tan(pntanx*dgtrad)
bmctry  = np.tan(pntany*dgtrad)
ML = M[(x > (bmctrx-slope)*z) & (x< (bmctrx+slope)*z) & (y>(bmctry-slope)*z) & (y< (bmctry+slope)*z) & (r < (max(z))) & (r > (min(z)*np.sqrt((max(z)**2+max(x)**2))/max(z)))]
XL = x[(x > (bmctrx-slope)*z) & (x< (bmctrx+slope)*z) & (y>(bmctry-slope)*z) & (y< (bmctry+slope)*z) & (r < (max(z))) & (r > (min(z)*np.sqrt((max(z)**2+max(x)**2))/max(z)))]
YL = y[(x > (bmctrx-slope)*z) & (x< (bmctrx+slope)*z) & (y>(bmctry-slope)*z) & (y< (bmctry+slope)*z) & (r < (max(z))) & (r > (min(z)*np.sqrt((max(z)**2+max(x)**2))/max(z)))]
ZL = z[(x > (bmctrx-slope)*z) & (x< (bmctrx+slope)*z) & (y>(bmctry-slope)*z) & (y< (bmctry+slope)*z) & (r < (max(z))) & (r > (min(z)*np.sqrt((max(z)**2+max(x)**2))/max(z)))]
RL = r[(x > (bmctrx-slope)*z) & (x< (bmctrx+slope)*z) & (y>(bmctry-slope)*z) & (y< (bmctry+slope)*z) & (r < (max(z))) & (r > (min(z)*np.sqrt((max(z)**2+max(x)**2))/max(z)))]
bins   = np.linspace(min(RL), max(RL), nzbins+1)
zspect = np.histogram(RL, weights=ML, bins=bins)[0]
binctr = (bins[:-1] + bins[1:]) / 2

v  = h*10.*r         
zd = np.sqrt((1+(h*100.*r/c))/(1-(h*100.*r/c)))-1


#for i in xrange(len(x)):
#    if (x[i]>-590.*z[i]) and (x[i]<590.*z[i]) and (y[i]>-590.*z[i]) and (y[i]>-590.*z[i]):
#        print i
#        ML.append(M[i])
        

