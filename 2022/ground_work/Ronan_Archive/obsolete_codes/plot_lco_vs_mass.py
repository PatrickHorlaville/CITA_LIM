import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

delta_mf=1.0
alpha=1.37
beta=-1.74
z   = 3.1

#LOAD SFR TABLE
dat_zp1, dat_logm, dat_logsfr, _ = np.loadtxt("sfr_release.dat", unpack=True) # Columns are: z+1, logmass, logsfr, logstellarmass
# Intermediate processing of tabulated data                                                                         
dat_logzp1 = np.log10(dat_zp1)
dat_sfr    = 10.**dat_logsfr

# Reshape arrays                                                                                                    
dat_logzp1  = np.unique(dat_logzp1)  # log(z+1), 1D                                                                 
dat_logm    = np.unique(dat_logm)  # log(Mhalo), 1D                                                                 
dat_sfr     = np.reshape(dat_sfr, (dat_logm.size, dat_logzp1.size))

# Get interpolated SFR value(s)                                                                                     
rbv         = scipy.interpolate.RectBivariateSpline(dat_logm, dat_logzp1, dat_sfr, kx=1, ky=1)

Mco = np.logspace(9,13.6,1000)      
sfr = rbv.ev(np.log10(Mco), np.log10(z+1))

####SFR#to#LCO##################################################################                                                                                                  
lir      = sfr * 1e10 / delta_mf
alphainv = 1./alpha
lcop     = lir**alphainv * 10**(-beta * alphainv)
Lco      =  4.9e-5 * lcop
Lco = Lco/Lco.max()
##PLOT#TRENDLINE################################################################
fig = plt.figure(figsize=(6,5))
ax = plt.subplot(111)
plt.plot(Mco,Lco,color='k',linewidth=4)
ax.set_xlabel(r"Mass [M$_\odot$]")
ax.set_ylabel(r"$L_{CO}$")         
ax.set_xlim(1e9,3.3e13)
ax.set_ylim(2e-3,2)
plt.show()


