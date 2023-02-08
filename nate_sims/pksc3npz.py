import numpy as np
import sys

file_in  = open( sys.argv[1] , 'rb' )
file_out = sys.argv[2]

Non      = np.fromfile(file_in, dtype=np.int32, count=1)[0]
RTHmax   = np.fromfile(file_in, dtype=np.float32, count=1)[0]
redshift = np.fromfile(file_in, dtype=np.float32, count=1)[0]
peakdata = np.fromfile(file_in, dtype=np.float32, count=int(Non*11)).reshape((Non,11))

# print(len('Length of "Redshift" object=', len(redshift)))



# Eulerian position of halos
xpk = peakdata[:,0]
ypk = peakdata[:,1]
zpk = peakdata[:,2]

# Halo velocities
vxpk = peakdata[:,3]
vypk = peakdata[:,4]
vzpk = peakdata[:,5]

# Halos found at filter scalae
RTH = peakdata[:,6]

# Halo comoving distance
xLpk = peakdata[:,7]
yLpk = peakdata[:,8]
zLpk = peakdata[:,9]

# Field density at halo peak
Fpk  = peakdata[:,10]

# Halo mass
h      = .6735
omegam = .3138
rho    = 2.775e11*omegam*h**2
M      = RTH**3 *4/3 * np.pi * rho

# Cosmo parameters
cosmo_header = {'Omega_M': 0.286, 'Omega_B': 0.047, 'Omega_L': 0.714, 'h': 0.7, 'ns': 0.96, 'sigma8': 0.82}

# Formation redshift
zform = 999


np.savez( file_out, Nhalo=Non, M=M, x=xpk, y=ypk, z=zpk, xL=xLpk, yL=yLpk, zL=zLpk, vx=vxpk, vy=vypk, vz=vzpk, cosmo_header = cosmo_header, zform = zform, peak = peakdata )



