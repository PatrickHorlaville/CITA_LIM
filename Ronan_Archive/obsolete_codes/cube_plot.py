#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import plotly as ply
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objs as go
from matplotlib import cm

lem = 115.27

#this code loads the maps in a datacube

#Get number of pixels and field of view
#f        = sys.argv[1]
f = '1160Mpc_CO_nmap_256_novel_fs.map'
infile   = open(f,"rb")
npix_x   = np.fromfile(infile,dtype=np.int32,count=1)[0]
npix_y   = np.fromfile(infile,dtype=np.int32,count=1)[0]
fov_radx = np.fromfile(infile,dtype=np.float32,count=1)
fov_rady = np.fromfile(infile,dtype=np.float32,count=1)

fov_x    = float(fov_radx/2/np.pi*360)
fov_y    = float(fov_rady/2/np.pi*360)
print "\nnpix_x, npix_y          ", npix_x, npix_y
print "fov_x, fov_y:           ", fov_x, fov_y


#find number of maps in cube
infile.seek(0,2) # 2 corresponds to end of file
infile_size = infile.tell()
nmaps       = (infile_size-16)/(4*npix_y*npix_x)
print "\nnumber of maps in cube: ", nmaps
infile.seek(16);


#load in full data cube
data_cube = np.array([np.reshape(np.fromfile(infile,dtype=np.float32,count=npix_x*npix_y),
                                       (npix_x,npix_y)) for i in range(nmaps)])

nnu = data_cube.shape[0]
print "shape of data cube: ", data_cube.shape[0], data_cube.shape[1], data_cube.shape[2]

def plotslice(cube, ax):
    return np.sum(cube, axis=ax)

###SKY#PROJECTION################################################################
#arrrr=plotslice(data_cube, 0)
#
#print arrrr.shape
#
#nx=256
#ny=256
#x, y = np.linspace(0,9.300000019073, nx), np.linspace(0,9.300000019073, ny)
#
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#
#X, Y = np.meshgrid(x, y)  # `plot_surface` expects `x` and `y` data to be 2D
#ax.plot_surface(X, Y, arrrr, cmap=cm.coolwarm)
#
#ax.set_xlabel('degrees')
#ax.set_ylabel('degrees')
#ax.set_zlabel('CO(1-0) Line Intensity (K)')
#
#plt.show()
##COLLAPSE#Y####################################################################
#arrrr=plotslice(data_cube, 1)
#
arrrr=data_cube[:,:,250]

nx=256
nz=1000

x, z = np.linspace(0,9.300000019073, nx),np.linspace(26.0, 33.8125, nz)

rs=(lem/z)-1

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

X, Z = np.meshgrid(x, rs)  # `plot_surface` expects `x` and `y` data to be 2D
fgr=ax.plot_surface(X+1E-7, Z, np.log10(arrrr+1E-7), cmap=cm.coolwarm)
fig.colorbar(fgr, aspect=10)
ax.set_xlabel('degrees')
ax.set_ylabel('Z')
ax.set_zlabel('CO(1-0) Line Intensity (K)')

plt.show()
##FOR#1D#PLOTTING################################################################
##arrrr=data_cube[0,11,250]
#plt.xlabel('Z')
#plt.ylabel('CO(1-0) Line Intensity (K)')
#arrrr=data_cube[:,11,250]
#novel = plt.semilogy(z,arrrr+1E-8, label='No Velocity Correction')
#plt.legend()
#plt.show()

####COLLAPSE#X#################################################################
#arrrr=plotslice(data_cube, 2)
#
#ny=256
#nz=1000
#
#y, z = np.linspace(0,9.300000019073, ny), np.linspace(26.0, 33.8125, nz)
#
#rs=(lem/z)-1
#
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#
#Y, Z = np.meshgrid(y, rs)
#fgr=ax.plot_surface(Y, Z, arrrr, cmap=cm.coolwarm)
#fig.colorbar(fgr, aspect=10)
#
#ax.set_xlabel('degrees')
#ax.set_ylabel('Z')
#ax.set_zlabel('CO(1-0) Line Intensity (K)')
#
#plt.show()
################################################################################

#plt.imshow(plotslice(data_cube, 2))

#for i in xrange(40):
#    print i
#    for k in xrange(40):
#        print k 
#        for m in xrange(40):
#            ax.scatter(i,k,m, alpha=(data_cube[i][k][m]*1331.), color='black')           
        
        
        
    
