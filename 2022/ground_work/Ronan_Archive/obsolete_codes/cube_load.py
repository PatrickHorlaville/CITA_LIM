#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

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
