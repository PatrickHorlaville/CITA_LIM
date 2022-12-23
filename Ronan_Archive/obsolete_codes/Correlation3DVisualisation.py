#######################################################
#         CORRELATION 3D SPATIAL VISUALISATION        #
#######################################################

# Essential so that IPython kernel will not crash
import os
os.environ ['ETS_TOOLKIT'] = 'wx'

import sys
import numpy as np
import pandas as pds
from mayavi import mlab

print ("\n========== 3D Correlation Spatial Visualisation ==========\n")

############## Parameters ##############

# Correlation parameters
rMax = 50.
rRes = 0.5
L = 2 * int (np.ceil(rMax/rRes))
csvFilename = "/home/bruno/Bureau/CITA/Python/Correlation plots/3D/Improved/Eulerian/1/dGrid.csv"
densitySelection = 3.

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
############## Grid initialization ##############

print ("Loading CSV...")
sys.stdout.flush ()

xP = []
yP = []
zP = []
for i in range (L):
    for j in range (L):
        for k in range (L):
            xP.append (-rMax + rRes/2. + i * rRes)
            yP.append (-rMax + rRes/2. + j * rRes)
            zP.append (-rMax + rRes/2. + k * rRes)
sP_tmp = data_cube.as_matrix ()

# We create numpy objects
xP = np.array (xP)
yP = np.array (yP)
zP = np.array (zP)

print ("Loaded!")
sys.stdout.flush ()

############## Data selection and normalization ##############

# We want to get the density contrast delta

rhoMean = sP_tmp.sum () / float (L ** 3)
sP = sP_tmp/rhoMean - 1.

print ("\nRho mean = {0}\n".format (rhoMean))
sys.stdout.flush ()

# Data density selection 

print ("\nDensity selection : delta > {0}\n".format (densitySelection))
sys.stdout.flush ()

sPMask = sP >= densitySelection
xP = xP [sPMask]
yP = yP [sPMask]
zP = zP [sPMask]
sP = sP [sPMask]

############## 3D plot ##############

print ("3D Plot...")
sys.stdout.flush ()

# Figure creation
fig = mlab.figure (bgcolor = (0,0,0))
fig.scene.disable_render = True # We disable temporarily the rendering for performance

# Mayavi Points3d plot
pts = mlab.points3d (xP, yP, zP, sP, colormap = 'jet', mode = 'cube')
pts.glyph.scale_mode = 'data_scaling_off'               # Disable default scaling with sP
pts.glyph.glyph.scale_factor = rRes                     # rRes scaling instead
pts.glyph.glyph.clamping = False                        # Not sure of what it is, but disable it for performance
pts.actor.mapper.resolve_coincident_topology = 'off'    # Not sure of what it is, but disable it for performance
pts.actor.property.opacity = 0.5                        # Needed to avoid strange disappearing of points
pts.actor.property.representation = 'wireframe'         # Wireframe representation to avoid overloading

# We add a transparency effect following scalar value of objects to highlight overdensities
lut = pts.module_manager.scalar_lut_manager.lut.table.to_array()
lut[:, -1] = np.linspace(0, 255, num = 256) # Lut alpha value modification
pts.module_manager.scalar_lut_manager.lut.table = lut

# Now we are ready for the actual display
mlab.colorbar ()                    # Colorbar display
fig.scene.disable_render = False    # We can finally enable rendering
mlab.show ()                        # Mayavi loop
