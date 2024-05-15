#!/usr/bin/env python
from mods     import *
from hod      import *
from fluxmap  import *
from utils    import *

Non    = get_catalog_Non()
Nreads = int(max(1,params.size/2))

report('Reading, Shuffling, Culling Catalogue',2)
#To load balance without reading in entire catalogue on each processor
for i in range(Nreads):  
    ceni, Mmax = get_catalog(i,Nreads)
    ceni = cull_catalog(ceni)
    ceni = shuffle_catalog(ceni)
    if i==0:
        cen = distribute_catalog(ceni)
    else:
        cen  = np.concatenate((cen,distribute_catalog(ceni)))

del ceni

#Center map on most massive halo
#cen = center_catalog(cen,Mmax)

#Get number of satellites for each halo
report('cen2ns starting',2)
ns  = cen2ns(cen)
report('cen2ns done',2)

#Populate halos with satellites
report('populate_centrals starting',2)
sat = populate_centrals(cen,ns)
report('populate_centrals done',2)

#Write timing
write_time("HOD completed", params.rank)

#np.savetxt('cen.txt', cen)

#nu_planck_GHz = ([217]) #([100,143,217,353,545,857])
print('make arrays')
#halos2map(cen,ns,sat)
C = np.zeros((np.shape(cen)[0],1),dtype='float32')
S = np.zeros((np.shape(sat)[0],2),dtype='float32')
print('make more arrays')
C[:,0] = cen[:,3] # Mass
S[:,:] = halocyno.cen2sat(np.column_stack((C[:,0],ns)),ns) # Parent Mass and Nsat

print("start saving")
np.save("centrals.npy",np.c_[cen,ns])
print('centrals saved')
np.save("satellites.npy",np.c_[sat[:,0],sat[:,1],sat[:,2],S[:,0],S[:,1]])
print('done')
'''
for nu_obs_GHz in nu_planck_GHz:

    #Put halos in map
    params.nu_obs = nu_obs_GHz * 1.e9
    cib_tot, flux_cen, flux_sat = halos2map(cen,ns,sat)
    params.comm.Barrier()

    #SAVE MAP TO DISK
    ns_str   = str(params.nside)
    nu_str   = str(nu_obs_GHz)
    zmin_str = str(params.min_redshift)
    zmax_str = str(params.max_redshift)

    #creating output file names
    dirout = "maps/" 
    if params.numdens==0: base = dirout+"cib_"
    if params.numdens==1: base = dirout+"opt_"

    if params.flat==0: base += "fullsky_"
    if params.flat==1: base += "flatsky_"

    if params.numdens==0: base += ('ns'+ns_str+'_zmin'+zmin_str+
                                '_zmax'+zmax_str+'_nu'+nu_str)
    if params.numdens==1: base += 'ns'+ns_str+'_zmin'+zmin_str+'_zmax'+zmax_str

#    if(params.rank==0): 
#        if params.flat==0:
#            file_tot = base+'_tot.fits'
#            hp.write_map(file_tot,cib_tot)
#            print '\nWrote Total Map' 

#        elif params.flat==1:
#            file_tot = base+'_tot.fits'
#            np.save(file_tot,cib_tot)
#    #Write timing
#    write_time("Halo projection completed", params.rank)
'''

