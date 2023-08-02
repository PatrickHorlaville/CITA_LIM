from stacking_params import *
from datetime import datetime
plt.rcParams["mathtext.fontset"] = "dejavuserif"


# Calculating luminosities

# 270-lightcone average:

lc_paths = '/home/dongwooc/scratchspace/pprun_hiz_npz/'
from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir(lc_paths) if isfile(join(lc_paths, f))]
onlyfiles.remove('pksc2npz_5313591.out')
onlyfiles.remove('pksc2npz.sh')
for i in range(len(onlyfiles)):
    onlyfiles[i] = lc_paths+onlyfiles[i]
  
 
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Beginning time =", current_time)
print("------")

for i in range(len(onlyfiles)):
    
    lim_sim.update(catalogue_file = f"{onlyfiles[i]}")
    
    print('Loading', i, 'th lightcone...')
    
    thresh = lim_sim.halos.M > mass_cut
    
    print('Finished loading', i, 'th lightcone! Now starting the stack computation...')

    map_zs = (lim_sim.mapinst.nu_rest/lim_sim.mapinst.nu_bincents) - 1

    halo_zs = lim_sim.halos.redshift[thresh]
    good_halo_zs = np.where(np.logical_and(halo_zs >= map_zs[ind] - err, halo_zs <= map_zs[ind] + err))
    
    halo_xs = lim_sim.halos.ra[thresh][good_halo_zs]
    halo_ys = lim_sim.halos.dec[thresh][good_halo_zs]
    halo_zs = halo_zs[good_halo_zs]

    print('------------------------')
    print(' - The total forecast observing time has been set to', t_obs, '-')
    print(' - Redshift of selected slice is', round(map_zs[ind], 2), 'and accepted halos are in the redshift range [', round(map_zs[ind] - err, 2), ',', round(map_zs[ind] + err, 2), '], which accounts for', len(halo_xs), 'halos -')
    print(' - The lightcone is in the redshift range z = [', round(np.min(map_zs), 3), ', ', round(np.max(map_zs), 3), '] -') 
    print(' - Stacked map is', n, 'by', n, 'which covers', stack_dim, 'deg by', stack_dim, 'deg -')
    print(' - We beam the stacked map with a width of', beam_width, ', which corresponds to a Gaussian filter of radius', beam_res, 'pixels - ')
    print('------------------------')
    
    pure_map, noisy_map = lum(lim_sim, n, halo_xs, halo_ys, halo_zs)
    pure_stack, noisy_stack = np.nanmean(pure_map, axis = 0), np.nanmean(noisy_map, axis = 0)
    
    np.save('/mnt/scratch-lustre/horlaville/nuObs270/zdex04/alpha_cii_0-024/alpha_mhi_0-74/Mmin_2-10e10/stacks_z6_6/sig/sig'+str(i)+'.npy', pure_stack)
    np.save('/mnt/scratch-lustre/horlaville/nuObs270/zdex04/alpha_cii_0-024/alpha_mhi_0-74/Mmin_2-10e10/stacks_z6_6/for/for'+str(i)+'.npy', noisy_stack)
    
    print('Finished loading', i, 'th stack!')
    
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("------")
print("Finished at =", current_time)



