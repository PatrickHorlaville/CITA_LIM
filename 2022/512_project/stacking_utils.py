import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy.units as u

from lim import lim

matplotlib.rcParams.update({'font.size': 18,'figure.figsize':[8,7]}) 
from scipy.ndimage import gaussian_filter



def halo_centpix(lim_obj, halo_xpos, halo_ypos, halo_zpos):

    '''
    Parameters
    ----------

    halo_xpos: array_like
               list of RA coordinates of halos, from the original distribution of halos

    halo_ypos: array_like
               list of DEC coordinates of halos, from the original distribution of halos


    Returns
    -------

    halo_centpix_x: array_like
                    list of the x-pixel positions of halos, mapped on the lim map


    halo_centpix_x: array_like
                    list of the y-pixel positions of halos, mapped on the lim map

    '''

    map_xs = lim_obj.mapinst.pix_bincents_x
    map_ys = lim_obj.mapinst.pix_bincents_y
    map_zs = (lim_obj.mapinst.nu_rest/lim_obj.mapinst.nu_bincents) - 1
    
    pixcents_x_mesh, halo_xs_mesh = np.meshgrid(map_xs, halo_xpos)
    halo_centpix_x = np.argmin(np.abs(halo_xs_mesh - pixcents_x_mesh), axis=1)
    
    pixcents_y_mesh, halo_ys_mesh = np.meshgrid(map_ys, halo_ypos)
    halo_centpix_y = np.argmin(np.abs(halo_ys_mesh - pixcents_y_mesh), axis=1)
    
    pixcents_z_mesh, halo_zs_mesh = np.meshgrid(map_zs, halo_zpos)
    halo_centpix_z = np.argmin(np.abs(halo_zs_mesh - pixcents_z_mesh), axis=1)
    
    return halo_centpix_x, halo_centpix_y, halo_centpix_z







def halo_map(lim_obj, n, halo_xpos, halo_ypos):

    '''
    Parameters
    ----------

    n: int
       Size of the desired stacked map

    halo_xpos: array_like
               List of RA positions of halos

    halo_ypos: array_like
               List of DEC positions of halos


    Returns
    -------

    halo_mapx: array_like
               Range of x-pixels for the stacked map

    halo_mapy: array_like
               Range of y-pixels for the stacked map
    '''

    halo_mapx = [0 for i in range(len(halo_xpos))]
    halo_mapy = [0 for i in range(len(halo_ypos))]

    halo_centpix_x, halo_centpix_y, _ = halo_centpix(lim_obj, halo_xpos, halo_ypos)

    for i in range(len(halo_xpos)):
        halo_mapx[i] = np.linspace(halo_centpix_x[i] - (n - 1)/2, halo_centpix_x[i] + (n - 1)/2, n)
    for i in range(len(halo_ypos)):
        halo_mapy[i] = np.linspace(halo_centpix_y[i] - (n - 1)/2, halo_centpix_y[i] + (n - 1)/2, n)

    return halo_mapx, halo_mapy




def lum(lim_obj, ind, n, halo_xpos, halo_ypos, halo_zpos):
    
    '''
    Parameters
    ----------
    
    lim_obj: lim_object
             The `lim` object of the generated line simulation
    
    n: int
       Size of the stacked map
       
    halo_xpos: array_like
               List of RA positions of halos
    
    halo_ypos: array_like
               List of DEC positions of halos
    
    Returns
    -------
    
    lum: array_like
         n by n by n_halo list of the luminosities of all the halo maps to be averaged in order to retrieve the stacked map    
    '''
    
    halo_mapx, halo_mapy = halo_map(lim_obj, n, halo_xpos, halo_ypos, halo_zpos)
    
    npix_x, npix_y = lim_obj.mapinst.npix_x + 1, lim_obj.mapinst.npix_y + 1
    n_halos = len(halo_xpos)

    lum_pure = [[[[] for i in range(n)] for j in range(n)] for k in range(n_halos)]
    lum_noisy= [[[[] for i in range(n)] for j in range(n)] for k in range(n_halos)]

    pure_map = lim_obj.maps.value
    noisy_map = lim_obj.noise_added_map
    
    for i in range(n_halos):
        for j in range(len(halo_mapx[i])):
            for k in range(len(halo_mapy[i])):
                if int(halo_mapx[i][j]) < npix_x:
                    if int(halo_mapy[i][k]) < npix_y:
                            lum_pure[i][j][k] = pure_map[int(halo_mapx[i][j]), int(halo_mapy[i][k]), ind]
                            lum_noisy[i][j][k]= noisy_map[int(halo_mapx[i][j]), int(halo_mapy[i][k]), ind]
                    else:
                        lum_noisy[i][j][k]   = np.nan
                        lum_nonoise[i][j][k] = np.nan
                else:
                    lum_noisy[i][j][k]   = np.nan
                    lum_nonoise[i][j][k] = np.nan
        
    return lum_pure, lum_noisy



