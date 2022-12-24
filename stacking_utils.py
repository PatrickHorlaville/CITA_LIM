import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy.units as u

from lim import lim

matplotlib.rcParams.update({'font.size': 18,'figure.figsize':[8,7]}) 
from scipy.ndimage import gaussian_filter



def halo_centpix(CII_map, halo_xpos, halo_ypos):

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

    halo_centpix_x = [0 for i in range(len(halo_xpos))]
    halo_centpix_y = [0 for i in range(len(halo_ypos))]

    for i in range(len(halo_xpos)):
        halo_centpix_x[i] = np.argmin(np.abs(halo_xpos[i] - CII_map.mapinst.pix_bincents_x))
    for i in range(len(halo_ypos)):
        halo_centpix_y[i] = np.argmin(np.abs(halo_ypos[i] - CII_map.mapinst.pix_bincents_y))

    return halo_centpix_x, halo_centpix_y







def halo_map(CII_map, n, halo_xpos, halo_ypos):

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

    halo_centpix_x, halo_centpix_y = halo_centpix(CII_map, halo_xpos, halo_ypos)

    for i in range(len(halo_xpos)):
        halo_mapx[i] = np.linspace(halo_centpix_x[i] - (n - 1)/2, halo_centpix_x[i] + (n - 1)/2, n)
    for i in range(len(halo_ypos)):
        halo_mapy[i] = np.linspace(halo_centpix_y[i] - (n - 1)/2, halo_centpix_y[i] + (n - 1)/2, n)

    return halo_mapx, halo_mapy




def lum(CII_map, ind, n, halo_xpos, halo_ypos):
    
    '''
    Parameters
    ----------
    
    CII_map: lim_object
             The `lim` object of the generated [CII] simulation
    
    n: int
       Size of the stacked map
       
    halo_xpos: array_like
               List of RA positions of halos
    
    halo_ypos: array_like
               List of DEC positions of halos
    
    Returns
    -------
    
    lum: array_like
         n by n by n_halo list of the [CII] luminosities of all the halo maps to be averaged in order to retrieve the stacked map    
    '''
    
    halo_mapx, halo_mapy = halo_map(CII_map, n, halo_xpos, halo_ypos)
    
    npix_x, npix_y = CII_map.mapinst.npix_x + 1, CII_map.mapinst.npix_y + 1
    n_halos = len(halo_xpos)

    lum_noisy   = [[[[] for i in range(n)] for j in range(n)] for k in range(n_halos)]
    lum_nonoise = [[[[] for i in range(n)] for j in range(n)] for k in range(n_halos)]

    for i in range(n_halos):
        for j in range(len(halo_mapx[i])):
            for k in range(len(halo_mapy[i])):
                if int(halo_mapx[i][j]) < npix_x:
                    if int(halo_mapy[i][k]) < npix_y:
                            lum_noisy[i][j][k]   = CII_map.noise_added_map[int(halo_mapx[i][j]), int(halo_mapy[i][k]), ind]
                            lum_nonoise[i][j][k] = (CII_map.maps[int(halo_mapx[i][j]), int(halo_mapy[i][k]), ind]).value
                    else:
                        lum_noisy[i][j][k]   = np.nan
                        lum_nonoise[i][j][k] = np.nan
                else:
                    lum_noisy[i][j][k]   = np.nan
                    lum_nonoise[i][j][k] = np.nan
        
    return lum_noisy, lum_nonoise



