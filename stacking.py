from stacking_params import *


thresh = m_noise.halos.M > mass_cut

map_zs = (m_noise.mapinst.nu_rest/m_noise.mapinst.nu_bincents) - 1

halo_zs = m_noise.halos.redshift[thresh]

good_halo_zs = np.where(np.logical_and(halo_zs >= map_zs[ind] - err, halo_zs <= map_zs[ind] + err))

halo_xs = m_noise.halos.ra[thresh][good_halo_zs]
halo_ys = m_noise.halos.dec[thresh][good_halo_zs]

print('The total forecast observing time has been set to', t_obs)

CII_total_noise, CII_total_pure = lum(m_noise, ind, n, halo_xs, halo_ys)

CII_stack_noise, CII_stack_pure = np.nanmean(CII_total_noise, axis = 0), np.nanmean(CII_total_pure, axis = 0)


fig , axes = plt.subplots(nrows = 1, ncols = 2, figsize = (20, 10))

plt.subplot(121)
plt.imshow(CII_stack_pure, cmap = 'CMRmap')
plt.title(r'$Pure\ \rm{[CII]}\ Stack$', math_fontfamily = 'dejavuserif')
plt.xlabel(r'$RA-pixels$', math_fontfamily='dejavuserif')
plt.ylabel(r'$DEC-pixels$', math_fontfamily='dejavuserif')
plt.colorbar()



plt.subplot(122)
plt.imshow(CII_stack_noise, cmap = 'CMRmap')
plt.title(r'$FYST\ \rm{[CII]}\ Forecast\ Stack$', math_fontfamily = 'dejavuserif')
plt.xlabel(r'$RA-pixels$', math_fontfamily='dejavuserif')
plt.ylabel(r'$DEC-pixels$', math_fontfamily='dejavuserif')
cb = plt.colorbar()
cb.set_label(label = r'$Signal\ [Jy/sr]$', math_fontfamily = 'dejavuserif')

plt.show()


