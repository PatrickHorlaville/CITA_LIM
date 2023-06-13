from stacking_params import *



# Calculating luminosities

thresh = lim_sim.halos.M > mass_cut

map_zs = (lim_sim.mapinst.nu_rest/lim_sim.mapinst.nu_bincents) - 1

halo_zs = lim_sim.halos.redshift[thresh]
good_halo_zs = np.where(np.logical_and(halo_zs >= map_zs[ind] - err, halo_zs <= map_zs[ind] + err))

halo_xs = lim_sim.halos.ra[thresh][good_halo_zs]
halo_ys = lim_sim.halos.dec[thresh][good_halo_zs]
halo_zs = halo_zs[good_halo_zs]

print('------------------------')
print(' - The total forecast observing time has been set to', t_obs, '.')
print(' - Redshift of selected slice is', round(map_zs[ind], 2), 'and accepted halos are in the redshift range [', round(map_zs[ind] - err, 2), ',', round(map_zs[ind] + err, 2), '].')
print(' - Stacked map is', n, 'by', n, 'which covers', stack_dim, 'deg by', stack_dim, 'deg.')
print('------------------------')


pure_map, noisy_map = lum(lim_sim, n, halo_xs, halo_ys, halo_zs)

pure_stack, noisy_stack = np.nanmean(pure_map, axis = 0), np.nanmean(noisy_map, axis = 0)


# Plotting

fig , axes = plt.subplots(nrows = 2, ncols = 1, figsize = (8, 16))

plt.subplot(211)
plt.imshow(pure_stack, cmap = 'CMRmap', extent = [-stack_dim/2, stack_dim/2, -stack_dim/2, stack_dim/2])
plt.title(r'$Pure\ Signal\ Stacked\ Map$', math_fontfamily = 'dejavuserif')
plt.xlabel(r'$RA\ [Degrees]$', math_fontfamily='dejavuserif')
plt.ylabel(r'$Dec\ [Degrees]$', math_fontfamily='dejavuserif')
cb1 = plt.colorbar(label = r'$Jy/sr$')
ax1 = cb1.ax
text1 = ax1.yaxis.label
font = matplotlib.font_manager.FontProperties(math_fontfamily='dejavuserif', style='italic')
text1.set_font_properties(font)


plt.subplot(212)
plt.imshow(noisy_stack, cmap = 'CMRmap', extent = [-stack_dim/2, stack_dim/2, -stack_dim/2, stack_dim/2])
plt.title(r'$Forecast\ Stacked\ Map$', math_fontfamily = 'dejavuserif')
plt.xlabel(r'$RA\ [Degrees]$', math_fontfamily='dejavuserif')
plt.ylabel(r'$Dec\ [Degrees]$', math_fontfamily='dejavuserif')
cb2 = plt.colorbar(label = r'$Jy/sr$')
ax2 = cb2.ax
text2 = ax2.yaxis.label
font = matplotlib.font_manager.FontProperties(math_fontfamily='dejavuserif', style='italic')
text2.set_font_properties(font)

#fig.tight_layout()
plt.savefig('Stacking/stacking_ngauss_june12.png', bbox_inches = "tight")
plt.show()


