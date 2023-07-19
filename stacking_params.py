from stacking_utils import *


# Set the limlam_mocker object
lim_sim = lim('Lichen_v4', doSim=True)

# Update parameters as desired: tobs for noise, catalogue_file for input lightcone
# Be aware: make sure to adapt nuObs and Delta_nu to your lightcone
# dnu for the amount of redshift slices

t_obs = 40000 * u.hr

# [C II] example:
lim_sim.update(model_par = {'zdex': 0.3,
                            'M0': 1900000000.0,
                            'Mmin': 20000000000,
                            'alpha_MH1': 0.74,
                            'alpha_LCII': 0.017,
                            'BehrooziFile': 'sfr_reinterp.dat'},
               tobs = t_obs,
               nuObs = 270*u.GHz,
               Delta_nu = 40*u.GHz,
               Omega_field = 4*(u.deg**2),
               dnu = 1.11*u.GHz,
               catalogue_file = '/home/dongwooc/scratchspace/pprun_hiz_npz/COMAP_z5.8-7.9_960Mpc_seed_13819.npz'
            )


# Lower bound of halo masses to be considered
mass_cut = 2*(10**10) # in Solar masses

# Error of acceptance on the redshift of the 2D map for surrounding halos
err = 0.03

# Index of the redshift map to be selected.
ind = 17

# Size of the stacked map to be produced (n by n)
n = 50

# Size of the stacked map, in angular dimension
ang_side = np.sqrt(lim_sim.Omega_field)
nside = lim_sim.Nside
stack_dim = round(n*(ang_side/nside).value, 2)
