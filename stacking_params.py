from stacking_utils import *


# Set the limlam_mocker object
lim_sim = lim('Lichen_v4', doSim=True)

# Update parameters as desired: tobs for noise, catalogue_file for input lightcone
# Be aware: make sure to adapt nuObs and Delta_nu to your lightcone
# dnu for the amount of redshift slices

t_obs = 2000 * u.hr

lim_sim.update(model_par = {'zdex': 0.3,
 'M0': 1900000000.0,
 'Mmin': 20000000000,
 'alpha_MH1': 0.74,
 'alpha_LCII': 0.017,
 'BehrooziFile': 'sfr_reinterp.dat'}, 
             tobs = t_obs,
             nuObs = 400*u.GHz,
             Delta_nu = 28*u.GHz,
             Omega_field = 16*(u.deg**2),
               catalogue_file = '/mnt/scratch-lustre/horlaville/nate_sims/june_9/ngaussian/ngauss_mlambda25_0.npz'
            )

# Lower bound of halo masses to be considered
mass_cut = 2*(10**10) # in Solar masses

# Error of acceptance on the redshift of the 2D map for surrounding halos
err = 0.03

# Index of the redshift map to be selected. There are 21 redshift slices by default, so ind can be any integer between 0 and 20
ind = 3

# Size of the stacked map to be produced (n by n)
n = 75

# Size of the stacked map, in angular dimension
ang_side = np.sqrt(lim_sim.Omega_field)
nside = lim_sim.Nside
stack_dim = round(n*(ang_side/nside).value, 2)
