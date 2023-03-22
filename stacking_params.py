from stacking_utils import *


# Set the limlam_mocker object
lim_sim = lim('Lichen_v2', doSim=True)

# Default tobs is 2000h. Update if desired
t_obs = 40000 * u.hr
lim_sim.update(tobs = t_obs) 

# Lower bound of halo masses to be considered
mass_cut = 2*(10**10) # in Solar masses

# Error of acceptance on the redshift of the 2D map for surrounding halos
err = 0.03

# Index of the redshift map to be selected. There are 21 redshift slices by default, so ind can be any integer between 0 and 20
ind = 3

# Size of the stacked map to be produced (n by n)
n = 50
