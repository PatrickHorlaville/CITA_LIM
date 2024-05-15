import gc
import psutil

# import pyximport; pyximport.install()
import halo_t as halocyno

import subprocess
import sys
import os 
import scipy

import numpy              as np
import healpy             as hp
import matplotlib.pylab   as plt
import scipy              as sp
import params             as params
import cosmolopy.distance as cd

from   pykTools           import *
from   progressbar        import *

from   iohod          import *
from   utils          import *
from   ana            import *
from   getparams      import *

from   cosmolopy          import *

from   scipy              import special
from   scipy.integrate    import *
from   scipy.interpolate  import *

from   mpi4py             import MPI

getparams('params.py')
params.proc     = psutil.Process(os.getpid())
params.overhead = params.proc.memory_info().rss
params.maxmem   = 0

params.comm = MPI.COMM_WORLD
params.rank = params.comm.Get_rank()
params.size = params.comm.Get_size()

fmt       = '%H:%M:%S on %m/%d/%Y'
timestamp = datetime.datetime.now().strftime(fmt)
if(params.rank==0):
    bar = 72*'-'
    print('')
    print(bar)
    print('Halo running on',params.size,'processors')
    print('Time:      '+timestamp)
    print('Directory: '+os.getcwd())
    print(bar)
    print('')

params.justify = 25*' '
