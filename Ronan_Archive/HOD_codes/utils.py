from mods import *

def r2m(r):
    rho = 2.78e11 * params.omegam * (params.h**2)
    return 4./3.*np.pi*rho*r**3

def m2r(m):
    rho = 2.78e11 * params.omegam * (params.h**2)
    return (3.*m/4./np.pi/rho)**(1./3.)

def random_phi(N):
    return (np.random.uniform(size=N)*2*np.pi).astype('float32')

def random_theta(N):
    return (np.arccos(2*np.random.uniform(size=N)-1)).astype('float32')

def resetoverhead():
    params.overhead = params.proc.get_memory_info().rss

def mz2c(m,z):
    # concentration factor from Duffy et al. 2008
    return 7.85 * (m / (2e+12/params.h))**-0.081 / (1+z)**0.71
               
def r2z(r):
	zrange   = np.linspace(0,4,1000)
	r_to_z   = sp.interpolate.interp1d(
		cd.comoving_distance_transverse(
			zrange, **params.cosmo), zrange)

	return r_to_z(r).astype('float32')

def report(description,verbosity):
    if(params.rank==0 and params.verbose>=verbosity): 
        print(params.justify,description)

def check_memory(description,N):
  
    memt = params.proc.memory_info().rss
    mem = memt - params.overhead
    mem_per_object = mem/(N*1.0)
    flt_per_object = mem_per_object / 4

    mem /= 1024.**3
    memt /= 1024.**3

    params.maxmem=max(params.maxmem,memt)

    print('                     ',description,memt,'GB total')
    print('                         ',mem_per_object,'bytes per array element')
    print('                         ',flt_per_object,'floats per array element')
  
    return

def distribute_catalog(data):

#    report('Distributing...',2)

    N = np.shape(data)[0]

    Nextra = N % params.size

    Nlocal = (N-Nextra) / (params.size)
        
    start = Nlocal*(params.rank  )
    end   = Nlocal*(params.rank+1)

    if(params.rank==params.size-1): end=N

    Nlocal = end - start

    return data[int(start):int(end),:]


def shuffle_catalog(data):

#    report('Shuffling...',2)
    np.random.seed(13579) #added seed or each process does its own permutation
                          # this was an error before and halos were double counted
    p = np.random.permutation(np.shape(data)[0])

    return data[p,:]

def cull_catalog(data):

#    report('Culling...',2)

    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    M = data[:,3]

    r = np.sqrt((x**2)+(y**2)+(z**2))
    redshift = r2z(r)
    
#    print 'total number of halos before culling: ',np.shape(x)[0]
#    print 'min and max redshift of halos: ',redshift.min(),redshift.max()
#    print 'max redshift:',params.max_redshift

    # filtering halos in the sphere of r<box/2 and z<max_redshift
    dm = ([
            (redshift > params.min_redshift ) & 
            (redshift < params.max_redshift ) & 
            (  abs(r) < (params.box_size)/2 ) & 
            (       M > params.min_mass     )
            ])
    
    x = x[dm]
    y = y[dm]
    z = z[dm]
    M = M[dm]
    
    r = np.sqrt((x**2)+(y**2)+(z**2))
    redshift = r2z(r)

#    print "r max min mean = ", np.max(r), np.min(r), np.mean(r)   
#    print "redshift max min mean= ", np.max(redshift), np.min(redshift), np.mean(redshift)
#    print 'total number of halos after culling: ',np.shape(x)[0]
#    print 'min and max redshift of halos: ',redshift.min(),redshift.max()
#    print 'max redshift:',params.max_redshift

    return np.column_stack((x,y,z,M))

def jiang_shmf(m,M_halo):	
    gamma_1    = 0.13
    alpha_1    = -0.83
    gamma_2    = 1.33
    alpha_2    = -0.02
    beta_2     = 5.67
    zeta       = 1.19 	
    
    dndm = (((gamma_1*((m/M_halo)**alpha_1))+
             (gamma_2*((m/M_halo)**alpha_2)))*
            (np.exp(-(beta_2)*((m/M_halo)**zeta))))
    
    return dndm
'''
def center_catalog(data,Mmax):
   #centre on largest halo
   # find x,y,z position amd angles of largest halo     
    xmmaxl = np.array([-1.0e5],dtype='float64')
    ymmaxl = np.array([-1.0e5],dtype='float64')
    zmmaxl = np.array([-1.0e5],dtype='float64')
    xmmax  = np.array([0.0],dtype='float64')
    ymmax  = np.array([0.0],dtype='float64')
    zmmax  = np.array([0.0],dtype='float64')

    indi=np.where(data[:,3] > Mmax*0.99)[0]
    if len(indi)>0:
        indi=int(indi)
        print("Max halo mass ", data[indi,3])
        print("index and pos of largest halo", indi, data[indi,2], data[indi,1], data[indi,0])
        xmmax = np.array([data[indi,2]],dtype='float64')
        ymmax = np.array([data[indi,1]],dtype='float64')
        zmmax = np.array([data[indi,0]],dtype='float64')
        xmmaxl = xmmax
        ymmaxl = ymmax
        zmmaxl = zmmax
    else:
        indi=-1
    #rotate so largest halo is in first quadrant to make angles easier
    params.comm.Allreduce(xmmaxl,xmmax, op=MPI.MAX)
    params.comm.Allreduce(ymmaxl,ymmax, op=MPI.MAX)
    params.comm.Allreduce(zmmaxl,zmmax, op=MPI.MAX)

    if(xmmax <0.0):   data[:,2] = -data[:,2]
    if(ymmax <0.0):   data[:,1] = -data[:,1]
    if(zmmax <0.0):   data[:,0] = -data[:,0]
    xmmax = abs(xmmax)
    ymmax = abs(ymmax)
    zmmax = abs(zmmax)

    rmaxtheta = np.arcsin(zmmax/np.sqrt(xmmax**2+zmmax**2)) # rotate in z-x plane 
    xmmaxl    = np.sqrt(xmmax**2+zmmax**2)
    rmaxphi   = np.arcsin(ymmax/np.sqrt(ymmax**2+xmmaxl**2))  # rotate in z-y plane 

    xh       =   data[:,2]*np.cos(rmaxtheta) + data[:,0]*np.sin(rmaxtheta)
    data[:,0] =  -data[:,2]*np.sin(rmaxtheta) + data[:,0]*np.cos(rmaxtheta)

    data[:,2] =   xh*np.cos(rmaxphi) + data[:,1]*np.sin(rmaxphi)
    data[:,1] =  -xh*np.sin(rmaxphi) + data[:,1]*np.cos(rmaxphi)

    if(indi>0):
        print("after rotation: x, y, z, RTH = ",data[indi,2],data[indi,1],data[indi,0],data[indi,3])
        print("at a distance of ", (data[indi,0]**2+data[indi,1]**2+data[indi,2]**2)**(1./2))

    return data
'''
