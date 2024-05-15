from   mods   import * 
from   utils  import *
import params as     params

def get_catalog(i,Nreads):
    '''
    Gets the file. Checks what the format given in the parameters file 
    is and returns x,y,z,M values as column stack data accordingly.
    '''
    
#    report('Reading catalog...',2)

    extension  = params.format 
    if extension=='pksc':
        peakdata  = readpks_lean(params.inputfile,i,Nreads)
        data      = peakdata.data
        rho       = 2.775e11*params.omegam*params.h**2
        data[:,3] = (4*(np.pi)/3)*(data[:,3]**3)*rho        
        Mmax      = (4*(np.pi/3))*(peakdata.RTHmax)**3*rho
            
    elif extension=='fits':
        peaks    = pykTools.fits2pyk(params.inputfile)    
        data = np.column_stack((peaks.x,peaks.y,peaks.z,peaks.M))

    elif extension=='asci':
        data = np.loadtxt(params.inputfile)
        data = np.column_stack((data[:,0],data[:,1],data[:,2],data[:,3]))

    else:
        print(" extension "+extension+" not recognized")

    return data, Mmax

def get_catalog_Non():

        # Non  = read_header(params.folder+params.inputfile)
        
        Non  = read_header(params.inputfile)

        return Non
