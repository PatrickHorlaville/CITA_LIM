import numpy as np
import pylab as pl
import scipy as sp
import os

from astropy.io import fits
import datetime
from cosmolopy import fidcosmo
import cosmolopy.distance as cosdist

class pykTable():
    """
    @brief Class describing a PeakPatch Table.
    """
    def __init__(self):
        pass

    def copy(self):
        """
        @brief Creates a copy of the  PeakPatch Table.
        """
        return copy.copy(self)

    def set_cosmo(self):
        cosmo = fidcosmo
        
        cosmo['h'] =  self.cosmo.h
        cosmo['n'] = self.cosmo.ns
        cosmo['omega_M_0'] = self.cosmo.omegam
        cosmo['omega_b_0'] = self.cosmo.omegab
        cosmo['omega_lambda_0'] = self.cosmo.omegal
        cosmo['sigma_8'] =  self.cosmo.sigma8
        
        return cosmo
        
    def set_redshifts(self):
        cosmo = self.set_cosmo() 
        
        #range of redshifts
        zrange=np.linspace(2,4,100)
        
        #distance of each point to the origin
        r=np.sqrt((self.x**2)+(self.y**2)+(self.z**2))
        
        r_to_z = sp.interpolate.interp1d(cosdist.comoving_distance_transverse(zrange, **cosmo), zrange)
        redshift = r_to_z(r)
        
        self.redshift = redshift
        self.distance = r



class cosmology():
    """
    @brief Class describing a cosmology.
    """
    def __init__(self, omegab=0.045, omegac=0.255, omegal=0.7, h=0.7,ns=0.96, sigma8=0.8):
        self.omegab= omegab#density of baryons    
        self.omegac=omegac #density of dark matter
        self.omegam=omegab+omegac
        self.omegal=omegal
        self.h=h
        self.ns=ns
        self.sigma8=sigma8


def readpks(filename):
    """
    @brief Reads peakpatch output in binary format into an object.
    """
    pkfile=open(filename,"rb")
    peakdata = pykTable()

    peakdata.Non = np.fromfile(pkfile, dtype=np.int32, count=1).astype('int')[0]
    peakdata.RTHmax = np.fromfile(pkfile, dtype=np.float32, count=1).astype('float')[0]
    npkdata=7*peakdata.Non
    
    print(peakdata.Non)

    peakdatain = np.fromfile(pkfile, dtype=np.float32, count=npkdata)
    peakdata.data = np.reshape(peakdatain,(peakdata.Non,7)).astype('float32')
    #xon(i),yon(i),zon(i),vxon(i),vyon(i),vzon(i),Fcollv(i),&
    #RTHLv(i),vTHvir(i),iwas(i),lagrange(i),hp(i),F_ev(i),F_pv(i),(strain_bar(k,i),k=1,6)
    return peakdata

def read_header(filename):
    """
    @brief Reads peakpatch output in binary format into an object.
    """
    pkfile   = open(filename,"rb")
    peakdata = pykTable()

    peakdata.Non    = np.fromfile(pkfile, dtype=np.int32, count=1).astype('int')[0]

    pkfile.close()

    return peakdata.Non


def readpks_lean(filename,i,Nreads):
    """
    @brief Reads peakpatch output in binary format into an object.
    """
    
    peakdata = pykTable()

    pkfile   = open(filename,"rb")

    peakdata.Non        = np.fromfile(pkfile, dtype=np.int32, count=1).astype('int')[0]
    peakdata.RTHmax     = np.fromfile(pkfile, dtype=np.float32, count=1).astype('float')[0]
    peakdata.redshiftin = np.fromfile(pkfile, dtype=np.float32, count=1).astype('float')[0]

    start = int(np.ceil(float(peakdata.Non)/Nreads)*i)
    end   = int(min( peakdata.Non , np.ceil(float(peakdata.Non)/Nreads)*(i+1) ))

    Nread = end-start

    outnum = 10 #10 floats per halo
    npkdata = Nread*outnum

    pkfile.seek( (3+start*outnum)*4 ,os.SEEK_SET)  #12 byte header plus array
    peakdata.data = np.fromfile(pkfile, dtype=np.float32, count=npkdata)
    peakdata.data = np.reshape(peakdata.data,(Nread,outnum)).astype('float32')

    pkfile.close()
    peakdata.data = np.column_stack((
                       peakdata.data[:,2],
                       peakdata.data[:,1],
                       peakdata.data[:,0],
                       peakdata.data[:,6]))            

    return peakdata


def pyk2fits(filename, peakdata, cosmo=cosmology()):
    """
    @brief Writes a peakpatch object into a fits file.
    """
    tbhdu = fits.new_table(
        fits.ColDefs([fits.Column(name='xon', format='E', array=peakdata.data[:,0]),
                      fits.Column(name='yon', format='E', array=peakdata.data[:,1]),
                      fits.Column(name='zon', format='E', array=peakdata.data[:,2]),
                      fits.Column(name='vxon', format='E', array=peakdata.data[:,3]),
                      fits.Column(name='vyon', format='E', array=peakdata.data[:,4]),
                      fits.Column(name='vzon', format='E', array=peakdata.data[:,5]),
                      fits.Column(name='Fcollv', format='E', array=peakdata.data[:,6]),
                      fits.Column(name='RTHLv', format='E', array=peakdata.data[:,7]),
                      fits.Column(name='vTHvir', format='E', array=peakdata.data[:,8]),
                      fits.Column(name='iwas', format='E', array=peakdata.data[:,9]),
                      fits.Column(name='lagrange', format='E', array=peakdata.data[:,10]),
                      fits.Column(name='hp', format='E', array=peakdata.data[:,11]),
                      fits.Column(name='F_ev', format='E', array=peakdata.data[:,12]),
                      fits.Column(name='F_pv', format='E', array=peakdata.data[:,13]),
                      fits.Column(name='strain_bar1', format='E', array=peakdata.data[:,14]),
                      fits.Column(name='strain_bar2', format='E', array=peakdata.data[:,15]),
                      fits.Column(name='strain_bar3', format='E', array=peakdata.data[:,16]),
                      fits.Column(name='strain_bar4', format='E', array=peakdata.data[:,17]),
                      fits.Column(name='strain_bar5', format='E', array=peakdata.data[:,18]),
                      fits.Column(name='strain_bar6', format='E', array=peakdata.data[:,19]),
                      ]),
        )
  


    prihdr = fits.Header()
    prihdr['Card1'] = 'Cosmological parameters'     
    prihdr['sep1'] = '-----------------------------------------------'
    prihdr['omegab'] =  cosmo.omegab
    prihdr['omegab'] = cosmo.omegab
    prihdr['omegac'] = cosmo.omegac 
    prihdr['omegam'] = cosmo.omegam
    prihdr['omegal'] = cosmo.omegal
    prihdr['h'] = cosmo.h
    prihdr['ns'] = cosmo.ns
    prihdr['sigma8'] = cosmo.sigma8

    prihdr['Card2'] = "PeakPatch Parameters"
    prihdr['sep2'] = '-----------------------------------------------'
    prihdr['nrec'] = peakdata.nrec 
    prihdr['nic'] = peakdata.nic
    prihdr['Non'] = peakdata.Non
    prihdr['biasold'] = peakdata.biasold 
    prihdr['biasnew'] = peakdata.biasnew
    prihdr['zzt'] = peakdata.zzt
    prihdr['h0'] = peakdata.h0 
    prihdr['omb'] = peakdata.omb 

    prihdr['Card3'] = "Date and Other Info"
    prihdr['sep3'] = '-----------------------------------------------'
    date=datetime.date.today()
    prihdr['Date'] = '%d-%d-%d'%(date.year, date.month, date.day)

    prihdu = fits.PrimaryHDU(header=prihdr)
    
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(filename)

    
    
    

def fits2pyk(filename):
    hdr = fits.getheader(filename)
    data = fits.getdata(filename)

    pt = pykTable()
    pt.cosmo = cosmology(omegab=hdr['OMEGAB'], omegac=hdr['OMEGAC'], omegal=hdr['OMEGAL'], h=hdr['H'],ns=hdr['NS'], sigma8=hdr['SIGMA8'])

    pt.x = data['xon'].astype('float')
    pt.y = data['yon'].astype('float')
    pt.z = data['zon'].astype('float')
    
    rho =  2.78e11 * pt.cosmo.omegam * pt.cosmo.h**2 # mean comoving total matter 
                                                     # density in Msun / Mpc^3
    
    RTH      = data['RTHLv'].astype('float')                            #radius
    Mhalo    = 4.0/3*np.pi*RTH**3*rho                #mass of halo in units of solar mass
    pt.M = Mhalo 
   
   
    pt.set_redshifts()                                # compute redshifts
    
  
    return pt


def write_time(string_in, rank):
    if rank==0:
        fmt       = '%H:%M:%S on %m/%d/%Y'
        timestamp = datetime.datetime.now().strftime(fmt)
        bar = 72*'-'
        print('')
        print(bar)
        print(string_in)
        print('Time:      '+timestamp)
        print(bar)
        print('')

    return
