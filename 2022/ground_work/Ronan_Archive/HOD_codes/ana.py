from mods import *

def tSZxCIB(tsz,cib):

    report('Correlating...',2)

    cib = hp.ud_grade(cib,1024)

    cl_tsz = hp.anafast(tsz)
    cl_cib = hp.anafast(cib)
    cl_txc = hp.anafast(tsz,cib)

    ell_cib = np.arange(len(cl_cib))+1
    ell_tsz = np.arange(len(cl_tsz))+1
    ell_txc = np.arange(len(cl_txc))+1

    cl_tsz *= (ell_tsz+1)*ell_tsz/(2*np.pi)
    cl_cib *= (ell_cib+1)*ell_cib/(2*np.pi)
    cl_txc *= (ell_txc+1)*ell_txc/(2*np.pi)

    bins = np.linspace(np.log(10.),np.log(len(ell_cib)),50)

    clb_tsz, bins = np.histogram(np.log(ell_cib),bins,weights=cl_tsz)
    clb_cib, bins = np.histogram(np.log(ell_cib),bins,weights=cl_cib)
    clb_txc, bins = np.histogram(np.log(ell_cib),bins,weights=cl_txc)

    num, bins = np.histogram(np.log(ell_cib),bins)

    center = (bins[:-1]+bins[1:])/2

    clb_tsz /= num
    clb_cib /= num
    clb_txc /= num

    ell = np.exp(center)
    
    data = np.column_stack((ell,clb_tsz,clb_cib,clb_txc))
    report('Writing power spectra to '+params.clfile,1)
    np.savetxt(params.clfile,data,fmt="%.5e")


