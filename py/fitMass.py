import copy
import numpy
from galpy.potential import MiyamotoNagaiPotential, \
    NFWPotential,\
    HernquistPotential,\
    DoubleExponentialDiskPotential, \
    evaluateRforces    
from scipy import optimize, integrate
from fitvc import _REFR0, _REFV0, cb
_convertrho= (_REFV0**2./_REFR0**2.*0.0002321707838637446)**-1.
_convertmass= (_REFV0**2.*_REFR0*10.**9.*0.0002321707838637446)/10.**(12)
_rhodm= 0.008*_convertrho
_rhodmerr= 0.003*_convertrho
_diskscale=3.25/_REFR0
_diskscaleerr=.25/_REFR0
_dexp= False
def obj(x,data,bp,dp,hp):
    #x=[fd,fh,vc,rs,hdisk]
    if x[0] > 1. or x[0] < 0.: return numpy.finfo(numpy.dtype(numpy.float64)).max
    if x[1] > 1. or x[1] < 0.: return numpy.finfo(numpy.dtype(numpy.float64)).max
    if (1.-x[0]-x[1]) > 1. or (1.-x[0]-x[1]) < 0.: return numpy.finfo(numpy.dtype(numpy.float64)).max
    if x[2] < 0. or x[3] < 0. or x[4] < 0.: return numpy.finfo(numpy.dtype(numpy.float64)).max
    #Renormalize potentials, intially normalized to 1./3. each
    if False:
        bp= copy.deepcopy(bp)
        dp= copy.deepcopy(dp)
        hp= copy.deepcopy(hp)
        bp._amp*= (1.-x[0]-x[1])**2.*9.
        dp._amp*= x[0]**2.*9.
        hp._amp*= x[1]**2.*9.
        #Re-define disk scale length and halo scale length
        if _dexp:
            dp._hr= x[4]
        else:
            dp._a= x[4]
        hp._a= x[3]
    else:
        #Set-up
        bp= HernquistPotential(a=0.6/_REFR0,normalize=1.-x[0]-x[1])
        if _dexp:
            dp= DoubleExponentialDiskPotential(normalize=x[0],
                                               hr=x[4],
                                               hz=0.3/_REFR0)
        else:
            dp= MiyamotoNagaiPotential(normalize=x[0],
                                       a=x[4],
                                       b=0.3/_REFR0)
        hp= NFWPotential(normalize=x[1],
                         a=x[3])
    #Re-normalize data
    vcircdata= copy.copy(data)
    vcircdata[:,1]/= x[2]
    vcircdata[:,2]/= x[2]
    #Vcirc chi2
    vcmodel= numpy.zeros(vcircdata.shape[0])
    for ii in range(vcircdata.shape[0]):
        vcmodel[ii]= numpy.sqrt(vcircdata[ii,0]\
                                    *numpy.fabs(evaluateRforces(vcircdata[ii,0],
                                                                0.,[bp,dp,hp])))
    #print vcircdata[:,0], vcmodel
    chi2= numpy.sum((vcircdata[:,1]-vcmodel)**2./vcircdata[:,2]**2.)
    #Add scale length measurement
    chi2+= (x[4]-_diskscale)**2./_diskscaleerr**2.
    #Add dark matter density at the Solar radius
    print hp.dens(1.,0.),_rhodm*x[2]**2.
    chi2+= (hp.dens(1.,0.)-_rhodm*x[2]**2.)**2./_rhodmerr**2./x[2]**4.
    return chi2

def fitMass():
    #Read data
    vcircdata= numpy.loadtxt('vcirc.txt',comments='#',delimiter='|')
    vcircdata[:,0]/= _REFR0
    vcircdata[:,1]/= _REFV0
    vcircdata[:,2]/= _REFV0
    #Set-up
    bp= HernquistPotential(a=0.6/_REFR0,normalize=1./3.)
    if _dexp:
        dp= DoubleExponentialDiskPotential(normalize=1./3.,
                                           hr=3.25/_REFR0,
                                           hz=0.3/_REFR0)
    else:
        dp= MiyamotoNagaiPotential(normalize=1./3.,
                                   a=3.25/_REFR0,
                                   b=0.3/_REFR0)
    hp= NFWPotential(normalize=1./3.,
                     a=5./_REFR0)
    init= numpy.array([0.6,0.35,218./_REFV0,15./_REFR0,3.25/_REFR0])
    #print _convertrho
    out= optimize.fmin_powell(obj,init,args=(vcircdata,bp,dp,hp),
                              callback=cb)
    print out
    #Calculate mass
    #Set-up best-fit
    bp= HernquistPotential(a=0.6/_REFR0,normalize=1.-out[0]-out[1])
    if _dexp:
        dp= DoubleExponentialDiskPotential(normalize=out[0],
                                           hr=out[4],
                                           hz=0.3/_REFR0)
    else:
        dp= MiyamotoNagaiPotential(normalize=out[0],
                                   a=out[4],
                                   b=0.3/_REFR0)
    hp= NFWPotential(normalize=out[1],
                     a=out[3])
    #Halo mass
    halo= integrate.quad((lambda x: hp.dens(x,0.)*x**2.),0.,250./_REFR0)[0]*4.*numpy.pi*_convertmass/out[2]**2.
    bulge= integrate.quad((lambda x: bp.dens(x,0.)*x**2.),0.,250./_REFR0)[0]*4.*numpy.pi*_convertmass/out[2]**2.
    disk= integrate.dblquad((lambda x,y: dp.dens(y,x)*y),0.,250./_REFR0,lambda x: 0.,lambda x: 20./_REFR0)[0]*2.*numpy.pi*_convertmass/out[2]**2.
    print halo, disk, bulge, halo+disk+bulge
    return None

if __name__ == '__main__':
    fitMass()
