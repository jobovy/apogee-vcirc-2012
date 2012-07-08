import copy
import numpy
from galpy.util import bovy_plot
from galpy.potential import MiyamotoNagaiPotential, \
    NFWPotential,\
    HernquistPotential,\
    DoubleExponentialDiskPotential, \
    evaluateRforces    
from scipy import optimize, integrate
from fitvc import _REFR0, _REFV0, cb
import bovy_mcmc
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
    if x[2] > 2. or x[3] > 10. or x[4] > 1.: return numpy.finfo(numpy.dtype(numpy.float64)).max
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
    #print hp.dens(1.,0.),_rhodm*x[2]**2.
    chi2+= (hp.dens(1.,0.)-_rhodm*x[2]**2.)**2./_rhodmerr**2./x[2]**4.
    return chi2

def fitMass():
    numpy.random.seed(1)
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
    #Halo mass
    halo= halomass(out)
    bulge= halomass(out)
    disk= diskmass(out)
    print halo, disk, bulge, totalmass(out)
    #Sample
    samples= bovy_mcmc.markovpy(out,0.05,
                                (lambda x: -obj(x,vcircdata,bp,dp,hp)),
                                (),
                                isDomainFinite=[[True,True],
                                                [True,True],
                                                [True,True],
                                                [True,True],
                                                [True,True]],
                                domain=[[0.,1.],
                                        [0.,1.],
                                        [0.,2.],
                                        [0.,10.],
                                        [0.,1.]],
                                nwalkers=10,
                                nsamples=10000)
    print "Done with sampling ..."
    print numpy.mean(numpy.array(samples),axis=0)
    print numpy.std(numpy.array(samples),axis=0)
    samples= numpy.random.permutation(samples)[0:500]
    #total
    totalmasssamples= []
    for s in samples:
        totalmasssamples.append(totalmass(s))
    totalmasssamples= numpy.array(totalmasssamples)
    print "total mass: ", numpy.mean(totalmasssamples), numpy.std(totalmasssamples)
    bovy_plot.bovy_print()
    bovy_plot.bovy_hist(totalmasssamples,bins=16,range=[0.,2.])
    bovy_plot.bovy_end_print('totalmass.png')
    #halo
    totalmasssamples= []
    for s in samples:
        totalmasssamples.append(halomass(s))
    totalmasssamples= numpy.array(totalmasssamples)
    print "halo mass: ", numpy.mean(totalmasssamples), numpy.std(totalmasssamples)
    return None
    #disk
    totalmasssamples= []
    for s in samples:
        totalmasssamples.append(diskmass(s))
    totalmasssamples= numpy.array(totalmasssamples)
    print "disk mass: ", numpy.mean(totalmasssamples), numpy.std(totalmasssamples)
    #bulge
    totalmasssamples= []
    for s in samples:
        totalmasssamples.append(bulgemass(s))
    totalmasssamples= numpy.array(totalmasssamples)
    print "bulge mass: ", numpy.mean(totalmasssamples), numpy.std(totalmasssamples)
    return None

def totalmass(out):
    return halomass(out)+bulgemass(out)+diskmass(out)

def bulgemass(out):
    bp= HernquistPotential(a=0.6/_REFR0,normalize=1.-out[0]-out[1])
    return integrate.quad((lambda x: bp.dens(x,0.)*x**2.),0.,250./_REFR0)[0]*4.*numpy.pi*_convertmass/out[2]**2.
def halomass(out):
    hp= NFWPotential(normalize=out[1],
                     a=out[3])
    return integrate.quad((lambda x: hp.dens(x,0.)*x**2.),0.,250./_REFR0)[0]*4.*numpy.pi*_convertmass/out[2]**2.

def diskmass(out):
    if _dexp:
        dp= DoubleExponentialDiskPotential(normalize=out[0],
                                           hr=out[4],
                                           hz=0.3/_REFR0)
    else:
        dp= MiyamotoNagaiPotential(normalize=out[0],
                                   a=out[4],
                                   b=0.3/_REFR0)
    return integrate.dblquad((lambda x,y: dp.dens(y,x)*y),0.,25./_REFR0,lambda x: 0.,lambda x: 10./_REFR0)[0]*2.*numpy.pi*_convertmass/out[2]**2.

if __name__ == '__main__':
    fitMass()
