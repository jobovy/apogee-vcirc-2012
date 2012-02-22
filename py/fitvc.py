###############################################################################
# fitvc.py: fit the MW rotation curve from APOGEE data
###############################################################################
#
# REMEMBER: Everything is wrt a _REFV0 and a _REFR0, likelihood is calculated
#           in terms of current v0, R0
#
# To Do: 
#        - Correct Ro prior
#
# BaseModel: Flat rotation curve, Ro, 
#            Gaussian velocity distribution + outlier model (frac + mean, 
#            fixed 150 km/s dispersion)
#            --dwarf adds same for dwarfs
#
###############################################################################
import sys
import os, os.path
import cPickle as pickle
from optparse import OptionParser
import math
import numpy
from scipy import integrate, optimize
from scipy.maxentropy import logsumexp
from scipy.stats import norm
import bovy_mcmc
from galpy.util import save_pickles
from readVclosData import readVclosData
import isomodel
_PMSGRA= 30.24 #km/s/kpc
_PMSGRA_ERR= 0.12 #km/s/kpc
_VRSUN= -11.1
#Reference values of parameters, such that fit should be around 1.
_REFV0= 235. #km/s
_REFR0= 8. #kpc
#Integration parameters when bintegrating
_BINTEGRATENBINS= 1001
_BINTEGRATEDMIN= 0.001 #kpc
_BINTEGRATEDMAX= 30. #kpc
_BINTEGRATEDMIN_DWARF= 0.001 #kpc
_BINTEGRATEDMAX_DWARF= .1 #kpc
_DEGTORAD= math.pi/180.
_ERASESTR= "                                                                                "
def fitvc(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    #Check whether the savefile already exists
    if os.path.exists(args[0]):
        savefile= open(args[0],'rb')
        params= pickle.load(savefile)
        savefile.close()
        if options.mcsample:
            for kk in range(len(params[0])):
                xs= numpy.array([s[kk] for s in params])
                print numpy.mean(xs), numpy.std(xs)
        else:
            print params
        print "Savefile already exists, not re-fitting and overwriting"
        return None
    #Read the data
    print "Reading the data ..."
    data= readVclosData(postshutdown=options.postshutdown,
                        fehcut=options.fehcut,
                        cohort=options.cohort,
                        lmin=options.lmin,
                        bmax=options.bmax,
                        ak=True,
                        jkmax=options.jkmax)
    if not options.location is None:
        data= data[(data['LOCATION'] == options.location)]
    #data= data[(data['J0MAG']-data['K0MAG'] > 0.82)]
    #data= data[(data['GLON'] > 200.)*(data['GLON'] < 360.)*(data['LOGG'] < 3.)]
    print "Using %i data points ..." % len(data)
    #Pre-calculate stuff
    l= data['GLON']*_DEGTORAD
    b= data['GLAT']*_DEGTORAD
    sinl= numpy.sin(l)
    cosl= numpy.cos(l)
    sinb= numpy.sin(b)
    cosb= numpy.cos(b)
    jk= data['J0MAG']-data['K0MAG']
    jk[(jk < 0.5)]= 0.5 #BOVY: FIX THIS HACK BY EMAILING GAIL
    h= data['H0MAG']
    #Set up the isochrone
    print "Setting up the isochrone model ..."
    iso= isomodel.isomodel(imfmodel=options.imfmodel,Z=options.Z)
    if options.dwarf:
        iso= [iso,
              isomodel.isomodel(imfmodel=options.imfmodel,Z=options.Z,dwarf=True)]
    else:
        iso= [iso]
    df= None
    #Initial condition for fit/sample
    init_params, isDomainFinite, domain= _initialize_params(options)
    if not options.init is None:
        #Load initial parameters from file
        savefile= open(options.init,'rb')
        init_params= pickle.load(savefile)
        savefile.close()
    #Pre-calculate p(J-K,Hs|...)[isochrone] (in fact does not depend on the parameters
    print "Pre-calculating isochrone distance prior ..."
    logpiso= numpy.zeros((len(data),_BINTEGRATENBINS))
    ds= numpy.linspace(_BINTEGRATEDMIN,_BINTEGRATEDMAX,
                       _BINTEGRATENBINS)
    dm= _dm(ds)
    for ii in range(len(data)):
        mh= h[ii]-dm
        logpiso[ii,:]= iso[0](numpy.zeros(_BINTEGRATENBINS)+jk[ii],mh)
    if options.dwarf:
        logpisodwarf= numpy.zeros((len(data),_BINTEGRATENBINS))
        dwarfds= numpy.linspace(_BINTEGRATEDMIN_DWARF,_BINTEGRATEDMAX_DWARF,
                                    _BINTEGRATENBINS)
        dm= _dm(dwarfds)
        for ii in range(len(data)):
            mh= h[ii]-dm
            logpisodwarf[ii,:]= iso[1](numpy.zeros(_BINTEGRATENBINS)+jk[ii],mh)
    else:
        logpisodwarf= None
    #Fit/sample
    if not options.mcsample:
        print "Fitting ..."
        params= optimize.fmin_powell(mloglike,init_params,
                                     args=(data['VHELIO'],
                                           l,b,jk,h,df,options,
                                           sinl,cosl,cosb,sinb,
                                           logpiso,logpisodwarf),
                                     callback=cb)
        print params
        save_pickles(args[0],params)       
    else:
        samples= bovy_mcmc.markovpy(init_params,
                                    0.01,
                                    loglike,
                                    (data['VHELIO'],
                                     l,b,jk,h,df,options,
                                     sinl,cosl,cosb,sinb,
                                     logpiso,logpisodwarf),
                                    isDomainFinite=isDomainFinite,
                                    domain=domain,
                                    nsamples=options.nsamples,
                                    nwalkers=4*len(init_params),
                                    threads=options.multi)
        for kk in range(len(init_params)):
            xs= numpy.array([s[kk] for s in samples])
            print numpy.mean(xs), numpy.std(xs)
        save_pickles(args[0],samples)
    return None

def _initialize_params(options):
    if options.rotcurve.lower() == 'flat' and options.dfmodel.lower() == 'simplegaussian':
        if options.dwarf:
            return ([235./_REFV0,8./_REFR0,numpy.log(35./_REFV0),0.1,0.,0.2],
                    [[True,False],[True,True],[False,False],
                    [True,True],[False,False],[True,True]],
                    [[0.,0.],[5./_REFR0,11./_REFR0],
                     [0.,0.],[0.,1.],[0.,0.],
                     [0.,1.]])
        else:
            return ([235./_REFV0,8./_REFR0,numpy.log(35./_REFV0),0.1,0.],
                    [[True,False],[True,True],[False,False],
                     [True,True],[False,False]],
                    [[0.,0.],[5./_REFR0,11./_REFR0],
                     [0.,0.],[0.,1.],[0.,0.]])

def cb(x): print x

def loglike(params,vhelio,l,b,jk,h,df,options,sinl,cosl,cosb,sinb,
            logpiso,logpisodwarf):
    return -mloglike(params,vhelio,l,b,jk,h,df,options,sinl,cosl,cosb,sinb,
                     logpiso,logpisodwarf)
def mloglike(params,vhelio,l,b,jk,h,df,options,sinl,cosl,cosb,sinb,
             logpiso,logpisodwarf):
    """minus log likelihood Eqn (1),
    params= [vc(ro),ro,sigmar]"""
    #Boundaries
    if params[0] <= 0.: #Don't allow counter-rotation
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    elif params[1] <= 5./_REFR0 or params[1] > 11./_REFR0: #So it does not go nuts
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if params[3] < 0. or params[3] > 1.:
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if options.dwarf and (params[5] < 0. or params[5] > 1.):
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    #For each star, marginalize over distance
    if options.dontbintegrate:
        out= 0.
        for ii in range(len(vhelio)):
            sys.stdout.write('\r'+"Working on %i / %i ...\r" % (ii,len(vhelio)))
            sys.stdout.flush()
            raise NotImplementedError("Direct scipy.integrate integration does not work bc of returning the p(d)")
            thisout, thispd= integrate.quad(_mloglikedIntegrand,
                                            0.,numpy.inf,
                                            args=(params,vhelio[ii]/params[0]/_REFV0,
                                                  l[ii],b[ii],jk[ii],h[ii],
                                                  iso[0],df,options,
                                                  sinl[ii],cosl[ii],cosb[ii],sinb[ii],
                                                  False))[0]
            if options.dwarf:
                thisout= (1.-params[5])*thisout\
                    +params[5]*integrate.quad(_mloglikedIntegrand,
                                              0.,numpy.inf,
                                              args=(params,vhelio[ii]/params[0]/_REFV0,
                                                    l[ii],b[ii],jk[ii],h[ii],
                                                    iso[1],df,options,
                                                    sinl[ii],cosl[ii],cosb[ii],sinb[ii],
                                                    False))[0]
            if thisout == 0.:
                print vhelio[ii], l[ii]/_DEGTORAD, jk[ii], h[ii]
            out+= -numpy.log(thisout)
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
        print out, params
    else: #Calculate the integral by binning
        thisout= numpy.zeros((len(vhelio),_BINTEGRATENBINS))
        logpd= numpy.zeros((len(vhelio),_BINTEGRATENBINS))
        ds= numpy.linspace(_BINTEGRATEDMIN,_BINTEGRATEDMAX,
                           _BINTEGRATENBINS)/params[1]/_REFR0
        if options.dwarf:
            thisxtraout= numpy.zeros((len(vhelio),_BINTEGRATENBINS))
            logpddwarf= numpy.zeros((len(vhelio),_BINTEGRATENBINS))
            dwarfds= numpy.linspace(_BINTEGRATEDMIN_DWARF,_BINTEGRATEDMAX_DWARF,
                                    _BINTEGRATENBINS)/params[1]/_REFR0
        for ii in range(_BINTEGRATENBINS):
            thisout[:,ii], logpd[:,ii]= _mloglikedIntegrand(ds[ii],
                                                            params,vhelio/params[0]/_REFV0,
                                                            l,b,jk,h,df,options,
                                                            sinl,cosl,cosb,sinb,
                                                            True,logpiso[:,ii])
            if options.dwarf:
                thisxtraout[:,ii], logpddwarf[:,ii]= _mloglikedIntegrand(dwarfds[ii],
                                                                         params,vhelio/params[0]/_REFV0,
                                                                         l,b,jk,h,df,options,
                                                                         sinl,cosl,cosb,sinb,
                                                                         True,logpisodwarf[:,ii])
        #Sum each one
        out= 0.
        for ii in range(len(vhelio)):
            if options.dwarf: #Combine
                tmpthisout= logsumexp(thisout[ii,:])+numpy.log(1.-params[5])-logsumexp(logpd[ii,:])
                thislogpdwarf= logsumexp(logpddwarf[ii,:])
                if thislogpdwarf == -numpy.finfo(numpy.dtype(numpy.float64)).max:
                    out-= tmpthisout
                else:
                    tmpthisxtraout= logsumexp(thisxtraout[ii,:])+numpy.log(params[5])-thislogpdwarf
                    #print tmpthisxtraout, -logsumexp(logpddwarf[ii,:])
                    c= numpy.amax([tmpthisout,tmpthisxtraout])
                    out-= numpy.log(numpy.exp(tmpthisout-c)+numpy.exp(tmpthisxtraout-c))+c
            else:
                out+= -logsumexp(thisout[ii,:])+logsumexp(logpd[ii,:]) #so we can compare to the +dwarf case
        print out, params
    #BOVY:Apply Ro prior correcly
    return out+(params[1]*_REFR0-8.2)**2./0.5#+(params[0]*_REFV0+2.25*math.exp(2.*params[2])*_REFV0/params[0]+12.24-_PMSGRA*params[1]*_REFR0)**2./200. #params[1]=Ro, SBD10 Solar motion

def _mloglikedIntegrand(d,params,vhelio,l,b,jk,h,
                        df,options,sinl,cosl,cosb,sinb,returnlog,logpiso):
    #All positions are /ro, all velocities are /vo (d and vhelio have this already
    #Calculate coordinates, distances are /Ro (/params[1])
    #
    #Returns the loglike and the p(d) factor
    R= numpy.sqrt(1.+d**2.-2.*d*cosl)
    if isinstance(vhelio,numpy.ndarray):
        indx= (R == 0.)
        R[indx]+= 0.0001
    else:
        if R == 0.:
            R+= 0.0001
            d+= 0.0001
    if isinstance(vhelio,numpy.ndarray):
        theta= numpy.arcsin(d/R*sinl)
        indx= (1./cosl < d)*(cosl > 0.)
        theta[indx]= numpy.pi-theta[indx]
    else:
        if 1./cosl < d and cosl > 0.:
            theta= numpy.pi-numpy.arcsin(d/R*sinl)
        else:
            theta= numpy.arcsin(d/R*sinl)
    vgal= _vgal(params,vhelio,l,b,options,sinl,cosl)
    vpec= _vpec(params,vgal,R,options,l,theta)
    #Calculate probabilities
    logpvlos= _logdf(params,vpec,R,options,df,l,theta)+numpy.log(1.-params[3])
    logpvlos_outlier= _logoutlierdf(params,vgal)+numpy.log(params[3])
    c= numpy.amax(numpy.array([logpvlos,logpvlos_outlier]),axis=0)
    logpvlos= numpy.log(numpy.exp(logpvlos-c)
                        +numpy.exp(logpvlos_outlier-c))+c
    logpd= _logpd(params,d,l,b,jk,h,df,options,R,theta,cosb,sinb,logpiso)
    if returnlog: return (logpvlos+logpd,logpd)
    else: return (numpy.exp(logpvlos+logpd),numpy.exp(logpd))

def _logpd(params,d,l,b,jk,h,df,options,R,theta,cosb,sinb,logpiso):
    logpddf= _logpddf(params,d,l,b,R,theta,cosb,sinb,options)
    return logpiso+logpddf

def _logpddf(params,d,l,b,R,theta,cosb,sinb,options):
    """Prior on distance distribution, including Jacobian"""
    absZ= numpy.fabs(d*sinb)
    if options.densmodel.lower() == 'expdisk':
        logdensRZ= (-(R-1.)/options.hr-absZ/options.hz)*params[1]*_REFR0
    return logdensRZ+2.*numpy.log(d*params[1])+numpy.log(cosb)

def _dm(d):
    """Distance modulus w/ d in kpc"""
    return 5.*numpy.log10(d)+10.

def _logdf(params,vpec,R,options,df,l,theta):
    if options.dfmodel.lower() == 'simplegaussian':
        sinlt= numpy.sin(l+theta)
        slos= numpy.exp(params[2])/params[0]\
            *numpy.sqrt(1.-0.5*sinlt**2.)*numpy.exp(-(R-1.)/options.hs*params[1]*_REFR0)
        t= vpec/slos
        return norm.logpdf(t)-numpy.log(slos*params[0]*_REFV0)

def _logoutlierdf(params,vgal):
    t= (vgal*params[0]-params[4])*_REFV0/200.
    return norm.logpdf(t)-numpy.log(200.)

def _vc(params,R,options):
    """Circular velocity at R for different models"""
    if options.rotcurve.lower() == 'flat':
        return 1. #vc/vo

def _vpec(params,vgal,R,options,l,theta):
    return vgal-_vc(params,R,options)*numpy.sin(l+theta)

def _vgal(params,vhelio,l,b,options,sinl,cosl):
    return vhelio-cosl*_VRSUN/params[0]/_REFV0+sinl*_PMSGRA*params[1]*_REFR0/params[0]/_REFV0 #params[1]=Ro

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the fit/samples will be saved to"
    parser = OptionParser(usage=usage)
    #Initial conditions file
    parser.add_option("--init",dest='init',default=None,
                      help="Initial parameters")
    #Rotation curve parameters/model
    parser.add_option("--rotcurve",dest='rotcurve',default='flat',
                      help="Rotation curve model to fit")
    #Velocity distribution model
    parser.add_option("--dfmodel",dest='dfmodel',default='simplegaussian',
                      help="DF model to use")
    #Density model
    parser.add_option("--densmodel",dest='densmodel',default='expdisk',
                      help="Density model to use")
    parser.add_option("--hr",dest='hr',default=3.,type='float',
                      help="scale length in kpc")
    parser.add_option("--hz",dest='hz',default=0.3,type='float',
                      help="scale height in kpc")
    parser.add_option("--hs",dest='hs',default=8.,type='float',
                      help="dispersion scale length in kpc")
    #Data options
    parser.add_option("--lmin",dest='lmin',default=35.,type='float',
                      help="readVclosData 'lmin'")
    parser.add_option("--bmax",dest='bmax',default=2.,type='float',
                      help="readVclosData 'bmax'")
    parser.add_option("--allshutdown",action="store_false", 
                      dest="postshutdown",
                      default=True,
                      help="setting this sets postshutdown to False in data")
    parser.add_option("--fehcut",action="store_true", 
                      dest="fehcut",
                      default=False,
                      help="readVclosData 'fehcut'")
    parser.add_option("--cohort",dest='cohort',default=None,
                      help="readVclosData 'cohort'")
    parser.add_option("--jkmax",dest='jkmax',default=1.2,type='float',
                      help="readVclosData 'jkmax'")
    parser.add_option("--location",dest='location',default=None,type='int',
                      help="location id when looking at single los")
    #Isochrone IMF
    parser.add_option("--imfmodel",dest='imfmodel',default='lognormalChabrier2001',
                      help="imfmodel for isochrone model")
    parser.add_option("--Z",dest='Z',default=.019,type='float',
                      help="Metallicity of isochrone")
    #Add dwarf part?
    parser.add_option("--dwarf",action="store_true", 
                      dest="dwarf",
                      default=False,
                      help="setting this adds dwarf contamination")
    #Distance marginalization by binning or not
    parser.add_option("--dontbintegrate",action="store_true",
                      dest="dontbintegrate",
                      default=False,
                      help="If set, *don't* calculate the integral by binning in d")
    #Sample?
    parser.add_option("--mcsample",action="store_true", dest="mcsample",
                      default=False,
                      help="If set, sample around the best fit, save in args[1]")
    parser.add_option("--nsamples",dest='nsamples',default=1000,type='int',
                      help="Number of MCMC samples to obtain")
    #Multiprocessing?
    parser.add_option("-m","--multi",dest='multi',default=1,type='int',
                      help="number of cpus to use for sampling")
    return parser

if __name__ == '__main__':
    fitvc(get_options())
