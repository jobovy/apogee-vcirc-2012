###############################################################################
# fitvc.py: fit the MW rotation curve from APOGEE data
###############################################################################
import sys
import os, os.path
import cPickle as pickle
from optparse import OptionParser
import math
import numpy
from scipy import integrate, optimize
from scipy.stats import norm
from readVclosData import readVclosData
import isomodel
_PMSGRA= 30.24 #km/s/kpc
_PMSGRA_ERR= 0.12 #km/s/kpc
_VRSUN= -11.1
_DEGTORAD= math.pi/180.
_ERASESTR= "                                                                                "
def fitvc(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return
    #Check whether the savefile already exists
    if os.path.exists(args[0]):
        raise IOError("Savefile already exists, not re-fitting and overwriting")
    #Read the data
    print "Reading the data ..."
    data= readVclosData(postshutdown=options.postshutdown,
                        fehcut=options.fehcut,
                        cohort=options.cohort,
                        lmin=options.lmin,
                        bmax=options.bmax,
                        ak=True,
                        jkmax=options.jkmax)
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
    df= None
    #Initial condition for fit/sample
    init_params= _initialize_params(options)
    #Fit/sample
    if not options.mcsample:
        print "Fitting ..."
        params= optimize.fmin_powell(mloglike,init_params,
                                     args=(data['VHELIO'],
                                           l,b,jk,h,iso,df,options,
                                           sinl,cosl,cosb,cosl),
                                     callback=cb)
        print params
    #Save
    return None

def _initialize_params(options):
    if options.rotcurve.lower() == 'flat' and options.dfmodel.lower() == 'simplegaussian':
        return [235.,8.,numpy.log(35.),20.]

def cb(x): print numpy.exp(x)

def mloglike(params,vhelio,l,b,jk,h,iso,df,options,sinl,cosl,cosb,sinb):
    """minus log likelihood Eqn (1),
    params= [vc(ro),ro,sigmar,meanvt]"""
    #For each star, marginalize over distance, first divide out vo
    out= 0.
    for ii in range(len(vhelio)):
        sys.stdout.write('\r'+"Working on %i / %i ...\r" % (ii,len(vhelio)))
        sys.stdout.flush()
        thisout= integrate.quad(_mloglikedIntegrand,
                                0.,numpy.inf,
                                args=(params,vhelio[ii]/params[1],
                                      l[ii],b[ii],jk[ii],h[ii],
                                      iso,df,options,
                                      sinl[ii],cosl[ii],cosb[ii],sinb[ii]))[0]
        if thisout == 0.:
            print vhelio[ii], l[ii]/_DEGTORAD, jk[ii], h[ii]
        out+= -numpy.log(thisout)
    #BOVY:Apply Ro prior correcly
    sys.stdout.write('\r'+_ERASESTR+'\r')
    sys.stdout.flush()
    print out, params
    return out+(params[1]-8.2)**2./0.5

def _mloglikedIntegrand(d,params,vhelio,l,b,jk,h,
                        iso,df,options,sinl,cosl,cosb,sinb):
    #All positions are /ro, all velocities are /vo (d and vhelio have this already
    #Calculate coordinates, distances are /Ro (/params[1])
    R= numpy.sqrt(1.+d**2.-2.*d*cosl)
    if R == 0.:
        R+= 0.0001
        d+= 0.0001
    if isinstance(d,numpy.ndarray):
        theta= numpy.arcsin(d/R*sinl)
        indx= (1./cosl < d)*(cosl > 0.)
        theta[indx]= numpy.pi-theta[indx]
    else:
        if 1./cosl < d and cosl > 0.:
            theta= numpy.pi-numpy.arcsin(d/R*sinl)
        else:
            theta= numpy.arcsin(d/R*sinl)
    vgal= _vgal(params,vhelio,l,b,options,sinl,cosl,cosb)
    vpec= _vpec(params,vgal,R,options,l,theta)
    #Calculate probabilities
    logpvlos= _logdf(params,vpec,R,options,df,l,theta)
    logpd= _logpd(params,d,l,b,jk,h,df,iso,options,R,theta,cosb,sinb)
    return numpy.exp(logpvlos+logpd)

def _logpd(params,d,l,b,jk,h,df,iso,options,R,theta,cosb,sinb):
    dm= _dm(params,d)
    mh= h-dm
    logpiso= iso(jk,mh)
    logpddf= _logpddf(params,d,l,b,R,theta,cosb,sinb,options)
    return logpiso+logpddf

def _logpddf(params,d,l,b,R,theta,cosb,sinb,options):
    """Prior on distance distribution, including Jacobian"""
    absZ= numpy.fabs(d*sinb)
    if options.densmodel.lower() == 'expdisk':
        logdensRZ= (-(R-1.)/options.hr-absZ/options.hz)*params[1]
    return logdensRZ+2.*numpy.log(d)+numpy.log(cosb)

def _dm(params,d):
    """Distance modulus w/ d in ro"""
    return 5.*numpy.log10(d*params[1])+10.

def _logdf(params,vpec,R,options,df,l,theta):
    if options.dfmodel.lower() == 'simplegaussian':
        sinlt= numpy.sin(l+theta)
        meanvlos= params[3]*sinlt
        slos= numpy.exp(params[2]*2.)*numpy.sqrt(1.-0.5*sinlt**2.)
        t= (vpec-meanvlos)/slos
        return norm.logpdf(t)-numpy.log(slos)

def _vc(params,R,options):
    """Circular velocity at R for different models"""
    if options.rotcurve.lower() == 'flat':
        return 1. #vc/vo

def _vpec(params,vgal,R,options,l,theta):
    return vgal-_vc(params,R,options)*numpy.sin(l+theta)

def _vgal(params,vhelio,l,b,options,sinl,cosl,cosb):
    return vhelio-cosl*_VRSUN/params[0]+sinl*_PMSGRA*params[1]/params[0] #params[1]=Ro

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the fit/samples will be saved to"
    parser = OptionParser(usage=usage)
    parser.add_option("--rotcurve",dest='rotcurve',default='flat',
                      help="Rotation curve model to fit")
    parser.add_option("--dfmodel",dest='dfmodel',default='simplegaussian',
                      help="DF model to use")
    parser.add_option("--densmodel",dest='densmodel',default='expdisk',
                      help="Density model to use")
    parser.add_option("--hr",dest='hr',default=3.,type='float',
                      help="scale length in kpc")
    parser.add_option("--hz",dest='hz',default=0.3,type='float',
                      help="scale height in kpc")
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
    #Isochrone IMF
    parser.add_option("--imfmodel",dest='imfmodel',default='lognormalChabrier2001',
                      help="imfmodel for isochrone model")
    parser.add_option("--Z",dest='Z',default=.019,type='float',
                      help="Metallicity of isochrone")
    #Sample?
    parser.add_option("--mcsample",action="store_true", dest="mcsample",
                      default=False,
                      help="If set, sample around the best fit, save in args[1]")
    parser.add_option("--nsamples",dest='nsamples',default=1000,type='int',
                      help="Number of MCMC samples to obtain")
    return parser

if __name__ == '__main__':
    fitvc(get_options())
