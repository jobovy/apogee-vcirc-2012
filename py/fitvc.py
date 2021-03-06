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
import copy
from optparse import OptionParser
import math
import numpy
from scipy import integrate, optimize, interpolate
from scipy.maxentropy import logsumexp
from scipy.stats import norm
import multi
import multiprocessing
import acor
import bovy_mcmc, bovy_mcmc.elliptical_slice
import flexgp.fast_cholesky
from galpy.util import save_pickles
from galpy.df import dehnendf
from readVclosData import readVclosData
import isomodel
import isodist
import asymmetricDriftModel
import skewnormal
_PMSGRA= 30.24 #km/s/kpc
_PMSGRA_ERR= 0.12 #km/s/kpc
_VRSUN= -11.1 #SBD10
_VTSUN= 12.24 #SBD10
#Reference values of parameters, such that fit should be around 1.
_REFV0= 235. #km/s
_REFR0= 8. #kpc
#Integration parameters when bintegrating
_BINTEGRATENBINS= 501
_BINTEGRATEVNBINS= 1001
_BINTEGRATEDMIN= 0.001 #kpc
_BINTEGRATEDMAX= 30. #kpc
_BINTEGRATEDMIN_DWARF= 0.001 #kpc
_BINTEGRATEDMAX_DWARF= .1 #kpc
_NMULTIPLEPOPS= 20
_TAU1= 0.1
_TAUM= 10.
_TAUBETA= 0.38
_TAU0= 8.
_DEGTORAD= math.pi/180.
_ERASESTR= "                                                                                "
def dummyIso(jk,h):
    return 0. #Log
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
            print_samples_qa(params)
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
                        cutmultiples=options.cutmultiples,
                        correctak=options.correctak,
                        loggcut=options.loggcut,
                        validfeh=options.indivfeh, #if indivfeh, we need validfeh
                        jkmax=options.jkmax,
                        datafilename=options.fakedata)
    if options.justinner:
        data= data[(data['GLON'] < 97.)]
    elif options.justouter:
        data= data[(data['GLON'] >= 97.)]
    if not options.location is None:
        data= data[(data['LOCATION'] == options.location)]
    if not options.removelocation is None:
        data= data[(data['LOCATION'] != options.removelocation)]
    if not options.downsample is None:
        indx= numpy.random.permutation(len(data))[0:int(math.floor(len(data)/options.downsample))]
        data= data[indx]
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
    if options.basti: jk= 1./0.961*jk-0.047/0.961 #Carpenter 2001
    h= data['H0MAG']
    #Set up the isochrone
    if not options.isofile is None and os.path.exists(options.isofile):
        print "Loading the isochrone model ..."
        isofile= open(options.isofile,'rb')
        iso= pickle.load(isofile)
        if options.indivfeh:
            zs= pickle.load(isofile)
        elif options.varfeh:
            locl= pickle.load(isofile)
        isofile.close()
    else:
        print "Setting up the isochrone model ..."
        if options.nods:
            iso= dummyIso
            global _BINTEGRATEDMAX
            _BINTEGRATEDMAX= 10. #kpc
        elif options.indivfeh:
            #Load all isochrones
            iso= []
            if options.basti:
                zs= numpy.array([0.0001,0.0003,0.0006,0.001,0.002,0.004,
                                 0.008,0.01,0.0198,0.03,0.04])
            else:
                zs= numpy.arange(0.0005,0.03005,0.0005)
            for ii in range(len(zs)):
                iso.append(isomodel.isomodel(imfmodel=options.imfmodel,
                                             expsfh=options.expsfh,
                                             Z=zs[ii],
                                             basti=options.basti,
                                             loggmax=options.loggcut))
        elif options.varfeh:
            locs= list(set(data['LOCATION']))
            iso= []
            for ii in range(len(locs)):
                indx= (data['LOCATION'] == locs[ii])
                locl= numpy.mean(data['GLON'][indx]*_DEGTORAD)
                iso.append(isomodel.isomodel(imfmodel=options.imfmodel,
                                             expsfh=options.expsfh,
                                             marginalizefeh=True,
                                             loggmax=options.loggcut,
                                             glon=locl))
        else:
            iso= isomodel.isomodel(imfmodel=options.imfmodel,Z=options.Z,
                                   expsfh=options.expsfh,
                                   loggmax=options.loggcut)
        if options.dwarf:
            iso= [iso,
                  isomodel.isomodel(imfmodel=options.imfmodel,Z=options.Z,
                                    dwarf=True,expsfh=options.expsfh,
                                    loggmax=options.loggcut)]
        else:
            iso= [iso]
        if not options.isofile is None:
            isofile= open(options.isofile,'wb')
            pickle.dump(iso,isofile)
            if options.indivfeh:
                pickle.dump(zs,isofile)
            elif options.varfeh:
                pickle.dump(locl,isofile)
            isofile.close()
    df= None
    #Initial condition for fit/sample
    init_params, isDomainFinite, domain= _initialize_params(options)
    if not options.init is None:
        #Load initial parameters from file
        savefile= open(options.init,'rb')
        init_params= pickle.load(savefile)
        savefile.close()
    if options.sbdvpec:
        init_params= [ 0.96025887 , 1.00661491, -1.9503565 ,  0.01060685,  0.54155172 , 0.09690153,  0.32055182,0.22680105]
    #Pre-calculate p(J-K,Hs|...)[isochrone] (in fact does not depend on the parameters if dm or ah are not fitted
    print "Pre-calculating isochrone distance prior ..."
    logpiso= numpy.zeros((len(data),_BINTEGRATENBINS))
    ds= numpy.linspace(_BINTEGRATEDMIN,_BINTEGRATEDMAX,
                       _BINTEGRATENBINS)
    dm= _dm(ds)
    for ii in range(len(data)):
        mh= h[ii]-dm
        if options.indivfeh:
            #Find closest Z
            thisZ= isodist.FEH2Z(data[ii]['FEH'])
            indx= numpy.argmin(numpy.fabs(thisZ-zs))
            logpiso[ii,:]= iso[0][indx](numpy.zeros(_BINTEGRATENBINS)+jk[ii],mh)
        elif options.varfeh:
            #Find correct iso
            indx= (locl == data[ii]['LOCATION'])
            logpiso[ii,:]= iso[0][indx](numpy.zeros(_BINTEGRATENBINS)+jk[ii],mh)
        else:
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
                                           logpiso,logpisodwarf,False,None,
                                           iso,data['FEH']),
                                     callback=cb)
        print params
        save_pickles(args[0],params)       
    elif options.rotcurve.lower() == 'gp':
        if options.gpemcee:
            #Set up GP
            gprs= numpy.linspace(options.gprmin,options.gprmax,options.gpnr)
            if options.gpfixtau is None:
                init_hyper= [numpy.log(20./_REFV0),numpy.log(1./_REFR0)]
            else:
                init_hyper= [numpy.log(20./_REFV0)]
            init_f= (gprs/_REFR0)**-0.1-1.
            #Add the hyper-parameters and the GP parameters
            init_params= list(init_params)
            init_params.extend(init_hyper)
            init_params.extend(init_f)
            init_params= numpy.array(init_params)
            if options.gpfixtau is None:
                isDomainFinite.extend([[True,True],[True,True]])
                domain.extend([[-4.,2.],[-4.,2.5]])
            else:
                isDomainFinite.extend([[True,True]])
                domain.extend([[-4.,2.]])
            for ii in range(options.gpnr):
                isDomainFinite.append([False,False])
                domain.append([0.,0.])
            #Now sample
            samples= bovy_mcmc.markovpy(init_params,
                                        0.05,
                                        loglike_gpemcee,
                                        (data['VHELIO'],
                                         l,b,jk,h,df,options,
                                         sinl,cosl,cosb,sinb,
                                         logpiso,logpisodwarf,gprs),
                                        isDomainFinite=isDomainFinite,
                                        domain=domain,
                                        nsamples=options.nsamples,
                                        nwalkers=2*len(init_params),
                                        threads=options.multi)
            print_samples_qa(samples)
            save_pickles(args[0],samples)           
        else:
            #In this case we only sample, and we do it in a special way
            #Set up GP
            gprs= numpy.linspace(options.gprmin,options.gprmax,options.gpnr)
            if options.gpfixtau is None:
                init_hyper= [numpy.log(20./_REFV0),numpy.log(1./_REFR0)]
            else:
                init_hyper= [numpy.log(20./_REFV0)]
            #cov= _calc_covar(gprs,init_hyper,options)
            #chol= flexgp.fast_cholesky.fast_cholesky(cov)[0]
            #init_f= numpy.dot(chol,numpy.random.normal(size=options.gpnr))
            init_f= (gprs/_REFR0)**-0.1-1.
            #Slice steps
            if options.dwarf:
                raise NotImplementedError("'-dwarf' w/ rotcurve=GP not implemented")
            else:
                one_step= [0.03,0.03,0.03,0.02,0.05]
            if options.fitvpec:
                one_step.extend([0.03,0.03])
            if options.gpfixtau is None:
                hyper_step= [0.05,0.05]
            else:
                hyper_step= 0.05
            #Set up samples
            samples= []
            f_samples= []
            hyper_samples= []
            lnps= []
            these_params= numpy.array(init_params)
            these_f= init_f
            these_hyper= numpy.array(init_hyper)
            totalnsamples= int(math.ceil(1.15*options.nsamples))
            for ii in range(totalnsamples):
                print ii, totalnsamples
                #Sample in three times
                #Slice sample the hyper-parameters
                if options.gpfixtau is None:
                    new_hyper= bovy_mcmc.slice(these_hyper,
                                               hyper_step,
                                               hyperloglike,
                                               (these_f,gprs,options),
                                               isDomainFinite=[[True,True],
                                                               [True,True]],
                                               domain=[[-4.,2.],[-4.,2.5]],
                                               nsamples=1)
                else:
                    new_hyper= bovy_mcmc.slice(these_hyper,
                                               hyper_step,
                                               hyperloglike_fixtau,
                                               (these_f,gprs,options,
                                                numpy.log(options.gpfixtau/_REFR0)),
                                               isDomainFinite=[True,True],
                                               domain=[-4.,2.],
                                               nsamples=1)
                hyper_samples.append(new_hyper)
                these_hyper= copy.copy(new_hyper)
                print these_hyper
                if ii > 0.0*options.nsamples: #Let the GP burn-in on its own NOT SET
                    #regular parameters
                    #Spline interpolation of these_f
                    vcf= interpolate.InterpolatedUnivariateSpline(gprs,these_f)
                    thesesamples= bovy_mcmc.slice(these_params,
                                                  one_step,
                                                  loglike,
                                                  (data['VHELIO'],
                                                   l,b,jk,h,df,options,
                                                   sinl,cosl,cosb,sinb,
                                                   logpiso,logpisodwarf,vcf),
                                                  isDomainFinite=isDomainFinite,
                                                  domain=domain,
                                                  nsamples=1)
                    samples.append(thesesamples)
                    these_params= copy.copy(thesesamples)
                #Elliptical slice sample the node-values, first calculate the cov
                if options.gpfixtau is None:
                    cov= _calc_covar(gprs,these_hyper,options)
                else:
                    cov= _calc_covar(gprs,[these_hyper,numpy.log(options.gpfixtau/_REFR0)],
                                     options)
                chol= flexgp.fast_cholesky.fast_cholesky(cov)[0]
                new_f, lnp= bovy_mcmc.elliptical_slice.elliptical_slice(these_f,
                                                                        chol,
                                                                        floglike,
                                                                        pdf_params=(gprs,these_params,data['VHELIO'],
                                                                                    l,b,jk,h,df,options,
                                                                                    sinl,cosl,cosb,sinb,
                                                                                    logpiso,logpisodwarf),)
                f_samples.append(new_f)
                these_f= copy.copy(new_f)
                lnps.append(lnp)
                print numpy.mean(these_f), these_f
            #Save
            samples= samples[len(samples)-options.nsamples:len(samples)]
            f_samples= f_samples[len(f_samples)-options.nsamples:len(f_samples)]
            hyper_samples= hyper_samples[len(hyper_samples)-options.nsamples:len(hyper_samples)]
            save_pickles(args[0],samples,f_samples,hyper_samples,lnps)
            print_samples_qa(f_samples)
            print_samples_qa(samples)
            print_samples_qa(hyper_samples)
    else:
        samples= bovy_mcmc.markovpy(init_params,
                                    0.01,
                                    loglike,
                                    (data['VHELIO'],
                                     l,b,jk,h,df,options,
                                     sinl,cosl,cosb,sinb,
                                     logpiso,logpisodwarf,None,iso,
                                     data['FEH']),
                                    isDomainFinite=isDomainFinite,
                                    domain=domain,
                                    nsamples=options.nsamples,
                                    nwalkers=4*len(init_params),
                                    threads=options.multi)
        print_samples_qa(samples)
        save_pickles(args[0],samples)
    return None

def loglike_gpemcee(params,vhelio,l,b,jk,h,df,options,sinl,cosl,cosb,sinb,
                    logpiso,logpisodwarf,gprs):
    #Break up params into 1) basic params, 2) hyper-params, 3) f
    f= params[len(params)-options.gpnr:len(params)]
    if options.gpfixtau is None:
        hyper_params= params[len(params)-options.gpnr-2:len(params)-options.gpnr]
        hyperlnl= hyperloglike(hyper_params,f,gprs,options)
        basicparams= params[0:len(params)-options.gpnr-2]
    else:
        hyper_params= params[len(params)-options.gpnr-1:len(params)-options.gpnr]
        hyperlnl= hyperloglike_fixtau(hyper_params,f,gprs,options,
                                      numpy.log(options.gpfixtau/_REFR0))
        basicparams= params[0:len(params)-options.gpnr-1]
    print hyper_params
    print numpy.mean(f), f
    return floglike(f,gprs,
                    basicparams,vhelio,l,b,jk,h,df,options,sinl,cosl,cosb,sinb,
                    logpiso,logpisodwarf)+hyperlnl

def hyperloglike(hyper_params,f,gprs,options):
    #Just the Gaussian
    #print hyper_params
    cov= _calc_covar(gprs,hyper_params,options)
    covinv, detcov= flexgp.fast_cholesky.fast_cholesky_invert(cov,logdet=True)
    return -0.5*numpy.dot(f,numpy.dot(covinv,f))-0.5*detcov

def hyperloglike_fixtau(hyper_params,f,gprs,options,logtau):
    #Just the Gaussian
    #print hyper_params
    cov= _calc_covar(gprs,[hyper_params,logtau],options)
    covinv, detcov= flexgp.fast_cholesky.fast_cholesky_invert(cov,logdet=True)
    return -0.5*numpy.dot(f,numpy.dot(covinv,f))-0.5*detcov

def floglike(f,gprs,params,vhelio,l,b,jk,h,df,options,sinl,cosl,cosb,sinb,
            logpiso,logpisodwarf):
    vcf= interpolate.InterpolatedUnivariateSpline(gprs,f)
    return loglike(params,vhelio,l,b,jk,h,df,options,sinl,cosl,cosb,sinb,
                   logpiso,logpisodwarf,vcf)

def _calc_covar(rs,hyper_params,options):
    """GP covariance function"""
    l2= numpy.exp(2.*hyper_params[1])*_REFR0**2.
    out= numpy.zeros((len(rs),len(rs)))
    for ii in range(len(rs)):
        out[ii,:]= numpy.exp(-0.5*(rs-rs[ii])**2./l2)
    return numpy.exp(2.*hyper_params[0])*out

def _initialize_params(options):
    init_params= []
    isDomainFinite= []
    domain= []
    if options.fixvo is None:
        init_params.append(235./_REFV0)
        isDomainFinite.append([True,False])
        domain.append([0.,0.])
    if options.fixro is None:
        init_params.append(8./_REFR0)
        isDomainFinite.append([True,True])
        domain.append([5./_REFR0,11./_REFR0])
    init_params.append(numpy.log(35./_REFV0))
    init_params.append(0.1)
    isDomainFinite.append([False,False])
    isDomainFinite.append([True,True])
    domain.append([0.,0.])
    domain.append([0.,1.])
    if not options.nooutliermean:
        init_params.append(0.1)
        isDomainFinite.append([False,False])
        domain.append([0.,0.])      
    if options.dwarf:
        init_params.append(0.1)
        isDomainFinite.append([True,True])
        domain.append([0.,1.])
    if options.rotcurve.lower() == 'flatplusuniform':
        init_params.append(0.05)
        isDomainFinite.append([False,False])
        domain.append([0.,0.])
    if options.rotcurve.lower() == 'linear':
        init_params.append(0.)
        isDomainFinite.append([False,False])
        domain.append([0.,0.])
    if options.rotcurve.lower() == 'powerlaw':
        init_params.append(0.)
        isDomainFinite.append([False,False])
        domain.append([0.,0.])
    if options.rotcurve.lower() == 'quadratic':
        init_params.append(0.)
        init_params.append(0.)
        isDomainFinite.append([False,False])
        isDomainFinite.append([False,False])
        domain.append([0.,0.])
        domain.append([0.,0.])
    if options.rotcurve.lower() == 'cubic':
        init_params.append(0.)
        init_params.append(0.)
        init_params.append(0.)
        isDomainFinite.append([False,False])
        isDomainFinite.append([False,False])
        isDomainFinite.append([False,False])
        domain.append([0.,0.])
        domain.append([0.,0.])
        domain.append([0.,0.])
    if options.fitvpec:
        init_params.append(1.)
        init_params.append(1.)
        isDomainFinite.append([False,False])
        isDomainFinite.append([False,False])
        domain.append([0.,0.])
        domain.append([0.,0.])
    if options.fitsratio:
        init_params.append(0.5)
        isDomainFinite.append([True,False])
        domain.append([0.,0.])
    if options.fitsratioinnerouter:
        init_params.append(0.5)
        init_params.append(0.5)
        isDomainFinite.append([True,False])
        isDomainFinite.append([True,False])
        domain.append([0.,0.])
        domain.append([0.,0.])
    if options.fitdm:
        init_params.append(-0.1)
        isDomainFinite.append([False,False])
        domain.append([0.,0.])
    if options.fitah:
        init_params.append(-0.05)
        isDomainFinite.append([False,False])
        domain.append([0.,0.])
    if options.fitfeh:
        init_params.append(-0.1)
        isDomainFinite.append([False,False])
        domain.append([0.,0.])
    if options.fiths:
        init_params.append(1.)
        isDomainFinite.append([False,False])
        domain.append([0.,0.])
    if options.fitm2:
        init_params.append(0.)
        init_params.append(0.)
        isDomainFinite.append([True,False])
        isDomainFinite.append([True,True])
        domain.append([0.,0.])
        domain.append([-0.5,.5])
    if options.fitsrinnerouter:
        init_params.append(0.)
        isDomainFinite.append([False,False])
        domain.append([0.,0.])
    if options.dwarfinnerouter:
        init_params.append(0.1)
        isDomainFinite.append([True,True])
        domain.append([0.,1.])
    if options.fitdminnerouter:
        init_params.append(-0.1)
        isDomainFinite.append([False,False])
        domain.append([0.,0.])
    if options.fitahinnerouter:
        init_params.append(-0.05)
        isDomainFinite.append([False,False])
        domain.append([0.,0.])       
    if options.fitfehinnerouter:
        init_params.append(-0.1)
        isDomainFinite.append([False,False])
        domain.append([0.,0.])       
    if options.dfmodel.lower() == 'twopops':
        init_params.append(0.5)
        init_params.append(numpy.log(35./_REFV0)+1.)
        isDomainFinite.append([True,True])
        isDomainFinite.append([False,False])
        domain.append([0.,1.])       
        domain.append([0.,0.])
    if options.fitdl:
        init_params.append(1.*_DEGTORAD)
        isDomainFinite.append([False,False])
        domain.append([0.,0.])
    return (init_params,isDomainFinite,domain)

def cb(x): print x

def loglike(params,vhelio,l,b,jk,h,df,options,sinl,cosl,cosb,sinb,
            logpiso,logpisodwarf,vcf,iso,feh):
    return -mloglike(params,vhelio,l,b,jk,h,df,options,sinl,cosl,cosb,sinb,
                     logpiso,logpisodwarf,False,vcf,iso,feh)
def mloglike(params,vhelio,l,b,jk,h,df,options,sinl,cosl,cosb,sinb,
             logpiso,logpisodwarf,arrout,vcf,iso,feh):
    """minus log likelihood Eqn (1),
    params= [vc(ro),ro,sigmar]"""
    #If vo is fixed, insert it before proceeding
    if not options.fixvo is None:
        params= list(copy.copy(params))
        params.insert(0,options.fixvo/_REFV0)
        params= numpy.array(params)
    #If vo is fixed, insert it before proceeding
    if not options.fixro is None:
        params= list(copy.copy(params))
        params.insert(1,options.fixro/_REFR0)
        params= numpy.array(params)
    if options.fitdl:
        l= copy.copy(l)
        sinl= copy.copy(sinl)
        cosl= copy.copy(cosl)
        l-= params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitdm+options.fitah+options.fitfeh]
        sinl= numpy.sin(l)
        cosl= numpy.cos(l)
    #Boundaries
    if params[0] <= 0.: #Don't allow counter-rotation
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    elif params[1] <= 5./_REFR0 or params[1] > 11./_REFR0: #So it does not go nuts
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if params[3] < 0. or params[3] > 1.:
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if options.dwarf and (params[5-options.nooutliermean] < 0. or params[5-options.nooutliermean] > 1.):
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if options.fitsratio and params[5+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf-options.nooutliermean] < 0.:
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if options.fitsratioinnerouter and (params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf] < 0. or params[6-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf] < 0.):
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if options.fitdm and (params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] > 1. or params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] < -1.):
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if options.fitah and (params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] < -0.2 or params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] > 0.4):
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if options.fitfeh and (params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] > 1. or params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] < -1.):
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if options.fitdminnerouter and (params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitah+options.fitdm+options.fitfeh+options.fiths] > 1. or params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitah+options.fitdm+options.fitfeh+options.fiths] < -1.):
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if options.fitahinnerouter and (params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter] < -0.2 or params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter] > 0.4):
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if options.fitfehinnerouter and (params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitah+options.fitdm+options.fitfeh+options.fiths] > 1. or params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitah+options.fitdm+options.fitfeh+options.fiths] < -1.):
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if options.fitm2 and (params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fiths] < 0. or params[6-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fiths] < -0.5 or params[6-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fiths] > 0.5):
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if options.dwarfinnerouter and (params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fiths+options.fitsrinnerouter] < 0. or params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fiths+options.fitsrinnerouter] > 1.):
        return numpy.finfo(numpy.dtype(numpy.float64)).max
    if options.dfmodel.lower() == 'twopops' and (params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter] <= 0. or params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter] >= 1.):
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
                                                  sinl[ii],cosl[ii],cosb[ii],sinb[ii],iso,
                                                  False))[0]
            if options.dwarf:
                thisout= (1.-params[5-options.nooutliermean])*thisout\
                    +params[5-options.nooutliermean]*integrate.quad(_mloglikedIntegrand,
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
        #If we are fitting a dehnendf, set it up
        if options.dfmodel.lower() == 'dehnen':
            #BOVY: Only works for flat rotation curve currently
            if options.fiths: thishs= options.hs/params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter]
            else: thishs= options.hs
            df= dehnendf(correct=False,profileParams=(options.hr/_REFR0/params[1],
                                                      thishs/_REFR0/params[1],
                                                      numpy.exp(params[2])/params[0]))
        thisout= numpy.zeros((len(vhelio),_BINTEGRATENBINS))
        logpd= numpy.zeros((len(vhelio),_BINTEGRATENBINS))
        ds= numpy.linspace(_BINTEGRATEDMIN,_BINTEGRATEDMAX,
                           _BINTEGRATENBINS)/params[1]/_REFR0
        if options.dwarf:
            thisxtraout= numpy.zeros((len(vhelio),_BINTEGRATENBINS))
            logpddwarf= numpy.zeros((len(vhelio),_BINTEGRATENBINS))
            dwarfds= numpy.linspace(_BINTEGRATEDMIN_DWARF,_BINTEGRATEDMAX_DWARF,
                                    _BINTEGRATENBINS)/params[1]/_REFR0
        if options.multi > 1: #BOVY: MAY BE BROKEN
            raise NotImplementedError("'multi' evaluation of the log likelihood in all likelihood broken")
            thisout= numpy.zeros((len(vhelio),_BINTEGRATENBINS,2))
            thisout= multi.parallel_map((lambda x: _mloglikedIntegrand(ds[x],
                                                                       params,vhelio/params[0]/_REFV0,
                                                                       l,b,jk,h,df,options,
                                                                       sinl,cosl,cosb,sinb,
                                                                       True,logpiso[:,x],vcf,iso,feh,True)),
                                               range(len(ds)),
                                        numcores=numpy.amin([len(ds),multiprocessing.cpu_count(),options.multi]))
            #thisout is list of len(vhelio)
            out= numpy.zeros((len(vhelio),_BINTEGRATENBINS))
            for ii in range(_BINTEGRATENBINS):
                out[:,ii]= thisout[ii][:,0]
                logpd[:,ii]= thisout[ii][:,1]
            thisout= out
        else:
            for ii in range(_BINTEGRATENBINS):
                thisout[:,ii], logpd[:,ii]= _mloglikedIntegrand(ds[ii],
                                                                params,vhelio/params[0]/_REFV0,
                                                                l,b,jk,h,df,options,
                                                                sinl,cosl,cosb,sinb,
                                                                True,logpiso[:,ii],vcf,iso,feh)
                if options.dwarf:
                    if ii > 0: continue #HACK
                    thisxtraout[:,ii], logpddwarf[:,ii]= _mloglikedIntegrand(dwarfds[ii],
                                                                             params,vhelio/params[0]/_REFV0,
                                                                             l,b,jk,h,df,options,
                                                                             sinl,cosl,cosb,sinb,
                                                                             True,logpisodwarf[:,ii],vcf,iso,feh,dwarf=True)
        #Sum each one
        if arrout:
            out= numpy.zeros(len(vhelio))
        else:
            out= 0.
        for ii in range(len(vhelio)):
            if options.dwarf: #Combine
                tmpthisout= logsumexp(thisout[ii,:])+numpy.log(1.-params[5-options.nooutliermean])-logsumexp(logpd[ii,:])
                if options.dwarfinnerouter and l[ii] < 97.*_DEGTORAD:
                    tmpthisout+= -numpy.log(1.-params[5-options.nooutliermean])+numpy.log(1.-params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.fitdm+options.fitah+options.fitfeh])
                thislogpdwarf= logpddwarf[ii,0]
                tmpthisxtraout= thisxtraout[ii,0]+numpy.log(params[5-options.nooutliermean])-thislogpdwarf
                if options.dwarfinnerouter and l[ii] < 97.*_DEGTORAD:
                    tmpthisxtraout+= -numpy.log(params[5-options.nooutliermean])+numpy.log(params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.fitdm+options.fitah+options.fitfeh])
                #print tmpthisout, tmpthisxtraout, thisxtraout[ii,0]
                c= numpy.amax([tmpthisout,tmpthisxtraout])
                if arrout:
                    out[ii]= -(numpy.log(numpy.exp(tmpthisout-c)+numpy.exp(tmpthisxtraout-c))+c)
                    if numpy.log10(numpy.fabs(out[ii])) <=-2. : out[ii]= numpy.finfo(numpy.dtype(numpy.float64)).max #means all distances were zero prob
                else:
                    out-= numpy.log(numpy.exp(tmpthisout-c)+numpy.exp(tmpthisxtraout-c))+c
            else:
                if arrout:
                    out[ii]= -logsumexp(thisout[ii,:])+logsumexp(logpd[ii,:]) #so we can compare to the +dwarf case
                    if numpy.log10(numpy.fabs(out[ii])) <=-2. : out[ii]= numpy.finfo(numpy.dtype(numpy.float64)).max #means all distances were zero prob
                else:
                    out+= -logsumexp(thisout[ii,:])+logsumexp(logpd[ii,:]) #so we can compare to the +dwarf case
        #print out, params
    #BOVY:Apply Ro prior correcly
    if options.noroprior:
        return out
    else:
        return out+(params[1]*_REFR0-8.2)**2./0.5#+(params[0]*_REFV0+2.25*math.exp(2.*params[2])*_REFV0/params[0]+12.24-_PMSGRA*params[1]*_REFR0)**2./200. #params[1]=Ro, SBD10 Solar motion

def _mloglikedIntegrand(d,params,vhelio,l,b,jk,h,
                        df,options,sinl,cosl,cosb,sinb,returnlog,logpiso,vcf,
                        iso,feh,
                        pdinout=False,dwarf=False):
    #All positions are /ro, all velocities are /vo (d and vhelio have this already
    #Calculate coordinates, distances are /Ro (/params[1])
    #
    #Returns the loglike and the p(d) factor
    R= numpy.sqrt(1.+d**2.-2.*d*cosl)
    if isinstance(vhelio,numpy.ndarray) and len(vhelio) > 1:
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
    vpec= _vpec(params,vgal,R,options,l,theta,vcf)
    #Calculate probabilities
    if dwarf: #assume R=1, theta=0
        vpec= _vpec(params,vgal,1.,options,l,0.,vcf)
        logpvlos= _logdf(params,vpec,numpy.ones(len(vpec)),
                         options,df,l,numpy.zeros(len(vpec)),
                         vcf)+numpy.log(1.-params[3])
    else:
        logpvlos= _logdf(params,vpec,R,options,df,l,theta,vcf)+numpy.log(1.-params[3])
    logpvlos_outlier= _logoutlierdf(params,vgal,options)+numpy.log(params[3])
    c= numpy.amax(numpy.array([logpvlos,logpvlos_outlier]),axis=0)
    logpvlos= numpy.log(numpy.exp(logpvlos-c)
                        +numpy.exp(logpvlos_outlier-c))+c
    if options.indivfeh and (options.fitdm or options.fitah or options.fitfeh):
        zs= numpy.arange(0.0005,0.03005,0.0005)
    if dwarf:
        logpd= 0.
    #Re-calculate logpiso if necessary
    elif options.fitdm: # and not (params[5+2*options.fitvpec+options.dwarf+options.fitsratio] == 0.):
        dm= _dm(d*params[1]*_REFR0)\
            -params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter]
        if options.fitdminnerouter: dminnerouter= _dm(d*params[1]*_REFR0)\
                -params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitdm+options.fitah+options.fitfeh]
        for ii in range(len(vhelio)):
            mh= h[ii]-dm
            if options.fitdminnerouter and l[ii] < 97.*_DEGTORAD:
                mh= h[ii]-dminnerouter
            if options.indivfeh:
                #Find closest Z
                thisZ= isodist.FEH2Z(feh[ii])
                indx= numpy.argmin(numpy.fabs(thisZ-zs))
                logpiso[ii]= iso[0][indx](jk[ii],mh)
            elif options.varfeh:
                raise NotImplementedError("changing dm with varfeh not implemented yet")
                for jj in range(len(iso[0])):
                    logpiso[ii,jj]= iso[0][jj](jk[ii],mh)
            else:
                logpiso[ii]= iso[0](jk[ii],mh)
        logpd= _logpd(params,d,l,b,jk,h,df,options,R,theta,cosb,sinb,logpiso)
    elif options.fitah:
        ah= params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter]
        if options.fitahinnerouter: ahinnerouter= params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitah+options.fitdm+options.fitfeh]
        dm= _dm(d*params[1]*_REFR0)
        for ii in range(len(vhelio)):
            mh= h[ii]-dm+ah
            if options.indivfeh:
                #Find closest Z
                thisZ= isodist.FEH2Z(feh[ii])
                indx= numpy.argmin(numpy.fabs(thisZ-zs))
                logpiso[ii]= iso[0][indx](jk[ii]+1.5/1.55*ah,mh) #Accounts for red. law
            elif options.varfeh:
                raise NotImplementedError("changing ah with varfeh not implemented yet")
                for jj in range(len(iso[0])):
                    logpiso[ii,jj]= iso[0][jj](jk[ii]+1.5/1.55*ah,mh) #Accounts for red. law
            else:
                logpiso[ii]= iso[0](jk[ii]+1.5/1.55*ah,mh) #Accounts for red. law
            if options.fitahinnerouter and l[ii] < 97.*_DEGTORAD:
                mh= h[ii]-dm+ahinnerouter
                if options.indivfeh:
                    #Find closest Z
                    thisZ= isodist.FEH2Z(feh[ii])
                    indx= numpy.argmin(numpy.fabs(thisZ-zs))
                    logpiso[ii]= iso[0][indx](jk[ii]+1.5/1.55*ah,mh) #Accounts for red. law
                elif options.varfeh:
                    raise NotImplementedError("changing ah with varfeh not implemented yet")
                    for jj in range(len(iso[0])):
                        logpiso[ii,jj]= iso[0][jj](jk[ii]+1.5/1.55*ahinnerouter,mh) #Accounts for red. law
                else:
                    logpiso[ii]= iso[0](jk[ii]+1.5/1.55*ahinnerouter,mh) #Accounts for red. law
        logpd= _logpd(params,d,l,b,jk,h,df,options,R,theta,cosb,sinb,logpiso)
    elif options.fitfeh:
        dm= _dm(d*params[1]*_REFR0)
        dfeh= params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter]
        if options.fitfehinnerouter: dfehinnerouter= params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitdm+options.fitah+options.fitfeh]
        for ii in range(len(vhelio)):
            mh= h[ii]-dm
            if options.indivfeh:
                #Find closest Z
                if options.fitfehinnerouter and l[ii] < 97.*_DEGTORAD:
                    thisZ= isodist.FEH2Z(feh[ii]+dfehinnerouter)
                else:
                    thisZ= isodist.FEH2Z(feh[ii]+dfeh)
                indx= numpy.argmin(numpy.fabs(thisZ-zs))
                logpiso[ii]= iso[0][indx](jk[ii],mh)
            elif options.varfeh:
                raise NotImplementedError("changing feh with varfeh not implemented yet")
                for jj in range(len(iso[0])):
                    logpiso[ii,jj]= iso[0][jj](jk[ii],mh)
            else:
                logpiso[ii]= iso[0](jk[ii],mh)
        logpd= _logpd(params,d,l,b,jk,h,df,options,R,theta,cosb,sinb,logpiso)
    else:
        logpd= _logpd(params,d,l,b,jk,h,df,options,R,theta,cosb,sinb,logpiso)
    if pdinout:
        out= numpy.empty((logpvlos.shape[0],2))
        if returnlog:
            out[:,0]= logpvlos+logpd
            out[:,1]= logpd
        else:
            out[:,0]= numpy.exp(logpvlos+logpd)
            out[:,1]= numpy.exp(logpd)
        return out
    else:
        if returnlog: return (logpvlos+logpd,logpd)
        else: return (numpy.exp(logpvlos+logpd),numpy.exp(logpd))

def _logpd(params,d,l,b,jk,h,df,options,R,theta,cosb,sinb,logpiso):
    logpddf= _logpddf(params,d,l,b,R,theta,cosb,sinb,options)
    return logpiso+logpddf

def _logpddf(params,d,l,b,R,theta,cosb,sinb,options):
    """Prior on distance distribution, including Jacobian"""
    absZ= numpy.fabs(d*sinb)
    if options.dfmodel.lower() == 'twopops':
        npops= 2
        hrpops= [options.hr,
                 options.hr2]
        ppops= numpy.log(numpy.array([1.-params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fitfehinnerouter],
                                      params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fitfehinnerouter]]))
        logdensRZ= numpy.zeros((len(R),npops))
        for ii in range(npops):
            logdensRZ[:,ii]= (-(R-1.)/hrpops[ii]\
#                                   -absZ/options.hz)*params[1]*_REFR0\
                                   )-ppops[ii]
        logdensRZ= mylogsumexp(logdensRZ,axis=1)
    elif options.densmodel.lower() == 'expdisk':
        logdensRZ= (-(R-1.)/options.hr-0.*absZ/options.hz)*params[1]*_REFR0
    return logdensRZ+2.*numpy.log(d*params[1])+numpy.log(cosb)

def _dm(d):
    """Distance modulus w/ d in kpc"""
    return 5.*numpy.log10(d)+10.

def _logdf(params,vpec,R,options,df,l,theta,vcf):
    sinlt= numpy.sin(l+theta)
    thishs= options.hs
    if options.fiths: thishs/= params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitdm+options.fitah+options.fitfeh]
    if options.dfmodel.lower() == 'simplegaussian':
        if options.fitsratio:
            slos= numpy.exp(params[2])/params[0]\
                *numpy.sqrt(1.+sinlt**2.*(params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf]-1.))*numpy.exp(-(R-1.)/thishs*params[1]*_REFR0)
        elif options.fitsratioinnerouter:
            innerl= (l < 97.*_DEGTORAD)
            slos= numpy.exp(params[2])/params[0]\
                *numpy.sqrt(1.+sinlt**2.*(params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf]-1.))*numpy.exp(-(R-1.)/thishs*params[1]*_REFR0)
            slos[innerl]= numpy.exp(params[2])/params[0]\
                *numpy.sqrt(1.+sinlt[innerl]**2.*(params[6-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf]-1.))*numpy.exp(-(R[innerl]-1.)/thishs*params[1]*_REFR0)
        else:
            slos= numpy.exp(params[2])/params[0]\
                *numpy.sqrt(1.-0.5*sinlt**2.)*numpy.exp(-(R-1.)/thishs*params[1]*_REFR0)
        t= vpec/slos
        return norm.logpdf(t)-numpy.log(slos*params[0]*_REFV0)
    elif options.dfmodel.lower() == 'simplegaussiandrift':
        va= asymmetricDriftModel.va(R,numpy.exp(params[2])/params[0],
                                    hR=options.hr/params[1]/_REFR0,
                                    hs=thishs/params[1]/_REFR0,
                                    vc=_vc(params,R,options,vcf))*sinlt 
        #va= vc- <v>
        if options.fitsratio:
            slos= numpy.exp(params[2])/params[0]\
                *numpy.sqrt(1.+sinlt**2.*(params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf]-1.))*numpy.exp(-(R-1.)/thishs*params[1]*_REFR0)          
            if options.fitsrinnerouter and not numpy.all(R == 1.): #Don't do this for dwarfs
                innerl= (l < 97.*_DEGTORAD)
                slos[innerl]*= numpy.exp(params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fiths])
                va[innerl]= asymmetricDriftModel.va(R[innerl],numpy.exp(params[2]+params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fiths])/params[0],
                                                    hR=options.hr/params[1]/_REFR0,
                                                    hs=thishs/params[1]/_REFR0,
                                                    vc=_vc(params,R[innerl],options,vcf))*sinlt[innerl]
        elif options.fitsratioinnerouter:
            innerl= (l < 97.*_DEGTORAD)
            slos= numpy.exp(params[2])/params[0]\
                *numpy.sqrt(1.+sinlt**2.*(params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf]-1.))*numpy.exp(-(R-1.)/thishs*params[1]*_REFR0)          
            slos[innerl]= numpy.exp(params[2])/params[0]\
                *numpy.sqrt(1.+sinlt[innerl]**2.*(params[6-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf]-1.))*numpy.exp(-(R[innerl]-1.)/thishs*params[1]*_REFR0)          
            if options.fitsrinnerouter and not numpy.all(R == 1.): #Don't do this for dwarfs
                innerl= (l < 97.*_DEGTORAD)
                slos[innerl]*= numpy.exp(params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fiths])
                va[innerl]= asymmetricDriftModel.va(R[innerl],numpy.exp(params[2]+params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fiths])/params[0],
                                                    hR=options.hr/params[1]/_REFR0,
                                                    hs=thishs/params[1]/_REFR0,
                                                    vc=_vc(params,R[innerl],options,vcf))*sinlt[innerl]
        else:
            slos= numpy.exp(params[2])/params[0]\
                *numpy.sqrt(1.-0.5*sinlt**2.)*numpy.exp(-(R-1.)/thishs*params[1]*_REFR0)
        t= (vpec+va)/slos
        return norm.logpdf(t)-numpy.log(slos*params[0]*_REFV0)
    elif options.dfmodel.lower() == 'multiplepops':
        out= numpy.zeros((len(vpec),_NMULTIPLEPOPS))
        taupops= numpy.linspace(1.,10.,_NMULTIPLEPOPS)
        lognorm= logsumexp(taupops/_TAU0)
        sigmapops= ((taupops+_TAU1)/(_TAUM+_TAU1))**_TAUBETA
        for ii in range(_NMULTIPLEPOPS):
            va= asymmetricDriftModel.va(R,
                                        numpy.exp(params[2])/params[0]*sigmapops[ii],
                                        hR=options.hr/params[1]/_REFR0,
                                        hs=thishs/params[1]/_REFR0,
                                        vc=_vc(params,R,options,vcf))*sinlt 
            #va= vc- <v>
            if options.fitsratio:
                slos= numpy.exp(params[2])/params[0]\
                    *numpy.sqrt(1.+sinlt**2.*(params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf]-1.))*numpy.exp(-(R-1.)/thishs*params[1]*_REFR0)*sigmapops[ii]
            elif options.fitsratioinnerouter:
                innerl= (l < 97.*_DEGTORAD)
                slos= numpy.exp(params[2])/params[0]\
                    *numpy.sqrt(1.+sinlt**2.*(params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf]-1.))*numpy.exp(-(R-1.)/thishs*params[1]*_REFR0)*sigmapops[ii]
                slos[innerl]= numpy.exp(params[2])/params[0]\
                    *numpy.sqrt(1.+sinlt[innerl]**2.*(params[6-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf]-1.))*numpy.exp(-(R[innerl]-1.)/thishs*params[1]*_REFR0)*sigmapops[ii]
            else:
                slos= numpy.exp(params[2])/params[0]\
                    *numpy.sqrt(1.-0.5*sinlt**2.)*numpy.exp(-(R-1.)/thishs*params[1]*_REFR0)*sigmapops[ii]
            t= (vpec+va)/slos
            out[:,ii]= taupops[ii]/_TAU0-lognorm\
                +norm.logpdf(t)-numpy.log(slos*params[0]*_REFV0)
        #Sum
        out= mylogsumexp(out,axis=1)
        return out
    elif options.dfmodel.lower() == 'twopops':
        npops= 2
        out= numpy.zeros((len(vpec),npops))
        sigmapops= [numpy.exp(params[2])/params[0],
                    numpy.exp(params[6-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fitfehinnerouter])/params[0]]
        ppops= numpy.log(numpy.array([1.-params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fitfehinnerouter],
                                      params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fitfehinnerouter]]))
        hrpops= [options.hr/params[1]/_REFR0,
                 options.hr2/params[1]/_REFR0]
        logpk= numpy.zeros((len(vpec),npops))
        for ii in range(npops):
            va= asymmetricDriftModel.va(R,
                                        sigmapops[ii],
                                        hR=hrpops[ii],
                                        hs=thishs/params[1]/_REFR0,
                                        vc=_vc(params,R,options,vcf))*sinlt
            #va= vc- <v>
            if options.fitsratio:
                slos= numpy.sqrt(1.+sinlt**2.*(params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf]-1.))*numpy.exp(-(R-1.)/thishs*params[1]*_REFR0)*sigmapops[ii]
            elif options.fitsratioinnerouter:
                innerl= (l < 97.*_DEGTORAD)
                slos= numpy.sqrt(1.+sinlt**2.*(params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf]-1.))*numpy.exp(-(R-1.)/thishs*params[1]*_REFR0)*sigmapops[ii]
                slos[innerl]= numpy.sqrt(1.+sinlt[innerl]**2.*(params[6-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf]-1.))*numpy.exp(-(R[innerl]-1.)/thishs*params[1]*_REFR0)*sigmapops[ii]
            else:
                slos= numpy.sqrt(1.-0.5*sinlt**2.)*numpy.exp(-(R-1.)/thishs*params[1]*_REFR0)*sigmapops[ii]
            t= (vpec+va)/slos
            logpk[:,ii]= ppops[ii]-(R-1.)/hrpops[ii]
            out[:,ii]= logpk[:,ii]\
                +norm.logpdf(t)-numpy.log(slos*params[0]*_REFV0)
        #Sum
        lognorm= mylogsumexp(logpk,axis=1)
        out= mylogsumexp(out,axis=1)
        return out-lognorm
    elif options.dfmodel.lower() == 'simpleskeweddrift':
        coslt= numpy.cos(l+theta)
        va= asymmetricDriftModel.va(R,numpy.exp(params[2])/params[0],
                                    hR=options.hr/params[1]/_REFR0,
                                    hs=thishs/params[1]/_REFR0,
                                    vc=_vc(params,R,options,vcf))#no sinlt
        #va= vc- <v>
        alphaskew= params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter]
        delta= alphaskew/math.sqrt(1.+alphaskew**2.)
        sigmaR2= numpy.exp(2.*params[2])/params[0]**2.*numpy.exp(-2.*(R-1.)/thishs*params[1]*_REFR0)
        sigmaR= numpy.sqrt(sigmaR2)
        if options.fitsratio:
            sigmaT2= params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf]*sigmaR2
        elif options.fitsratio:
            innerl= (l < 97.*_DEGTORAD)
            sigmaT2= params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf]*sigmaR2
            sigmaT2[innerl]= params[6-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf]*sigmaR2[innerl]
        else:
            sigmaT2= 0.5*sigmaR2
        omega= numpy.sqrt(sigmaT2/(1.-2.*delta**2./math.pi))
        xi= -va-omega*delta*math.sqrt(2./math.pi)
        out= numpy.log(evalSkeweddf(vpec,sinlt,coslt,sigmaR,xi,omega,alphaskew))
        return out-numpy.log(_REFV0) #For comparison
    elif options.dfmodel.lower() == 'dehnen':
        coslt= numpy.cos(l+theta)
        sigmaR= numpy.exp(params[2])/params[0]*numpy.exp(-(R-1.)/thishs*params[1]*_REFR0)
        out= numpy.log(evaldehnendf(vpec,sinlt,coslt,df,R,sigmaR,options,
                                    params))
        return out-numpy.log(_REFV0) #For comparison

def _logoutlierdf(params,vgal,options):
    if options.nooutliermean:
        t= vgal*params[0]*_REFV0/100.
    else:
        t= (vgal*params[0]-params[4])*_REFV0/100.
    return norm.logpdf(t)-numpy.log(100.)

def _vc(params,R,options,vcf):
    """Circular velocity at R for different models"""
    if options.rotcurve.lower() == 'flat':
        return 1. #vc/vo
    elif options.rotcurve.lower() == 'flatplusuniform': #Useless?
        return 1.+R*params[5-options.nooutliermean+options.dwarf] #vc/vo
    elif options.rotcurve.lower() == 'powerlaw':
        return R**params[5-options.nooutliermean+options.dwarf] #vc/vo
    elif options.rotcurve.lower() == 'linear':
        return 1.+(R-1.)*params[5-options.nooutliermean+options.dwarf] #vc/vo
    elif options.rotcurve.lower() == 'quadratic':
        return 1.+(R-1.)*params[5-options.nooutliermean+options.dwarf]\
            +(R-1.)**2.*params[6-options.nooutliermean+options.dwarf]
    elif options.rotcurve.lower() == 'cubic':
        return 1.+(R-1.)*params[5-options.nooutliermean+options.dwarf]\
            +(R-1.)**2.*params[6-options.nooutliermean+options.dwarf]\
            +(R-1.)**3.*params[7-options.nooutliermean+options.dwarf]
    elif options.rotcurve.lower() == 'gp':
        return 1.+vcf(R*params[1]*_REFR0)/params[0] #interpolation of GP f

def _vpec(params,vgal,R,options,l,theta,vcf):
    if options.fitm2:
        eps= params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fiths]
        phio= math.pi*params[6-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fitdm+options.fitah+options.fitfeh+options.fiths]
        #Defined as in Kuijken & Tremaine (1994)
        dvt= -eps*numpy.cos(2.*(theta-phio))
        dvr= -eps*numpy.sin(2.*(theta-phio))
        return vgal-(_vc(params,R,options,vcf)+dvt)*numpy.sin(l+theta)\
            +dvr*numpy.cos(l+theta)
    else:
        return vgal-_vc(params,R,options,vcf)*numpy.sin(l+theta)

def _vgal(params,vhelio,l,b,options,sinl,cosl):
    if options.fitvpec:
        return vhelio-params[5-options.nooutliermean+options.dwarf+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')]*cosl*_VRSUN/params[0]/_REFV0\
            +params[6-options.nooutliermean+(options.rotcurve.lower() == 'linear')+(options.rotcurve.lower() == 'flatplusuniform') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+options.dwarf]*sinl*_PMSGRA*params[1]*_REFR0/params[0]/_REFV0 #params[1]=Ro
    elif options.sbdvpec:
        return vhelio-cosl*_VRSUN/params[0]/_REFV0+sinl*(_VTSUN/params[0]/_REFV0+1.) #params[1]=Ro
    else:
        return vhelio-cosl*_VRSUN/params[0]/_REFV0+sinl*_PMSGRA*params[1]*_REFR0/params[0]/_REFV0 #params[1]=Ro

def evaldehnendf(vpec,sinlt,coslt,df,R,sigmaR,options,params):
    """Evaluate the dehnen df by integrating over the tangential velocity"""
    if isinstance(vpec,numpy.ndarray):
        #Integrate by binning
        vs= numpy.linspace(-5.,5.,_BINTEGRATEVNBINS)
        out= numpy.zeros((len(vpec),_BINTEGRATEVNBINS))
        cotlt= coslt/sinlt
        tanlt= sinlt/coslt
        for ii in range(len(vpec)):
            if sinlt[ii]**2. > 0.5:
                out[ii,:]= sigmaR[ii]/numpy.fabs(sinlt[ii])\
                    *evaldehnendfIntegrandLargelt(vs,cotlt[ii],
                                              sinlt[ii],
                                              vpec[ii],df,R[ii],sigmaR[ii])
            else:
                out[ii,:]= sigmaR[ii]/numpy.fabs(coslt[ii])\
                    *evaldehnendfIntegrandSmalllt(vs,tanlt[ii],
                                                  coslt[ii],
                                                  vpec[ii],df,R[ii],sigmaR[ii])
        return numpy.sum(out,axis=1)*(vs[1]-vs[0])/numpy.exp(-R/options.hr*params[1]*_REFR0) #This last factor corrects for the normalization in galpy (-ish since the uncorrected surface-mass density isn't quite exponential)

def evaldehnendfIntegrandLargelt(vR,cotlt,sinlt,vlos,df,R,sigmaR):
    #Create input for evaluation
    dfin= numpy.zeros((3,len(vR)))
    dfin[0,:]= R
    dfin[1,:]= vR*sigmaR
    dfin[2,:]= 1.+cotlt*vR*sigmaR+vlos/sinlt
    return df(dfin)
              
def evaldehnendfIntegrandSmalllt(vT,tanlt,coslt,vlos,df,R,sigmaR):
    #Create input for evaluation
    dfin= numpy.zeros((3,len(vT)))
    dfin[0,:]= R
    dfin[1,:]= tanlt*vT*sigmaR-vlos/coslt
    dfin[2,:]= vT*sigmaR+1.
    return df(dfin)

def evalSkeweddf(vpec,sinlt,coslt,sigmaR,xi,omega,alphaskew):
    """Evaluate the skewed df by integrating over the tangential velocity"""
    if isinstance(vpec,numpy.ndarray):
        #Integrate by binning
        vs= numpy.linspace(-5.,5.,_BINTEGRATEVNBINS)
        out= numpy.zeros((len(vpec),_BINTEGRATEVNBINS))
        cotlt= coslt/sinlt
        tanlt= sinlt/coslt
        for ii in range(len(vpec)):
            if sinlt[ii]**2. > 0.5:
                out[ii,:]= sigmaR[ii]/numpy.fabs(sinlt[ii])\
                *evalSkeweddfIntegrandLargelt(vs,cotlt[ii],
                                              sinlt[ii],
                                              vpec[ii],sigmaR[ii],
                                              xi[ii],omega[ii],alphaskew)
            else:
                out[ii,:]= sigmaR[ii]/numpy.fabs(coslt[ii])\
                    *evalSkeweddfIntegrandSmalllt(vs,tanlt[ii],
                                              coslt[ii],
                                              vpec[ii],sigmaR[ii],
                                              xi[ii],omega[ii],
                                              alphaskew)
        return numpy.sum(out,axis=1)*(vs[1]-vs[0])
    if sinlt**2. > 0.5:
        cotlt= coslt/sinlt
        return sigmaR/numpy.fabs(sinlt)*integrate.quadrature(evalSkeweddfIntegrandLargelt,
                                                             -5.,5., #5 sigma
                                                        args=(cotlt,
                                                              sinlt,vpec,sigmaR,
                                                              xi,omega,alphaskew))[0]
    else:
        tanlt= sinlt/coslt
        return sigmaR/numpy.fabs(coslt)*integrate.quadrature(evalSkeweddfIntegrandSmalllt,
                                                             -5.,5., #5 sigma
                                                             args=(tanlt,
                                                                   coslt,
                                                                   vpec,sigmaR,
                                                              xi,omega,
                                                              alphaskew))[0]

def evalSkeweddfIntegrandLargelt(vR,cotlt,sinlt,vlos,
                                 sigmaR,xi,omega,alphaskew):
    return skeweddf(vR*sigmaR,cotlt*vR*sigmaR+vlos/sinlt,
                    sigmaR,xi,omega,alphaskew)
def evalSkeweddfIntegrandSmalllt(vT,tanlt,coslt,vlos,
                                 sigmaR,xi,omega,alphaskew):
    #Definitely okay
    return skeweddf(tanlt*vT*sigmaR-vlos/coslt,
                    vT*sigmaR,
                    sigmaR,xi,omega,alphaskew)

def skeweddf(vR,vT,sigmaR,xi,omega,alphaskew):
    """skewed df= Gaussian in vR, Skew normal in vT
    vR,vT
    sigmaR = radial dispersion
    xi - xi parameter of skew normal 'mean'
    omega - omega parameter of skew normal 'variance'
    alphaskew= skew normal alpha parameter"""
    #First calculate the skew normal's parameters
    return norm.pdf(vR/sigmaR)/sigmaR\
        *skewnormal.skewnormal(vT,m=xi,s=omega,a=alphaskew)

def mylogsumexp(arr,axis=0):
    """Faster logsumexp?"""
    minarr= numpy.amax(arr,axis=axis)
    if axis == 1:
        minarr= numpy.reshape(minarr,(arr.shape[0],1))
    if axis == 0:
        minminarr= numpy.tile(minarr,(arr.shape[0],1))
    elif axis == 1:
        minminarr= numpy.tile(minarr,(1,arr.shape[1]))
    elif axis == None:
        minminarr= numpy.tile(minarr,arr.shape)
    else:
        raise NotImplementedError("'mylogsumexp' not implemented for axis > 2")
    if axis == 1:
        minarr= numpy.reshape(minarr,(arr.shape[0]))
    return minarr+numpy.log(numpy.sum(numpy.exp(arr-minminarr),axis=axis))

def print_samples_qa(samples):
    print "Mean, standard devs, acor tau, acor mean, acor s ..."
    for kk in range(len(samples[0])):
        xs= numpy.array([s[kk] for s in samples])
        #Auto-correlation time
        tau, m, s= acor.acor(xs)
        print numpy.mean(xs), numpy.std(xs), tau, m, s

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the fit/samples will be saved to"
    parser = OptionParser(usage=usage)
    #Initial conditions file
    parser.add_option("--init",dest='init',default=None,
                      help="Initial parameters")
    #Rotation curve parameters/model
    parser.add_option("--rotcurve",dest='rotcurve',default='flat',
                      help="Rotation curve model to fit")
    parser.add_option("--gpnr",dest='gpnr',default=101,type='int',
                      help="If using a GP for the rotcurve, this is the number of nodes")
    parser.add_option("--gprmin",dest='gprmin',default=5.,type='float',
                      help="If using a GP for the rotcurve, this is minimum R")
    parser.add_option("--gprmax",dest='gprmax',default=18.,type='float',
                      help="If using a GP for the rotcurve, this is maximum R")
    parser.add_option("--gpfixtau",dest='gpfixtau',default=None,type='float',
                      help="If using a GP for the rotcurve and this is set, fix tau (correlation length) to this value in kpc")
    parser.add_option("--gpemcee",action="store_true", dest="gpemcee",
                      default=False,
                      help="If set, use emcee for GP indeed")
    parser.add_option("--fitm2",action="store_true", dest="fitm2",
                      default=False,
                      help="If set, fit for an m=2 component")
    #Fix vo? CRAZY!!
    parser.add_option("--fixvo",dest='fixvo',default=None,type='float',
                      help="If set, fix vo to this value, and optimize other parameters")
    #Fix ro
    parser.add_option("--fixro",dest='fixro',default=None,type='float',
                      help="If set, fix ro to this value, and optimize other parameters")
    #Ro prior
    parser.add_option("--noroprior",action="store_true", dest="noroprior",
                      default=False,
                      help="If set, do not apply an Ro prior")
    #Sun's peculiar velocity
    parser.add_option("--fitvpec",action="store_true", dest="fitvpec",
                      default=False,
                      help="If set, fit for the peculiar velocity of the Sun as well")
    parser.add_option("--sbdvpec",action="store_true", dest="sbdvpec",
                      default=False,
                      help="If set, use the SBD10 value for the Solar motion + vo as vpec")
    #Velocity distribution model
    parser.add_option("--dfmodel",dest='dfmodel',default='simplegaussian',
                      help="DF model to use")
    parser.add_option("--fitsratio",action="store_true", dest="fitsratio",
                      default=False,
                      help="If set, fit for the ration squared of tangential to radial dispersion")
    parser.add_option("--fitsratioinnerouter",
                      action="store_true", dest="fitsratioinnerouter",
                      default=False,
                      help="If set, fit for the ratio squared of tangential to radial dispersion")
    parser.add_option("--nooutliermean",action="store_true", 
                      dest="nooutliermean",
                      default=False,
                      help="If set, use a zero mean for the outlier model (in Galactocentric coordinates)")
    parser.add_option("--fitsrinnerouter",
                      action="store_true", dest="fitsrinnerouter",
                      default=False,
                      help="If set, fit for the sigma_r separately for the inner disk")
    parser.add_option("--dwarfinnerouter",
                      action="store_true", dest="dwarfinnerouter",
                      default=False,
                      help="If set, fit for the dwarf contamination separately for the inner disk")
    parser.add_option("--fitahinnerouter",
                      action="store_true", dest="fitahinnerouter",
                      default=False,
                      help="If set, fit for the ah separately for the inner disk")
    parser.add_option("--fitdminnerouter",
                      action="store_true", dest="fitdminnerouter",
                      default=False,
                      help="If set, fit for the dm separately for the inner disk")
    parser.add_option("--fitfehinnerouter",
                      action="store_true", dest="fitfehinnerouter",
                      default=False,
                      help="If set, fit for the feh offset separately for the inner disk")
    #Density model
    parser.add_option("--densmodel",dest='densmodel',default='expdisk',
                      help="Density model to use")
    parser.add_option("--hr",dest='hr',default=3.,type='float',
                      help="scale length in kpc")
    parser.add_option("--hz",dest='hz',default=0.25,type='float',
                      help="scale height in kpc")
    parser.add_option("--hs",dest='hs',default=8.,type='float',
                      help="dispersion scale length in kpc")
    parser.add_option("--fiths",action="store_true", dest="fiths",
                      default=False,
                      help="If set, fit for a dispersion scale length offsett")
    #Data options
    parser.add_option("--lmin",dest='lmin',default=25.,type='float',
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
    parser.add_option("--jkmax",dest='jkmax',default=1.1,type='float',
                      help="readVclosData 'jkmax'")
    parser.add_option("--location",dest='location',default=None,type='int',
                      help="location id when looking at single los")
    parser.add_option("--removelocation",
                      dest='removelocation',default=None,type='int',
                      help="location id to remove")
    parser.add_option("--downsample",dest='downsample',default=None,
                      type='float',
                      help="Factor with which to downsample the data")
    parser.add_option("--cutmultiples",action="store_true", 
                      dest="cutmultiples",
                      default=False,
                      help="readVclosData 'cutmultiples'")
    parser.add_option("--correctak",action="store_true", 
                      dest="correctak",
                      default=False,
                      help="readVclosData 'correctak'")
    parser.add_option("--loggcut",dest='loggcut',type='float',
                      default=None,
                      help="readVclosData 'loggcut'")
    parser.add_option("--justinner",action="store_true", 
                      dest="justinner",
                      default=False,
                      help="Only fit inner (l < 97) fields")
    parser.add_option("--justouter",action="store_true", 
                      dest="justouter",
                      default=False,
                      help="Only fit outer (l > 97) fields")
    parser.add_option("-f",dest='fakedata',default=None,
                      help="Name of the fake data filename")
    parser.add_option("--fitdl",action="store_true", dest="fitdl",
                      default=False,
                      help="If set, fit for an offset in l for the GC")
    #Isochrone IMF
    parser.add_option("--imfmodel",dest='imfmodel',default='lognormalChabrier2001',
                      help="imfmodel for isochrone model")
    parser.add_option("--Z",dest='Z',default=.019,type='float',
                      help="Metallicity of isochrone")
    parser.add_option("--nods",action="store_true", 
                      dest="nods",
                      default=False,
                      help="setting this assumes no distance information from isochrones, NOT SUPPORTED WITH DWARF")
    parser.add_option("--fitdm",action="store_true", dest="fitdm",
                      default=False,
                      help="If set, fit for a distance modulus offset")
    parser.add_option("--fitah",action="store_true", dest="fitah",
                      default=False,
                      help="If set, fit for an extinction-correction offset")
    parser.add_option("--fitfeh",action="store_true", dest="fitfeh",
                      default=False,
                      help="If set, fit for a [Fe/H] offset (with indivfeh)")
    parser.add_option("--expsfh",action="store_true", dest="expsfh",
                      default=False,
                      help="If set, use an exponentially declining SFH")
    parser.add_option("--varfeh",action="store_true", dest="varfeh",
                      default=False,
                      help="If set, don't use a varying [Fe/H] distribution as a function of l")
    parser.add_option("--indivfeh",action="store_false", dest="indivfeh",
                      default=True,
                      help="If set, use distances calculated based on the individual [Fe/H] of the objects")
    parser.add_option("--isofile",dest="isofile",default=None,
                      help="if set, store or restore the isochrone model(s) in this file")
    parser.add_option("--basti",action="store_true", 
                      dest="basti",
                      default=False,
                      help="Use Basti isochrones")
    #Add dwarf part?
    parser.add_option("--dwarf",action="store_true", 
                      dest="dwarf",
                      default=False,
                      help="setting this adds dwarf contamination")
    #Second population
    parser.add_option("--hr2",dest='hr2',default=5.,type='float',
                      help="Second population scale length")
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
                      help="number of cpus to use for sampling or evaluation")
    return parser

if __name__ == '__main__':
    numpy.random.seed(1) #We need to seed to get, e.g., the same permutation when downsampling
    fitvc(get_options())
