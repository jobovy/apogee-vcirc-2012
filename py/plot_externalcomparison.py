import sys
import copy
import os, os.path
import cPickle as pickle
from optparse import OptionParser
import math
import numpy
from scipy.maxentropy import logsumexp
from scipy import special, stats
from galpy.util import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
import multi
import multiprocessing
import isomodel
import isodist
from readVclosData import readVclosData
from fitvc import mloglike, _dm, \
    _DEGTORAD, _REFV0, _REFR0, _PMSGRA, \
    _BINTEGRATEDMAX, _BINTEGRATEDMIN, _BINTEGRATENBINS, \
    _BINTEGRATEDMAX_DWARF, _BINTEGRATEDMIN_DWARF
from compareDataModel import pvlosplate
from plot_internalcomparison import calc_model
import logl
def plot_externalcomparison(parser):
    (options,args)= parser.parse_args()
    if len(args) == 0 or options.plotfilename is None:
        parser.print_help()
        return
    #Read the data
    print "Reading the data ..."
    data= readVclosData(postshutdown=options.postshutdown,
                        fehcut=options.fehcut,
                        cohort=options.cohort,
                        lmin=options.lmin,
                        bmax=options.bmax,
                        ak=True,
                        cutmultiples=options.cutmultiples,
                        validfeh=options.indivfeh, #if indivfeh, we need validfeh
                        jkmax=options.jkmax,
                        datafilename=options.fakedata)
    #data= data[0:20]
    #HACK
    indx= (data['J0MAG']-data['K0MAG'] < 0.5)
    data['J0MAG'][indx]= 0.5+data['K0MAG'][indx]
    #Cut outliers
    #data= data[(data['VHELIO'] < 200.)*(data['VHELIO'] > -200.)]
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
        if options.indivfeh:
            #Load all isochrones
            iso= []
            zs= numpy.arange(0.0005,0.03005,0.0005)
            for ii in range(len(zs)):
                iso.append(isomodel.isomodel(imfmodel=options.imfmodel,
                                             expsfh=options.expsfh,
                                             Z=zs[ii]))
        elif options.varfeh:
            locs= list(set(data['LOCATION']))
            iso= []
            for ii in range(len(locs)):
                indx= (data['LOCATION'] == locs[ii])
                locl= numpy.mean(data['GLON'][indx]*_DEGTORAD)
                iso.append(isomodel.isomodel(imfmodel=options.imfmodel,
                                             expsfh=options.expsfh,
                                             marginalizefeh=True,
                                             glon=locl))
        else:
            iso= isomodel.isomodel(imfmodel=options.imfmodel,Z=options.Z,
                                   expsfh=options.expsfh)
        if options.dwarf:
            iso= [iso, 
                  isomodel.isomodel(imfmodel=options.imfmodel,Z=options.Z,
                                    dwarf=True,expsfh=options.expsfh)]
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
    print "Pre-calculating isochrone distance prior ..."
    logpiso= numpy.zeros((len(data),_BINTEGRATENBINS))
    ds= numpy.linspace(_BINTEGRATEDMIN,_BINTEGRATEDMAX,
                       _BINTEGRATENBINS)
    dm= _dm(ds)
    for ii in range(len(data)):
        mh= data['H0MAG'][ii]-dm
        if options.indivfeh:
            #Find closest Z
            thisZ= isodist.FEH2Z(data[ii]['FEH'])
            indx= numpy.argmin((thisZ-zs))
            logpiso[ii,:]= iso[0][indx](numpy.zeros(_BINTEGRATENBINS)+(data['J0MAG']-data['K0MAG'])[ii],mh)
        elif options.varfeh:
            #Find correct iso
            indx= (locl == data[ii]['LOCATION'])
            logpiso[ii,:]= iso[0][indx](numpy.zeros(_BINTEGRATENBINS)+(data['J0MAG']-data['K0MAG'])[ii],mh)
        else:
            logpiso[ii,:]= iso[0](numpy.zeros(_BINTEGRATENBINS)
                                  +(data['J0MAG']-data['K0MAG'])[ii],mh)
    if options.dwarf:
        logpisodwarf= numpy.zeros((len(data),_BINTEGRATENBINS))
        dwarfds= numpy.linspace(_BINTEGRATEDMIN_DWARF,_BINTEGRATEDMAX_DWARF,
                                    _BINTEGRATENBINS)
        dm= _dm(dwarfds)
        for ii in range(len(data)):
            mh= data['H0MAG'][ii]-dm
            logpisodwarf[ii,:]= iso[1](numpy.zeros(_BINTEGRATENBINS)
                                       +(data['J0MAG']-data['K0MAG'])[ii],mh)
    else:
        logpisodwarf= None
    #Calculate data means etc.
    #Calculate means
    locations= list(set(data['LOCATION']))
    nlocs= len(locations)
    l_plate= numpy.zeros(nlocs)
    avg_plate= numpy.zeros(nlocs)
    sig_plate= numpy.zeros(nlocs)
    siga_plate= numpy.zeros(nlocs)
    for ii in range(nlocs):
        indx= (data['LOCATION'] == locations[ii])
        l_plate[ii]= numpy.mean(data['GLON'][indx])
        avg_plate[ii]= numpy.mean(data['VHELIO'][indx])
        sig_plate[ii]= numpy.std(data['VHELIO'][indx])
        siga_plate[ii]= numpy.std(data['VHELIO'][indx])/numpy.sqrt(numpy.sum(indx))
    #Calculate plate means and variances from the model
    #Load initial parameters from file
    savefile= open(args[0],'rb')
    params= pickle.load(savefile)
    if not options.index is None:
        params= params[options.index]
    savefile.close()
    #First calculate fiducial model
    if not options.dwarf:
        logpisodwarf= None
    fid_logl= numpy.sum(logl.logl(init=params,data=data,options=copy.copy(options))) #Copy options, since logl sets multi=1
    avg_plate_model_fid= calc_model(params,options,data,
                                logpiso,logpisodwarf,
                                df,nlocs,locations,iso)
    #Plot everything
    bovy_plot.bovy_print(fig_height=6.,fig_width=7.)
    dx= 0.8/5.
    left, bottom, width, height= 0.1, 0.9-dx, 0.8, dx
    axTop= pyplot.axes([left,bottom,width,height])
    allaxes= [axTop]
    fig= pyplot.gcf()
    fig.sca(axTop)
    bovy_plot.bovy_plot([0.,360.],[0.,0.],'-',color='0.5',overplot=True)
    bovy_plot.bovy_plot(l_plate,
                        avg_plate-avg_plate_model_fid,
                        'ko',overplot=True)
    pyplot.errorbar(l_plate,avg_plate-avg_plate_model_fid,
                    yerr=siga_plate,marker='o',color='k',
                    linestyle='none',elinestyle='-')
    bovy_plot.bovy_text(r'$\mathrm{fiducial}$',top_right=True,size=14.)
    thisax= pyplot.gca()
    thisax.set_ylim(-19.5,19.5)
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    nullfmt   = NullFormatter()         # no labels
    axTop.xaxis.set_major_formatter(nullfmt)
    #pyplot.xlabel(r'$\mathrm{Galactic\ longitude}\ [\mathrm{deg}]$')
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    #Second is powerlaw with beta=0.15
    #fid_slope= params[5-options.nooutliermean+options.dwarf]
    #params[5-options.nooutliermean+options.dwarf]= 0.
    fid_params= copy.copy(params)
    if options.rotcurve.lower() == 'flat':
        params= list(params)
        params.insert(6,0.15)
        params= numpy.array(params)
        fid_rotcurve= options.rotcurve
        options.rotcurve= 'powerlaw'
    elif options.rotcurve.lower() == 'powerlaw':
        params[6]= 0.15
        fid_rotcurve= options.rotcurve
    avg_plate_model= calc_model(params,options,data,
                                logpiso,logpisodwarf,
                                df,nlocs,locations,iso)
    #params[5-options.nooutliermean+options.dwarf]= fid_slope
    #Undo changes
    params= fid_params
    options.rotcurve= fid_rotcurve
    left, bottom, width, height= 0.1, 0.9-2.*dx, 0.8, dx
    thisax= pyplot.axes([left,bottom,width,height])
    allaxes.append(thisax)
    fig.sca(thisax)
    bovy_plot.bovy_plot([0.,360.],[0.,0.],'-',color='0.5',overplot=True,
                        zorder=-1)
    bovy_plot.bovy_plot(l_plate,
                        avg_plate-avg_plate_model_fid,
                        'o',overplot=True,color='0.6')
    pyplot.errorbar(l_plate,avg_plate-avg_plate_model_fid,
                    yerr=siga_plate,marker='o',color='0.6',
                    linestyle='none',elinestyle='-')
    bovy_plot.bovy_plot(l_plate,
                        avg_plate-avg_plate_model,
                        'o',overplot=True,color='k')
    pyplot.errorbar(l_plate,avg_plate-avg_plate_model,
                    yerr=siga_plate,marker='o',color='k',
                    linestyle='none',elinestyle='-')
    bovy_plot.bovy_text(r'$\mathrm{power\!\!-\!\!law\ rotation\ curve}, \beta = 0.15$',
                        top_right=True,size=14.)
    #pyplot.ylabel(r'$\langle V_{\mathrm{los}}\rangle_{\mathrm{data}}-\langle V_{\mathrm{los}}\rangle_{\mathrm{model}}$')
    thisax.set_ylim(-19.5,19.5)
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    thisax.xaxis.set_major_formatter(nullfmt)
    #pyplot.xlabel(r'$\mathrm{Galactic\ longitude}\ [\mathrm{deg}]$')
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    #Third is powerlaw with beta=-0.15
    #fid_slope= params[5-options.nooutliermean+options.dwarf]
    #params[5-options.nooutliermean+options.dwarf]= 0.
    fid_params= copy.copy(params)
    if options.rotcurve.lower() == 'flat':
        params= list(params)
        params.insert(6,-0.15)
        params= numpy.array(params)
        fid_rotcurve= options.rotcurve
        options.rotcurve= 'powerlaw'
    elif options.rotcurve.lower() == 'powerlaw':
        params[6]= -0.15
        fid_rotcurve= options.rotcurve
    avg_plate_model= calc_model(params,options,data,
                                logpiso,logpisodwarf,
                                df,nlocs,locations,iso)
    #params[5-options.nooutliermean+options.dwarf]= fid_slope
    #Undo changes
    params= fid_params
    options.rotcurve= fid_rotcurve
    left, bottom, width, height= 0.1, 0.9-3.*dx, 0.8, dx
    thisax= pyplot.axes([left,bottom,width,height])
    allaxes.append(thisax)
    fig.sca(thisax)
    bovy_plot.bovy_plot([0.,360.],[0.,0.],'-',color='0.5',overplot=True,
                        zorder=-1)
    bovy_plot.bovy_plot(l_plate,
                        avg_plate-avg_plate_model_fid,
                        'o',overplot=True,color='0.6')
    pyplot.errorbar(l_plate,avg_plate-avg_plate_model_fid,
                    yerr=siga_plate,marker='o',color='0.6',
                    linestyle='none',elinestyle='-')
    bovy_plot.bovy_plot(l_plate,
                        avg_plate-avg_plate_model,
                        'o',overplot=True,color='k')
    pyplot.errorbar(l_plate,avg_plate-avg_plate_model,
                    yerr=siga_plate,marker='o',color='k',
                    linestyle='none',elinestyle='-')
    bovy_plot.bovy_text(r'$\mathrm{power\!\!-\!\!law\ rotation\ curve}, \beta = -0.15$',
                        top_right=True,size=14.)
    pyplot.ylabel(r'$\bar{V}_{\mathrm{data}}-\bar{V}_{\mathrm{model}}$')
    thisax.set_ylim(-19.5,19.5)
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    thisax.xaxis.set_major_formatter(nullfmt)
    #pyplot.xlabel(r'$\mathrm{Galactic\ longitude}\ [\mathrm{deg}]$')
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    #Second is vc=250
    #Load initial parameters from file
    savefile= open(options.vo250file,'rb')
    params= pickle.load(savefile)
    savefile.close()
    options.fixvo= 250.
    options.fixro= 8.4
    vo250_logl= numpy.sum(logl.logl(init=params,data=data,options=copy.copy(options)))
    avg_plate_modelr84= calc_model(params,options,data,
                                logpiso,logpisodwarf,
                                df,nlocs,locations,iso)
    #Load initial parameters from file
    savefile= open(options.vo250r8file,'rb')
    paramsr8= pickle.load(savefile)
    savefile.close()
    options.fixvo= 250.
    options.fixro= 8.
    vo250r8_logl= numpy.sum(logl.logl(init=paramsr8,data=data,options=copy.copy(options)))
    avg_plate_modelr8= calc_model(paramsr8,options,data,
                                logpiso,logpisodwarf,
                                df,nlocs,locations,iso)
    options.fixro= None
    options.fixvo= None
    left, bottom, width, height= 0.1, 0.9-4.*dx, 0.8, dx
    thisax= pyplot.axes([left,bottom,width,height])
    allaxes.append(thisax)
    fig.sca(thisax)
    bovy_plot.bovy_plot([0.,360.],[0.,0.],'-',color='0.5',overplot=True,
                        zorder=-1)
    bovy_plot.bovy_plot(l_plate,
                        avg_plate-avg_plate_model_fid,
                        'o',overplot=True,color='0.6')
    pyplot.errorbar(l_plate,avg_plate-avg_plate_model_fid,
                    yerr=siga_plate,marker='o',color='0.6',
                    linestyle='none',elinestyle='-')
    bovy_plot.bovy_plot(l_plate,
                        avg_plate-avg_plate_modelr84,
                        'x',overplot=True,color='k')
    pyplot.errorbar(l_plate,avg_plate-avg_plate_modelr84,
                    yerr=siga_plate,marker='x',color='k',
                    linestyle='none',elinestyle='-')
    bovy_plot.bovy_plot(l_plate,
                        avg_plate-avg_plate_modelr8,
                        'o',overplot=True,color='k')
    pyplot.errorbar(l_plate,avg_plate-avg_plate_modelr8,
                    yerr=siga_plate,marker='o',color='k',
                    linestyle='none',elinestyle='-')
    #top_right=True,size=14.)
    bovy_plot.bovy_text(r'$\mathrm{best\!\!-\!\!fit}\ V_c(R_0) = 250\ \mathrm{km\ s}^{-1}, R_0 = %.1f\,\mathrm{kpc}, \Delta \chi^2 = %.0f\ (%.1f\sigma)$' % (8.,2.*(fid_logl-vo250r8_logl),special.erfinv(stats.chi2.cdf(2.*(fid_logl-vo250r8_logl),10.))*numpy.sqrt(2.)),
                        top_right=True,size=13.)
    bovy_plot.bovy_text(r'$\mathrm{best\!\!-\!\!fit}\ V_c(R_0) = 250\ \mathrm{km\ s}^{-1}, R_0 = %.1f\,\mathrm{kpc}, \Delta \chi^2 = %.0f\ (%.1f\sigma)$' % (8.4,2.*(fid_logl-vo250_logl),special.erfinv(stats.chi2.cdf(2.*(fid_logl-vo250_logl),10.))*numpy.sqrt(2.)),
                        bottom_right=True,size=13.)
    bovy_plot.bovy_plot([40.],[-15.],'kx',overplot=True)
    bovy_plot.bovy_plot([40.],[12.],'ko',overplot=True)
    thisax.set_ylim(-19.5,19.5)
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    thisax.xaxis.set_major_formatter(nullfmt)
    #pyplot.xlabel(r'$\mathrm{Galactic\ longitude}\ [\mathrm{deg}]$')
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    #Fourth = vpec=SBD10
    #Load sbd fit
    #Load initial parameters from file
    savefile= open(options.sbdfile,'rb')
    params= pickle.load(savefile)
    savefile.close()
    options.sbdvpec= True
    options.fitvpec= False
    sbd_logl= numpy.nansum(logl.logl(init=params,data=data,options=copy.copy(options)))
    avg_plate_model= calc_model(params,options,data,
                                logpiso,logpisodwarf,
                                df,nlocs,locations,iso)
    left, bottom, width, height= 0.1, 0.9-5.*dx, 0.8, dx
    thisax= pyplot.axes([left,bottom,width,height])
    allaxes.append(thisax)
    fig.sca(thisax)
    bovy_plot.bovy_plot([0.,360.],[0.,0.],'-',color='0.5',overplot=True,zorder=-1)
    bovy_plot.bovy_plot(l_plate,
                        avg_plate-avg_plate_model_fid,
                        'o',overplot=True,color='0.6')
    pyplot.errorbar(l_plate,avg_plate-avg_plate_model_fid,
                    yerr=siga_plate,marker='o',color='0.6',
                    linestyle='none',elinestyle='-')
    bovy_plot.bovy_plot(l_plate,
                        avg_plate-avg_plate_model,
                        'o',overplot=True,color='k')
    pyplot.errorbar(l_plate,avg_plate-avg_plate_model,
                    yerr=siga_plate,marker='o',color='k',
                    linestyle='none',elinestyle='-')
    bovy_plot.bovy_text(r'$\mathrm{best\!\!-\!\!fit}\ \vec{V}_\odot = \vec{V}_c(R_0) + \mathrm{SBD10}, \Delta \chi^2 = %.0f$' % (2.*(fid_logl-sbd_logl)),#,special.erfinv(stats.chi2.cdf(2.*(fid_logl-sbd_logl),10.))*numpy.sqrt(2.)),
                        top_right=True,size=14.)
    #pyplot.ylabel(r'$\langle V_{\mathrm{los}}\rangle_{\mathrm{data}}-\langle V_{\mathrm{los}}\rangle_{\mathrm{model}}$')
    thisax.set_ylim(-19.5,19.5)
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    #thisax.xaxis.set_major_formatter(nullfmt)
    pyplot.xlabel(r'$\mathrm{Galactic\ longitude}\ [\mathrm{deg}]$')
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    #Save
    bovy_plot.bovy_end_print(options.plotfilename)
    print fid_logl, sbd_logl, 2.*(fid_logl-sbd_logl),special.erfinv(stats.chi2.cdf(2.*(fid_logl-sbd_logl),10.))*numpy.sqrt(2.)
    return None
        
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
                      help="If set, fit for the peculiar velocity of the Sun as well, CURRENTLY ASSUMES flat rotation curve")
    parser.add_option("--sbdvpec",action="store_true", dest="sbdvpec",
                      default=False,
                      help="If set, use the SBD10 value for the Solar motion + vo as vpec")
    #Velocity distribution model
    parser.add_option("--dfmodel",dest='dfmodel',default='simplegaussian',
                      help="DF model to use")
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
    parser.add_option("--fitdl",action="store_true", dest="fitdl",
                      default=False,
                      help="If set, fit for an offset in l for the GC")
    parser.add_option("--fitfehinnerouter",
                      action="store_true", dest="fitfehinnerouter",
                      default=False,
                      help="If set, fit for the feh offset separately for the inner disk")
    parser.add_option("--fitsratio",action="store_true", dest="fitsratio",
                      default=False,
                      help="If set, fit for the ration squared of tangential to radial dispersion")
    parser.add_option("--fitsratioinnerouter",
                      action="store_true", dest="fitsratioinnerouter",
                      default=False,
                      help="If set, fit for the ration squared of tangential to radial dispersion")
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
    parser.add_option("--downsample",dest='downsample',default=None,
                      type='float',
                      help="Factor with which to downsample the data")
    parser.add_option("--cutmultiples",action="store_true", 
                      dest="cutmultiples",
                      default=False,
                      help="readVclosData 'cutmultiples'")
    parser.add_option("-f",dest='fakedata',default=None,
                      help="Name of the fake data filename")
    #Isochrone IMF
    parser.add_option("--imfmodel",dest='imfmodel',default='lognormalChabrier2001',
                      help="imfmodel for isochrone model")
    parser.add_option("--Z",dest='Z',default=.019,type='float',
                      help="Metallicity of isochrone")
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
    parser.add_option("--varfeh",action="store_false", dest="varfeh",
                      default=True,
                      help="If set, don't use a varying [Fe/H] distribution as a function of l")
    parser.add_option("--indivfeh",action="store_false", dest="indivfeh",
                      default=True,
                      help="If set, use distances calculated based on the individual [Fe/H] of the objects")
    parser.add_option("--isofile",dest="isofile",default=None,
                      help="if set, store or restore the isochrone model(s) in this file")
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
    #Comparison options
    parser.add_option("--nvlos",dest='nvlos',default=11,type='int',
                      help="number of vlos")
    #Multiprocessing?
    parser.add_option("-m","--multi",dest='multi',default=None,type='int',
                      help="number of cpus to use for sampling")
    #Output file
    parser.add_option("-o",dest='plotfilename',default=None,
                      help="Name of the file for the plot")
    #
    parser.add_option("-i","--indx",dest='index',default=None,type='int',
                      help="If samples are given, use this index")
    #File with fits for vc=250
    parser.add_option("--vo250file",dest='vo250file',default=None,
                      help="Name of the file that has the vc=250. fit")
    #File with fits for vc=250, ro=8.
    parser.add_option("--vo250r8file",dest='vo250r8file',default=None,
                      help="Name of the file that has the vc=250., ro=8. fit")
    #File with fits for sbdvpec
    parser.add_option("--sbdfile",dest='sbdfile',default=None,
                      help="Name of the file that has the sbdvpec fit")
    return parser

if __name__ == '__main__':
    numpy.random.seed(1) #We need to seed to get, e.g., the same permutation when downsampling
    plot_externalcomparison(get_options())

