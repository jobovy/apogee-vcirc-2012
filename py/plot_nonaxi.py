import sys
import copy
import os, os.path
import cPickle as pickle
from optparse import OptionParser
import math
import numpy
from scipy.maxentropy import logsumexp
from galpy.util import bovy_plot, bovy_coords
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
import multi
import multiprocessing
import isomodel
import isodist
from readVclosData import readVclosData
from fitvc import mloglike, _dm, _logpd, \
    _DEGTORAD, _REFV0, _REFR0, _PMSGRA, \
    _BINTEGRATEDMAX, _BINTEGRATEDMIN, _BINTEGRATENBINS, \
    _BINTEGRATEDMAX_DWARF, _BINTEGRATEDMIN_DWARF
from compareDataModel import pvlosplate
from plot_internalcomparison import calc_model
import marginalize_phasedist
from velocity_field import read_output
def calc_nonaxi(params,barfile,options,data,logpiso,
                nlocs,locations,phio):
    #Average over data and integrate over distance
    surfmass,meanvr,meanvt,sigmar2,sigmat2,sigmart,vertexdev,surfmass_init,meanvt_init,sigmar2_init,sigmat2_init,ii,jj,grid = read_output(barfile)
    xgrid, ygrid= marginalize_phasedist.phiR_input(51,0.5,2.)
    interpObj= marginalize_phasedist.interp_velocity_field(xgrid,ygrid,surfmass,meanvr,meanvt)
    interpObjAxi= marginalize_phasedist.interp_velocity_field(xgrid,ygrid,surfmass_init,0.,meanvt_init)
    avg_plate_model= numpy.zeros(nlocs)
    for ii in range(nlocs):
        #Calculate vlos | los
        indx= (data['LOCATION'] == locations[ii])
        thesedata= data[indx]
        thislogpiso= logpiso[indx,:]
        if not options.multi is None:
            mvlos= multi.parallel_map((lambda x: mvlosnonaxi(params,
                                                             interpObj,
                                                             interpObjAxi,
                                                             thesedata[x],
                                                             options,
                                                             thislogpiso[x,:],
                                                             phio)),
                                      range(len(thesedata)),
                                      numcores=numpy.amin([len(thesedata),multiprocessing.cpu_count(),options.multi]))
        else:
            mvlos= numpy.zeros(len(data))
            for jj in range(len(thesedata)):
                print jj
                mvlos[jj]= mvlosnonaxi(params,
                                       interpObj,
                                       interpObjAxi,
                                       thesedata[jj],
                                       options,
                                       thislogpiso[jj,:],
                                       phio)
        #Calculate mean and velocity dispersion
        avg_plate_model[ii]= numpy.mean(mvlos)
    return avg_plate_model*_REFV0*params[0]

def mvlosnonaxi(params,
                interpObj,
                interpObjAxi,
                thesedata,options,logpiso,
                phio):
    l= thesedata['GLON']*_DEGTORAD
    b= thesedata['GLAT']*_DEGTORAD
    cosb= math.cos(b)
    sinb= math.sin(b)
    out= 0.
    norm= 0.
    ds= numpy.linspace(_BINTEGRATEDMAX,_BINTEGRATEDMIN,_BINTEGRATENBINS)/params[1]/_REFR0
    for ii in range(_BINTEGRATENBINS):
        #Calculate R and phi
        R, phi= bovy_coords.dl_to_rphi_2d(ds[ii],l,
                                          degree=False,phio=phio)
        if R < 0.5 or R > 2.: continue #Not in the model
        logpd= _logpd(params,ds[ii],l,b,0.,0.,None,options,R,phi,cosb,sinb,logpiso[ii])
        out+= math.exp(logpd)*((interpObj.meanvt(R,phi)*numpy.sin(phi+l-phio)-interpObj.meanvr(R,phi)*numpy.cos(phi+l-phio))-(interpObjAxi.meanvt(R,phi)*numpy.sin(phi+l-phio)-interpObjAxi.meanvr(R,phi)*numpy.cos(phi+l-phio)))
        norm+= math.exp(logpd)
    if norm == 0.: return 0.
    out/= norm
    #if l < 35.*_DEGTORAD: print out
    if options.dwarf:
        out*= (1.-params[5-options.nooutliermean])
        R, phi= 1., phio
        out+= params[5-options.nooutliermean]*((interpObj.meanvt(R,phi)*numpy.sin(phi+l-phio)-interpObj.meanvr(R,phi)*numpy.cos(phi+l-phio))-(interpObjAxi.meanvt(R,phi)*numpy.sin(phi+l-phio)-interpObjAxi.meanvr(R,phi)*numpy.cos(phi+l-phio)))
    return out[0][0]
def plot_nonaxi(parser):
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
                        validfeh=options.indivfeh, #if indivfeh, we need validfeh
                        cutmultiples=options.cutmultiples,
                        jkmax=options.jkmax)
    #data= data[0:20]
    #HACK
    indx= (data['J0MAG']-data['K0MAG'] < 0.5)
    data['J0MAG'][indx]= 0.5+data['K0MAG'][indx]
    indx= (data['GLON'] <= 30.)
    data['GLON'][indx]= 30.05 #Hack because the non-axi models don't go below Ro/2
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
            indx= numpy.argmin(numpy.fabs(thisZ-zs))
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
    #clean logpiso
    dataindx= []
    for ii in range(len(data)):
        thislogpiso= logpiso[ii,:]-logsumexp(logpiso[ii,:])
        if numpy.all(thislogpiso == 0.): dataindx.append(False)
        else: dataindx.append(True)
    dataindx= numpy.array(dataindx,dtype='bool')
    data= data[dataindx]
    logpiso= logpiso[dataindx,:]
    logpisodwarf= logpisodwarf[dataindx,:]
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
    savefile.close()
    #First calculate fiducial model
    if not options.dwarf:
        logpisodwarf= None
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
    thisax.set_ylim(-14.5,14.5)
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    nullfmt   = NullFormatter()         # no labels
    axTop.xaxis.set_major_formatter(nullfmt)
    #pyplot.xlabel(r'$\mathrm{Galactic\ longitude}\ [\mathrm{deg}]$')
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    #Second is bar
    print "Reading bar file ..."
    barfile= '/work/bovy/data/bovy/nonaximw/bar/bar_rect_vr_tform-150_tsteady-4pi_51_101.sav'
    avg_plate_bar= calc_nonaxi(params,barfile,options,
                               data,logpiso,
                               nlocs,locations,0.)
    print avg_plate_bar
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
                        avg_plate-avg_plate_model_fid-avg_plate_bar,
                        'o',overplot=True,color='k')
    pyplot.errorbar(l_plate,avg_plate-avg_plate_model_fid-avg_plate_bar,
                    yerr=siga_plate,marker='o',color='k',
                    linestyle='none',elinestyle='-')
    bovy_plot.bovy_text(r'$\mathrm{bar}$',
                        top_right=True,size=14.)
    #pyplot.ylabel(r'$\langle v_{\mathrm{los}}\rangle_{\mathrm{data}}-\langle v_{\mathrm{los}}\rangle_{\mathrm{model}}$')
    thisax.set_ylim(-14.5,14.5)
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    thisax.xaxis.set_major_formatter(nullfmt)
    #pyplot.xlabel(r'$\mathrm{Galactic\ longitude}\ [\mathrm{deg}]$')
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    #Third is spiral
    print "Reading spiral file ..."
    barfile= '/work/bovy/data/bovy/nonaximw/spiral/spiral_rect_vr_51_101.sav'
    avg_plate_bar= calc_nonaxi(params,barfile,options,
                               data,logpiso,
                               nlocs,locations,0.)
    print avg_plate_bar
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
                        avg_plate-avg_plate_model_fid-avg_plate_bar,
                        'o',overplot=True,color='k')
    pyplot.errorbar(l_plate,avg_plate-avg_plate_model_fid-avg_plate_bar,
                    yerr=siga_plate,marker='o',color='k',
                    linestyle='none',elinestyle='-')
    bovy_plot.bovy_text(r'$\mathrm{spiral}$',
                        top_right=True,size=14.)
    pyplot.ylabel(r'$\langle v_{\mathrm{los}}\rangle_{\mathrm{data}}-\langle v_{\mathrm{los}}\rangle_{\mathrm{model}}$')
    thisax.set_ylim(-14.5,14.5)
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    thisax.xaxis.set_major_formatter(nullfmt)
    #pyplot.xlabel(r'$\mathrm{Galactic\ longitude}\ [\mathrm{deg}]$')
    pyplot.xlim(0.,360.)
    #Fourth is first elliptical
    print "Reading elliptical file ..."
    barfile= '/work/bovy/data/bovy/nonaximw/elliptical/el_rect_so_0.2_res_51_grid_101_tform_-150._tsteady_125._cp_0.05_nsigma_4.sav'
    avg_plate_bar= calc_nonaxi(params,barfile,options,
                               data,logpiso,
                               nlocs,locations,0.)
    print avg_plate_bar
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
                        avg_plate-avg_plate_model_fid-avg_plate_bar,
                        'o',overplot=True,color='k')
    pyplot.errorbar(l_plate,avg_plate-avg_plate_model_fid-avg_plate_bar,
                    yerr=siga_plate,marker='o',color='k',
                    linestyle='none',elinestyle='-')
    bovy_plot.bovy_text(r'$\mathrm{elliptical}$',
                        top_right=True,size=14.)
    #pyplot.ylabel(r'$\langle v_{\mathrm{los}}\rangle_{\mathrm{data}}-\langle v_{\mathrm{los}}\rangle_{\mathrm{model}}$')
    thisax.set_ylim(-14.5,14.5)
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    thisax.xaxis.set_major_formatter(nullfmt)
    #pyplot.xlabel(r'$\mathrm{Galactic\ longitude}\ [\mathrm{deg}]$')
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    #Fourth is first elliptical
    print "Reading elliptical file ..."
    barfile= '/work/bovy/data/bovy/nonaximw/elliptical/el_rect_so_0.2_res_51_grid_101_tform_-150._tsteady_125._cp_0.05_nsigma_4.sav'
    avg_plate_bar= calc_nonaxi(params,barfile,options,
                               data,logpiso,
                               nlocs,locations,-45.*_DEGTORAD)
    print avg_plate_bar
    left, bottom, width, height= 0.1, 0.9-5.*dx, 0.8, dx
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
                        avg_plate-avg_plate_model_fid-avg_plate_bar,
                        'o',overplot=True,color='k')
    pyplot.errorbar(l_plate,avg_plate-avg_plate_model_fid-avg_plate_bar,
                    yerr=siga_plate,marker='o',color='k',
                    linestyle='none',elinestyle='-')
    bovy_plot.bovy_text(r'$\mathrm{elliptical}$',
                        top_right=True,size=14.)
    #pyplot.ylabel(r'$\langle v_{\mathrm{los}}\rangle_{\mathrm{data}}-\langle v_{\mathrm{los}}\rangle_{\mathrm{model}}$')
    thisax.set_ylim(-14.5,14.5)
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    #thisax.xaxis.set_major_formatter(nullfmt)
    pyplot.xlabel(r'$\mathrm{Galactic\ longitude}\ [\mathrm{deg}]$')
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    #Save
    bovy_plot.bovy_end_print(options.plotfilename)
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
    parser.add_option("--fitsratio",action="store_true", dest="fitsratio",
                      default=False,
                      help="If set, fit for the ration squared of tangential to radial dispersion")
    parser.add_option("--fitsratioinnerouter",
                      action="store_true", dest="fitsratioinnerouter",
                      default=False,
                      help="If set, fit for the ration squared of tangential to radial dispersion")
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
    parser.add_option("--jkmax",dest='jkmax',default=1.2,type='float',
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
    parser.add_option("--varfeh",action="store_true", dest="varfeh",
                      default=False,
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
    return parser

if __name__ == '__main__':
    numpy.random.seed(1) #We need to seed to get, e.g., the same permutation when downsampling
    plot_nonaxi(get_options())

