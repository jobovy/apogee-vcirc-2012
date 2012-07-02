import sys
import os, os.path
import cPickle as pickle
from optparse import OptionParser
import math
import numpy
from scipy.maxentropy import logsumexp
from galpy.util import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
import multi
import multiprocessing
import isomodel
import isodist
from readVclosData import readVclosData
from fitvc import mloglike, _dm, \
    _DEGTORAD, _REFV0, _REFR0, \
    _BINTEGRATEDMAX, _BINTEGRATEDMIN, _BINTEGRATENBINS, \
    _BINTEGRATEDMAX_DWARF, _BINTEGRATEDMIN_DWARF
from compareDataModel import pvlosplate
from plot_bestfit import get_options, bootstrap_sigerr
import logl
def plot_chi2(parser):
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
    #HACK
    indx= (data['J0MAG']-data['K0MAG'] < 0.5)
    data['J0MAG'][indx]= 0.5+data['K0MAG'][indx]
    #Cut inner disk locations
    #data= data[(data['GLON'] > 75.)]
    #Cut outliers
    #data= data[(data['VHELIO'] < 200.)*(data['VHELIO'] > -200.)]
    print "Using %i data points ..." % len(data)
    #Set up the isochrone
    if not options.isofile is None and os.path.exists(options.isofile):
        print "Loading the isochrone model ..."
        isofile= open(options.isofile,'rb')
        iso= pickle.load(isofile)
        if options.indivfeh:
            zs= pickle.load(isofile)
        if options.varfeh:
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
    #Load initial parameters from file
    savefile= open(args[0],'rb')
    params= pickle.load(savefile)
    if not options.index is None:
        params= params[options.index]
    savefile.close()
    #params[0]= 245./235.
    #params[1]= 8.5/8.
    #Calculate data means etc.
    #Calculate means
    locations= list(set(data['LOCATION']))
    nlocs= len(locations)
    l_plate= numpy.zeros(nlocs)
    avg_plate= numpy.zeros(nlocs)
    sig_plate= numpy.zeros(nlocs)
    siga_plate= numpy.zeros(nlocs)
    sigerr_plate= numpy.zeros(nlocs)
    fidlogl= logl.logl(init=params,data=data,options=options)
    logl_plate= numpy.zeros(nlocs)
    for ii in range(nlocs):
        indx= (data['LOCATION'] == locations[ii])
        l_plate[ii]= numpy.mean(data['GLON'][indx])
        avg_plate[ii]= numpy.mean(data['VHELIO'][indx])
        sig_plate[ii]= numpy.std(data['VHELIO'][indx])
        siga_plate[ii]= numpy.std(data['VHELIO'][indx])/numpy.sqrt(numpy.sum(indx))
        sigerr_plate[ii]= bootstrap_sigerr(data['VHELIO'][indx])
        #Logl
        logl_plate[ii]= -2.*(numpy.sum(fidlogl[indx])-numpy.sum(fidlogl)/len(indx)*numpy.sum(indx))
    #Calculate plate means and variances from the model
    avg_plate_model= numpy.zeros(nlocs)
    sig_plate_model= numpy.zeros(nlocs)
    for ii in range(nlocs):
        #Calculate vlos | los
        indx= (data['LOCATION'] == locations[ii])
        thesedata= data[indx]
        thislogpiso= logpiso[indx,:]
        if options.dwarf:
            thislogpisodwarf= logpisodwarf[indx,:]
        else:
            thislogpisodwarf= None
        vlos= numpy.linspace(-200.,200.,options.nvlos)
        pvlos= numpy.zeros(options.nvlos)
        if not options.multi is None:
            pvlos= multi.parallel_map((lambda x: pvlosplate(params,vlos[x],
                                                            thesedata,
                                                            df,options,
                                                            thislogpiso,
                                                            thislogpisodwarf,iso)),
                                      range(options.nvlos),
                                      numcores=numpy.amin([len(vlos),multiprocessing.cpu_count(),options.multi]))
        else:
            for jj in range(options.nvlos):
                print jj
                pvlos[jj]= pvlosplate(params,vlos[jj],thesedata,df,options,
                                      thislogpiso,thislogpisodwarf,iso)
        pvlos-= logsumexp(pvlos)
        pvlos= numpy.exp(pvlos)
        #Calculate mean and velocity dispersion
        avg_plate_model[ii]= numpy.sum(vlos*pvlos)
        sig_plate_model[ii]= numpy.sqrt(numpy.sum(vlos**2.*pvlos)\
                                            -avg_plate_model[ii]**2.)
    #Plot everything
    left, bottom, width, height= 0.1, 0.4, 0.8, 0.3
    axTop= pyplot.axes([left,bottom,width,height])
    left, bottom, width, height= 0.1, 0.1, 0.8, 0.3
    axChi2= pyplot.axes([left,bottom,width,height])
    #left, bottom, width, height= 0.1, 0.1, 0.8, 0.2
    #axSig= pyplot.axes([left,bottom,width,height])
    fig= pyplot.gcf()
    #Plot the difference
    fig.sca(axTop)
    bovy_plot.bovy_plot([0.,360.],[0.,0.],'-',color='0.5',overplot=True)
    bovy_plot.bovy_plot(l_plate,
                        avg_plate-avg_plate_model,
                        'ko',overplot=True)
    pyplot.errorbar(l_plate,avg_plate-avg_plate_model,
                    yerr=siga_plate,marker='o',color='k',linestyle='none',elinestyle='-')
    pyplot.ylabel(r'$\langle v_{\mathrm{los}}\rangle_{\mathrm{data}}-\langle v_{\mathrm{los}}\rangle_{\mathrm{model}}$')
    pyplot.ylim(-14.5,14.5)
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    nullfmt   = NullFormatter()         # no labels
    axTop.xaxis.set_major_formatter(nullfmt)
    #pyplot.xlabel(r'$\mathrm{Galactic\ longitude}\ [\mathrm{deg}]$')
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    #Plot the chi2
    fig.sca(axChi2)
    bovy_plot.bovy_plot([0.,360.],[0.,0.],'-',color='0.5',overplot=True)
    bovy_plot.bovy_plot(l_plate,
                        logl_plate,
                        'ko',overplot=True)
    pyplot.ylabel(r'$\Delta \chi^2$')
    #pyplot.ylim(numpy.amin(logl_plate),numpy.amax(logl_plate))
    pyplot.ylim(-150.,150.)
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    pyplot.xlabel(r'$\mathrm{Galactic\ longitude}\ [\mathrm{deg}]$')
    pyplot.xlim(0.,360.)
    bovy_plot._add_ticks()
    #Save
    bovy_plot.bovy_end_print(options.plotfilename)
    return None

if __name__ == '__main__':
    numpy.random.seed(1) #We need to seed to get, e.g., the same permutation when downsampling
    plot_chi2(get_options())

