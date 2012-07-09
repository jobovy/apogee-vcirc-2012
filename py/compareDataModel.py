###################################################################################
#  compareDataModel.py: module for comparing the data to the model
###################################################################################
import os, os.path
import copy
from optparse import OptionParser
import cPickle as pickle
import numpy
from scipy.maxentropy import logsumexp
from galpy.util import bovy_plot
import multi
import multiprocessing
from fitvc import mloglike, _dm, \
    _DEGTORAD, _REFV0, _REFR0, \
    _BINTEGRATEDMAX, _BINTEGRATEDMIN, _BINTEGRATENBINS, \
    _BINTEGRATEDMAX_DWARF, _BINTEGRATEDMIN_DWARF
from readVclosData import readVclosData
import isomodel
import isodist
_PLOTZERO= False
_PLOTALLJHK= False
def pvlosplate(params,vhelio,data,df,options,logpiso,logpisodwarf,iso):
    """
    NAME:
       pvlosplate
    PURPOSE:
       calculate the vlos probability for a given location
    INPUT:
       params - parameters of the model
       vhelio - heliocentric los velocity to evaluate
       data - data array for this location
       df - df object(s) (?)
       options - options
       logpiso, logpisodwarf - precalculated isochrones
    OUTPUT:
       log of the probability
    HISTORY:
       2012-02-20 - Written - Bovy (IAS)
    """
    #Output is sum over data l,b,jk,h
    l= data['GLON']*_DEGTORAD
    b= data['GLAT']*_DEGTORAD
    sinl= numpy.sin(l)
    cosl= numpy.cos(l)
    sinb= numpy.sin(b)
    cosb= numpy.cos(b)
    jk= data['J0MAG']-data['K0MAG']
    try:
        jk[(jk < 0.5)]= 0.5 #BOVY: FIX THIS HACK BY EMAILING GAIL
    except TypeError:
        pass #HACK
    h= data['H0MAG']
    options.multi= 1 #To avoid conflict
    out= -mloglike(params,numpy.zeros(len(data))+vhelio,
                   l,
                   b,
                   jk,
                   h,
                   df,options,
                   sinl,
                   cosl,
                   cosb,
                   sinb,
                   logpiso,
                   logpisodwarf,True,None,iso,data['FEH']) #None iso for now
    #indx= (out >= -0.1)*(out <= 0.1)
    #print out[indx], jk[indx], h[indx]
    return logsumexp(out)

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the fit/samples will be saved to"
    parser = OptionParser(usage=usage)
    #Initial conditions file
    parser.add_option("--init",dest='init',default=None,
                      help="Initial parameters")
    #Rotation curve parameters/model
    parser.add_option("--rotcurve",dest='rotcurve',default='flat',
                      help="Rotation curve model to fit")
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
    parser.add_option("--fitsratio",action="store_true", dest="fitsratio",
                      default=False,
                      help="If set, fit for the ration squared of tangential to radial dispersion")
    parser.add_option("--fitsratioinnerouter",
                      action="store_true", dest="fitsratioinnerouter",
                      default=False,
                      help="If set, fit for the ration squared of tangential to radial dispersion")
    parser.add_option("--fitdl",action="store_true", dest="fitdl",
                      default=False,
                      help="If set, fit for an offset in l for the GC")
    parser.add_option("--fitfehinnerouter",
                      action="store_true", dest="fitfehinnerouter",
                      default=False,
                      help="If set, fit for the feh offset separately for the inner disk")
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
    parser.add_option("--location",dest='location',default=4318,type='int',
                      help="location id when looking at single los (if 0, all los)")
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
    #Comparison options
    parser.add_option("--plottype",dest='plottype',default='pvloslos',
                       help="Type of plot to make")
    parser.add_option("-o","--plotfile",dest='plotfile',
                      help="name of file for plot")
    parser.add_option("--nvlos",dest='nvlos',default=11,type='int',
                      help="number of vlos")
    #Multiprocessing?
    parser.add_option("-m","--multi",dest='multi',default=None,type='int',
                      help="number of cpus to use")
    parser.add_option("--seed",dest='seed',default=1,type='int',
                      help="seed for random number generator")
    #
    parser.add_option("-i","--indx",dest='index',default=None,type='int',
                      help="If samples are given, use this index")
    #A_K quantiles
    parser.add_option("--aklow",dest='aklow',type='float',default=None,
                      help="lower AK quantile of using AK-quantiles")
    parser.add_option("--akhigh",dest='akhigh',type='float',default=None,
                      help="upper AK quantile of using AK-quantiles")
    return parser

if __name__ == '__main__':
    #Get options
    parser= get_options()
    options,args= parser.parse_args()
    #Read data
    #AK quantiles?
    if not options.aklow is None:
        akquantiles= [options.aklow,options.akhigh]
    else:
        akquantiles= None
    data= readVclosData(postshutdown=options.postshutdown,
                        fehcut=options.fehcut,
                        cohort=options.cohort,
                        lmin=options.lmin,
                        bmax=options.bmax,
                        ak=True,
                        jkmax=options.jkmax,
                        datafilename=options.fakedata,
                        validfeh=options.indivfeh, #if indivfeh, we need validfeh
                        akquantiles=akquantiles)
    if options.location == 0:
        locations= list(set(data['LOCATION']))
    else:
        locations= [options.location]
    #HACK
    indx= (data['J0MAG']-data['K0MAG'] < 0.5)
    data['J0MAG'][indx]= 0.5+data['K0MAG'][indx]
    #Fit parameters
    if options.dwarf:
        params= [270./_REFV0,8./_REFR0,numpy.log(35./_REFV0),0.025,0.5,1.]
    else:
        params= [270./_REFV0,8./_REFR0,numpy.log(35./_REFV0),0.025,0.5]
    if not options.init is None:
        #Load initial parameters from file
        savefile= open(options.init,'rb')
        params= pickle.load(savefile)
        if not options.index is None:
            params= params[options.index]
        savefile.close()
    #params[2]= -2#numpy.log(20./235.)
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
            iso= isomodel.isomodel(imfmodel=options.imfmodel,Z=options.Z,expsfh=options.expsfh)
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
        if options.fitah:
            ah= params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter]
            if options.fitahinnerouter and data[ii]['GLON'] < 35.: ah= params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitah+options.fitdm]
            mh= data['H0MAG'][ii]-dm+ah
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
                                      +(data['J0MAG']-data['K0MAG'])[ii]
                                      +1.5/1.55*ah,mh)
        elif options.fitdm:
            if not options.fitdminnerouter or data[ii]['GLON'] >= 35.: xtra= params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter]
            if options.fitdminnerouter and data[ii]['GLON'] < 35.: xtra= params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitdm+options.fitah]
            mh= data[ii]['H0MAG']-dm+xtra
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
        else:
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
            if options.fitah:
                ah= params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter]
                if options.fitahinnerouter and data[ii]['GLON'] < 35.: ah= params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitah+options.fitdm]
                mh= data['H0MAG'][ii]-dm+ah
                logpisodwarf[ii,:]= iso[1](numpy.zeros(_BINTEGRATENBINS)
                                           +(data['J0MAG']-data['K0MAG'])[ii]
                                           +1.5/1.55*ah,mh)
            elif options.fitdm:
                pass
            else:
                mh= data['H0MAG'][ii]-dm
                logpisodwarf[ii,:]= iso[1](numpy.zeros(_BINTEGRATENBINS)
                                           +(data['J0MAG']-data['K0MAG'])[ii],mh)
    else:
        logpisodwarf= None
    if options.fitah:
        params= list(params)
        params.pop(5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter)
        options.fitah= False
        if options.fitahinnerouter:
            params.pop(5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitah+options.fitdm)
            options.fitahinnerouter= False
        params= numpy.array(params)
    if options.fitdm:
        params= list(params)
        params.pop(5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter)
        options.fitdm= False
        if options.fitdminnerouter:
            params.pop(5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitah+options.fitdm)
            options.fitdminnerouter= False
        params= numpy.array(params)
    #test
    if options.plottype.lower() == 'pvloslos':
        for location in locations:
            indx= (data['LOCATION'] == location)
            thesedata= data[indx]
            thislogpiso= logpiso[indx,:]
            if options.dwarf:
                thislogpisodwarf= logpisodwarf[indx,:]
            else:
                thislogpisodwarf= None
            #Calculate vlos | los
            vlos= numpy.linspace(-200.,200.,options.nvlos)
            pvlos= numpy.zeros(options.nvlos)
            if not options.multi is None:
                pvlos= multi.parallel_map((lambda x: pvlosplate(params,vlos[x],
                                                                thesedata,df,options,
                                                                thislogpiso,
                                                                thislogpisodwarf,iso)),
                                          range(options.nvlos),
                                          numcores=numpy.amin([len(vlos),multiprocessing.cpu_count(),options.multi]))
            else:
                for ii in range(options.nvlos):
                    print ii
                    pvlos[ii]= pvlosplate(params,vlos[ii],thesedata,df,options,
                                          thislogpiso,thislogpisodwarf,iso)
            pvlos-= logsumexp(pvlos)
            pvlos= numpy.exp(pvlos)
            if _PLOTZERO:
                pvloszero= numpy.zeros(options.nvlos)
                params[2]= -3.8
                if not options.multi is None:
                    pvloszero= multi.parallel_map((lambda x: pvlosplate(params,vlos[x],
                                                                        thesedata,df,options,
                                                                        thislogpiso,
                                                                        thislogpisodwarf,iso)),
                                                  range(options.nvlos),
                                                  numcores=numpy.amin([len(vlos),multiprocessing.cpu_count(),options.multi]))
                else:
                    for ii in range(options.nvlos):
                        print ii
                        pvloszero[ii]= pvlosplate(params,vlos[ii],thesedata,df,options,
                                                  thislogpiso,thislogpisodwarf,iso)
                pvloszero-= logsumexp(pvloszero)
                pvloszero= numpy.exp(pvloszero)
            if _PLOTALLJHK:
                pvlosalljhk= numpy.zeros(options.nvlos)
                alljhkdata= copy.copy(data)
                alljhkdata['GLON']= numpy.mean(thesedata['GLON'])
                alljhkdata['GLAT']= numpy.mean(thesedata['GLAT'])
                if not options.multi is None:
                    pvlosalljhk= multi.parallel_map((lambda x: pvlosplate(params,vlos[x],
                                                                          alljhkdata,df,options,
                                                                          logpiso,
                                                                        logpisodwarf,iso)),
                                                  range(options.nvlos),
                                                  numcores=numpy.amin([len(vlos),multiprocessing.cpu_count(),options.multi]))
                else:
                    for ii in range(options.nvlos):
                        print ii
                        pvlosalljhk[ii]= pvlosplate(params,vlos[ii],alljhkdata,df,options,
                                                    logpiso,logpisodwarf,iso)
                pvlosalljhk-= logsumexp(pvlosalljhk)
                pvlosalljhk= numpy.exp(pvlosalljhk)
            #Plot data
            bovy_plot.bovy_print()
            hist, xvec, p= bovy_plot.bovy_hist(thesedata['VHELIO'],
                                               range=[-200.,200.],
                                               bins=31,
                                               histtype='step',color='k',
                                               xlabel=r'$\mathrm{heliocentric}\ V_{\mathrm{los}}\ [\mathrm{km\ s}^{-1}]$')
            #Normalize prediction
            data_int= numpy.sum(hist)*(xvec[1]-xvec[0])
            pvlos*= data_int/numpy.sum(pvlos)/(vlos[1]-vlos[0])
            bovy_plot.bovy_plot(vlos,pvlos,'-',color='0.65',overplot=True,lw=2.)
            if _PLOTZERO:
                pvloszero*= data_int/numpy.sum(pvloszero)/(vlos[1]-vlos[0])
                bovy_plot.bovy_plot(vlos,pvloszero,'--',color='0.65',overplot=True,lw=2.)
            if _PLOTALLJHK:
                pvlosalljhk*= data_int/numpy.sum(pvlosalljhk)/(vlos[1]-vlos[0])
                bovy_plot.bovy_plot(vlos,pvlosalljhk,'--',color='0.65',overplot=True,lw=2.)
            #Add text
            bovy_plot.bovy_text(#r'$\mathrm{location}\ =\ %i$' % location
                                #+'\n'
                                r'$l\ \approx\ %.0f^\circ$' % numpy.mean(thesedata['GLON']),
                                top_right=True,size=14.)
            bovy_plot.bovy_end_print(options.plotfile+'_%i.ps' % location)

