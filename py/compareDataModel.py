###################################################################################
#  compareDataModel.py: module for comparing the data to the model
###################################################################################
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
def pvlosplate(params,vhelio,data,df,options,logpiso,logpisodwarf):
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
                   logpisodwarf,True,None)
    #indx= (out == 0.)
    #print l[indx], b[indx], jk[indx], h[indx]
    return logsumexp(out)
    

    out= numpy.zeros(len(data))
    for ii in range(len(out)):
        if options.dwarf:
            pisodwarf= logpisodwarf[ii,:].reshape((1,logpiso.shape[1])),
        else:
            pisodwarf= None
        out[ii]= -mloglike(params,numpy.array([vhelio]),
                           numpy.array([l[ii]]),
                           numpy.array([b[ii]]),
                           numpy.array([jk[ii]]),
                           numpy.array([h[ii]]),
                           df,options,
                           numpy.array([sinl[ii]]),
                           numpy.array([cosl[ii]]),
                           numpy.array([cosb[ii]]),
                           numpy.array([sinb[ii]]),
                           logpiso[ii,:].reshape((1,logpiso.shape[1])),
                           pisodwarf)
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
    #Ro prior
    parser.add_option("--noroprior",action="store_true", dest="noroprior",
                      default=False,
                      help="If set, do not apply an Ro prior")
    #Sun's peculiar velocity
    parser.add_option("--fitvpec",action="store_true", dest="fitvpec",
                      default=False,
                      help="If set, fit for the peculiar velocity of the Sun as well, CURRENTLY ASSUMES flat rotation curve")
    #Velocity distribution model
    parser.add_option("--dfmodel",dest='dfmodel',default='simplegaussian',
                      help="DF model to use")
    #Density model
    parser.add_option("--densmodel",dest='densmodel',default='expdisk',
                      help="Density model to use")
    parser.add_option("--hr",dest='hr',default=3.,type='float',
                      help="scale length in kpc")
    parser.add_option("--hz",dest='hz',default=0.25,type='float',
                      help="scale height in kpc")
    parser.add_option("--hs",dest='hs',default=8.,type='float',
                      help="dispersion scale length in kpc")
    parser.add_option("--fitsratio",action="store_true", dest="fitsratio",
                      default=False,
                      help="If set, fit for the ration squared of tangential to radial dispersion")
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
    parser.add_option("--location",dest='location',default=4318,type='int',
                      help="location id when looking at single los")
    parser.add_option("--cutmultiples",action="store_true", 
                      dest="cutmultiples",
                      default=False,
                      help="readVclosData 'cutmultiples'")
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
    return parser

if __name__ == '__main__':
    #Get options
    parser= get_options()
    options,args= parser.parse_args()
    #Read data
    data= readVclosData(postshutdown=options.postshutdown,
                        fehcut=options.fehcut,
                        cohort=options.cohort,
                        lmin=options.lmin,
                        bmax=options.bmax,
                        ak=True,
                        jkmax=options.jkmax)
    #HACK
    indx= (data['J0MAG']-data['K0MAG'] < 0.5)
    data['J0MAG'][indx]= 0.5+data['K0MAG'][indx]
    #Set up the isochrone
    print "Setting up the isochrone model ..."
    iso= isomodel.isomodel(imfmodel=options.imfmodel,Z=options.Z)
    if options.dwarf:
        iso= [iso, 
              isomodel.isomodel(imfmodel=options.imfmodel,Z=options.Z,dwarf=True)]
    else:
        iso= [iso]
    df= None
    print "Pre-calculating isochrone distance prior ..."
    logpiso= numpy.zeros((len(data),_BINTEGRATENBINS))
    ds= numpy.linspace(_BINTEGRATEDMIN,_BINTEGRATEDMAX,
                       _BINTEGRATENBINS)
    dm= _dm(ds)
    for ii in range(len(data)):
        mh= data['H0MAG'][ii]-dm
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
    #test
    if options.plottype.lower() == 'pvloslos':
        print set(list(data['LOCATION']))
        indx= (data['LOCATION'] == options.location)
        data= data[indx]
        logpiso= logpiso[indx,:]
        if options.dwarf:
            logpisodwarf= logpisodwarf[indx,:]
        #Fit parameters
        if options.dwarf:
            params= [270./_REFV0,8./_REFR0,numpy.log(35./_REFV0),0.025,0.5,1.]
        else:
            params= [270./_REFV0,8./_REFR0,numpy.log(35./_REFV0),0.025,0.5]
        if not options.init is None:
            #Load initial parameters from file
            savefile= open(options.init,'rb')
            params= pickle.load(savefile)
            savefile.close()
        #Calculate vlos | los
        vlos= numpy.linspace(-200.,200.,options.nvlos)
        pvlos= numpy.zeros(options.nvlos)
        if not options.multi is None:
            pvlos= multi.parallel_map((lambda x: pvlosplate(params,vlos[x],
                                                            data,df,options,
                                                            logpiso,
                                                            logpisodwarf)),
                                      range(options.nvlos),
                                      numcores=numpy.amin([len(vlos),multiprocessing.cpu_count(),options.multi]))
        else:
            for ii in range(options.nvlos):
                print ii
                pvlos[ii]= pvlosplate(params,vlos[ii],data,df,options,
                                      logpiso,logpisodwarf)
        pvlos-= logsumexp(pvlos)
        pvlos= numpy.exp(pvlos)
        #Plot data
        bovy_plot.bovy_print()
        hist, xvec, p= bovy_plot.bovy_hist(data['VHELIO'],range=[-200.,200.],
                                           bins=31,
                                           histtype='step',color='k',
                                           xlabel=r'$\mathrm{heliocentric}\ v_{\mathrm{los}}\ [\mathrm{km\ s}^{-1}]$')
        #Normalize prediction
        data_int= numpy.sum(hist)*(xvec[1]-xvec[0])
        pvlos*= data_int/numpy.sum(pvlos)/(vlos[1]-vlos[0])
        bovy_plot.bovy_plot(vlos,pvlos,'-',color='0.65',overplot=True)
        #Add text
        bovy_plot.bovy_text(r'$\mathrm{location}\ =\ %i$' % options.location
                            +'\n'
                            +r'$l\ \approx\ %.0f^\circ$' % numpy.mean(data['GLON']),
                            top_right=True,size=14.)
        bovy_plot.bovy_end_print(options.plotfile)

