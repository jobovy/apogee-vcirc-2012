import numpy
import cPickle as pickle
from scipy.maxentropy import logsumexp
import multi
import multiprocessing
import fitsio
from readVclosData import readVclosData
import isomodel
from fitvc import mloglike, _dm, \
    _DEGTORAD, _REFV0, _REFR0, \
    _BINTEGRATEDMAX, _BINTEGRATEDMIN, _BINTEGRATENBINS, \
    _BINTEGRATEDMAX_DWARF, _BINTEGRATEDMIN_DWARF
from compareDataModel import get_options, pvlosplate
def createFakeData(parser):
    options, args= parser.parse_args()
    if len(args) == 0:
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
    #Load initial parameters from file
    savefile= open(args[0],'rb')
    params= pickle.load(savefile)
    savefile.close()
    #Prep data
    l= data['GLON']*_DEGTORAD
    b= data['GLAT']*_DEGTORAD
    sinl= numpy.sin(l)
    cosl= numpy.cos(l)
    sinb= numpy.sin(b)
    cosb= numpy.cos(b)
    jk= data['J0MAG']-data['K0MAG']
    jk[(jk < 0.5)]= 0.5 #BOVY: FIX THIS HACK BY EMAILING GAIL
    h= data['H0MAG']
    #Re-sample
    vlos= numpy.linspace(-200.,200.,options.nvlos)
    pvlos= numpy.zeros((len(data),options.nvlos))
    if options.dwarf:
        thislogpisodwarf= logpisodwarf
    else:
        thislogpisodwarf= None
    for jj in range(options.nvlos):
        pvlos[:,jj]= -mloglike(params,numpy.zeros(len(data))+vlos[jj],
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
                               thislogpisodwarf,True,None,None)
    for ii in range(len(data)):
        pvlos[ii,:]-= logsumexp(pvlos[ii,:])
        pvlos[ii,:]= numpy.exp(pvlos[ii,:])
        pvlos[ii,:]= numpy.cumsum(pvlos[ii,:])
        pvlos[ii,:]/= pvlos[ii,-1]
        #Draw
        randindx= numpy.random.uniform()
        kk= 0
        while pvlos[ii,kk] < randindx:
            kk+= 1
        data['VHELIO'][ii]= vlos[kk]
    #Dump raw
    fitsio.write(options.plotfile,data,clobber=True)

if __name__ == '__main__':
    createFakeData(get_options())
