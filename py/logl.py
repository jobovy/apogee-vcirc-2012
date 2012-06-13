import os, os.path
import numpy
import cPickle as pickle
import isomodel
from readVclosData import readVclosData
from fitvc import _BINTEGRATENBINS, _BINTEGRATEDMAX, _BINTEGRATEDMIN, \
    _BINTEGRATEDMIN_DWARF, _BINTEGRATEDMAX_DWARF, _dm, get_options, \
    _DEGTORAD, mloglike
def logl(init=None,data=None,options=None):
    if options is None:
        parser= get_options()
        options, args= parser.parse_args([])
    if data is None:
        #Read data
        data= readVclosData(lmin=25.,
                            bmax=2.,
                            ak=True,
                            jkmax=1.1)
    #HACK
    indx= (data['J0MAG']-data['K0MAG'] < 0.5)
    data['J0MAG'][indx]= 0.5+data['K0MAG'][indx]
    #Set up the isochrone
    if not options.isofile is None and os.path.exists(options.isofile):
        print "Loading the isochrone model ..."
        isofile= open(options.isofile,'rb')
        iso= pickle.load(isofile)
        locl= pickle.load(isofile)
        isofile.close()
    else:
        if options.varfeh:
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
            iso= isomodel.isomodel(Z=0.019)
        if options.dwarf:
            iso= [iso, 
                  isomodel.isomodel(Z=0.019,
                                    dwarf=True)]
        else:
            iso= [iso]
        if not options.isofile is None:
            isofile= open(options.isofile,'wb')
            pickle.dump(iso,isofile)
            pickle.dump(locl,isofile)
            isofile.close()
    df= None
    #Pre-calculate distance prior
    logpiso= numpy.zeros((len(data),_BINTEGRATENBINS))
    ds= numpy.linspace(_BINTEGRATEDMIN,_BINTEGRATEDMAX,
                       _BINTEGRATENBINS)
    dm= _dm(ds)
    for ii in range(len(data)):
        mh= data['H0MAG'][ii]-dm
        if options.varfeh:
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
    if isinstance(init,str): #FILE
        #Load initial parameters from file
        savefile= open(init,'rb')
        params= pickle.load(savefile)
        savefile.close()
    else: #Array
        params= init
    #Prep data
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
    options.multi= 1
    out= -mloglike(params,data['VHELIO'],
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
                   logpisodwarf,True,None,iso) #None iso for now
    return out
