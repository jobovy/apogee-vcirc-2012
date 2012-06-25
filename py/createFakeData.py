import os, os.path
import numpy
import cPickle as pickle
from scipy.maxentropy import logsumexp
import multi
import multiprocessing
import fitsio
from readVclosData import readVclosData
import isomodel
import isodist
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
    if os.path.exists(options.plotfile):
        print "Outfile "+options.plotfile+" exists ..."
        print "Returning ..."
        return None
    #Read the data
    numpy.random.seed(options.seed)
    print "Reading the data ..."
    data= readVclosData(postshutdown=options.postshutdown,
                        fehcut=options.fehcut,
                        cohort=options.cohort,
                        lmin=options.lmin,
                        bmax=options.bmax,
                        validfeh=options.indivfeh, #if indivfeh, we need validfeh
                        ak=True,
                        cutmultiples=options.cutmultiples,
                        jkmax=options.jkmax)
    #HACK
    indx= (data['J0MAG']-data['K0MAG'] < 0.5)
    data['J0MAG'][indx]= 0.5+data['K0MAG'][indx]
    #Set up the isochrone
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
    if not options.multi is None and options.multi > 1:
        thismulti= options.multi
        options.multi= 1 #To avoid conflict
        thispvlos= multi.parallel_map((lambda x: -mloglike(params,
                                                           numpy.zeros(len(data))+vlos[x],
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
                                                           thislogpisodwarf,
                                                           True,
                                                           None,None,None)),
                                      range(options.nvlos),
                                      numcores=numpy.amin([len(vlos),multiprocessing.cpu_count(),thismulti]))
        for jj in range(options.nvlos):
            pvlos[:,jj]= thispvlos[jj]
    else:
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
                                   thislogpisodwarf,True,None,None,None)
    """
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
                               thislogpisodwarf,True,None,None,None)
    """
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
