import sys
import os, os.path
import cPickle as pickle
import numpy
from scipy import special
from fitvc import get_options, _DEGTORAD, _REFR0, _REFV0, _VRSUN, _PMSGRA
from galpy.util import bovy_plot
from matplotlib import pyplot
_vcrange=[180.,250.]
_vclabel= r'$v_0\ [\mathrm{km\ s}^{-1}]$'
_JACKFILES=['../fits/allwoloc4151_simpledrift_noro_dwarf_vpec_sratio_hs.sav',
            '../fits/allwoloc4154_simpledrift_noro_dwarf_vpec_sratio_hs.sav',
            '../fits/allwoloc4157_simpledrift_noro_dwarf_vpec_sratio_hs.sav',
            '../fits/allwoloc4240_simpledrift_noro_dwarf_vpec_sratio_hs.sav',
            '../fits/allwoloc4241_simpledrift_noro_dwarf_vpec_sratio_hs.sav',
            '../fits/allwoloc4242_simpledrift_noro_dwarf_vpec_sratio_hs.sav',
            '../fits/allwoloc4243_simpledrift_noro_dwarf_vpec_sratio_hs.sav',
            '../fits/allwoloc4270_simpledrift_noro_dwarf_vpec_sratio_hs.sav',
            '../fits/allwoloc4271_simpledrift_noro_dwarf_vpec_sratio_hs.sav',
            '../fits/allwoloc4272_simpledrift_noro_dwarf_vpec_sratio_hs.sav',
            '../fits/allwoloc4273_simpledrift_noro_dwarf_vpec_sratio_hs.sav',
            '../fits/allwoloc4318_simpledrift_noro_dwarf_vpec_sratio_hs.sav',
            '../fits/allwoloc4319_simpledrift_noro_dwarf_vpec_sratio_hs.sav',
            '../fits/allwoloc4321_simpledrift_noro_dwarf_vpec_sratio_hs.sav']
_JACKFILES_PL=['../fits/allwoloc4151_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav',
            '../fits/allwoloc4154_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav',
            '../fits/allwoloc4157_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav',
            '../fits/allwoloc4240_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav',
            '../fits/allwoloc4241_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav',
            '../fits/allwoloc4242_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav',
            '../fits/allwoloc4243_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav',
            '../fits/allwoloc4270_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav',
            '../fits/allwoloc4271_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav',
            '../fits/allwoloc4272_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav',
            '../fits/allwoloc4273_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav',
            '../fits/allwoloc4318_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav',
            '../fits/allwoloc4319_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav',
            '../fits/allwoloc4321_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav']

_JACKFILES_LINEAR=['../fits/allwoloc4151_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav',
            '../fits/allwoloc4154_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav',
            '../fits/allwoloc4157_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav',
            '../fits/allwoloc4240_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav',
            '../fits/allwoloc4241_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav',
            '../fits/allwoloc4242_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav',
            '../fits/allwoloc4243_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav',
            '../fits/allwoloc4270_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav',
            '../fits/allwoloc4271_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav',
            '../fits/allwoloc4272_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav',
            '../fits/allwoloc4273_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav',
            '../fits/allwoloc4318_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav',
            '../fits/allwoloc4319_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav',
            '../fits/allwoloc4321_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav']
_JACKCOLOR='w'
_JACKSYMBOL='x'
def load_jack_params(linear=False,powerlaw=False):
    out= []
    njacks= len(_JACKFILES)
    for ii in range(njacks):
        if linear:
            savefile= open(_JACKFILES_LINEAR[ii],'rb')
        elif powerlaw:
            savefile= open(_JACKFILES_PL[ii],'rb')
        else:
            savefile= open(_JACKFILES[ii],'rb')
        out.append(pickle.load(savefile))
        savefile.close()
    return out
_PRELOADJACK= True
if _PRELOADJACK:
    _JACKPARAMS= load_jack_params()
_PLOTJACK= True
def set_options(options):
    if options is None:
        parser= get_options()
        options, args= parser.parse_args([])
    #Set up to default fit
    options.dwarf= True
    options.noroprior= True
    options.fitvpec= True
    options.dfmodel= 'simplegaussiandrift'
    options.fitsratio= True
    options.fiths= True
    return options
def load_samples(filename):
    if not os.path.exists(filename):
        raise IOError("given filename does not exist")
    savefile= open(filename,'rb')
    params= pickle.load(savefile)
    savefile.close()
    return params
def rovo(filename=None,options=None,bins=31):
    """Confusingly, this plots ro Omegao"""
    options= set_options(options)
    params= load_samples(filename)
    vos= numpy.array([s[0] for s in params])*_REFV0
    ros= numpy.array([s[1] for s in params])*_REFR0
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(vos/ros,ros,'k,',levels=levels,
                          xlabel=r'$\Omega_0\ [\mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}]$',
                          ylabel=r'$R_0\ [\mathrm{kpc}]$',
                          bins=bins,
                          xrange=[150./8.,250./8.],
                          yrange=[6.,12.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def rovc(filename=None,options=None,bins=31,multipops=False):
    options= set_options(options)
    if isinstance(filename,str):
        params= load_samples(filename)
        vcs= numpy.array([s[0] for s in params])*_REFV0
        ros= numpy.array([s[1] for s in params])*_REFR0
        bovy_plot.bovy_print()
        levels= list(special.erf(0.5*numpy.arange(1,4)))
        levels.append(1.01) #HACK to not plot outliers
        bovy_plot.scatterplot(ros,vcs,'k,',levels=levels,
                              ylabel=_vclabel,
                              xlabel=r'$R_0\ [\mathrm{kpc}]$',
                              bins=bins,
                              yrange=_vcrange,
                              xrange=[6.,12.],
                              contours=True,
                              cntrcolors='k',
                              cmap='gist_yarg',onedhistx=True,
                              onedhists=multipops)
        if multipops:
            bovy_plot.bovy_text(r'$\mathrm{Multiple\ populations,\ SFR} = \exp\left( -t/8\ \mathrm{Gyr}\right)$',top_right=True,size=14.)
        if _PLOTJACK:
            if _PRELOADJACK:
                for ii in range(len(_JACKPARAMS)):
                    rovc(filename=_JACKPARAMS[ii])
    else:
        params= filename
        bovy_plot.bovy_plot(params[1]*_REFR0,params[0]*_REFV0,
                            color=_JACKCOLOR,marker=_JACKSYMBOL,
                            overplot=True)        
    return None
def rosr(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    ros= numpy.array([s[1] for s in params])*_REFR0
    srs= numpy.exp(numpy.array([s[2] for s in params]))*_REFV0
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(srs,ros,'k,',levels=levels,
                          xlabel=r'$\sigma_R\ [\mathrm{km\ s}^{-1}]$',
                          ylabel=r'$R_0\ [\mathrm{kpc}]$',
                          bins=bins,
                          xrange=[20.,50.],
                          yrange=[6.,12.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def vcsr(filename=None,options=None,bins=31):
    options= set_options(options)
    if isinstance(filename,str):
        params= load_samples(filename)
        vcs= numpy.array([s[0] for s in params])*_REFV0
        srs= numpy.exp(numpy.array([s[2] for s in params]))*_REFV0
        bovy_plot.bovy_print()
        levels= list(special.erf(0.5*numpy.arange(1,4)))
        levels.append(1.01) #HACK to not plot outliers
        bovy_plot.scatterplot(srs,vcs,'k,',levels=levels,
                              xlabel=r'$\sigma_R\ [\mathrm{km\ s}^{-1}]$',
                              ylabel=r'$v_0\ [\mathrm{km\ s}^{-1}]$',
                              bins=bins,
                              xrange=[20.,50.],
                              yrange=_vcrange,
                              contours=True,
                              cntrcolors='k',
                              onedhistx=True,
                              cmap='gist_yarg')
        if _PLOTJACK:
            if _PRELOADJACK:
                for ii in range(len(_JACKPARAMS)):
                    vcsr(filename=_JACKPARAMS[ii])
    else:
        params= filename
        bovy_plot.bovy_plot(numpy.exp(params[2])*_REFV0,
                            params[0]*_REFV0,
                            color=_JACKCOLOR,marker=_JACKSYMBOL,
                            overplot=True)        
    return None
def rohs(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    ros= numpy.array([s[1] for s in params])*_REFR0
    hss= options.hs/numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] for s in params])
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(hss,ros,'k,',levels=levels,
                          xlabel=r'$h_\sigma [\mathrm{kpc}]$',
                          ylabel=r'$R_0\ [\mathrm{kpc}]$',
                          bins=bins,
                          xrange=[0.,50.],
                          yrange=[6.,12.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def vchs(filename=None,options=None,bins=31):
    options= set_options(options)
    if isinstance(filename,str):
        params= load_samples(filename)
        vos= numpy.array([s[0] for s in params])*_REFV0
        ros= numpy.array([s[1] for s in params])*_REFR0
        hss= options.hs/numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] for s in params])
        bovy_plot.bovy_print()
        levels= list(special.erf(0.5*numpy.arange(1,4)))
        levels.append(1.01) #HACK to not plot outliers
        bovy_plot.scatterplot(ros/hss,vos,'k,',levels=levels,
                              xlabel=r'$R_0/h_\sigma$',
                              ylabel=_vclabel,
                              bins=bins,
                              xrange=[-1.,1.],
                              yrange=_vcrange,
                              contours=True,
                              cntrcolors='k',
                              onedhistx=True,
                              cmap='gist_yarg')
        if _PLOTJACK:
            if _PRELOADJACK:
                for ii in range(len(_JACKPARAMS)):
                    vchs(filename=_JACKPARAMS[ii])
    else:
        params= filename
        bovy_plot.bovy_plot(params[1]*_REFR0/(options.hs/params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter]),
                            params[0]*_REFV0,
                            color=_JACKCOLOR,marker=_JACKSYMBOL,
                            overplot=True)        
    return None
def vcx2(filename=None,options=None,bins=31):
    options= set_options(options)
    if isinstance(filename,str):
        params= load_samples(filename)
        vos= numpy.array([s[0] for s in params])*_REFV0
        x2s= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf] for s in params])
        bovy_plot.bovy_print()
        levels= list(special.erf(0.5*numpy.arange(1,4)))
        levels.append(1.01) #HACK to not plot outliers
        bovy_plot.scatterplot(x2s,vos,'k,',levels=levels,
                              xlabel=r'$\sigma_\phi^2/\sigma_R^2$',
                              ylabel=_vclabel,
                              bins=bins,
                              xrange=[0.,1.5],
                              yrange=_vcrange,
                              contours=True,
                              cntrcolors='k',
                              onedhists=True,
                              cmap='gist_yarg')
        if _PLOTJACK:
            if _PRELOADJACK:
                for ii in range(len(_JACKPARAMS)):
                    vcx2(filename=_JACKPARAMS[ii])
    else:
        params= filename
        bovy_plot.bovy_plot(params[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf],
                            params[0]*_REFV0,
                            color=_JACKCOLOR,marker=_JACKSYMBOL,
                            overplot=True)        
    return None
def roah(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    ros= numpy.array([s[1] for s in params])*_REFR0
    #ahs= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitah+options.fitdm] for s in params])
    ahs= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] for s in params])
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(ahs,ros,'k,',levels=levels,
                          xlabel=r'$\Delta A_H [\mathrm{mag}]$',
                          ylabel=r'$R_0\ [\mathrm{kpc}]$',
                          bins=bins,
                          xrange=[-.1,.1],
                          yrange=[6.,12.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def ahah(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    ahinners= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitah+options.fitdm] for s in params])
    ahs= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] for s in params])
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(ahs,ahinners,'k,',levels=levels,
                          xlabel=r'$\Delta A_H\ (\mathrm{outer}) [\mathrm{mag}]$',
                          ylabel=r'$\Delta A_H\ (\mathrm{inner}) [\mathrm{mag}]$',
                          bins=bins,
                          xrange=[-0.1,.1],
                          yrange=[-0.1,.1],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def vcfeh(filename=None,filename2=None,options=None,bins=31):
    options= set_options(options)
    options.fitfeh= True
    params= load_samples(filename)
    vcs= numpy.array([s[0] for s in params])*_REFV0
    dms= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] for s in params])
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    dmrange= [-1.,1.]
    axScatter, axHistx, axHisty= bovy_plot.scatterplot(dms,vcs,'k,',levels=levels,
                                                       xlabel=r'$\Delta [\mathrm{Fe/H}]\ [\mathrm{dex}]$',
                          ylabel=_vclabel,
                          bins=bins,
                          xrange=dmrange,
                          yrange=_vcrange,
                          contours=True,
                          cntrcolors='k',
                          onedhistx=True,
                                                       retAxes=True,
                          cmap='gist_yarg')
    #Now plot inner-outer
    params= load_samples(filename2)
    vcs= numpy.array([s[0] for s in params])*_REFV0
    dms= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] for s in params])
    dmsinner= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitdm+options.fitah+options.fitfeh] for s in params])
    #Just plot contours
    data= numpy.array([dms,vcs]).T
    hist, edges= numpy.histogramdd(data,bins=bins,range=[dmrange,_vcrange])
    X= hist.T
    X[numpy.isnan(X)]= 0.
    sortindx= numpy.argsort(X.flatten())[::-1]
    cumul= numpy.cumsum(numpy.sort(X.flatten())[::-1])/numpy.sum(X.flatten())
    cntrThis= numpy.zeros(numpy.prod(X.shape))
    cntrThis[sortindx]= cumul
    cntrThis= numpy.reshape(cntrThis,X.shape)
    extent=[dmrange[0],dmrange[1],_vcrange[0],_vcrange[1]]
    cont= pyplot.contour(cntrThis,levels,colors='k',
                         extent=extent,
                         linestyles='dashed',
                         origin='lower')
    histx= numpy.nansum(X.T,axis=1)*numpy.fabs(_vcrange[1]-_vcrange[0])/X.shape[1]
    histx[numpy.isnan(histx)]= 0.
    dx= (extent[1]-extent[0])/float(len(histx))
    axHistx.plot(numpy.linspace(extent[0]+dx,extent[1]-dx,len(histx)),histx/3800.,
                 drawstyle='steps-mid',color='k',ls='--')
    #Just plot contours
    data= numpy.array([dmsinner,vcs]).T
    hist, edges= numpy.histogramdd(data,bins=bins,range=[dmrange,_vcrange])
    X= hist.T
    X[numpy.isnan(X)]= 0.
    sortindx= numpy.argsort(X.flatten())[::-1]
    cumul= numpy.cumsum(numpy.sort(X.flatten())[::-1])/numpy.sum(X.flatten())
    cntrThis= numpy.zeros(numpy.prod(X.shape))
    cntrThis[sortindx]= cumul
    cntrThis= numpy.reshape(cntrThis,X.shape)
    cont= pyplot.contour(cntrThis,levels,colors='k',
                         extent=extent,
                         linestyles='dashed',
                         origin='lower')
    histx= numpy.nansum(X.T,axis=1)*numpy.fabs(_vcrange[1]-_vcrange[0])/X.shape[1]
    histx[numpy.isnan(histx)]= 0.
    dx= (extent[1]-extent[0])/float(len(histx))
    axHistx.plot(numpy.linspace(extent[0]+dx,extent[1]-dx,len(histx)),histx/3800.,
                 drawstyle='steps-mid',color='k',ls='--')
    return None
def vcah(filename=None,filename2=None,options=None,bins=31):
    options= set_options(options)
    options.fitah= True
    params= load_samples(filename)
    vcs= numpy.array([s[0] for s in params])*_REFV0
    dms= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] for s in params])
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    dmrange= [-.25,.25]
    axScatter, axHistx, axHisty= bovy_plot.scatterplot(dms,vcs,'k,',levels=levels,
                                                       xlabel=r'$\Delta A_H\ [\mathrm{mag}]$',
                          ylabel=_vclabel,
                          bins=bins,
                          xrange=dmrange,
                          yrange=_vcrange,
                          contours=True,
                          cntrcolors='k',
                          onedhistx=True,
                                                       retAxes=True,
                          cmap='gist_yarg')
    #Now plot inner-outer
    params= load_samples(filename2)
    vcs= numpy.array([s[0] for s in params])*_REFV0
    dms= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] for s in params])
    dmsinner= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitdm+options.fitah] for s in params])
    #Just plot contours
    data= numpy.array([dms,vcs]).T
    hist, edges= numpy.histogramdd(data,bins=bins,range=[dmrange,_vcrange])
    X= hist.T
    X[numpy.isnan(X)]= 0.
    sortindx= numpy.argsort(X.flatten())[::-1]
    cumul= numpy.cumsum(numpy.sort(X.flatten())[::-1])/numpy.sum(X.flatten())
    cntrThis= numpy.zeros(numpy.prod(X.shape))
    cntrThis[sortindx]= cumul
    cntrThis= numpy.reshape(cntrThis,X.shape)
    extent=[dmrange[0],dmrange[1],_vcrange[0],_vcrange[1]]
    cont= pyplot.contour(cntrThis,levels,colors='k',
                         extent=extent,
                         linestyles='dashed',
                         origin='lower')
    histx= numpy.nansum(X.T,axis=1)*numpy.fabs(_vcrange[1]-_vcrange[0])/X.shape[1]
    histx[numpy.isnan(histx)]= 0.
    dx= (extent[1]-extent[0])/float(len(histx))
    axHistx.plot(numpy.linspace(extent[0]+dx,extent[1]-dx,len(histx)),histx/800.,
                 drawstyle='steps-mid',color='k',ls='--')
    #Just plot contours
    data= numpy.array([dmsinner,vcs]).T
    hist, edges= numpy.histogramdd(data,bins=bins,range=[dmrange,_vcrange])
    X= hist.T
    X[numpy.isnan(X)]= 0.
    sortindx= numpy.argsort(X.flatten())[::-1]
    cumul= numpy.cumsum(numpy.sort(X.flatten())[::-1])/numpy.sum(X.flatten())
    cntrThis= numpy.zeros(numpy.prod(X.shape))
    cntrThis[sortindx]= cumul
    cntrThis= numpy.reshape(cntrThis,X.shape)
    cont= pyplot.contour(cntrThis,levels,colors='k',
                         extent=extent,
                         linestyles='dashed',
                         origin='lower')
    histx= numpy.nansum(X.T,axis=1)*numpy.fabs(_vcrange[1]-_vcrange[0])/X.shape[1]
    histx[numpy.isnan(histx)]= 0.
    dx= (extent[1]-extent[0])/float(len(histx))
    axHistx.plot(numpy.linspace(extent[0]+dx,extent[1]-dx,len(histx)),histx/800.,
                 drawstyle='steps-mid',color='k',ls='--')
    return None
def vcdm(filename=None,filename2=None,options=None,bins=31):
    options= set_options(options)
    options.fitdm= True
    params= load_samples(filename)
    vcs= numpy.array([s[0] for s in params])*_REFV0
    dms= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] for s in params])
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    dmrange= [-1.25,1.25]
    axScatter, axHistx, axHisty= bovy_plot.scatterplot(dms,vcs,'k,',levels=levels,
                          xlabel=r'$\Delta \mu\ [\mathrm{mag}]$',
                          ylabel=_vclabel,
                          bins=bins,
                          xrange=dmrange,
                          yrange=_vcrange,
                          contours=True,
                          cntrcolors='k',
                          onedhistx=True,
                                                       retAxes=True,
                          cmap='gist_yarg')
    #Now plot inner-outer
    params= load_samples(filename2)
    vcs= numpy.array([s[0] for s in params])*_REFV0
    dms= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] for s in params])
    dmsinner= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitdm+options.fitah] for s in params])
    #Just plot contours
    data= numpy.array([dms,vcs]).T
    hist, edges= numpy.histogramdd(data,bins=bins,range=[dmrange,_vcrange])
    X= hist.T
    X[numpy.isnan(X)]= 0.
    sortindx= numpy.argsort(X.flatten())[::-1]
    cumul= numpy.cumsum(numpy.sort(X.flatten())[::-1])/numpy.sum(X.flatten())
    cntrThis= numpy.zeros(numpy.prod(X.shape))
    cntrThis[sortindx]= cumul
    cntrThis= numpy.reshape(cntrThis,X.shape)
    extent=[dmrange[0],dmrange[1],_vcrange[0],_vcrange[1]]
    cont= pyplot.contour(cntrThis,levels,colors='k',
                         extent=extent,
                         linestyles='dashed',
                         origin='lower')
    histx= numpy.nansum(X.T,axis=1)*numpy.fabs(_vcrange[1]-_vcrange[0])/X.shape[1]
    histx[numpy.isnan(histx)]= 0.
    dx= (extent[1]-extent[0])/float(len(histx))
    axHistx.plot(numpy.linspace(extent[0]+dx,extent[1]-dx,len(histx)),histx/3000.,
                 drawstyle='steps-mid',color='k',ls='--')
    #Just plot contours
    data= numpy.array([dmsinner,vcs]).T
    hist, edges= numpy.histogramdd(data,bins=bins,range=[dmrange,_vcrange])
    X= hist.T
    X[numpy.isnan(X)]= 0.
    sortindx= numpy.argsort(X.flatten())[::-1]
    cumul= numpy.cumsum(numpy.sort(X.flatten())[::-1])/numpy.sum(X.flatten())
    cntrThis= numpy.zeros(numpy.prod(X.shape))
    cntrThis[sortindx]= cumul
    cntrThis= numpy.reshape(cntrThis,X.shape)
    cont= pyplot.contour(cntrThis,levels,colors='k',
                         extent=extent,
                         linestyles='dashed',
                         origin='lower')
    histx= numpy.nansum(X.T,axis=1)*numpy.fabs(_vcrange[1]-_vcrange[0])/X.shape[1]
    histx[numpy.isnan(histx)]= 0.
    dx= (extent[1]-extent[0])/float(len(histx))
    axHistx.plot(numpy.linspace(extent[0]+dx,extent[1]-dx,len(histx)),histx/5000.,
                 drawstyle='steps-mid',color='k',ls='--')
    return None
def rodm(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    ros= numpy.array([s[1] for s in params])*_REFR0
    dms= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter] for s in params])
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(dms,ros,'k,',levels=levels,
                          xlabel=r'$\Delta m [\mathrm{mag}]$',
                          ylabel=r'$R_0\ [\mathrm{kpc}]$',
                          bins=bins,
                          xrange=[-0.5,.5],
                          yrange=[6.,12.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def rodminner(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    ros= numpy.array([s[1] for s in params])*_REFR0
    dms= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf+options.fitsratio+2*options.fitsratioinnerouter+options.fiths+options.fitsrinnerouter+options.dwarfinnerouter+options.fitdm+options.fitah] for s in params])
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(dms,ros,'k,',levels=levels,
                          xlabel=r'$\Delta m\ (\mathrm{inner})\ [\mathrm{mag}]$',
                          ylabel=r'$R_0\ [\mathrm{kpc}]$',
                          bins=bins,
                          xrange=[-1.,1.],
                          yrange=[6.,12.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def rovpecr(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    vpecrs= numpy.array([s[5-options.nooutliermean+options.dwarf+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')]*_VRSUN/s[0] for s in params])
    ros= numpy.array([s[1] for s in params])*_REFR0
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(vpecrs,ros,'k,',levels=levels,
                          xlabel=r'$v_{R,\odot}\ [\mathrm{km\ s}^{-1}]$',
                          ylabel=r'$R_0\ [\mathrm{kpc}]$',
                          bins=bins,
                          xrange=[-20.,10.],
                          yrange=[6.,12.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def rovpect(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    vpects= numpy.array([s[6-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+options.dwarf]*_PMSGRA*s[1]*_REFR0 for s in params])
    ros= numpy.array([s[1] for s in params])*_REFR0
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(vpects,ros,'k,',levels=levels,
                          xlabel=r'$v_{\phi,\odot}\ [\mathrm{km\ s}^{-1}]$',
                          ylabel=r'$R_0\ [\mathrm{kpc}]$',
                          bins=bins,
                          xrange=[220.,290.],
                          yrange=[6.,12.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def vovpect(filename=None,options=None,bins=31):
    options= set_options(options)
    if isinstance(filename,str):
        params= load_samples(filename)
        vpects= numpy.array([s[6-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+options.dwarf]*_PMSGRA*s[1]*_REFR0 for s in params])
        vos= numpy.array([s[0] for s in params])*_REFV0
        bovy_plot.bovy_print()
        levels= list(special.erf(0.5*numpy.arange(1,4)))
        levels.append(1.01) #HACK to not plot outliers
        bovy_plot.scatterplot(vpects-vos,vos,'k,',levels=levels,
                              xlabel=r'$v_{\phi,\odot}-v_0\ [\mathrm{km\ s}^{-1}]$',
                              ylabel=_vclabel,
                              bins=bins,
                              xrange=[0.,40.],
                              yrange=_vcrange,
                              contours=True,
                              cntrcolors='k',
                              onedhistx=True,
                              cmap='gist_yarg')
        if _PLOTJACK:
            if _PRELOADJACK:
                for ii in range(len(_JACKPARAMS)):
                    vovpect(filename=_JACKPARAMS[ii])
    else:
        params= filename
        bovy_plot.bovy_plot(params[6-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+options.dwarf]*_PMSGRA*params[1]*_REFR0-params[0]*_REFV0,
                            params[0]*_REFV0,
                            color=_JACKCOLOR,marker=_JACKSYMBOL,
                            overplot=True)        
    return None
def vopdwarf(filename=None,options=None,bins=31):
    options= set_options(options)
    if isinstance(filename,str):
        params= load_samples(filename)
        pdwarfs= numpy.array([s[5-options.nooutliermean] for s in params])
        vos= numpy.array([s[0] for s in params])*_REFV0
        bovy_plot.bovy_print()
        levels= list(special.erf(0.5*numpy.arange(1,4)))
        levels.append(1.01) #HACK to not plot outliers
        bovy_plot.scatterplot(pdwarfs,vos,'k,',levels=levels,
                              xlabel=r'$P(\mathrm{dwarf})$',
                              ylabel=_vclabel,
                              bins=bins,
                              xrange=[0.,1.],
                              yrange=_vcrange,
                              contours=True,
                              cntrcolors='k',
                              onedhists=True,
                              cmap='gist_yarg')
        if _PLOTJACK:
            if _PRELOADJACK:
                for ii in range(len(_JACKPARAMS)):
                    vopdwarf(filename=_JACKPARAMS[ii])
    else:
        params= filename
        bovy_plot.bovy_plot(params[5-options.nooutliermean],
                            params[0]*_REFV0,
                            color=_JACKCOLOR,marker=_JACKSYMBOL,
                            overplot=True)        
    return None
def roosun(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    vpects= numpy.array([s[6-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+options.dwarf] for s in params])
    ros= numpy.array([s[1] for s in params])*_REFR0
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(vpects,ros,'k,',levels=levels,
                          xlabel=r'$\Omega_{\odot}/\mu_{\mathrm{Sgr\ A}^*}\ [\mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}]$',
                          ylabel=r'$R_0\ [\mathrm{kpc}]$',
                          bins=bins,
                          xrange=[0.75,1.25],
                          yrange=[6.,12.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def voosun(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    vpects= numpy.array([s[6-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+options.dwarf] for s in params])
    vos= numpy.array([s[0] for s in params])*_REFV0
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(vpects,vos,'k,',levels=levels,
                          xlabel=r'$\Omega_{\odot}/\mu_{\mathrm{Sgr\ A}^*}\ [\mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}]$',
                          ylabel=r'$v_0\ [\mathrm{km\ s}^{-1}]$',
                          bins=bins,
                          xrange=[0.75,1.25],
                          yrange=[180.,240.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def rooutmean(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    outmeans= numpy.array([s[4]*_REFV0 for s in params])
    ros= numpy.array([s[1] for s in params])*_REFR0
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(outmeans,ros,'k,',levels=levels,
                          xlabel=r'$\bar{v}_{\mathrm{outlier}\ [\mathrm{km\ s}^{-1}]$',
                          ylabel=r'$R_0\ [\mathrm{kpc}]$',
                          bins=bins,
                          xrange=[0.,250.],
                          yrange=[6.,12.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def rooutfrac(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    outfracs= numpy.array([s[3] for s in params])
    ros= numpy.array([s[1] for s in params])*_REFR0
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(outfracs,ros,'k,',levels=levels,
                          xlabel=r'$P(\mathrm{bad})$',
                          ylabel=r'$R_0\ [\mathrm{kpc}]$',
                          bins=bins,
                          xrange=[0.,1.],
                          yrange=[6.,12.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def rodwarffrac(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    dwarffracs= numpy.array([s[5-options.nooutliermean] for s in params])
    ros= numpy.array([s[1] for s in params])*_REFR0
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(dwarffracs,ros,'k,',levels=levels,
                          xlabel=r'$P(\mathrm{dwarf})$',
                          ylabel=r'$R_0\ [\mathrm{kpc}]$',
                          bins=bins,
                          xrange=[0.,.2],
                          yrange=[6.,12.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def rosratio(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    sratios= numpy.array([s[5-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+2*options.fitvpec+options.dwarf] for s in params])
    ros= numpy.array([s[1] for s in params])*_REFR0
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(sratios,ros,'k,',levels=levels,
                          xlabel=r'$\sigma_T^2/\sigma_R^2$',
                          ylabel=r'$R_0\ [\mathrm{kpc}]$',
                          bins=bins,
                          xrange=[0.,1.],
                          yrange=[6.,12.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def rodvcdr(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    dvcdrs= numpy.array([s[5-options.nooutliermean+options.dwarf] for s in params])
    ros= numpy.array([s[1] for s in params])*_REFR0
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(dvcdrs,ros,'k,',levels=levels,
                          xlabel=r'$\beta$',
                          ylabel=r'$R_0\ [\mathrm{kpc}]$',
                          bins=bins,
                          xrange=[-0.5,.5],
                          yrange=[6.,12.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def vcdvcdr(filename=None,options=None,bins=31):
    options= set_options(options)
    if isinstance(filename,str):
        params= load_samples(filename)
        dvcdrs= numpy.array([s[5-options.nooutliermean+options.dwarf] for s in params])*_REFV0/_REFR0
        vcs= numpy.array([s[0] for s in params])*_REFV0
        bovy_plot.bovy_print()
        levels= list(special.erf(0.5*numpy.arange(1,4)))
        levels.append(1.01) #HACK to not plot outliers
        bovy_plot.scatterplot(dvcdrs,vcs,'k,',levels=levels,
                              ylabel=_vclabel,
                              xlabel=r"$\mathrm{d} v_0 / \mathrm{d} R\ [\mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}]$",
                              bins=bins,
                              xrange=[-0.5*_REFV0/_REFR0,0.5*_REFV0/_REFR0],
                              yrange=_vcrange,
                              contours=True,
                              cntrcolors='k',
                              onedhistx=True,
                              cmap='gist_yarg')
        if _PLOTJACK:
            #Load powerlaw jacks
            jackparams= load_jack_params(linear=True)
            for ii in range(len(jackparams)):
                vcdvcdr(filename=jackparams[ii])
    else:
        params= filename
        bovy_plot.bovy_plot(params[5-options.nooutliermean+options.dwarf]*_REFV0/_REFR0,
                            params[0]*_REFV0,                            
                            color=_JACKCOLOR,marker=_JACKSYMBOL,
                            overplot=True)        
    return None
def vcbeta(filename=None,options=None,bins=31):
    options= set_options(options)
    if isinstance(filename,str):
        params= load_samples(filename)
        dvcdrs= numpy.array([s[5-options.nooutliermean+options.dwarf] for s in params])
        vcs= numpy.array([s[0] for s in params])*_REFV0
        bovy_plot.bovy_print()
        levels= list(special.erf(0.5*numpy.arange(1,4)))
        levels.append(1.01) #HACK to not plot outliers
        bovy_plot.scatterplot(dvcdrs,vcs,'k,',levels=levels,
                              ylabel=_vclabel,
                              xlabel=r'$\beta$',
                              bins=bins,
                              xrange=[-0.5,0.5],
                              yrange=_vcrange,
                              contours=True,
                              cntrcolors='k',
                              onedhistx=True,
                              cmap='gist_yarg')
        if _PLOTJACK:
            #Load powerlaw jacks
            jackparams= load_jack_params(powerlaw=True)            
            for ii in range(len(jackparams)):
                vcbeta(filename=jackparams[ii])
    else:
        params= filename
        bovy_plot.bovy_plot(params[5-options.nooutliermean+options.dwarf],
                            params[0]*_REFV0,                            
                            color=_JACKCOLOR,marker=_JACKSYMBOL,
                            overplot=True)        
    return None
def vcd2vcdr2(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    dvcdrs= numpy.array([s[6-options.nooutliermean+options.dwarf] for s in params])*_REFV0/_REFR0**2.
    vcs= numpy.array([s[0] for s in params])*_REFV0
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(dvcdrs,vcs,'k,',levels=levels,
                          ylabel=_vclabel,
                          xlabel=r'$\mathrm{second\ derivative}$',
                          bins=bins,
                          xrange=[-0.5*_REFV0/_REFR0**2.,0.5*_REFV0/_REFR0**2.],
                          yrange=_vcrange,
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def vcd3vcdr3(filename=None,options=None,bins=31):
    options= set_options(options)
    params= numpy.array(load_samples(filename))
    indx= numpy.zeros(len(params),dtype='bool')
    indx[:]= True
    for ii in range(len(params)):
        if params[ii][1] > (9./_REFR0): indx[ii]= False
    params= params[indx]
    dvcdrs= numpy.array([s[7-options.nooutliermean+options.dwarf] for s in params])*_REFV0/_REFR0**3.
    vcs= numpy.array([s[0] for s in params])*_REFV0
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(dvcdrs,vcs,'k,',levels=levels,
                          ylabel=_vclabel,
                          xlabel=r'$\mathrm{third\ derivative}$',
                          bins=bins,
                          xrange=[-1.*_REFV0/_REFR0**3.,1.*_REFV0/_REFR0**3.],
                          yrange=_vcrange,
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    bovy_plot.bovy_text(r'$R_0 < 9\,\mathrm{kpc}$',
                        top_right=True,size=16.)      
    return None
def dvcdrd2vcdr2(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    dvcdrs= numpy.array([s[6-options.nooutliermean+options.dwarf] for s in params])
    d2vcdr2s= numpy.array([s[5-options.nooutliermean+options.dwarf] for s in params])
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(dvcdrs,d2vcdr2s,'k,',levels=levels,
                          xlabel=r'$\mathrm{first\ derivative}$',
                          ylabel=r'$\mathrm{second\ derivative}$',
                          bins=bins,
                          xrange=[-0.5,0.5],
                          yrange=[-0.5,0.5],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None
def dvcdrd3vcdr3(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    dvcdrs= numpy.array([s[5-options.nooutliermean+options.dwarf] for s in params])
    d2vcdr2s= numpy.array([s[7-options.nooutliermean+options.dwarf] for s in params])
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(dvcdrs,d2vcdr2s,'k,',levels=levels,
                          xlabel=r'$\mathrm{first\ derivative}$',
                          ylabel=r'$\mathrm{third\ derivative}$',
                          bins=bins,
                          xrange=[-0.5,0.5],
                          yrange=[-1.,1.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    return None

if __name__ == '__main__':
    if len(sys.argv) < 2:
        outdir= '../tex/'
    else:
        outdir= sys.argv[1]
    filename= '../fits/all_simpledrift_noro_dwarf_vpec_sratio_hs_10000samples.sav'
    bins= 16
    options= set_options(None)
    ext= 'ps'
    #First plot vc,ro
    rovc(filename=filename,bins=bins,options=options)
    bovy_plot.bovy_end_print(os.path.join(outdir,'pdf_vc_ro.'+ext))
    #vc,pdwarf
    vopdwarf(filename=filename,bins=bins,options=options)
    bovy_plot.bovy_end_print(os.path.join(outdir,'pdf_vc_pdwarf.'+ext))
    #vc,vpect
    vovpect(filename=filename,bins=bins,options=options)
    bovy_plot.bovy_end_print(os.path.join(outdir,'pdf_vc_vphisun.'+ext))
    #vc,hs
    vchs(filename=filename,bins=bins,options=options)
    bovy_plot.bovy_end_print(os.path.join(outdir,'pdf_vc_hs.'+ext))
    #vc,X2
    vcx2(filename=filename,bins=bins,options=options)
    bovy_plot.bovy_end_print(os.path.join(outdir,'pdf_vc_x2.'+ext))
    #vc,sr
    vcsr(filename=filename,bins=bins,options=options)
    bovy_plot.bovy_end_print(os.path.join(outdir,'pdf_vc_sr.'+ext))
    ###OTHER FITS
    #vc,beta
    filename= '../fits/all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_10000samples.sav'
    vcbeta(filename=filename,bins=bins,options=options)
    bovy_plot.bovy_end_print(os.path.join(outdir,'pdf_vc_beta.'+ext))
    #vc,linear dvcdr
    filename= '../fits/all_simpledrift_noro_dwarf_linear_vpec_sratio_hs_10000samples.sav'
    vcdvcdr(filename=filename,bins=bins,options=options)
    bovy_plot.bovy_end_print(os.path.join(outdir,'pdf_vc_linear.'+ext))
    #vc, dm
    filename= '../fits/all_simpledrift_noro_dwarf_vpec_sratio_dm_hs_10000samples.sav'
    filename2= '../fits/all_simpledrift_noro_dwarf_vpec_sratio_dm_hs_innerouterdm_10000samples.sav'
    vcdm(filename=filename,filename2=filename2,bins=bins,options=options)
    bovy_plot.bovy_end_print(os.path.join(outdir,'pdf_vc_dm.'+ext))
    options.fitdm= False
    #vc, ah
    filename= '../fits/all_simpledrift_noro_dwarf_vpec_sratio_ah_hs_10000samples.sav'
    filename2= '../fits/all_simpledrift_noro_dwarf_vpec_sratio_ah_hs_innerouterah_10000samples.sav'
    vcah(filename=filename,filename2=filename2,bins=bins,options=options)
    bovy_plot.bovy_end_print(os.path.join(outdir,'pdf_vc_ah.'+ext))
    options.fitah= False
    #vc, feh
    filename= '../fits/all_simpledrift_noro_dwarf_vpec_sratio_feh_hs_10000samples.sav'
    filename2= '../fits/all_simpledrift_noro_dwarf_vpec_sratio_feh_hs_innerouterfeh_10000samples.sav'
    vcfeh(filename=filename,filename2=filename2,bins=bins,options=options)
    bovy_plot.bovy_end_print(os.path.join(outdir,'pdf_vc_feh.'+ext))
    options.fitfeh= False
    #vc, ro, multiple pops
    filename= '../fits/all_multiplepops_noro_dwarf_vpec_sratio_hs_10000samples.sav'
    rovc(filename=filename,bins=bins,options=options,multipops=True)
    bovy_plot.bovy_end_print(os.path.join(outdir,'pdf_vc_ro_multipops.'+ext))
    #vc,quadratic dvcdr
    #filename= '../fits/all_simpledrift_noro_dwarf_quadratic_vpec_sratio_hs_10000samples.sav'
    #vcd2vcdr2(filename=filename,bins=bins,options=options)
    #bovy_plot.bovy_end_print(os.path.join(outdir,'pdf_vc_quadratic.'+ext))
    #vc,cubic d3vcdr3
    #filename= '../fits/all_simpledrift_noro_dwarf_cubic_vpec_sratio_hs_10000samples.sav'
    #vcd3vcdr3(filename=filename,bins=bins,options=options)
    #bovy_plot.bovy_end_print(os.path.join(outdir,'pdf_vc_cubic.'+ext))
