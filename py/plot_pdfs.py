import os, os.path
import cPickle as pickle
import numpy
from scipy import special
from fitvc import get_options, _DEGTORAD, _REFR0, _REFV0, _VRSUN, _PMSGRA
from galpy.util import bovy_plot
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
def rovc(filename=None,options=None,bins=31):
    options= set_options(options)
    params= load_samples(filename)
    vcs= numpy.array([s[0] for s in params])*_REFV0
    ros= numpy.array([s[1] for s in params])*_REFR0
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(vcs,ros,'k,',levels=levels,
                          xlabel=r'$v_0\ [\mathrm{km\ s}^{-1}]$',
                          ylabel=r'$R_0\ [\mathrm{kpc}]$',
                          bins=bins,
                          xrange=[180.,250.],
                          yrange=[6.,12.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
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
    vpects= numpy.array([s[6-options.nooutliermean+(options.rotcurve.lower() == 'linear') +(options.rotcurve.lower() == 'powerlaw') + 2*(options.rotcurve.lower() == 'quadratic')+3*(options.rotcurve.lower() == 'cubic')+options.dwarf]*_PMSGRA*s[1]*_REFR0/s[0] for s in params])
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
