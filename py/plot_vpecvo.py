import sys
import os, os.path
import cPickle as pickle
import math
import numpy
from scipy import special
from galpy.util import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
from galpy.orbit import Orbit
from galpy.potential import EllipticalDiskPotential, LogarithmicHaloPotential
from fitvc import _REFV0, _REFR0, _PMSGRA, _VRSUN
def plot_vpecvo(filename,plotfilename):
    if not os.path.exists(filename):
        raise IOError("given filename does not exist")
    savefile= open(filename,'rb')
    params= pickle.load(savefile)
    savefile.close()
    vos= numpy.array([s[0] for s in params])*_REFV0
    ros= numpy.array([s[1] for s in params])*_REFR0
    vpec= numpy.array([s[6] for s in params])*_PMSGRA*ros -vos#7 w/ dwarf
    vpecR= numpy.array([s[5] for s in params])*_VRSUN#6 w/ dwarf
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    axScatter, axHistx,axHisty= bovy_plot.scatterplot(vpecR,#vos/ros+vpec/ros,
                                                      vpec,'k,',levels=levels,
                          xlabel=r'$v_{R,\odot}\ [\mathrm{km\ s}^{-1}]$',
                          ylabel=r'$v_{\phi,\odot}-v_0\ [\mathrm{km\ s}^{-1}]$',
                          bins=31,
                          xrange=[-15.,0.],
                          yrange=[0.,35.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg',retAxes=True)
    #SBD10 value
    bovy_plot.bovy_plot([-20.,40.],[12.24,12.24],'--',color='0.5',overplot=True)
    bovy_plot.bovy_text(-4.,12.7,r'$\mathrm{SBD10}$')
    axHisty.plot([0.,100.],[12.24,12.24],'--',color='0.5')
    bovy_plot.bovy_plot([_VRSUN,_VRSUN],[-100.,100.],
                        '--',color='0.5',overplot=True)
    bovy_plot.bovy_text(_VRSUN-.75,7.,r'$\mathrm{SBD10}$',rotation=90.)
    axHistx.plot([_VRSUN,_VRSUN],[0.,100.],'--',color='0.5')
    #Reid / Brunthaler
    #bovy_plot.bovy_plot([_PMSGRA,_PMSGRA],[-10.,100.],
    #                    '--',color='0.5',overplot=True)
    #axHistx.plot([_PMSGRA,_PMSGRA],[0.,100.],'--',color='0.5')
    #bovy_plot.bovy_text(29.4,5.,r'$\mathrm{RB04}$')

    #Inset, closed orbit at the Sun
    lp= LogarithmicHaloPotential(normalize=1.)
    ep= EllipticalDiskPotential(phib=numpy.pi/2.,p=0.,tform=-100.,tsteady=-100.,twophio=0.045)
    o= Orbit([1.,0.,1.+10./220.,0.])
    oc= Orbit([1.,0.,1.,0.])
    ts= numpy.linspace(0.,4.*numpy.pi,1001)
    o.integrate(ts,[lp,ep])
    oc.integrate(ts,lp)
    left, bottom, width, height= 0.45, 0.45, 0.25, 0.25
    axInset= pyplot.axes([left,bottom,width,height])
    pyplot.sca(axInset)
    pyplot.plot(o._orb.orbit[:,0]*numpy.cos(o._orb.orbit[:,3]),
                o._orb.orbit[:,0]*numpy.sin(o._orb.orbit[:,3]),
                'k-')
    pyplot.plot(oc._orb.orbit[:,0]*numpy.cos(oc._orb.orbit[:,3]),
                oc._orb.orbit[:,0]*numpy.sin(oc._orb.orbit[:,3]),
                '--',color='0.6')
    pyplot.xlim(-1.5,1.5)
    pyplot.ylim(-1.5,1.5)
    #pyplot.xlabel(r'$x / R_0$')
    #pyplot.ylabel(r'$y / R_0$')
    nullfmt   = NullFormatter()         # no labels
    axInset.xaxis.set_major_formatter(nullfmt)
    axInset.yaxis.set_major_formatter(nullfmt)
    bovy_plot.bovy_end_print(plotfilename)

if __name__ == '__main__':
    plot_vpecvo(sys.argv[1],sys.argv[2])
