import sys
import os, os.path
import cPickle as pickle
import math
import numpy
from scipy import special
from galpy.util import bovy_plot
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
    bovy_plot.bovy_end_print(plotfilename)

if __name__ == '__main__':
    plot_vpecvo(sys.argv[1],sys.argv[2])
