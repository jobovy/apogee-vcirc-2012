import sys
import os, os.path
import cPickle as pickle
import math
import numpy
from scipy import special
from galpy.util import bovy_plot
from fitvc import _REFV0, _REFR0
def plot_rovo(filename,plotfilename):
    if not os.path.exists(filename):
        raise IOError("given filename does not exist")
    savefile= open(filename,'rb')
    params= pickle.load(savefile)
    savefile.close()
    vos= numpy.array([s[0] for s in params])*_REFV0
    ros= numpy.array([s[1] for s in params])*_REFR0
    bovy_plot.bovy_print()
    levels= list(special.erf(0.5*numpy.arange(1,4)))
    levels.append(1.01) #HACK to not plot outliers
    bovy_plot.scatterplot(vos/ros,ros,'k,',levels=levels,
                          xlabel=r'$\Omega_0\ [\mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}]$',
                          ylabel=r'$R_0\ [\mathrm{kpc}]$',
                          bins=31,
                          xrange=[200./8.,250./8.],
                          yrange=[7.,9.],
                          contours=True,
                          cntrcolors='k',
                          onedhists=True,
                          cmap='gist_yarg')
    bovy_plot.bovy_end_print(plotfilename)

if __name__ == '__main__':
    plot_rovo(sys.argv[1],sys.argv[2])
