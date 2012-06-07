import sys
import os, os.path
import cPickle as pickle
import copy
import multi
import multiprocessing
import math
import numpy
from scipy import special
from scipy.maxentropy import logsumexp
from galpy.util import bovy_plot
from fitvc import _REFV0, _REFR0
import logl, plot_pdfs
_ANALYTIC=False
_MULTI=15
def plot_rovo(filename,plotfilename):
    if not os.path.exists(filename):
        raise IOError("given filename does not exist")
    savefile= open(filename,'rb')
    params= pickle.load(savefile)
    savefile.close()
    if _ANALYTIC: #Calculate by fixing everything except for Ro anv vo
        options= plot_pdfs.set_options(None)
        nros= 15
        noos= 15
        ros= numpy.linspace(7.,13.,nros)
        oos= numpy.linspace(20.,30.,noos)
        ll= numpy.zeros((noos,nros))
        for ii in range(noos):
            if not _MULTI is None:
                theseparamss= []
                for jj in range(nros):
                    theseparams= copy.copy(params)
                    theseparams[0]= oos[ii]*ros[jj]/_REFV0
                    theseparams[1]= ros[jj]/_REFR0
                    theseparamss.append(theseparams)
                thisll= multi.parallel_map((lambda x: numpy.sum(logl.logl(init=theseparamss[x],options=options))),
                                           range(nros),
                                           numcores=numpy.amin([nros,_MULTI,multiprocessing.cpu_count()]))
                ll[ii,:]= thisll
            else:
                for jj in range(nros):
                    theseparams= copy.copy(params)
                    theseparams[0]= oos[ii]*ros[jj]/_REFV0
                    theseparams[1]= ros[jj]/_REFR0
                    ll[ii,jj]= numpy.sum(logl.logl(init=theseparams,
                                                   options=options))
        #Normalize
        ll-= logsumexp(ll)
        ll= numpy.exp(ll)
        levels= list(special.erf(0.5*numpy.arange(1,4)))
        bovy_plot.bovy_dens2d(ll.T,origin='lower',levels=levels,
                              xlabel=r'$\Omega_0\ [\mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}]$',
                              ylabel=r'$R_0\ [\mathrm{kpc}]$',
                              xrange=[20.,35.],
                              yrange=[7.,13.],
                              contours=True,
                              cntrcolors='k',
                              onedhists=True,
                              cmap='gist_yarg')
    else:
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
