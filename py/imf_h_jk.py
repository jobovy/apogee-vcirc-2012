#Plot the IMF vs. H and J-K, sensibly
import sys
import os, os.path
import numpy
import isodist
from galpy.util import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import FuncFormatter, MultipleLocator
OUTDIR= os.path.join(os.getenv('HOME'),'Desktop')
OUTDIR= '../tex/'
OUTEXT= 'ps'
def imf_h_jk(Z=None,h=12.):
    #Read isochrones
    zs= numpy.arange(0.0005,0.03005,0.0005)
    if Z is None:
        Zs= zs
    else:
        if Z < 0.01:
            Zs= [Z-0.001,Z-0.0005,Z,Z+0.0005,Z+0.001] #build up statistics
        else:
            Zs= [Z-0.0005,Z,Z+0.0005] #build up statistics
    p= isodist.PadovaIsochrone(Z=Zs)
    #Get relevant data
    sample= []
    weights= []
    for logage in p.logages():
        for z in Zs:
            thisiso= p(logage,z)
            dmpm= numpy.roll(thisiso['M_ini'],-1)-thisiso['M_ini']
            for ii in range(1,len(thisiso['M_ini'])-1):
                JK= thisiso['J'][ii]-thisiso['Ks'][ii]
                H= thisiso['H'][ii]
                if JK < 0.5 or H > 3.5:
                    continue
                if dmpm[ii] > 0.: 
                    sample.append([thisiso['J'][ii]-thisiso['Ks'][ii],
                                   thisiso['H'][ii]])
                    weights.append(dmpm[ii]*10**(logage-7.))
                else: 
                    continue #no use in continuing here
    #Form array
    sample= numpy.array(sample)
    weights= numpy.array(weights)
    #Histogram
    hist, edges= numpy.histogramdd(sample,weights=weights,bins=36,
                                   range=[[0.5,1.6],[-11.,2]])
    #Normalize each J-K
    for ii in range(len(hist[:,0])):
        hist[ii,:]/= numpy.nanmax(hist[ii,:])/numpy.nanmax(hist)
    hist[:,:]= hist[:,::-1]
    #Plot
    bovy_plot.bovy_print()
    bovy_plot.bovy_dens2d(hist.T,origin='lower',cmap='gist_yarg',
                          xrange=[edges[0][0],edges[0][-1]],
                          yrange=[edges[1][-1],edges[1][0]],
                          aspect=(edges[0][-1]-edges[0][0])/float(edges[1][-1]-edges[1][0]),
                          xlabel=r'$J-K_s\ [\mathrm{mag}]$',
                          ylabel=r'$M_H\ [\mathrm{mag}]$',
                          interpolation='nearest')
    #Add twin y axis
    ax= pyplot.gca()
    def my_formatter(x, pos):
        """distance in kpc for m=12.2"""
        xs= 10.**((h-x)/5.-2.)
        return r'$%.0f$' % xs
    ax2= pyplot.twinx()
    major_formatter = FuncFormatter(my_formatter)
    ax2.yaxis.set_major_formatter(major_formatter)
    ystep= ax.yaxis.get_majorticklocs()
    ystep= ystep[1]-ystep[0]
    ax2.yaxis.set_minor_locator(MultipleLocator(ystep/5.))
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    ymin, ymax= ax.yaxis.get_view_interval()
    ax2.yaxis.set_view_interval(ymin,ymax,ignore=True)
    ax2.set_ylabel('$\mathrm{distance\ for}\ H_0\ =\ %.1f\ [\mathrm{kpc}]$' % h)
    xstep= ax.xaxis.get_majorticklocs()
    xstep= xstep[1]-xstep[0]
    ax2.xaxis.set_minor_locator(MultipleLocator(xstep/5.))
    if Z is None:
        bovy_plot.bovy_end_print(os.path.join(OUTDIR,'imf_h_jk.'+OUTEXT))
    else:
        bovy_plot.bovy_text(r'$Z\ =\ %.3f$' % Z,top_right=True,size=14.)
        bovy_plot.bovy_end_print(os.path.join(OUTDIR,'imf_h_jk_%.4f.' % Z
                                              +OUTEXT))
    return None

if __name__ == '__main__':
    if len(sys.argv) > 1:
        imf_h_jk(Z=float(sys.argv[1]))
    else:
        imf_h_jk()
                          
