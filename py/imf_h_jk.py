#Plot the IMF vs. H and J-K, sensibly
import sys
import os, os.path
import copy
import math
import numpy
from optparse import OptionParser
import isodist, isodist.imf
from galpy.util import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import FuncFormatter, MultipleLocator
OUTDIR= os.path.join(os.getenv('HOME'),'Desktop')
OUTDIR= '../tex/'
OUTEXT= 'ps'
def imf_h_jk(plotfile,Z=None,dwarf=False,log=False,h=12.,basti=False,
             dartmouth=False,kroupa=False):
    #Read isochrones
    if basti:
        zs= numpy.array([0.0001,0.0003,0.0006,0.001,0.002,0.004,0.008,
                         0.01,0.0198,0.03,0.04])
    elif dartmouth:
        zs= isodist.FEH2Z(numpy.array([-2.5,-2.,-1.5,-1.,-0.5,0.,0.2,0.3,0.5]))
    else:
        zs= numpy.arange(0.0005,0.03005,0.0005)
    if Z is None:
        Zs= zs
    elif not basti and not dartmouth:
        if Z < 0.01:
            Zs= [Z-0.001,Z-0.0005,Z,Z+0.0005,Z+0.001] #build up statistics
        else:
            Zs= [Z-0.0005,Z,Z+0.0005] #build up statistics
    else:
        Zs= [Z]
    if basti:
        p= isodist.BastiIsochrone(Z=Zs)
    elif dartmouth:
        p= isodist.DartmouthIsochrone(feh=isodist.Z2FEH(Zs),onlyold=True)
    else:
        p= isodist.PadovaIsochrone(Z=Zs)
    #Get relevant data
    sample= []
    weights= []
    for logage in p.logages():
        for z in Zs:
            thisiso= p(logage,z)
            if basti: mini= thisiso['M_ini']
            elif dartmouth: mini= thisiso['M']
            else: mini= thisiso['M_ini']
            if basti:
                int_IMF= isodist.imf.lognormalChabrier2001(thisiso['M_ini'],
                                                           int=True)
                dmpm= numpy.roll(int_IMF,-1)-int_IMF
            elif dartmouth:
                int_IMF= isodist.imf.lognormalChabrier2001(thisiso['M'],
                                                           int=True)
                dmpm= numpy.roll(int_IMF,-1)-int_IMF
            else:
                if kroupa:
                    int_IMF= isodist.imf.kroupa2003(thisiso['M_ini'],int=True)
                    dmpm= numpy.roll(int_IMF,-1)-int_IMF
                else:
                    dmpm= numpy.roll(thisiso['int_IMF'],-1)-thisiso['int_IMF']
            for ii in range(1,len(mini)-1):
                if basti:
                    JK= thisiso['J'][ii]-thisiso['K'][ii]
                else:
                    JK= thisiso['J'][ii]-thisiso['Ks'][ii]
                H= thisiso['H'][ii]
                if JK < 0.: # or thisiso['logg'][ii] > 3.5:
                    continue
                if dmpm[ii] > 0.: 
                    if basti:
                        sample.append([thisiso['J'][ii]-thisiso['K'][ii],
                                       thisiso['H'][ii]])
                    else:
                        sample.append([thisiso['J'][ii]-thisiso['Ks'][ii],
                                       thisiso['H'][ii]])
                    if dartmouth:
                        if logage > numpy.log10(5.)+9.:
                            weights.append(2.*dmpm[ii]) #Dartmouth are linearly spaced, but spacing is bigger at > 5 Gyr
                        else:
                            weights.append(dmpm[ii]) #Dartmouth are linearly spaced?
                    else:
                        weights.append(dmpm[ii]*10**(logage-7.))
                    #weights.append(dmpm[ii]*10**(logage-7.)*numpy.exp((10.**(logage-7.))/800.))
                else: 
                    continue #no use in continuing here
    #Form array
    sample= numpy.array(sample)
    weights= numpy.array(weights)
    #Histogram
    if dwarf:
        hist, edges= numpy.histogramdd(sample,weights=weights,bins=51,
                                       range=[[0.,1.6],[2.,9.]])
    else:
        hist, edges= numpy.histogramdd(sample,weights=weights,bins=49,
                                       range=[[0.3,1.6],[-11.,2]])
    #Normalize each J-K
    for ii in range(len(hist[:,0])):
        hist[ii,:]/= numpy.nanmax(hist[ii,:])/numpy.nanmax(hist)
        rev= copy.copy(hist[ii,::-1]) #reverse, but in one go does not always work
        hist[ii,:]= rev
    #Plot
    bovy_plot.bovy_print()
    if log:
        hist= numpy.log(hist)
    bovy_plot.bovy_dens2d(hist.T,origin='lower',cmap='gist_yarg',
                          xrange=[edges[0][0],edges[0][-1]],
                          yrange=[edges[1][-1],edges[1][0]],
                          aspect=(edges[0][-1]-edges[0][0])/float(edges[1][-1]-edges[1][0]),
                          xlabel=r'$(J-K_s)_0\ [\mathrm{mag}]$',
                          ylabel=r'$M_H\ [\mathrm{mag}]$',
                          interpolation='nearest')
    #Add extinction arrow
    djk= 0.45
    dh= 1.55/1.5*djk
    from matplotlib.patches import FancyArrowPatch
    ax=pyplot.gca()
    ax.add_patch(FancyArrowPatch((1.,-2.),(1+djk,-2+dh),
                                 arrowstyle='->',mutation_scale=20,fill=True,
                                 lw=1.25))
    bovy_plot.bovy_text(1.05,-2.05,r'$\mathrm{extinction}$',
                        rotation=-math.atan(1.5/1.55*1.3/13.)/math.pi*180.,
                        size=14.)
    #Add color cut
    bovy_plot.bovy_plot([0.5,0.5],[-20.,20.],'--',color='0.6',overplot=True)
    ax.add_patch(FancyArrowPatch((0.5,-6.),(0.7,-6.),
                                 arrowstyle='->',mutation_scale=20,fill=True,
                                 lw=1.25,ls='dashed',color='0.6'))
    bovy_plot.bovy_text(0.43,-8.,r'$\mathrm{APOGEE\ color\ cut}$',rotation=90.,
                        size=14.)
    #Add twin y axis
    ax= pyplot.gca()
    def my_formatter(x, pos):
        """distance in kpc for m=h"""
        xs= 10.**((h-x)/5.-2.)
        return r'$%.0f$' % xs
    def my_formatter2(x, pos):
        """distance in kpc for m=h"""
        xs= 10.**((h-x)/5.+1.)
        return r'$%.0f$' % xs
    ax2= pyplot.twinx()
    if dwarf:
        major_formatter = FuncFormatter(my_formatter2)
    else:
        major_formatter = FuncFormatter(my_formatter)
    ax2.yaxis.set_major_formatter(major_formatter)
    ystep= ax.yaxis.get_majorticklocs()
    ystep= ystep[1]-ystep[0]
    ax2.yaxis.set_minor_locator(MultipleLocator(ystep/5.))
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    ymin, ymax= ax.yaxis.get_view_interval()
    ax2.yaxis.set_view_interval(ymin,ymax,ignore=True)
    if dwarf:
        ax2.set_ylabel('$\mathrm{distance\ for}\ H_0\ =\ %.1f\ [\mathrm{pc}]$' % h)
    else:
        ax2.set_ylabel('$\mathrm{distance\ for}\ H_0\ =\ %.1f\ [\mathrm{kpc}]$' % h)
    xstep= ax.xaxis.get_majorticklocs()
    xstep= xstep[1]-xstep[0]
    ax2.xaxis.set_minor_locator(MultipleLocator(xstep/5.))
    if Z is None:
        bovy_plot.bovy_end_print(plotfile)
    else:
        bovy_plot.bovy_text(r'$Z\ =\ %.3f$' % Z,top_right=True,size=14.)
        bovy_plot.bovy_end_print(plotfile)
    return None

def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-o",dest='plotfile',default=None,
                      help="Name of file for plot")
    parser.add_option("-Z",dest='Z',type='float',default=None,
                      help="Metallicity Z")
    parser.add_option("--ho",dest='ho',type='float',default=12.,
                      help="H_0 to use for distance scale")
    parser.add_option("--dwarf",action="store_true", 
                      dest="dwarf",
                      default=False,
                      help="Show the dwarf part of the isochrone")
    parser.add_option("--log",action="store_true", 
                      dest="log",
                      default=False,
                      help="Use a logarithmic grayscale")
    return parser

if __name__ == '__main__':
    parser= get_options()
    (options,args)= parser.parse_args()
    imf_h_jk(options.plotfile,Z=options.Z,dwarf=options.dwarf,log=options.log,
             h=options.ho)
