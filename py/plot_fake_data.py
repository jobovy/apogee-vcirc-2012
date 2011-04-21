import sys
import os, os.path
#import shutil
import cPickle as pickle
import math as m
import numpy as nu
from optparse import OptionParser
#import subprocess
#from galpy.orbit import Orbit
from galpy.df import dehnendf, shudf
from galpy.util import bovy_plot
_DEGTORAD= m.pi/180.
def plot_fake_data(parser):
    """
    NAME:
       plot_fake_data
    PURPOSE:
       plot the fake data created by create_fake_data
    INPUT:
       parser - from optparse
    OUTPUT:
       stuff as specified by the options
    HISTORY:
       2011-04-20 - Written - Bovy (NYU)
    """
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        sys.exit(-1)
    #Set up DF
    dfc= dehnendf(beta=options.beta,profileParams=(options.rd,options.rs,options.so),correct=True,niter=20)
    #Open pickle
    picklefile= open(args[0],'rb')
    out= pickle.load(picklefile)
    picklefile.close()
    if not options.vlosdname is None:
        #Plot distribution of vlos vs. d
        vloss=nu.array([o.vlos(obs=[1.,0.,0.,0.,1.,0.],ro=1.,vo=1.) \
                            for o in out]).flatten()
        ds= nu.array([o.dist(obs=[1.,0.,0.],ro=1.) \
                          for o in out]).flatten()
        #Also calculate the expected relation
        ntheory= 1001
        dx= nu.linspace(0.,2.,ntheory)
        vtheory= nu.zeros(ntheory)
        l= options.los*_DEGTORAD
        for ii in range(ntheory):
            R= nu.sqrt(1.+dx[ii]**2.-2.*dx[ii]*nu.cos(l))
            if 1./nu.cos(l) < dx[ii] and nu.cos(l) > 0.:
                phi= nu.pi-m.asin(dx[ii]/R*nu.sin(l))
            else:
                phi= m.asin(dx[ii]/R*nu.sin(l))
            vtheory[ii]= (R**options.beta-dfc.asymmetricdrift(R))*m.sin(phi+l)-m.sin(l) #last term is the LSR
        bovy_plot.bovy_print()
        bovy_plot.scatterplot(ds,vloss,'.',bins=options.scatterbins,
                              xlabel=r'$d / R_0$',
                              ylabel=r'$v_{los} / v_0$')
        bovy_plot.bovy_plot(dx,vtheory,overplot=True)
        bovy_plot.bovy_text(r'$l = %i^\circ$' % round(options.los),
                            top_left=True)
        bovy_plot.bovy_end_print(options.vlosdname)

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that holds the los (as a pickle)"
    parser = OptionParser(usage=usage)
    #Same options as create_fake_data and then some
    parser.add_option("--vlosdname",dest="vlosdname",
                      default=None,
                      help="Plot of the distribution of v_los and d")
    parser.add_option("--scatterbins",dest="scatterbins",type='int',
                      default=15,
                      help="Number of bins to use in scatterplot (in one dimension)")
    parser.add_option("-n","--nobjects",dest="nobjects",type='int',
                      default=500,
                      help="Number of objects to simulate")
    parser.add_option("--los",dest="los",type='float',
                      default=30.,
                      help="Galactic longitude of the LOS")
    parser.add_option("--beta",dest="beta",type='float',
                      default=0.,
                      help="logarithmic slope of rotation curve at R=1.")
    parser.add_option("--so",dest="so",type='float',
                      default=0.2,
                      help="Velocity dispersion at Ro")
    parser.add_option("--rd",dest="rd",type='float',
                      default=1./3.,
                      help="Disk scale-length")
    parser.add_option("--rs",dest="rs",type='float',
                      default=1,
                      help="Disk sigma_R scale-length")
    parser.add_option("--dmax",dest="dmax",type='float',
                      default=1.25,
                      help="Maximum distance")
    parser.add_option("--distuncertainty",dest="distuncertainty",type='float',
                      default=.1,
                      help="Relative distance uncertainty")
    parser.add_option("--local",action="store_true", 
                      default=False, dest="local",
                      help="Create a local sample")
    return parser

if __name__ == '__main__':
    plot_fake_data(get_options())

