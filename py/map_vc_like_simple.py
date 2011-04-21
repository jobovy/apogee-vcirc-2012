import sys
import os, os.path
#import shutil
import cPickle as pickle
import math as m
import numpy as nu
from scipy.maxentropy import logsumexp
from optparse import OptionParser
#import subprocess
from galpy.orbit import Orbit
from galpy.df import dehnendf, shudf
from galpy.util import bovy_plot
_DEGTORAD= nu.pi/180.
def map_vc_like_simple(parser):
    """
    NAME:
       map_vc_like_simple
    PURPOSE:
       map the vc likelihood assuming knowledge of the DF
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
    dfc= dehnendf(beta=0.,profileParams=(options.rd,options.rs,options.so),
                  correct=True,niter=20)
    #Load data
    picklefile= open(args[0],'rb')
    out= pickle.load(picklefile)
    picklefile.close()
    ndata= len(out)
    #Map likelihood
    vcirc= nu.linspace(options.vmin,options.vmax,options.nvcirc)
    if not options.nbeta is None:
        betas= nu.linspace(options.betamin,options.betamax,options.nbeta)
        like= nu.zeros((options.nvcirc,options.nbeta))
        for ii in range(options.nvcirc):
            for kk in range(options.nbeta):
                thislike= 0.
                for jj in range(ndata):
                    thislike+= single_vlos_loglike(vcirc[ii],out[jj],dfc,
                                                   options,
                                                   options.los*_DEGTORAD,
                                                   beta=betas[kk])
                like[ii,kk]= thislike
        like-= logsumexp(like.flatten())+m.log(vcirc[1]-vcirc[0])
        bovy_plot.bovy_print()
        bovy_plot.bovy_dens2d(nu.exp(like).T,
                              origin='lower',
                              xrange=[options.vmin,options.vmax],
                              yrange=[options.betamin,options.betamax],
                              aspect=(options.vmax-options.vmin)/\
                                  (options.betamax-options.betamin),
                              cmap='gist_yarg',
                              xlabel=r'$v_c / v_0$',
                              ylabel=r'$\beta$',
                              contours=True,cntrmass=True,
                              levels=[0.682,0.954,0.997])
        bovy_plot.bovy_text(r'$\sigma_R(R_0) = %4.2f \ v_0$' % options.so\
                                +'\n'+\
                                r'$l  = %i^\circ$' % round(options.los),
                            top_left=True)
        bovy_plot.bovy_end_print(options.plotfilename)
    else:
        like= nu.zeros(options.nvcirc)           
        for ii in range(options.nvcirc):
            thislike= 0.
            for jj in range(ndata):
                thislike+= single_vlos_loglike(vcirc[ii],out[jj],dfc,options,
                                               options.los*_DEGTORAD)
            like[ii]= thislike
        like-= logsumexp(like)+m.log(vcirc[1]-vcirc[0])
        #Calculate mean and sigma
        vcmean= nu.sum(vcirc*nu.exp(like)*(vcirc[1]-vcirc[0]))
        vc2mean= nu.sum(vcirc**2.*nu.exp(like)*(vcirc[1]-vcirc[0]))
        #Plot
        bovy_plot.bovy_print()
        bovy_plot.bovy_plot(vcirc,nu.exp(like),'k-',xlabel=r'$v_c / v_0$',
                            ylabel=r'$p(\mathrm{data} | v_c)$')
        bovy_plot.bovy_text(r'$\langle v_c \rangle = %4.2f \ v_0$' % vcmean +'\n'+
                            r'$\sqrt{\langle v_c^2 \rangle - \langle v_c \rangle^2} = %4.2f \ v_0$' % (m.sqrt(vc2mean-vcmean**2.)) +'\n'+\
                                r'$\sigma_R(R_0) = %4.2f \ v_0$' % options.so+'\n'+\
                                r'$l  = %i^\circ$' % round(options.los),
                            top_left=True)
        bovy_plot.bovy_end_print(options.plotfilename)

def single_vlos_loglike(vc,o,dfc,options,l,beta=0.):
    """Log likelihood of a single los velocity"""
    #Currently we only do distuncertain = 0.
    if options.distuncertainty == 0.:
        R= o.R()
        phi= o.phi()
        sigmaR2= dfc.targetSigma2(R)
        sigmaT2= sigmaR2/dfc._gamma
        vtmean= vc*R**beta-dfc.asymmetricdrift(R)
        vrmean= 0.
        AT= nu.zeros((2,2))
        AT[0,0]= m.sin(phi+l)
        AT[0,1]= m.cos(phi+l)
        AT[1,0]= -m.cos(phi+l)
        AT[1,1]= m.sin(phi+l)
        vlosmean= nu.dot(AT,nu.array([vrmean,vtmean]))[1]-vc*m.sin(l)
        V= nu.zeros((2,2))
        V[0,0]= sigmaR2
        V[1,1]= sigmaT2
        vlossigma2= nu.dot(AT,nu.dot(V,AT.T))[1,1]
        return -0.5*(o.vlos(obs=[1.,0.,0.,0.,1.,0.],ro=1.,vo=1.)-vlosmean)**2./vlossigma2

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that holds the los data (as a pickle)"
    parser = OptionParser(usage=usage)
    parser.add_option("-o",dest="plotfilename",
                      default=None,
                      help="Name of the file that will hold the plot")
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
    parser.add_option("--distuncertainty",dest="distuncertainty",type='float',
                      default=.1,
                      help="Relative distance uncertainty")
    parser.add_option("--vmin",dest="vmin",type='float',
                      default=0.5,
                      help="Minimum circular velocity to consider")
    parser.add_option("--vmax",dest="vmax",type='float',
                      default=1.5,
                      help="Maximum circular velocity to consider")   
    parser.add_option("--nvcirc",dest="nvcirc",type='int',
                      default=50.,
                      help="Number of circular velocities to consider")
    parser.add_option("--betamin",dest="betamin",type='float',
                      default=-0.2,
                      help="Minimum circular velocity to consider")
    parser.add_option("--betamax",dest="betamax",type='float',
                      default=0.2,
                      help="Maximum circular velocity to consider")   
    parser.add_option("--nbeta",dest="nbeta",type='int',
                      default=None,
                      help="Number of circular velocities to consider")
    return parser

if __name__ == '__main__':
    map_vc_like_simple(get_options())
