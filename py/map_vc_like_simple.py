_EPSREL=10.**-14.
import sys
import os, os.path
#import shutil
import cPickle as pickle
import math as m
import numpy as nu
from scipy import linalg, integrate
from scipy.maxentropy import logsumexp
from optparse import OptionParser
from matplotlib import pyplot
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
    if options.linearfit:
        plot_linear(out,options.los*_DEGTORAD,options,dfc)
        return None
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

def plot_linear(out,l,options,dfc):
    #Calculate sin(phi+l) for each
    ndata= len(out)
    sinphil= nu.zeros(ndata)
    vloss= nu.zeros(ndata)
    for ii in range(ndata):
        o= out[ii]
        sinphil[ii]= m.sin(o.phi()+l)
        vloss[ii]= o.vlos(obs=[1.,0.,0.,0.,1.,0.],ro=1.,vo=1.)
    #Perform linear fit
    vlosf= []
    sphilf= []
    syf= []
    nsphilbins, minobjs= 150, 5
    for ii in range(nsphilbins):
        if ii == 0.:
            sphilmin, sphilmax= -1., -1.+1./nsphilbins*2.
        elif ii == (nsphilbins - 1):
            sphilmin, sphilmax= 1.-1./nsphilbins*2., 1.
        else:
            sphilmin, sphilmax= -1.+float(ii)/nsphilbins*2., -1.+2.*(ii+1.)/nsphilbins
        indx= (sinphil <= sphilmax)*(sinphil > sphilmin)
        if len(set(indx)) == 1: continue
        if len(vloss[indx]) < minobjs: continue
        sphilf.append(nu.mean(sinphil[indx]))
        vlosf.append(nu.mean(vloss[indx]))
        syf.append(nu.var(vloss[indx])/len(vloss[indx]))
    pfit= _my_polyfit(sphilf,vlosf,sy=syf)
    #Get asymmetricdrift
    va= dfc.asymmetricdrift(1.)
    vsun= -m.sin(l)
    #Fit for sigma_R and sigma_T
    svlos= []
    sphil= []
    sy= []
    nsphilbins, minobjs= 50, 10
    for ii in range(nsphilbins):
        if ii == 0.:
            sphilmin, sphilmax= 0., 1./nsphilbins
        elif ii == (nsphilbins - 1):
            sphilmin, sphilmax= 1.-1./nsphilbins, 1.
        else:
            sphilmin, sphilmax= float(ii)/nsphilbins, (ii+1.)/nsphilbins
        indx= (sinphil**2. <= sphilmax)*(sinphil**2.> sphilmin)
        if len(set(indx)) == 1: continue
        if len(vloss[indx]) < minobjs: continue
        sphil.append(nu.mean(sinphil[indx]**2.))
        svlos.append(nu.var(vloss[indx]))
        sy.append(nu.var(vloss[indx])**2./(len(vloss[indx]-1)))
    #And fit that
    sfit= _my_polyfit(sphil,svlos,sy=sy)
    bovy_plot.bovy_print()
    fig= pyplot.figure()
    pyplot.subplot(211)
    bovy_plot.bovy_plot(sphilf,vlosf,'k.',overplot=True)
    pyplot.xlim(0.5,1.)
    pyplot.ylim(0.,.6)
    pyplot.xlabel(r'$sin(\phi + l)$')
    pyplot.ylabel(r'$v_{los} / v_0$')
    bovy_plot.bovy_plot(nu.array([-1.,1.]),
                        nu.array([pfit[1]-pfit[0],pfit[1]+pfit[0]]),'k-',overplot=True)
    bovy_plot.bovy_text(r'$v_c - v_a = %4.2f \ v_0$' % pfit[0]\
                            +'\n'+
                        r'$\mathrm{expected}\ v_c - v_a = %4.2f\ v_0$' % (1.-va)\
                            +'\n'+
                        r'$v_{local}\,\sin l = %4.2f\ v_0$' % (-pfit[1])\
                            +'\n'+
                        r'$\mathrm{expected}\ v_{local}\,\sin l = %4.2f$'\
                            % (-vsun)\
                            +'\n'+
                        r'$\sigma_R(R_0) = %4.2f \ v_0$' % options.so+'\n'+\
                            r'$l  = %i^\circ$' % round(options.los),
                        top_left=True)
    fig.subplots_adjust(hspace=0.3)
    pyplot.subplot(210)
    bovy_plot.bovy_plot(sphil,svlos,'k.',overplot=True)
    pyplot.xlabel(r'$sin^2(\phi + l)$')
    pyplot.ylabel(r'$\sigma^2_{v_{los}} / v_0^2$')
    bovy_plot.bovy_text(r'$\sigma_R = %4.2f \ v_0$' % (nu.sqrt(sfit[1]))\
                            +'\n'+
                        r'$\mathrm{expected}\ \sigma_R = %4.2f\ v_0$' % (options.so)\
                            +'\n'+
                        r'$\sigma_T = %4.2f \ v_0$' % (nu.sqrt(sfit[0]+sfit[1]))\
                            +'\n'+
                        r'$\mathrm{expected}\ \sigma_T = %4.2f\ v_0$' % (options.so/dfc._gamma),
                        top_left=True)
    pyplot.ylim(0.,0.1)
    bovy_plot.bovy_plot(nu.array([0.,1.]),
                        nu.array([sfit[1],sfit[1]+sfit[0]]),'k-',overplot=True)
    bovy_plot.bovy_end_print(options.plotfilename)
    return None

def _my_polyfit(x,y,sy=None):
    #Put the dat in the appropriate arrays and matrices
    nsample= len(x)
    Y= nu.zeros(nsample)
    A= nu.ones((nsample,2))
    C= nu.eye(nsample)
    for jj in range(nsample):
        Y[jj]= y[jj]
        A[jj,1]= x[jj]
        if not sy is None: C[jj,jj]= sy[jj]
    #Now compute the best fit and the uncertainties
    bestfit= nu.dot(linalg.inv(C),Y.T)
    bestfit= nu.dot(A.T,bestfit)
    bestfitvar= nu.dot(linalg.inv(C),A)
    bestfitvar= nu.dot(A.T,bestfitvar)
    bestfitvar= linalg.inv(bestfitvar)
    bestfit= nu.dot(bestfitvar,bestfit)
    return bestfit[::-1]

def single_vlos_loglike(vc,o,dfc,options,l,beta=0.):
    """Log likelihood of a single los velocity"""
    #Currently we only do distuncertain = 0.
    if options.distuncertainty == 0.:
        R= o.R()
        phi= o.phi()
        if options.fixdfmoments:
            sigmaR2= dfc.targetSigma2(1.)
            sigmaT2= sigmaR2/dfc._gamma
            vtmean= vc*R**beta-dfc.asymmetricdrift(1.)
        else:
            sigmaR2= dfc.targetSigma2(R)
            sigmaT2= sigmaR2/dfc._gamma
            vtmean= vc*R**beta-dfc.asymmetricdrift(R)
        vrmean= 0.
        vlosmean= m.sin(phi+l)*vtmean-vc*m.sin(l)
        vlossigma2= m.sin(phi+l)**2.*(sigmaT2-sigmaR2)+sigmaR2
        return -0.5*(o.vlos(obs=[1.,0.,0.,0.,1.,0.],ro=1.,vo=1.)-vlosmean)**2./vlossigma2
    elif options.distuncertainty == 10.: #infinity
        sigmaR2= dfc.targetSigma2(1.)
        sigmaT2= sigmaR2/dfc._gamma
        out= integrate.quad(_marginalizedIntegrand,
                              0.01,1.25,
                              (l,dfc,vc,beta,
                               o.vlos(obs=[1.,0.,0.,0.,1.,0.],ro=1.,vo=1.),
                               options),
                              epsrel=_EPSREL)[0]
        return m.log(out)

def _marginalizedIntegrand(d,l,dfc,vc,beta,vlos,options):
    pd= dfc.surfacemassLOS(d,l,deg=False)
    R= nu.sqrt(1.+d**2.-2.*d*m.cos(l))
    if 1./m.cos(l) < d and m.cos(l) > 0.:
        phi= m.pi-m.asin(d/R*m.sin(l))
    else:
        phi= m.asin(d/R*m.sin(l))
    if options.fixdfmoments:
        sigmaR2= dfc.targetSigma2(1.)
        sigmaT2= sigmaR2/dfc._gamma
        vtmean= vc*R**beta-dfc.asymmetricdrift(1.)
    else:
        sigmaR2= dfc.targetSigma2(R)
        sigmaT2= sigmaR2/dfc._gamma
        vtmean= vc*R**beta-dfc.asymmetricdrift(R)
    vlosmean= m.sin(phi+l)*vtmean-vc*m.sin(l)
    vlossigma2= m.sin(phi+l)**2.*(sigmaT2-sigmaR2)+sigmaR2
    return m.exp(-0.5*(vlos-vlosmean)**2./vlossigma2)/m.sqrt(vlossigma2)*pd*d

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
    parser.add_option("--linearfit",action="store_true", 
                      default=False, dest="linearfit",
                      help="Perform a simple linear fit assuming a flat rotation curve")
    parser.add_option("--fixdfmoments",action="store_true", 
                      default=False, dest="fixdfmoments",
                      help="Fix the moments of the DF to the values on the Solar circle")
    return parser

if __name__ == '__main__':
    map_vc_like_simple(get_options())
