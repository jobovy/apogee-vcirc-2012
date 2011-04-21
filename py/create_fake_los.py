import sys
import os, os.path
#import shutil
import cPickle as pickle
import math as m
import numpy as nu
from optparse import OptionParser
#import subprocess
from galpy.orbit import Orbit
from galpy.df import dehnendf, shudf
def create_fake_los(parser):
    """
    NAME:
       create_fake_los
    PURPOSE:
       create a fake APOGEE LOS
    INPUT:
       parser - from optparse
    OUTPUT:
       stuff as specified by the options
    HISTORY:
       2011-04-20 - Written - Bovy (NYU)
    """
    nu.random.seed(seed=1)
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        sys.exit(-1)
    #Set up DF
    print "Setting up DF ..."
    dfc= dehnendf(beta=options.beta,profileParams=(options.rd,options.rs,options.so),correct=True,niter=20)
    #Sample until we have nobjects
    if options.local:
        sample_local(options,args,dfc)
    else:
        sample_los(options,args,dfc)
    return None

def sample_local(options,args,dfc):
    print "Sampling locally ..."
    #We just sample the velocity distribution at (R,phi)= (1,0)
    out= dfc.sampleVRVT(1.,n=options.nobjects)
    #No uncertainty for now
    #Save
    picklefile= open(args[0],'wb')
    pickle.dump(out,picklefile)
    picklefile.close()
    return None

def sample_los(options,args,dfc):
    print "Sampling the LOS ..."
    n= 0
    out= []
    while n < options.nobjects:
        thisos= dfc.sampleLOS(los=options.los,n=options.nobjects-n)
        #add distance uncertainty
        thisout= []
        for o in thisos:
            basico= Orbit([o.R(),o.vR(),o.vT(),0.,0.,o.phi()])#HACK bc ll does not work in 2D currently
            thisd= basico.dist(obs=[1.,0.,0.],ro=1.)
            #Add uncertainty logarithmically, like photometric distance uncertainty TO DO
            thisd+= nu.random.normal()*thisd*options.distuncertainty
            thisout.append(Orbit([basico.ll(obs=[1.,0.,0.],ro=1.)[0],
                                  thisd[0],
                                  basico.pmll(obs=[1.,0.,0.,0.,1.,0.],
                                              ro=1.,vo=1.)[0],
                                  basico.vlos(obs=[1.,0.,0.,0.,1.,0.],
                                              ro=1.,vo=1.)[0]],
                                 lb=True,
                                 solarmotion=[0.,0.,0.],ro=1.,vo=1.))
        ds= nu.array([o.dist(obs=[1.,0.,0.],ro=1.) for o in thisout]).flatten()
        thisos= [thisout[ii] for ii in range(len(thisout)) \
                     if (ds[ii] < options.dmax)]
        out.extend(thisos)
        n+= len(thisos)
    if len(out) > options.nobjects:
        out= out[0:options.nobjects-1]
    #Save
    picklefile= open(args[0],'wb')
    pickle.dump(out,picklefile)
    picklefile.close()
    return None                   

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the los will be saved to (as a pickle)"
    parser = OptionParser(usage=usage)
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
    create_fake_los(get_options())
