###############################################################################
# plot_distanceprior: make a graphic representation of the distance prior
###############################################################################
import sys
import os, os.path
import cPickle as pickle
import copy
from optparse import OptionParser
import math
import numpy
from scipy import integrate, optimize, interpolate
from scipy.maxentropy import logsumexp
from matplotlib import pyplot
from galpy.util import bovy_plot
from readVclosData import readVclosData
import isomodel
from fitvc import _logpd, _dm, _REFR0, _DEGTORAD, \
    _BINTEGRATEDMIN, _BINTEGRATEDMAX, _BINTEGRATENBINS
def plot_distanceprior(parser):
    (options,args)= parser.parse_args()
    #Read the data
    print "Reading the data ..."
    data= readVclosData(postshutdown=options.postshutdown,
                        fehcut=options.fehcut,
                        cohort=options.cohort,
                        lmin=options.lmin,
                        bmax=options.bmax,
                        ak=True,
                        cutmultiples=options.cutmultiples,
                        jkmax=options.jkmax,
                        datafilename=options.fakedata)
    l= data['GLON']*_DEGTORAD
    b= data['GLAT']*_DEGTORAD
    sinl= numpy.sin(l)
    cosl= numpy.cos(l)
    sinb= numpy.sin(b)
    cosb= numpy.cos(b)
    jk= data['J0MAG']-data['K0MAG']
    jk[(jk < 0.5)]= 0.5 #BOVY: FIX THIS HACK BY EMAILING GAIL
    h= data['H0MAG']
    #Set up the isochrone
    print "Setting up the isochrone model ..."
    iso= isomodel.isomodel(imfmodel=options.imfmodel,Z=options.Z,
                           expsfh=options.expsfh)
    #Set up polar grid
    res= 31
    xgrid= numpy.linspace(0.,2.*math.pi*(1.-1./res/2.),
                       2*res)
    ygrid= numpy.linspace(0.5,2.2,res)
    plotxgrid= numpy.linspace(xgrid[0]-(xgrid[1]-xgrid[0])/2.,
                              xgrid[-1]+(xgrid[1]-xgrid[0])/2.,
                              len(xgrid)+1)
    plotygrid= numpy.linspace(ygrid[0]-(ygrid[1]-ygrid[0])/2.,
                              ygrid[-1]+(ygrid[1]-ygrid[0])/2.,
                              len(ygrid)+1)
    plotthis= numpy.zeros((2*res,res,len(data)))-numpy.finfo(numpy.dtype(numpy.float64)).max
    ds= numpy.linspace(_BINTEGRATEDMIN,_BINTEGRATEDMAX,_BINTEGRATENBINS)
    logpiso= numpy.zeros((len(data),_BINTEGRATENBINS))
    dm= _dm(ds)
    for ii in range(len(data)):
        mh= h[ii]-dm
        logpiso[ii,:]= iso(numpy.zeros(_BINTEGRATENBINS)+jk[ii],mh)
    for jj in range(_BINTEGRATENBINS):
        d= ds[jj]
        R= numpy.sqrt(1.+d**2.-2.*d*cosl)
        indx= (R == 0.)
        R[indx]+= 0.0001
        theta= numpy.arcsin(d/R*sinl)
        indx= (1./cosl < d)*(cosl > 0.)
        theta[indx]= numpy.pi-theta[indx]
        thisout= _logpd([0.,1.],d,None,None,
                        None,None,None,
                        options,R,theta,
                        1.,0.,logpiso[:,jj])
        #Find bin to which these contribute
        thetabin= numpy.floor((theta-plotxgrid[0])/(plotxgrid[1]-plotxgrid[0]))
        Rbin= numpy.floor((R-plotygrid[0])/(plotygrid[1]-plotygrid[0]))
        indx= (thetabin < 0)
        thetabin[indx]= 0
        Rbin[indx]= 0
        thisout[indx]= -numpy.finfo(numpy.dtype(numpy.float64)).max
        indx= (thetabin >= 2*res)
        thetabin[indx]= 0
        Rbin[indx]= 0
        thisout[indx]= -numpy.finfo(numpy.dtype(numpy.float64)).max
        indx= (Rbin < 0)
        thetabin[indx]= 0
        Rbin[indx]= 0
        thisout[indx]= -numpy.finfo(numpy.dtype(numpy.float64)).max
        indx= (Rbin >= res)
        thetabin[indx]= 0
        Rbin[indx]= 0
        thisout[indx]= -numpy.finfo(numpy.dtype(numpy.float64)).max
        thetabin= thetabin.astype('int')
        Rbin= Rbin.astype('int')
        for ii in range(len(data)):
            plotthis[thetabin,Rbin,ii]= thisout[ii]
    """
    for ii in range(2*res):
        for jj in range(res):
            #First calculate d
            R, theta= ygrid[jj],xgrid[ii]
            d= math.sqrt(R**2.+1.-2.*R*math.cos(theta))
            if 1./math.cos(theta) < R and math.cos(theta) > 0.:
                gridl= math.pi-math.asin(R/d*math.sin(theta))
            else:
                gridl= math.asin(R/d*math.sin(theta))
            dl= 0.2
            indx= (l > gridl -dl)*(l < gridl + dl)
            if numpy.sum(indx) == 0.: continue
            thesedata= data[indx]
            theseh= h[indx]
            thesejk= jk[indx]
            #Calculate isochrone probability
            logpiso= numpy.zeros(len(thesedata))
            if d == 0.: d= 0.0000001
            dm= _dm(d*_REFR0)
            for kk in range(len(thesedata)):
                mh= theseh[kk]-dm
                logpiso[kk]= iso(thesejk[kk],mh)
            plotthis[ii,jj]= logsumexp(_logpd([0.,1.],d,None,None,
                                              None,None,None,
                                              options,R,theta,
                                              1.,0.,logpiso))
    """                                         
    #Normalize
    for ii in range(2*res):
        for jj in range(res):
            plotthis[ii,jj,0]= logsumexp(plotthis[ii,jj,:])
    plotthis= plotthis[:,:,0]
    plotthis-= numpy.amax(plotthis)
    plotthis= numpy.exp(plotthis)
    bovy_plot.bovy_print()
    ax= pyplot.subplot(111,projection='galpolar')#galpolar is in bovy_plot
    vmin, vmax= 0., 1.
    out= ax.pcolor(plotxgrid,plotygrid,plotthis.T,cmap='gist_yarg',
                   vmin=vmin,vmax=vmax)
    bovy_plot.bovy_end_print(options.plotfile)
    

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the fit/samples will be saved to"
    parser = OptionParser(usage=usage)
    parser.add_option("-o","--plotfile",
                      dest='plotfile',default=None,
                      help="Plot filename")
    parser.add_option("--densmodel",dest='densmodel',default='expdisk',
                      help="Density model to use")
    parser.add_option("--hr",dest='hr',default=3.,type='float',
                      help="scale length in kpc")
    parser.add_option("--hz",dest='hz',default=0.25,type='float',
                      help="scale height in kpc")
    #Data options
    parser.add_option("--lmin",dest='lmin',default=35.,type='float',
                      help="readVclosData 'lmin'")
    parser.add_option("--bmax",dest='bmax',default=2.,type='float',
                      help="readVclosData 'bmax'")
    parser.add_option("--allshutdown",action="store_false", 
                      dest="postshutdown",
                      default=True,
                      help="setting this sets postshutdown to False in data")
    parser.add_option("--fehcut",action="store_true", 
                      dest="fehcut",
                      default=False,
                      help="readVclosData 'fehcut'")
    parser.add_option("--cohort",dest='cohort',default=None,
                      help="readVclosData 'cohort'")
    parser.add_option("--jkmax",dest='jkmax',default=1.2,type='float',
                      help="readVclosData 'jkmax'")
    parser.add_option("--location",dest='location',default=None,type='int',
                      help="location id when looking at single los")
    parser.add_option("--downsample",dest='downsample',default=None,
                      type='float',
                      help="Factor with which to downsample the data")
    parser.add_option("--cutmultiples",action="store_true", 
                      dest="cutmultiples",
                      default=False,
                      help="readVclosData 'cutmultiples'")
    parser.add_option("-f",dest='fakedata',default=None,
                      help="Name of the fake data filename")
    #Isochrone IMF
    parser.add_option("--imfmodel",dest='imfmodel',default='lognormalChabrier2001',
                      help="imfmodel for isochrone model")
    parser.add_option("--Z",dest='Z',default=.019,type='float',
                      help="Metallicity of isochrone")
    parser.add_option("--expsfh",action="store_true", dest="expsfh",
                      default=False,
                      help="If set, use an exponentially declining SFH")
    return parser

if __name__ == '__main__':
    plot_distanceprior(get_options())

