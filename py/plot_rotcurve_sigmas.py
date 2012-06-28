import os, os.path
import cPickle as pickle
from optparse import OptionParser
import math
import numpy
from scipy import special
from fitvc import get_options, _DEGTORAD, _REFR0, _REFV0, _VRSUN, _PMSGRA
from galpy.util import bovy_plot
from matplotlib import pyplot
from plot_pdfs import set_options, load_samples
def plot_rotcurve_samples(options=None,args=None):
    """Plot sample rotation curves"""
    options= set_options(options)
    params= numpy.random.permutation(load_samples(args[0]))
    if options.removerooutliers:
        indx= numpy.zeros(len(params),dtype='bool')
        indx[:]= True
        for ii in range(len(params)):
            if params[ii][1] > (9./_REFR0): indx[ii]= False
        params= params[indx]
        options.nsamples= numpy.amin([options.nsamples,len(params)])
        print len(params)
    bestfitparams= load_samples(args[1])
    rotcurves= []
    rs= numpy.linspace(3.5,14.,1001)
    fs= numpy.zeros((len(rs),options.nsamples))
    besfitrc= numpy.zeros_like(rs)
    add= 0
    if options.nooutliermean: add-= 1
    if not options.dwarf: add-= 1
    if options.rotcurve.lower() == 'flat':
        bestfitrc= bestfitparams[0]*rs**0.*_REFV0
    elif options.rotcurve.lower() == 'powerlaw':
        bestfitrc= bestfitparams[0]*(rs/_REFR0/bestfitparams[1])**bestfitparams[6+add]*_REFV0
    elif options.rotcurve.lower() == 'linear':
        bestfitrc= (bestfitparams[0]+bestfitparams[6+add]*(rs/_REFR0/bestfitparams[1]-1.))*_REFV0
    elif options.rotcurve.lower() == 'quadratic':
        bestfitrc= (bestfitparams[0]+bestfitparams[6+add]*(rs/_REFR0/bestfitparams[1]-1.)+bestfitparams[7+add]*(rs/_REFR0/bestfitparams[1]-1)**2.)*_REFV0
    elif options.rotcurve.lower() == 'cubic':
        bestfitrc= (bestfitparams[0]+bestfitparams[6+add]*(rs/_REFR0/bestfitparams[1]-1)+bestfitparams[7+add]*(rs/_REFR0/bestfitparams[1]-1.)**2.+bestfitparams[8+add]*(rs/_REFR0/bestfitparams[1]-1.)**3.)*_REFV0
    for ii in range(options.nsamples):
        if options.rotcurve.lower() == 'flat':
            thisrc= params[ii][0]*rs**0.*_REFV0
        elif options.rotcurve.lower() == 'powerlaw':
            thisrc= params[ii][0]*(rs/_REFR0/params[ii][1])**params[ii][6+add]*_REFV0
        elif options.rotcurve.lower() == 'linear':
            thisrc= (params[ii][0]+params[ii][6+add]*(rs/_REFR0/params[ii][1]-1.))*_REFV0
        elif options.rotcurve.lower() == 'quadratic':
            thisrc= (params[ii][0]+params[ii][6+add]*(rs/_REFR0/params[ii][1]-1.)+params[ii][7+add]*(rs/_REFR0/params[ii][1]-1)**2.)*_REFV0
        elif options.rotcurve.lower() == 'cubic':
            thisrc= (params[ii][0]+params[ii][6+add]*(rs/_REFR0/params[ii][1]-1)+params[ii][7+add]*(rs/_REFR0/params[ii][1]-1.)**2.+params[ii][8+add]*(rs/_REFR0/params[ii][1]-1.)**3.)*_REFV0
        fs[:,ii]= thisrc
    #Now plot the mean and std-dev from the posterior
    rcmean= numpy.zeros(len(rs))
    nsigs= 3
    rcsigs= numpy.zeros((len(rs),2*nsigs))
    #Record mean and std-devs
    rcmean[:]= numpy.mean(fs,axis=1)
    bovy_plot.bovy_print(fig_width=8.)
    bovy_plot.bovy_plot(rs,bestfitrc,'k-',
                        xlabel=r'$R\ [\mathrm{kpc}]$',
                        ylabel=r'$V_c\ [\mathrm{km\ s}^{-1}]$',
                        xrange=[0.,20.],
                        yrange=[150.,300.])
    for ii in range(nsigs):
        for jj in range(len(rs)):
            thisf= sorted(fs[jj,:])
            thiscut= 0.5*special.erfc((ii+1.)/math.sqrt(2.))
            rcsigs[jj,2*ii]= thisf[int(math.floor(thiscut*options.nsamples))]
            thiscut= 1.-thiscut
            rcsigs[jj,2*ii+1]= thisf[int(math.floor(thiscut*options.nsamples))]
    colord, cc= (1.-0.75)/(nsigs+1.), 1
    nsigma= nsigs
    pyplot.fill_between(rs,rcsigs[:,0],rcsigs[:,1],color='0.75')
    while nsigma > 1:
        pyplot.fill_between(rs,rcsigs[:,cc+1],rcsigs[:,cc-1],
                            color='%f' % (.75+colord*cc))
        pyplot.fill_between(rs,rcsigs[:,cc],rcsigs[:,cc+2],
                            color='%f' % (.75+colord*cc))
        cc+= 2
        nsigma-= 1
    bovy_plot.bovy_plot(rs,bestfitrc,'k-',overplot=True)
    if options.plotm31:
        #Read file
        m31data= numpy.loadtxt('../data/m31.dat',comments='#')
        rm31= m31data[:,0]
        vcm31= m31data[:,1]
        bovy_plot.bovy_plot(rm31,vcm31,'ks',overplot=True,mfc='none',mew=2.)
        bovy_plot.bovy_text(17.,260.,r'$\mathrm{M31}$',size=14.)
        indx= (rm31 > 15.2)*(rm31 <= 16.8)
        bovy_plot.bovy_plot([17.,rm31[indx]],[260.,vcm31[indx]],'k-',
                            overplot=True)
    #Labels
    if options.rotcurve.lower() == 'flat':
        print "This is silly"
    elif options.rotcurve.lower() == 'powerlaw':
        bovy_plot.bovy_text(r'$\mathrm{Power\!\!-\!\!law\ rotation\ curve\ model}$',
                            top_left=True,size=16.)
    elif options.rotcurve.lower() == 'linear':
        bovy_plot.bovy_text(r'$\mathrm{Linearly rising\ rotation\ curve\ model}$',
                            top_left=True,size=16.)
    elif options.rotcurve.lower() == 'quadratic':
        bovy_plot.bovy_text(r'$\mathrm{Quadratic polynomial\ rotation\ curve\ model}$',
                            top_left=True,size=16.)
    elif options.rotcurve.lower() == 'cubic':
        bovy_plot.bovy_text(r'$\mathrm{Cubic\ polynomial\ rotation\ curve\ model}$',
                            top_left=True,size=16.)
    if options.removerooutliers:
        bovy_plot.bovy_text(r'$R_0 < 9\,\mathrm{kpc}$',
                            top_right=True,size=16.)
        
    if options.addmock2:
        bovy_plot.bovy_text(r'$\mathrm{Mock\ 2}$',
                            bottom_right=True,size=16.)
        
    bovy_plot.bovy_end_print(options.plotfilename)
    return None

def get_options():
    usage = "usage: %prog [options] <savefilename>\n\nsavefilename= name of the file that the fit/samples will be saved to"
    parser = OptionParser(usage=usage)
    #Initial conditions file
    parser.add_option("--init",dest='init',default=None,
                      help="Initial parameters")
    #Rotation curve parameters/model
    parser.add_option("--rotcurve",dest='rotcurve',default='flat',
                      help="Rotation curve model to fit")
    #Ro prior
    parser.add_option("--noroprior",action="store_true", dest="noroprior",
                      default=False,
                      help="If set, do not apply an Ro prior")
    #Sun's peculiar velocity
    parser.add_option("--fitvpec",action="store_true", dest="fitvpec",
                      default=False,
                      help="If set, fit for the peculiar velocity of the Sun as well, CURRENTLY ASSUMES flat rotation curve")
    #Velocity distribution model
    parser.add_option("--dfmodel",dest='dfmodel',default='simplegaussian',
                      help="DF model to use")
    parser.add_option("--nooutliermean",action="store_true", 
                      dest="nooutliermean",
                      default=False,
                      help="If set, use a zero mean for the outlier model (in Galactocentric coordinates)")
    parser.add_option("--fitsrinnerouter",
                      action="store_true", dest="fitsrinnerouter",
                      default=False,
                      help="If set, fit for the sigma_r separately for the inner disk")
    parser.add_option("--fitfehinnerouter",
                      action="store_true", dest="fitfehinnerouter",
                      default=False,
                      help="If set, fit for the feh offset separately for the inner disk")
    parser.add_option("--dwarfinnerouter",
                      action="store_true", dest="dwarfinnerouter",
                      default=False,
                      help="If set, fit for the dwarf contamination separately for the inner disk")
    parser.add_option("--fitahinnerouter",
                      action="store_true", dest="fitahinnerouter",
                      default=False,
                      help="If set, fit for the ah separately for the inner disk")
    parser.add_option("--fitdminnerouter",
                      action="store_true", dest="fitdminnerouter",
                      default=False,
                      help="If set, fit for the dm separately for the inner disk")
    parser.add_option("--fitsratio",action="store_true", dest="fitsratio",
                      default=False,
                      help="If set, fit for the ration squared of tangential to radial dispersion")
    parser.add_option("--fitsratioinnerouter",
                      action="store_true", dest="fitsratioinnerouter",
                      default=False,
                      help="If set, fit for the ration squared of tangential to radial dispersion")
    parser.add_option("--fiths",action="store_true", dest="fiths",
                      default=False,
                      help="If set, fit for a dispersion scale length offsett")
    #Isochrone IMF
    parser.add_option("--fitdm",action="store_true", dest="fitdm",
                      default=False,
                      help="If set, fit for a distance modulus offset")
    parser.add_option("--fitah",action="store_true", dest="fitah",
                      default=False,
                      help="If set, fit for an extinction-correction offset")
    parser.add_option("--fitfeh",action="store_true", dest="fitfeh",
                      default=False,
                      help="If set, fit for a [Fe/H] offset (with indivfeh)")
    #Add dwarf part?
    parser.add_option("--dwarf",action="store_true", 
                      dest="dwarf",
                      default=False,
                      help="setting this adds dwarf contamination")
    parser.add_option("--nsamples",dest='nsamples',default=10,type='int',
                      help="Number ofsamples to plot")
    parser.add_option("--removerooutliers",action="store_true", 
                      dest="removerooutliers",
                      default=False,
                      help="setting this removes samples with R0 > 9 kpc")
    parser.add_option("--plotm31",action="store_true", 
                      dest="plotm31",
                      default=False,
                      help="setting this plots the rotation curve for M31")
    parser.add_option("--addmock2",action="store_true", 
                      dest="addmock2",
                      default=False,
                      help="setting this will add a label that this is a mock (2)")
    #Output file
    parser.add_option("-o",dest='plotfilename',default=None,
                      help="Name of the file for the plot")
    return parser

if __name__ == '__main__':
    numpy.random.seed(1) #We need to seed to get, e.g., the same permutation when downsampling
    parser= get_options()
    (options,args)= parser.parse_args()
    plot_rotcurve_samples(options,args)

