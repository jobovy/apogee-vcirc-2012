import numpy
#from readVclosData import readVclosData
from fitvc import get_options, _DEGTORAD, _REFR0, _REFV0
from galpy.util import bovy_plot
def set_options(options=options):
    if options is None:
        parser= get_options()
        options, args= parser.parse_args([])
def load_samples(filename):
    if not os.path.exists(filename):
        raise IOError("given filename does not exist")
    savefile= open(filename,'rb')
    params= pickle.load(savefile)
    savefile.close()
    return params
def rovo(filename=None,options=None):
    options= set_options(options=options)
    params= load_samples(filename)
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
    return None
