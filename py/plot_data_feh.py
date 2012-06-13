import sys
import os, os.path
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
from readVclosData import readVclosData
from isomodel import isomodel
OUTDIR= os.path.join(os.getenv('HOME'),'Desktop')
OUTDIR= os.path.join(os.getenv('HOME'),'Desktop','fehfigs')
#OUTDIR= '../tex/'
#OUTDIR= '../figs/'
OUTEXT= 'png'
def plot_data_feh(location=0,
                  plotfilename=os.path.join(OUTDIR,'data_h_jk.'+OUTEXT)):
    data= readVclosData()
    #Good feh
    data= data[(data['FEH']!= -9999.00)]
    if not location is None:
        if location == 0:
            locs= set(data['LOCATION'])
            for l in locs:
                plot_data_feh(location=l,
                              plotfilename=os.path.join(OUTDIR,'data_feh_%i.' % l +OUTEXT))
            return None
        data= data[(data['LOCATION'] == location)]
    meanfeh= numpy.mean(data['FEH'])
    sigfeh= numpy.std(data['FEH'])
    bovy_plot.bovy_print()
    bovy_plot.bovy_hist(data['FEH'],
                        xlabel=r'$[\mathrm{Fe/H}]$',
                        xrange=[meanfeh-1.,meanfeh+1],
                        bins=16)
    bovy_plot.bovy_text(r'$\sigma = %.2f$' % sigfeh,top_right=True)
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    if len(sys.argv) > 2:
        plot_data_feh(location=int(sys.argv[1]),plotfilename=sys.argv[2])
    elif len(sys.argv) > 1:
        plot_data_feh(location=int(sys.argv[1]))
    else:
        plot_data_feh()
