import sys
import os, os.path
from galpy.util import bovy_plot
from readVclosData import readVclosData
OUTDIR= os.path.join(os.getenv('HOME'),'Desktop')
OUTDIR= os.path.join(os.getenv('HOME'),'Desktop','hjkfigs')
#OUTDIR= '../tex/'
OUTEXT= 'png'
def plot_data_h_jk(location=None,
                   plotfilename=os.path.join(OUTDIR,'data_h_jk.'+OUTEXT)):
    data= readVclosData()
    if not location is None:
        if location == 0:
            locs= set(data['LOCATION'])
            for l in locs:
                plot_data_h_jk(location=l,
                               plotfilename=os.path.join(OUTDIR,'data_h_jk_%i.' % l +OUTEXT))
            return None
        data= data[(data['LOCATION'] == location)]
    bovy_plot.bovy_print()
    bovy_plot.bovy_plot(data['J0MAG']-data['K0MAG'],
                        data['H0MAG'],'k,',
                        xlabel=r'$(J-K_s)_0\ [\mathrm{mag}]$',
                        ylabel=r'$H_0\ [\mathrm{mag}]$',
                        xrange=[0.4,1.6],
                        yrange=[5.,14.],
                        onedhists=True,bins=31)
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    if len(sys.argv) > 2:
        plot_data_h_jk(location=int(sys.argv[1]),plotfilename=sys.argv[2])
    elif len(sys.argv) > 1:
        plot_data_h_jk(location=int(sys.argv[1]))
    else:
        plot_data_h_jk()
