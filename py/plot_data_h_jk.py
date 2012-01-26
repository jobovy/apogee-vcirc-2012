import os, os.path
from galpy.util import bovy_plot
from readVclosData import readVclosData
OUTDIR= os.path.join(os.getenv('HOME'),'Desktop')
OUTDIR= '../tex/'
OUTEXT= 'ps'
def plot_data_h_jk():
    data= readVclosData()
    bovy_plot.bovy_print()
    bovy_plot.bovy_plot(data['JMAG']-data['KMAG'],
                        data['HMAG'],'k,',
                        xlabel=r'$J-K_s\ [\mathrm{mag}]$',
                        ylabel=r'$H\ [\mathrm{mag}]$',
                        xrange=[0.5,1.8],
                        yrange=[5.,14.],
                        onedhists=True)
    bovy_plot.bovy_end_print(os.path.join(OUTDIR,'data_h_jk.'+OUTEXT))
    return None

if __name__ == '__main__':
    plot_data_h_jk()
