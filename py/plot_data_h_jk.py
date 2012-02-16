import os, os.path
from galpy.util import bovy_plot
from readVclosData import readVclosData
OUTDIR= os.path.join(os.getenv('HOME'),'Desktop')
OUTDIR= '../tex/'
OUTEXT= 'ps'
def plot_data_h_jk():
    data= readVclosData()
    bovy_plot.bovy_print()
    bovy_plot.bovy_plot(data['J0MAG']-data['K0MAG'],
                        data['H0MAG'],'k,',
                        xlabel=r'$(J-K_s)_0\ [\mathrm{mag}]$',
                        ylabel=r'$H_0\ [\mathrm{mag}]$',
                        xrange=[0.4,1.6],
                        yrange=[5.,14.],
                        onedhists=True)
    bovy_plot.bovy_end_print(os.path.join(OUTDIR,'data_h_jk.'+OUTEXT))
    return None

if __name__ == '__main__':
    plot_data_h_jk()
