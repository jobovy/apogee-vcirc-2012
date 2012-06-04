import sys
import os, os.path
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
from readVclosData import readVclosData
from isomodel import isomodel
OUTDIR= os.path.join(os.getenv('HOME'),'Desktop')
OUTDIR= os.path.join(os.getenv('HOME'),'Desktop','hjkfigs')
OUTDIR= '../tex/'
#OUTDIR= '../figs/'
OUTEXT= 'ps'
_PLOTDISTANCE= True
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
                        xrange=[0.4,1.4],
                        yrange=[5.,14.],
                        onedhists=True,bins=31)
    #Overplot distance if wanted
    if _PLOTDISTANCE:
        iso= isomodel(imfmodel='lognormalChabrier2001',
                      Z=0.019,
                      expsfh=True)
        nds= 101
        ds= numpy.zeros((nds,nds))
        jks= numpy.linspace(0.5,1.2,nds)
        hs= numpy.linspace(14.,5.,nds)
        for ii in range(nds):
            for jj in range(nds):
                ds[ii,jj]= iso.peak(jks[ii],hs[jj])
        #Now contour this
        levels=[1.,3.,10.,30.]
        colors='0.6'#['0.5','0.5','0.5','k','k','k']
        CS=pyplot.contour(jks,hs,ds.T,levels=levels,
                          colors=colors,zorder=10.,linewidths=2.)
        ys= [5.3,6.7,9.22,11.7]
        for ii in range(len(levels)):
            bovy_plot.bovy_text(1.21,ys[ii],r'$%.0f\ \mathrm{kpc}$' % levels[ii],
                                fontsize=14.,color='0.3')
        if False:
            pyplot.clabel(CS, levels,
                          inline=1,
                          fmt='%.0f',
                          fontsize=14,
                          colors=colors,zorder=10.)
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    if len(sys.argv) > 2:
        plot_data_h_jk(location=int(sys.argv[1]),plotfilename=sys.argv[2])
    elif len(sys.argv) > 1:
        plot_data_h_jk(location=int(sys.argv[1]))
    else:
        plot_data_h_jk()
