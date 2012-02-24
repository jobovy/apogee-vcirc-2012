import sys
import numpy
from galpy.util import bovy_plot
from fitvc import _REFV0, _REFR0
def plot_rotcurves(plotfilename):
    rs= numpy.linspace(5.,16.,1001)
    #1: flat
    vo,ro= 0.92848419, 1.01213729
    vflat= rs**0.*vo*_REFV0
    #2: power-law
    vo, ro, beta= 0.97316236, 1.00337476, -0.17313119
    vpl= vo*(rs/_REFR0/ro)**beta*_REFV0
    #3: linear
    vo, ro, dvdr= 0.96853544, 1.00568151, -0.08273548
    vlinear= (vo+dvdr*(rs/_REFR0-1.))*_REFV0
    #4: quadratic
    vo, ro, dvdr, d2vdr2= 1.00872965, 1.03703529, -0.3736882, 0.5263341
    vquadratic= (vo+dvdr*(rs/_REFR0-1.)+d2vdr2*(rs/_REFR0-1.)**2.)*_REFV0
    #5: cubic
    vo, ro, dvdr, d2vdr2, d3vdr3= 0.99909976, 1.02583349, -0.35153006, 0.37377148, 0.11687969
    vcubic= (vo+dvdr*(rs/_REFR0-1.)+d2vdr2*(rs/_REFR0-1.)**2.\
                 +d3vdr3*(rs/_REFR0-1.)**3.)*_REFV0
    #Plot all
    bovy_plot.bovy_print(fig_width=8.)
    bovy_plot.bovy_plot(rs,vflat,'k-',
                        xlabel=r'$R\ [\mathrm{kpc}]$',
                        ylabel=r'$V_c\ [\mathrm{km\ s}^{-1}]$',
                        xrange=[0.,20.],
                        yrange=[150.,300.])
    bovy_plot.bovy_plot(rs,vpl,'k--',overplot=True)
    bovy_plot.bovy_plot(rs,vlinear,'k-.',overplot=True)
    bovy_plot.bovy_plot(rs,vquadratic,'k:',overplot=True)
    bovy_plot.bovy_plot(rs,vcubic,'k-.',overplot=True,color='red')
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    plot_rotcurves(sys.argv[1])
