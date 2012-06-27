import sys
import numpy
from galpy.util import bovy_plot
from fitvc import _REFV0, _REFR0
_PLOTM31= True
def plot_rotcurves(plotfilename):
    rs= numpy.linspace(3.5,16.,1001)
    #1: flat
    vo,ro= 0.92772766, 1.00939298
    vflat= rs**0.*vo*_REFV0
    #2: power-law
    vo, ro, beta= 0.92881531,  1.0009845, 0.00750639
    vpl= vo*(rs/_REFR0/ro)**beta*_REFV0
    #3: linear
    vo, ro, dvdr= 0.93050113,  1.00377579, 0.01625223
    vlinear= (vo+dvdr*(rs/_REFR0-1.))*_REFV0
    #4: quadratic
    vo, ro, dvdr, d2vdr2= 0.92954738,  1.00440273, 0.02909439,  0.25129582
    vquadratic= (vo+dvdr*(rs/_REFR0-1.)+d2vdr2*(rs/_REFR0-1.)**2.)*_REFV0
    #5: cubic
    vo, ro, dvdr, d2vdr2, d3vdr3=  0.9416574, 1.01691786, -0.10708624, -0.07733835,  0.42241186
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
    if _PLOTM31:
        #Read file
        m31data= numpy.loadtxt('../data/m31.dat',comments='#')
        rm31= m31data[:,0]
        vcm31= m31data[:,1]
        bovy_plot.bovy_plot(rm31,vcm31,'ks',overplot=True,mfc='none',mew=2.)
        bovy_plot.bovy_text(17.,260.,r'$\mathrm{M31}$',size=14.)
        indx= (rm31 > 15.2)*(rm31 <= 16.8)
        bovy_plot.bovy_plot([17.,rm31[indx]],[260.,vcm31[indx]],'k-',
                            overplot=True)
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    plot_rotcurves(sys.argv[1])
