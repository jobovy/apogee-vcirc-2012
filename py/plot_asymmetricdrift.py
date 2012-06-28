import sys
import os, os.path
import cPickle as pickle
import numpy
from scipy import interpolate
from asymmetricDriftModel import va, dlnnusR2dlnR, hR_fid, hs_fid
from galpy.util import bovy_plot
from matplotlib import pyplot
def plot_asymmetricdrift(outdir='../tex'):
    #First plot the fiducial model for vcva/sr2
    Ro= 8.
    hR= 3./Ro
    hs= 8./Ro
    vc= 1.
    #Zero-ish
    kinematicsFile= '../data/axi_va-kinematics_0.0.sav'
    kinematicsFile= open(kinematicsFile,'rb')
    vas= pickle.load(kinematicsFile) #vas has va/sigmaR^2(R_0)
    xs= pickle.load(kinematicsFile)
    kinematicsFile.close()
    vas= vas[15] # dispersion = small
    nrs= 101
    rs= numpy.linspace(0.1,2.2,nrs)
    vaInterp= interpolate.InterpolatedUnivariateSpline(rs,vas)
    #Plot
    rs= numpy.linspace(4.,20.,1001)
    so= 0.2
    #Hot
    sigmaR_hot= 0.2*numpy.exp(-(rs/Ro-1.)/hs)
    vas_fid= va(rs/Ro,0.2,hR=hR,hs=hs)
    #Cold
    dlnnusR2dlnR_correct= dlnnusR2dlnR(hR_fid,hs_fid)-dlnnusR2dlnR(hR,hs)
    vas_cold= vaInterp(rs/Ro)
    sigmaR_cold= 0.1*numpy.exp(-(rs/Ro-1.)/hs)
    vas_cold= 0.1**2./vc*(vas_cold+0.5*dlnnusR2dlnR_correct*numpy.exp(-2.*(rs/Ro-1.)))
    bovy_plot.bovy_print(fig_height=3.)
    lines= []
    lines.append(bovy_plot.bovy_plot(rs,vas_fid/sigmaR_hot**2.,'k-',
                                     xlabel=r'$R\ [\mathrm{kpc}]$',
                                     ylabel=r'$v_c\,(v_c-\bar{v}_T) / \sigma_R^2$',
                                     xrange=[0.,25.],
                                     yrange=[0.,6.5]))
    lines.append(bovy_plot.bovy_plot(rs,vas_cold/sigmaR_cold**2.,'--',color='0.4',
                                     overplot=True))
    bovy_plot.bovy_text(r'$h_R = 3\, \mathrm{kpc}$'+'\n'+
                        r'$h_\sigma = 8\, \mathrm{kpc}$',top_left=True)
    labels=[r'$\sigma_R(R_0) = 0.2\,v_c(R_0)$',
            r'$\sigma_R(R_0) = 0.1\,v_c(R_0)$']
    l1= pyplot.legend(lines,labels,loc=4,
                      frameon=False,numpoints=1)
    bovy_plot.bovy_end_print(os.path.join(outdir,'vaR.ps'))

    #Now plot differences from fiducial
    bovy_plot.bovy_print(fig_height=3.)
    bovy_plot.bovy_plot([rs[0],rs[-1]],[0.,0.],'k-',
                        xlabel=r'$R\ [\mathrm{kpc}]$',
                        ylabel=r'$(v_a - v_a^{\mathrm{fid}}) / v_c(R_0)$',
                        xrange=[0.,25.],
                        yrange=[-0.1,0.1])
    #hR=2 kpc
    bovy_plot.bovy_plot(rs,va(rs/Ro,0.2,hR=2./Ro,hs=hs)-vas_fid,overplot=True,
                        ls='-.',color='k')
    bovy_plot.bovy_text(0.5,0.02,r'$h_R = 2\, \mathrm{kpc}$')
    #hR=4 kpc
    bovy_plot.bovy_plot(rs,va(rs/Ro,0.2,hR=4./Ro,hs=hs)-vas_fid,overplot=True,
                        ls='--',color='k')
    bovy_plot.bovy_text(0.5,-0.04,r'$h_R = 4\, \mathrm{kpc}$')
    #hs=6 kpc
    bovy_plot.bovy_plot(rs,va(rs/Ro,0.2,hR=hR,hs=6./Ro)-vas_fid,overplot=True,
                        ls=':',color='k')
    bovy_plot.bovy_text(.5,0.07,r'$h_\sigma = 6\, \mathrm{kpc}$')
    #hs=6 kpc
    bovy_plot.bovy_plot(rs,va(rs/Ro,0.14,hR=hR,hs=242./Ro)-vas_fid,overplot=True,
                        ls='-',color='0.65')
    bovy_plot.bovy_text(17.5,0.05,r'$R_0/h_\sigma = 0.034$')
    bovy_plot.bovy_end_print(os.path.join(outdir,'vaR_diffs.ps'))
    

if __name__ == '__main__':
    plot_asymmetricdrift(sys.argv[1])
