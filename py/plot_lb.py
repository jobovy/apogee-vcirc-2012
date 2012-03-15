#Plot the distribution of the sample in l and b
import sys
import os, os.path
import numpy
from galpy.df import dehnendf
from galpy.util import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
from readVclosData import readVclosData
from theoryDistributions import surfacemassLOSd
OUTDIR= os.path.join(os.getenv('HOME'),'Desktop')
OUTDIR= '../tex/'
OUTEXT= 'ps'
hz= 250.
ro= 8000.
dmax= 10000.
def plot_lb():
    #Read basic data
    data= readVclosData()
    ndata= len(data)
    yrange= [-6.,6.]
    #Plot
    if OUTEXT == 'ps':
        bovy_plot.bovy_print(fig_width=8.5,fig_height=6./16.*8.5)
    else:
        bovy_plot.bovy_print(fig_width=16.,fig_height=6.)
    fig= pyplot.figure()
    #Set up axes
    nullfmt   = NullFormatter()         # no labels
    # definitions for the axes
    left, width = 0.1, 0.7
    bottom, height = 0.1, 0.85
    bottom_h = left_h = left+width
    rect_scatter = [left, bottom, width, height]
    rect_histy = [left_h, bottom, 0.15, height]
    axScatter = pyplot.axes(rect_scatter)
    axHisty = pyplot.axes(rect_histy)
    # no labels
    axHisty.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    fig.sca(axScatter)  
    bovy_plot.bovy_plot(data['GLON'],data['GLAT'],'k.',ms=1.5,
                        overplot=True)
    pyplot.xlim(0.,250.)
    pyplot.ylim(yrange[0],yrange[1])
    pyplot.xlabel(r'$\mathrm{Galactic\ longitude\ [deg]}$')
    pyplot.ylabel(r'$\mathrm{Galactic\ latitude\ [deg]}$')
    bovy_plot._add_ticks()
    #bovy_plot.bovy_end_print(os.path.join(OUTDIR,'data_lb.'+OUTEXT))
    #Read other data
    otherdata= readVclosData(meanb=-4.,bmax=10.)
    nmdata= len(otherdata)
    otherl= list(otherdata['GLON'])
    otherb= list(otherdata['GLAT'])
    otherdata2= readVclosData(meanb=4.,bmax=10.)
    npdata= len(otherdata2)
    otherl.extend(list(otherdata2['GLON']))
    otherb.extend(list(otherdata2['GLAT']))
    bovy_plot.bovy_plot(otherl,otherb,'.',mec='0.5',mfc='0.5',ms=1.5,
                        overplot=True)
    #Add histogram of bs
    otherl.extend(list(data['GLON']))
    otherb.extend(list(data['GLAT']))
    histy, edges, patches= axHisty.hist(otherb,
                                        bins=49,
                                        orientation='horizontal',
                                        normed=True,
                                        histtype='step',
                                        range=yrange,
                                        color='k',zorder=2)
    #Overlay predictions
    #Center, constant
    ys= numpy.linspace(-1.5,1.5,1001)
    pdf= numpy.sqrt(1.5**2.-ys**2.)
    pdf*= float(ndata)/(ndata+nmdata+npdata)/numpy.sum(pdf)/(ys[1]-ys[0])
    axHisty.plot(pdf,ys,'k-',zorder=-2)
    constpdf= pdf
    #Center, exponential disk
    """
    dfc= dehnendf(beta=0.,correct=True,niter=20)
    pdf= numpy.zeros(len(ys))
    plates= list(set(data['PLATE']))
    for ii in range(len(plates)):
        print ii
        #For each plate, calculate pdf and add with data weights
        thisdata= data[(data['PLATE'] == plates[ii])]
        thisl= numpy.mean(data['GLON'])
        thisndata= len(thisdata)
        thispdf= numpy.zeros(len(ys))
        for jj in range(len(ys)):
            thispdf[jj]= surfacemassLOSd(ys[jj],thisl,0.,dmax/ro,0.,dfc,hz/ro)
        pdf+= thispdf/numpy.sum(thispdf)/(ys[1]-ys[0])*thisndata
    pdf*= float(ndata)/(ndata+nmdata+npdata)/numpy.sum(pdf)/(ys[1]-ys[0])
    axHisty.plot(pdf,ys,'-',color='0.',lw=1.,zorder=-1.)
    #axHisty.plot(.8*pdf+.2*constpdf,ys,'-',color='0.',lw=1.,zorder=-1.)
    """
    #positive b, constant
    constpdf*= float(npdata)/ndata
    axHisty.plot(constpdf,ys+4.,'k-',zorder=-2)
    """
    #Positive, exponential disk
    pdfp= numpy.zeros(len(ys))
    plates= list(set(otherdata2['PLATE']))
    for ii in range(len(plates)):
        print ii
        #For each plate, calculate pdf and add with data weights
        thisdata= otherdata2[(otherdata2['PLATE'] == plates[ii])]
        thisl= numpy.mean(otherdata2['GLON'])
        thisndata= len(thisdata)
        thispdf= numpy.zeros(len(ys))
        for jj in range(len(ys)):
            thispdf[jj]= surfacemassLOSd(4.+ys[jj],
                                          thisl,0.,dmax/ro,4.,dfc,hz/ro)
        pdfp+= thispdf/numpy.sum(thispdf)/(ys[1]-ys[0])*thisndata
    pdfp*= float(npdata)/(ndata+nmdata+npdata)/numpy.sum(pdfp)/(ys[1]-ys[0])
    axHisty.plot(pdfp,4.+ys,'-',color='0.',lw=1.,zorder=-1)
    #axHisty.plot(.8*pdfp+.2*constpdf,4.+ys,'-',color='0.',lw=1.,zorder=-1)
    """
    #negative b, constant
    constpdf*= nmdata/float(npdata)
    axHisty.plot(constpdf,ys-4.,'k-',zorder=-2.)
    """
    #Negative, exponential disk
    pdfm= numpy.zeros(len(ys))
    plates= list(set(otherdata['PLATE']))
    for ii in range(len(plates)):
        print ii
        #For each plate, calculate pdf and add with data weights
        thisdata= otherdata[(otherdata['PLATE'] == plates[ii])]
        thisl= numpy.mean(otherdata['GLON'])
        thisndata= len(thisdata)
        thispdf= numpy.zeros(len(ys))
        for jj in range(len(ys)):
            thispdf[jj]= surfacemassLOSd(-4.+ys[jj],
                                          thisl,0.,dmax/ro,-4.,dfc,hz/ro)
        pdfm+= thispdf/numpy.sum(thispdf)/(ys[1]-ys[0])*thisndata
    pdfm*= float(nmdata)/(ndata+nmdata+npdata)/numpy.sum(pdfm)/(ys[1]-ys[0])
    axHisty.plot(pdfm,-4.+ys,'-',color='0.',lw=1.,zorder=-1)
    #axHisty.plot(.8*pdfm+.2*constpdf,-4.+ys,'-',color='0.',lw=1.,zorder=-1)
    """
    axHisty.set_ylim(yrange[0],yrange[1])
    axHisty.set_xlim( 0, 1.2*numpy.amax(histy))
    pyplot.sca(axHisty)
    bovy_plot._add_ticks()
    bovy_plot.bovy_end_print(os.path.join(OUTDIR,'data_lb.'+OUTEXT))
    return None

if __name__ == '__main__':
    plot_lb()
