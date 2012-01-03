import os, os.path
import cPickle as pickle
import math
import numpy
import fitsio
from galpy.df import dehnendf
from galpy.util import bovy_plot
from matplotlib import pyplot
import apogee
#set up
nodups= False
datafile= apogee.tools.aprvallPath(nodups=nodups)
data= fitsio.read(datafile,1)
data=data[(numpy.fabs(data['GLAT']) < 2.)*(numpy.fabs(data['GLON']) > 15.)]
#Calculate means
plates= list(set(data['PLATE']))
nplates= len(plates)
l_plate= numpy.zeros(nplates)
avg_plate= numpy.zeros(nplates)
sig_plate= numpy.zeros(nplates)
for ii in range(nplates):
    indx= (data['PLATE'] == plates[ii])
    l_plate[ii]= numpy.mean(data['GLON'][indx])
    avg_plate[ii]= numpy.mean(data['VHELIO'][indx])
    sig_plate[ii]= numpy.std(data['VHELIO'][indx])
#Prediction
pred_file= 'predict_l_vhelio.sav'
if os.path.exists(pred_file):
    pred_file= open(pred_file,'rb')
    ls= pickle.load(pred_file)
    avg_pred= pickle.load(pred_file)
    pred_file.close()
else:
    dfc= dehnendf(beta=0.,correct=True,niter=20)
    nds= 51
    nls= 101
    ds= numpy.linspace(0.,10./8.,nds)
    ls= numpy.linspace(0.,360.,nls)
    avg_pred= numpy.zeros(nls)
    for ii in range(nls):
        meanvlos= 0.
        norm= 0.
        for jj in range(nds):
            d,l= ds[jj], ls[ii]/180.*numpy.pi
            R= math.sqrt(1.+d**2.-2.*d*math.cos(l))
            if R == 0.:
                R+= 0.0001
                d+= 0.0001
            if 1./math.cos(l) < d and math.cos(l) > 0.:
                theta= math.pi-math.asin(d/R*math.sin(l))
            else:
                theta= math.asin(d/R*math.sin(l))
            surfmass= dfc.surfacemassLOS(d,l,deg=False)
            meanvlos+= surfmass*(dfc.meanvT(R)*math.sin(theta+l)+dfc.meanvR(R)*math.cos(theta+l))
            norm+= surfmass
        avg_pred[ii]= meanvlos/norm
        print avg_pred[ii]
    pred_file= open(pred_file,'wb')
    pickle.dump(ls,pred_file)
    pickle.dump(avg_pred,pred_file)
    pred_file.close()
bovy_plot.bovy_plot(data['GLON'],data['VHELIO'],'k,',
                    xlabel=r'$\mathrm{Galactic\ longitude}\ [\mathrm{deg}]$',
                    ylabel=r'$\mathrm{Heliocentric\ velocity}\ [\mathrm{km\ s}^{-1}]$',
                    yrange=[-200.,200.],
                    xrange=[0.,360.])
bovy_plot.bovy_plot(l_plate,avg_plate,'bo',overplot=True)
pyplot.errorbar(l_plate,avg_plate,yerr=sig_plate,marker='o',color='b',ls='none')
#Solar motion
vsun= [-11.1,240.]
vsolar= numpy.zeros(len(ls))
for ii in range(len(ls)):
    l= ls[ii]/180.*math.pi
    vsolar[ii]= numpy.dot(vsun,numpy.array([math.cos(l),math.sin(l)]))
bovy_plot.bovy_plot(ls,235.*avg_pred-vsolar,'k-',overplot=True)
bovy_plot.bovy_plot(ls,220.*avg_pred-vsolar,'k--',overplot=True)
bovy_plot.bovy_plot(ls,250.*avg_pred-vsolar,'k-.',overplot=True)
bovy_plot.bovy_end_print('test.png')
