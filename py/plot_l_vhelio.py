import os, os.path
import cPickle as pickle
import math
import numpy
import fitsio
from scipy import interpolate
from galpy.df import dehnendf
from galpy.util import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
import apogee
#set up
nodups= True
postshutdown= True
betas, vcbetas= False, 210.
ext= 'png'
datafile= apogee.tools.apallPath(nodups=nodups)
data= fitsio.read(datafile,1)
#data cuts
data=data[(numpy.fabs(data['GLAT']) < 2.)*(numpy.fabs(data['GLON']) > 15.)]
data= data[((data['APOGEE_TARGET1'] & 2**9) == 0)] #no probable cluster members
indx= numpy.array(['STAR' in data['OBJTYPE'][ii] for ii in range(len(data))],dtype='bool')
data= data[indx]
if postshutdown:
    data= data[(data['MJD5'] > 55788)]
#Calculate means
plates= list(set(data['PLATE']))
nplates= len(plates)
l_plate= numpy.zeros(nplates)
avg_plate= numpy.zeros(nplates)
sig_plate= numpy.zeros(nplates)
siga_plate= numpy.zeros(nplates)
for ii in range(nplates):
    indx= (data['PLATE'] == plates[ii])
    l_plate[ii]= numpy.mean(data['GLON'][indx])
    avg_plate[ii]= numpy.mean(data['VHELIO'][indx])
    sig_plate[ii]= numpy.std(data['VHELIO'][indx])
    siga_plate[ii]= numpy.std(data['VHELIO'][indx])/numpy.sqrt(numpy.sum(indx))
#Prediction
if betas:
    pred_file= 'predict_l_vhelio_betas.sav'
else:
    pred_file= 'predict_l_vhelio.sav'
if os.path.exists(pred_file):
    pred_file= open(pred_file,'rb')
    ls= pickle.load(pred_file)
    avg_pred= pickle.load(pred_file)
    if betas:
        avg_pred2= pickle.load(pred_file)
        avg_pred3= pickle.load(pred_file)
    pred_file.close()
else:
    dfc= dehnendf(beta=0.,correct=True,niter=20)
    if betas:
        dfc2= dehnendf(beta=0.1,correct=True,niter=20)
        dfc3= dehnendf(beta=-0.1,correct=True,niter=20)
    nds= 51
    nls= 101
    ds= numpy.linspace(0.,10./8.,nds)
    ls= numpy.linspace(0.,360.,nls)
    avg_pred= numpy.zeros(nls)
    if betas:
        avg_pred2= numpy.zeros(nls)
        avg_pred3= numpy.zeros(nls)
    for ii in range(nls):
        meanvlos= 0.
        norm= 0.
        if betas:
            meanvlos2= 0.
            norm2= 0.
            meanvlos3= 0.
            norm3= 0.
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
            meanvlos+= surfmass*dfc.meanvT(R)*math.sin(theta+l)
            norm+= surfmass
            if betas:
                surfmass= dfc2.surfacemassLOS(d,l,deg=False)
                meanvlos2+= surfmass*dfc2.meanvT(R)*math.sin(theta+l)
                norm2+= surfmass
                surfmass= dfc3.surfacemassLOS(d,l,deg=False)
                meanvlos3+= surfmass*dfc3.meanvT(R)*math.sin(theta+l)
                norm3+= surfmass
        avg_pred[ii]= meanvlos/norm
        if betas:
            avg_pred2[ii]= meanvlos2/norm2
            avg_pred3[ii]= meanvlos3/norm3
            print avg_pred[ii], avg_pred2[ii], avg_pred3[ii]
        else:
            print avg_pred[ii]
    pred_file= open(pred_file,'wb')
    pickle.dump(ls,pred_file)
    pickle.dump(avg_pred,pred_file)
    if betas:
        pickle.dump(avg_pred2,pred_file)
        pickle.dump(avg_pred3,pred_file)
    pred_file.close()
left, bottom, width, height= 0.1, 0.3, 0.8, 0.6
axTop= pyplot.axes([left,bottom,width,height])
left, bottom, width, height= 0.1, 0.1, 0.8, 0.2
axBottom= pyplot.axes([left,bottom,width,height])
fig= pyplot.gcf()
fig.sca(axTop)
pyplot.ylabel(r'$\mathrm{Heliocentric\ velocity}\ [\mathrm{km\ s}^{-1}]$')
pyplot.xlim(0.,360.)
pyplot.ylim(-200.,200.)
nullfmt   = NullFormatter()         # no labels
axTop.xaxis.set_major_formatter(nullfmt)
bovy_plot.bovy_plot(data['GLON'],data['VHELIO'],'k,',
                    
                    yrange=[-200.,200.],
                    xrange=[0.,360.],overplot=True)
bovy_plot.bovy_plot(l_plate,avg_plate,'o',overplot=True,mfc='0.5',mec='none')
#pyplot.errorbar(l_plate,avg_plate,yerr=sig_plate,marker='o',mfc='k',mec='w',
#                ls='none',color='0.75')
#Solar motion
vsun= [-11.1,240.]
vsolar= numpy.zeros(len(ls))
for ii in range(len(ls)):
    l= ls[ii]/180.*math.pi
    vsolar[ii]= numpy.dot(vsun,numpy.array([-math.cos(l),math.sin(l)]))
if betas:
    line1= bovy_plot.bovy_plot(ls,vcbetas*avg_pred-vsolar,'k-',overplot=True)
    line2= bovy_plot.bovy_plot(ls,vcbetas*avg_pred2-vsolar,'k--',overplot=True)
    line3= bovy_plot.bovy_plot(ls,vcbetas*avg_pred3-vsolar,'k-.',overplot=True)
else:
    line1= bovy_plot.bovy_plot(ls,230.*avg_pred-vsolar,'k-',overplot=True)
    line2= bovy_plot.bovy_plot(ls,210.*avg_pred-vsolar,'k--',overplot=True)
    line3= bovy_plot.bovy_plot(ls,250.*avg_pred-vsolar,'k-.',overplot=True)
ndata_t= int(math.floor(len(data)/1000.))
ndata_h= len(data)-ndata_t*1000
if betas:
    bovy_plot.bovy_text(r'$|b|\ <\ 2^\circ,\ |l|\ >\ 15^\circ$'
                        +'\n'+r'$%i,%i\ \mathrm{stars}$' % (ndata_t,ndata_h)
                        +'\n'+r'$\mathrm{assuming}\ R_0\ =\ 8\ \mathrm{kpc}$'
                        +'\n'+r'$v_{\mathrm{circ}}\ = %i\ \mathrm{km\ s}^{-1}$' % (int(vcbetas)),
                        top_right=True)
    #Legend
    pyplot.legend((line1,line2,line3),(r'$\frac{\mathrm{d} v_{\mathrm{circ}}}{\mathrm{d}R}\ =\ \,\,0\,\,\,\,\ \mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}$',
                                       r'$\frac{\mathrm{d} v_{\mathrm{circ}}}{\mathrm{d}R}\ =\ \,\,2.75\ \mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}$',
                                       r'$\frac{\mathrm{d} v_{\mathrm{circ}}}{\mathrm{d}R}\ =\ -2.75\ \mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}$'),
                  loc='upper right',bbox_to_anchor=(.91,.375),
                  numpoints=2,
                  prop={'size':12},
                  frameon=False)
else:
    bovy_plot.bovy_text(r'$|b|\ <\ 2^\circ,\ |l|\ >\ 15^\circ$'
                        +'\n'+r'$%i,%03i\ \mathrm{stars}$' % (ndata_t,ndata_h)
                        +'\n'+r'$\mathrm{assuming}\ R_0\ =\ 8\ \mathrm{kpc}$'
                        +'\n'+r'$\frac{\mathrm{d} v_{\mathrm{circ}}}{\mathrm{d}R}\ =\ 0\ \mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}$',
                        top_right=True)
    #Legend
    pyplot.legend((line1,line2,line3),(r'$v_{\mathrm{circ}}\ =\ 230\ \mathrm{km\ s}^{-1}$',
                                       r'$v_{\mathrm{circ}}\ =\ 210\ \mathrm{km\ s}^{-1}$',
                                       r'$v_{\mathrm{circ}}\ =\ 250\ \mathrm{km\ s}^{-1}$'),
                  loc='upper right',bbox_to_anchor=(.9,.375),
                  numpoints=2,
                  prop={'size':12},
                  frameon=False)
bovy_plot._add_ticks()
fig.sca(axBottom)
#Interpolate prediction
interpolPred= interpolate.InterpolatedUnivariateSpline(ls,210.*avg_pred-vsolar)
bovy_plot.bovy_plot(l_plate,avg_plate-interpolPred(l_plate),'ko',overplot=True)
pyplot.errorbar(l_plate,avg_plate-interpolPred(l_plate),
                yerr=siga_plate,marker='o',color='k',ls='none')
bovy_plot.bovy_plot([0.,360.],[0.,0.],'k--',overplot=True)
if betas:
    bovy_plot.bovy_plot(ls,vcbetas*avg_pred2-vcbetas*avg_pred,'k-',overplot=True)
    bovy_plot.bovy_plot(ls,vcbetas*avg_pred3-vcbetas*avg_pred,'k-.',overplot=True)
else:
    bovy_plot.bovy_plot(ls,230.*avg_pred-210.*avg_pred,'k-',overplot=True)
    bovy_plot.bovy_plot(ls,250.*avg_pred-210.*avg_pred,'k-.',overplot=True)
pyplot.xlabel(r'$\mathrm{Galactic\ longitude}\ [\mathrm{deg}]$')
pyplot.ylabel(r'$\langle v_{\mathrm{los}}^{\mathrm{helio}}\rangle^{\mathrm{data}}-\langle v_{\mathrm{los}}^{\mathrm{helio}}\rangle^{\mathrm{model}}$')
if betas:
    bovy_plot.bovy_text(r'$\mathrm{flat\ for}\ \frac{\mathrm{d} v_{\mathrm{circ}}}{\mathrm{d}R}\ =\ 0\ \mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}$',
                    top_right=True)
else:
    bovy_plot.bovy_text(r'$\mathrm{flat\ for}\ v_{\mathrm{circ}}\ =\ 210\ \mathrm{km\ s}^{-1}$',
                        top_right=True)
pyplot.ylim(-14.5,14.5)
pyplot.xlim(0.,360.)
bovy_plot._add_ticks()
if betas:
    bovy_plot.bovy_end_print('apogee_vcirc_l_vhelio_betas.'+ext)
else:
    bovy_plot.bovy_end_print('apogee_vcirc_l_vhelio.'+ext)
