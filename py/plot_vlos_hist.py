import os, os.path
import cPickle as pickle
import math
import numpy
from galpy.orbit import Orbit
from galpy.df import dehnendf
from galpy.util import bovy_plot, save_pickles
from readVclosData import readVclosData
_DEGTORAD= numpy.pi/180.
_PREDICTVC= 220.
_PREDICTSIGMA= 0.2
_PREDICTCACHE= True
_PREDICTVSOLAR= [-11.1,240.]
_PREDICTNDS= 201
def simplePlot(location=4242,plotfile=None,predict=False,nv=11,dmax=10./8.,
               **kwargs):
    """
    NAME:
       simplePlot
    PURPOSE:
       make a simple histogram for a given location
    INPUT:
       location - location ID
       +readVclosData inputs
    OPTIONAL INPUT:
       plotfile= if set, save plot to this file
    OUTPUT:
       plot to display or file
    HISTORY:
       2012-01-25 - Written - Bovy (IAS)
    """
    #Read data
    data= readVclosData(**kwargs)
    data= data[(data['LOCATION'] == location)]
    if not plotfile is None:
        bovy_plot.bovy_print()
    range= [-200.,200.]
    hist, xvec, p= bovy_plot.bovy_hist(data['VHELIO'],range=range,bins=31,
                                       histtype='step',color='k',
                                       xlabel=r'$\mathrm{heliocentric}\ v_{\mathrm{los}}\ [\mathrm{km\ s}^{-1}]$')
    #Prediction
    if predict:
        pred_vs= numpy.linspace(range[0],range[1],nv)
        pred_dist= _calc_pred(pred_vs,location,numpy.mean(data['GLON']),dmax)
        data_int= numpy.sum(hist)*(xvec[1]-xvec[0])
        pred_dist*= data_int/numpy.sum(pred_dist)/(pred_vs[1]-pred_vs[0])
        print pred_dist
        bovy_plot.bovy_plot(pred_vs,pred_dist,'-',color='0.65',overplot=True)
    #Add text
    bovy_plot.bovy_text(r'$\mathrm{location}\ =\ %i$' % location
                        +'\n'
                        +r'$l\ \approx\ %.0f^\circ$' % numpy.mean(data['GLON']),
                        top_right=True,size=14.)
    if not plotfile is None:
        bovy_plot.bovy_end_print(plotfile)
    
def _calc_pred(vs,location,l,dmax):
    """Calculate the predicted distribution"""
    thisvs= vs
    if _PREDICTCACHE:
        savefilename= os.path.join(os.getenv('DATADIR'),'bovy','vclos',
                                   'predictedVclosDistribution_%i_%.0f_%i_%.0f_%.2f_%i_%.2f.sav' % (location,l,len(vs),_PREDICTVC,_PREDICTSIGMA,_PREDICTNDS,dmax))
    if _PREDICTCACHE and os.path.exists(savefilename):
        #Load prediction
        savefile= open(savefilename,'rb')
        pred_dist= pickle.load(savefile)
        savefile.close()
    else:
        ds= numpy.linspace(0.,dmax,_PREDICTNDS)
        dfc= dehnendf(beta=beta,profileParams=(1./3.,1.,_PREDICTSIGMA),correct=True,niter=20)
        thisl= l*_DEGTORAD
        thisvsolar= numpy.dot(_PREDICTVSOLAR,numpy.array([-numpy.cos(thisl),numpy.sin(thisl)]))
        thisvs+= thisvsolar
        thisvs/= _PREDICTVC
        ii= 0
        pred_dist= numpy.zeros(len(vs))
        while ii < len(vs):
            out= 0.
            norm= 0.
            for jj in range(_PREDICTNDS):
                d= ds[jj]
                R= math.sqrt(1.+d**2.-2.*d*math.cos(thisl))
                if R == 0.:
                    R+= 0.0001
                    d+= 0.0001
                if 1./math.cos(thisl) < d and math.cos(thisl) > 0.:
                    theta= math.pi-math.asin(d/R*math.sin(thisl))
                else:
                    theta= math.asin(d/R*math.sin(thisl))
                vR= -thisvs[ii]*math.cos(theta+thisl) #vperp=0 for marginalization
                vT= thisvs[ii]*math.sin(theta+thisl)
                o= Orbit([R,vR,vT,theta])
                surfmass= dfc.surfacemassLOS(d,thisl,deg=False)
                out+= surfmass*dfc(o,marginalizeVperp=True)
                norm+= surfmass
            pred_dist[ii]= out/norm
            ii+= 1
        if _PREDICTCACHE:
            savefile= open(savefilename,'wb')
            try:
                pickle.dump(pred_dist,savefile)
            finally:
                savefile.close()
    return pred_dist

if __name__ == '__main__':
    data= readVclosData()
    locations= list(set(data['LOCATION']))
    for l in locations:
        print l
        simplePlot(location=l,plotfile='/home/bovy/Desktop/locationPlots/vlosdist_%i.png' % l,predict=True,nv=201)
