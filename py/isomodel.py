###############################################################################
# isomodel.py: isochrone model for the distribution in (J-Ks,M_H)
###############################################################################
import math
import numpy
from scipy import maxentropy
import scipy.interpolate
import isodist, isodist.imf
try:
    from galpy.util import bovy_plot
    _BOVY_PLOT_LOADED= True
except ImportError:
    _BOVY_PLOT_LOADED= False   
class isomodel:
    """isomodel: isochrone model for the distribution in (J-Ks,M_H)"""
    def __init__(self,dwarf=False,imfmodel='lognormalChabrier2001',Z=None,
                 interpolate=False):
        """
        NAME:
           __init__
        PURPOSE:
           initialize isomodel
        INPUT:
           Z= metallicity (if not set, use flat prior in Z over all Z; can be list)
           dwarf= (default: False) if True, use dwarf part
           imfmodel= (default: 'lognormalChabrier2001') IMF model to use (see isodist.imf code for options)
           interpolate= (default: False) if True, interpolate the binned representation
        OUTPUT:
           object
        HISTORY:
           2012-02-17 - Written - Bovy (IAS)
        """
       #Read isochrones
        zs= numpy.arange(0.0005,0.03005,0.0005)
        if Z is None:
            Zs= zs
        elif isinstance(Z,float):
            if Z < 0.01:
                Zs= [Z-0.001,Z-0.0005,Z,Z+0.0005,Z+0.001] #build up statistics
            else:
                Zs= [Z-0.0005,Z,Z+0.0005] #build up statistics
        p= isodist.PadovaIsochrone(Z=Zs)
        #Get relevant data
        sample= []
        weights= []
        for logage in p.logages():
            for z in Zs:
                thisiso= p(logage,z,asrecarray=True)
                #Calculate int_IMF for this IMF model
                if not imfmodel == 'lognormalChabrier2001': #That would be the default
                    if imfmodel == 'exponentialChabrier2001':
                        thisiso.int_IMF= isodist.imf.exponentialChabrier2001(thisiso.M_ini,int=True)
                    elif imfmodel == 'kroupa2003':
                        thisiso.int_IMF= isodist.imf.kroupa2003(thisiso.M_ini,int=True)
                    elif imfmodel == 'chabrier2003':
                        thisiso.int_IMF= isodist.imf.chabrier2003(thisiso.M_ini,int=True)
                    else:
                        raise IOError("imfmodel option not understood (non-existing model)")
                dN= numpy.roll(thisiso.int_IMF,-1)-thisiso.int_IMF
                for ii in range(1,len(thisiso.int_IMF)-1):
                    JK= thisiso.J[ii]-thisiso.Ks[ii]
                    H= thisiso.H[ii]
                    if JK < 0.5: # or thisiso['logg'][ii] > 3.5:
                        continue
                    if dN[ii] > 0.: 
                        sample.append([JK,H])
                        weights.append(dN[ii]*10**(logage-7.))
                    else: 
                        continue #no use in continuing here   
        #Form array
        sample= numpy.array(sample)
        weights= numpy.array(weights)
        #Histogram
        self._jkmin, self._jkmax= 0.5,1.6
        if dwarf:
            self._hmin, self._hmax= 2.,9.
            self._nbins= 51
        else: 
            self._hmin, self._hmax= -11.,2.
            self._nbins= 41
        self._djk= (self._jkmax-self._jkmin)/float(self._nbins)
        self._dh= (self._hmax-self._hmin)/float(self._nbins)
        self._hist, self._edges= numpy.histogramdd(sample,weights=weights,
                                                   bins=self._nbins,
                                                   range=[[self._jkmin,self._jkmax],
                                                          [self._hmin,self._hmax]])
        #Save
        self._Zs= Zs
        self._interpolate= interpolate
        self._dwarf= dwarf
        self._loghist= numpy.log(self._hist)
        self._loghist[(self._hist == 0.)]= -numpy.finfo(numpy.dtype(numpy.float64)).max
        if interpolate:
            #Form histogram grid
            jks= numpy.linspace(self._jkmin+self._djk/2.,
                                self._jkmax-self._djk/2.,
                                self._nbins)
            hs= numpy.linspace(self._hmin+self._dh/2.,
                               self._hmax-self._dh/2.,
                               self._nbins)
            self._interpolatedhist= scipy.interpolate.RectBivariateSpline(jks,hs,self._hist,
                                      bbox=[self._jkmin,self._jkmax,
                                            self._hmin,self._hmax],
                                      s=0.)
            #raise NotImplementedError("'interpolate=True' option for isomodel not implemented yet")
        return None

    def __call__(self,jk,h):
        """
        NAME:
           __call__
        PURPOSE:
           calls logpjkh
        INPUT:
           see inputs forlogpjkh
        OUTPUT:
           see output forlogpjkh
        HISTORY:
           2012-02-17 - Written - Bovy (IAS)
        """
        return self.logpjkh(jk,h)

    def logpjkh(self,jk,h):
        """
        NAME:
           logpjkh
        PURPOSE:
           return the probability of the (J-Ks,M_H) pair
        INPUT:
           jk - J-Ks
           h - M_H (absolute magnitude)
        OUTPUT:
           log of the probability
        HISTORY:
           2012-02-17 - Written - Bovy (IAS)
        """
        if self._interpolate:
            return numpy.log(self._interpolatedhist(jk,h))
        else:
            jkbin= int(math.floor((jk-self._jkmin)/self._djk))
            hbin= int(math.floor((h-self._hmin)/self._dh))
            if jkbin < 0 or jkbin >= self._nbins:
                return -numpy.finfo(numpy.dtype(numpy.float64)).max
            if hbin < 0 or hbin >= self._nbins:
                return -numpy.finfo(numpy.dtype(numpy.float64)).max
            return self._loghist[jkbin,hbin]
    
    def plot(self,log=False,conditional=False,nbins=None):
        """
        NAME:
           plot
        PURPOSE:
           plot the resulting (J-Ks,H) distribution
        INPUT:
           log= (default: False) if True, plot log
           conditional= (default: False) if True, plot conditional distribution
                        of H given J-Ks
           nbins= if set, set the number of bins
        OUTPUT:
           plot to output device
        HISTORY:
           2012-02-17 - Written - Bovy (IAS)
        """
        if not _BOVY_PLOT_LOADED:
            raise ImportError("'galpy.util.bovy_plot' plotting package not found")
        #Form histogram grid
        if nbins is None:
            nbins= self._nbins
        jks= numpy.linspace(self._jkmin+self._djk/2.,
                            self._jkmax-self._djk/2.,
                            nbins)
        hs= numpy.linspace(self._hmax-self._dh/2.,#we reverse
                           self._hmin+self._dh/2.,
                           nbins)
        plotthis= numpy.zeros((nbins,nbins))
        for ii in range(nbins):
            for jj in range(nbins):
                plotthis[ii,jj]= self(jks[ii],hs[jj])
        if not log:
            plotthis= numpy.exp(plotthis)
        if conditional: #normalize further
            for ii in range(nbins):
                plotthis[ii,:]/= numpy.nanmax(plotthis[ii,:])/numpy.nanmax(plotthis)
        return bovy_plot.bovy_dens2d(plotthis.T,origin='lower',cmap='gist_yarg',
                                     xrange=[self._jkmin,self._jkmax],
                                     yrange=[self._hmax,self._hmin],
                                     aspect=(self._jkmax-self._jkmin)/(self._hmax-self._hmin),
                                     xlabel=r'$(J-K_s)_0\ [\mathrm{mag}]$',
                                     ylabel=r'$M_H\ [\mathrm{mag}]$',
                                     interpolation='nearest')
    
