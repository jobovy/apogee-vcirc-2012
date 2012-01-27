#Module with theoretical distributions
import math
from scipy import integrate
_DEGTORAD= math.pi/180.
def surfacemassLOSbd(d,b,l,bmean,diskdf,hz):
    """
    NAME:
       surfacemassLOSbd
    PURPOSE:
       calculate the number density on a plate as a function of d and b (integrated over l already)
    INPUT:
       d, b, l- distance in /ro units, and angles in degree
       bmean - mean b of the plate
       diskdf - galpy diskdf object
       hz - scale height in /ro units
    OUTPUT:
       density
    HISTORY:
       2012-01-27 - Written - Bovy (IAS)
    """
    return diskdf.surfacemassLOS(d,l,deg=True)*d*math.cos(b*_DEGTORAD)*math.exp(-d*math.fabs(math.sin(b*_DEGTORAD))/hz)*math.sqrt(1.5**2.-(b-bmean)**2.)

def surfacemassLOSb(d,l,bmin,bmax,diskdf,hz):
    """
    NAME:
       surfacemassLOSb
    PURPOSE:
       calculate the number density on a plate as a function of d (integrated over l and balready)
    INPUT:
       d, l- distance in /ro units, and angles in degree
       bmin, bmax - range in b in degree
       diskdf - galpy diskdf object
       hz - scale height in /ro units
    OUTPUT:
       density
    HISTORY:
       2012-01-27 - Written - Bovy (IAS)
    """
    return integrate.quad(_surfacemassLOSbIntegrand,bmin,bmax,
                          args=(d,l,(bmin+bmax)/2.,diskdf,hz))[0]

def _surfacemassLOSbIntegrand(b,d,l,bmean,diskdf,hz):
    return surfacemassLOSbd(d,b,l,bmean,diskdf,hz)

def surfacemassLOSd(b,l,dmin,dmax,bmean,diskdf,hz):
    """
    NAME:
       surfacemassLOSd
    PURPOSE:
       calculate the number density on a plate as a function of b (integrated over l and d already)
    INPUT:
       b, l- angles in degree
       dmin, dmax - range in d/ro
       bmean - center of the plate
       diskdf - galpy diskdf object
       hz - scale height in /ro units
    OUTPUT:
       density
    HISTORY:
       2012-01-27 - Written - Bovy (IAS)
    """
    return integrate.quad(_surfacemassLOSdIntegrand,dmin,dmax,
                          args=(b,l,bmean,diskdf,hz))[0]

def _surfacemassLOSdIntegrand(d,b,l,bmean,diskdf,hz):
    return surfacemassLOSbd(d,b,l,bmean,diskdf,hz)
