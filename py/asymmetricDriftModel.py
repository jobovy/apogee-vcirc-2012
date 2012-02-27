###############################################################################
# asymmetricDriftModel.py: a model for the asymmetric drift
#
# DO NOT IMPORT *
###############################################################################
import cPickle as pickle
import numpy
from scipy import interpolate
#We load this all gloabally, BAD
kinematicsFile= '../data/axi_va-kinematics_0.0.sav'
kinematicsFile= open(kinematicsFile,'rb')
vas= pickle.load(kinematicsFile) #vas has va/sigmaR^2(R_0)
xs= pickle.load(kinematicsFile)
kinematicsFile.close()
vas= vas[-6] # dispersion = 0.2 vo
nrs= 101
rs= numpy.linspace(0.1,2.2,nrs)
vaInterp= interpolate.InterpolatedUnivariateSpline(rs,vas)
#Fiducial
hR_fid= 1./3.
hs_fid= 1.
def dlnnusR2dlnR(hr,hs):
    return -1./hr-2./hs
dlnnusR2dlnR_fid= dlnnusR2dlnR(hR_fid,hs_fid)
def va(R,sigmaR,vc=1.,hR=1./3.,hs=1.):
    """
    NAME:
       va
    PURPOSE:
       return the asymmetric drift, in normalized coordinates (/vo)
    INPUT:
       R - radius (/Ro)
       sigmaR - radial velocity dispersion at Ro (/vo)
       vc= circular velocity at R (/vo)
       hR= (default: Ro/3.) scale length of the disk (/Ro)
       hs (default: Ro) dispersion scale length (/Ro)
    OUTPUT:
       asymmetric drift vc-<v> (/vo)
    HISTORY:
       2012-02-23 - Written - Bovy (IAS)
    """
    dlnnusR2dlnR_correct= dlnnusR2dlnR_fid-dlnnusR2dlnR(hR,hs)
    va= vaInterp(R)
    return sigmaR**2./vc*numpy.exp(-2.*(R-1.)/hs+2.*(R-1)/hs_fid)\
        *(va+0.5*R*dlnnusR2dlnR_correct*numpy.exp(-2.*(R-1.)))
