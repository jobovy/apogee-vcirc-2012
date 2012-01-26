import numpy
import fitsio
import apogee
def readVclosData(lmin=35.,bmax=2.,postshutdown=True,fehcut=False):
    """
    NAME:
       readVclosData
    PURPOSE:
       single routine to read the data
    INPUT:
       lmin - minimal Galactic longitude
       bmax - maximal Galactic latitude
       postshutdown= if True, only use post-shutdown data (default: True)
       fehcut= if True, cut to rough FeH > -0.5 (default: False)
    OUTPUT:
    HISTORY:
       2012-01-25 - Written - Bovy (IAS)
    """
    datafile= apogee.tools.apallPath()
    data= fitsio.read(datafile,1)
    #Primary data only
    data= data[(data['SPECPRIMARY'] == 1)]
    #data cuts
    data=data[(numpy.fabs(data['GLAT']) < bmax)*(data['GLON'] > lmin)\
                  *(data['GLON'] < (360.-lmin))]
    data= data[((data['APOGEE_TARGET1'] & 2**9) == 0)] #no probable cluster members
    indx= numpy.array(['STAR' in data['OBJTYPE'][ii] for ii in range(len(data))],dtype='bool')
    data= data[indx]
    if postshutdown:
        data= data[(data['MJD5'] > 55788)]
    if fehcut:
        data= data[(data['FEH'] > -0.5)]
    return data
