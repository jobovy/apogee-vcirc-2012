import numpy
import fitsio
import apogee
def readVclosData(lmin=35.,bmax=2.,postshutdown=True,fehcut=False,cohort=None,
                  meanb=0.,meanb_tol=0.5):
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
       meanb= (default: 0.) require the mean b to be meanb within meanb_tol
              [useful to make sure one has plates in the plane]
       meanb_tol= (default:0.5) tolerance on meanb
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
    #Remove anything with vraderr= 0.
    data= data[(data['VRADERR'] != 0.)]
    if postshutdown:
        data= data[(data['MJD5'] > 55788)]
    if fehcut:
        data= data[(data['FEH'] > -0.5)]
    if not cohort is None:
        if cohort.lower() == 'short':
            data= data[((data['APOGEE_TARGET1'] & 2L**11) != 0)]
        elif cohort.lower() == 'medium' or cohort.lower() == 'intermediate':
            data= data[((data['APOGEE_TARGET1'] & 2L**12) != 0)]
        elif cohort.lower() == 'long':
            data= data[((data['APOGEE_TARGET1'] & 2L**13) != 0)]           
    #For every plate, calculate meanb, then cut on it
    plates= list(set(data['PLATE']))
    nplates= len(plates)
    for ii in range(nplates):
        #Calculate meanb
        mb= numpy.mean(data[(data['PLATE'] == plates[ii])]['GLAT'])
        if (mb-meanb)**2./meanb_tol**2. > 1.:
            data= data[(data['PLATE'] != plates[ii])]
    return data
