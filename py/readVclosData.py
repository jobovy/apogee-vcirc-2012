import numpy
import fitsio
import apogee
def readVclosData(lmin=35.,bmax=2.,postshutdown=True,fehcut=False,cohort=None,
                  meanb=0.,meanb_tol=0.5,jkmax=1.2,ak=True,
                  specprimary=True,
                  cutmultiples=False):
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
       cohort= (default: None) if set to 'short', 'medium', or 'long', cut to those cohorts
       meanb= (default: 0.) require the mean b to be meanb within meanb_tol
              [useful to make sure one has plates in the plane]
       meanb_tol= (default:0.5) tolerance on meanb
       ak= (default: True) only use objects for which dereddened mags exist
       jkmax - maximum (J-K)_0 (only in conjunction with ak=True)
       specprimary= (default= True) select only primary objects
       cutmultiples= (default: False) cut objects suspected to be in multiples (repeated Vlos std dev > 1 km/s) ONLY WORKS WITH SPECPRIMARY=True
    OUTPUT:
    HISTORY:
       2012-01-25 - Written - Bovy (IAS)
    """
    datafile= apogee.tools.apallPath()
    data= fitsio.read(datafile,1)
    if specprimary and not cutmultiples: #If cutmultiples we will cut to primary later
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
    if ak:
        data= data[(data['AK'] != -9999.9999)]
        if not jkmax is None:
            data= data[(data['J0MAG']-data['K0MAG'] < jkmax)]
    #For every plate, calculate meanb, then cut on it
    plates= list(set(data['PLATE']))
    nplates= len(plates)
    for ii in range(nplates):
        #Calculate meanb
        mb= numpy.mean(data[(data['PLATE'] == plates[ii])]['GLAT'])
        if (mb-meanb)**2./meanb_tol**2. > 1.:
            data= data[(data['PLATE'] != plates[ii])]
    if cutmultiples:
        #Build repeats
        primarydata= data[(data['SPECPRIMARY'] ==  1)]
        ndata= len(primarydata)
        keepindx= numpy.zeros(ndata,dtype='bool')
        keepindx[:]= False
        for ii in range(ndata):
            thesedata= data[(data['UNIQID'] == primarydata['SPECID'][ii])]
            indx= (thesedata['VRADERR'] != 0.)
            if numpy.sum(indx) < 2:
                continue
            thesedata= thesedata[indx]
            if numpy.std(thesedata['VRAD']) <= 1.: keepindx[ii]= True
        data= primarydata[keepindx]
    return data
