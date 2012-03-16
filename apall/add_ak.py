###############################################################################
#  add_ak: add the Ak values
###############################################################################
import sys
import os, os.path
import numpy
import esutil
import fitsio
_APOGEE_DATA= os.path.join(os.getenv('DATADIR'),'bovy','apogee')
def add_ak(infilename,outfilename):
    """
    NAME:
       add_ak
    PURPOSE:
       add the Ak values
    INPUT:
       infilename - filename of catalog file
       outfilename - filename of the output file (catalog+ak)
    OUTPUT:
       (none)
    HISTORY:
       2012-02-14 - Written - Bovy (IAS)
    """
    #Read input file
    cat= fitsio.read(infilename)
    #Read aK
    ak= _load_ak()
    #Now add this to the catalog
    cat= esutil.numpy_util.add_fields(cat,[('AK', float),
                                           ('AK_METHOD','|S18')])
    for ii in range(len(cat)):
        try:
            cat['AK'][ii]= float(ak[cat['OBJID'][ii].lower()][1])
            cat['AK_METHOD'][ii]= ak[cat['OBJID'][ii].lower()][2]
        except KeyError:
            cat['AK'][ii]= -9999.9999
            cat['AK_METHOD'][ii]= 'none'
    #Calculate AH and AJ, and apply to get H0MAG, etc.
    #Assume AH/AK= 1.55, AJ/AK=2.5 (A_l \propto l^-1.66)
    aj= cat['AK']*2.5
    ah= cat['AK']*1.55
    cat= esutil.numpy_util.add_fields(cat,[('J0MAG', float),
                                           ('H0MAG', float),
                                           ('K0MAG', float)])
    cat['J0MAG']= cat['JMAG']-aj
    cat['H0MAG']= cat['HMAG']-ah
    cat['K0MAG']= cat['KMAG']-cat['AK']
    cat['J0MAG'][(cat['AK'] == -9999.9999)]= -9999.9999
    cat['H0MAG'][(cat['AK'] == -9999.9999)]= -9999.9999
    cat['K0MAG'][(cat['AK'] == -9999.9999)]= -9999.9999
    #Save
    fitsio.write(outfilename,cat,clobber=True)
    
def _load_ak():
    akfile= os.path.join(_APOGEE_DATA,'ak.dat')
    input= numpy.genfromtxt(akfile,dtype='string',
                            delimiter='|',
                            skip_header=2,
                            skip_footer=2,
                            autostrip=True)
    #Create hash
    out= {}
    for ii in range(input.shape[0]):
        out[input[ii,0].lower()]= input[ii,:]
    return out
    

if __name__ == '__main__':
    add_ak(sys.argv[1],sys.argv[2])
