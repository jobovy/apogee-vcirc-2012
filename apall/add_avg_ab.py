###############################################################################
#  add_ak: add the Ak values
###############################################################################
import sys
import os, os.path
import numpy
import esutil
import fitsio
def add_avg_ab(infilename,outfilename):
    """
    NAME:
       add_avg_ab
    PURPOSE:
       add the average abundances
    INPUT:
       infilename - filename of catalog file
       outfilename - filename of the output file (catalog+ak)
    OUTPUT:
       (none)
    HISTORY:
       2012-11-15 - Written - Bovy (IAS)
    """
    #Read input file
    cat= fitsio.read(infilename)
    #Now add this to the catalog
    cat= esutil.numpy_util.add_fields(cat,[('FEH_AVG', float),
                                           ('LOGG_AVG', float),
                                           ('TEFF_AVG', float),
                                           ('VHELIO_AVG', float)])
    cat['FEH_AVG']= -9999.00
    cat['LOGG_AVG']= -9999.00
    cat['TEFF_AVG']= -9999.00
    cat['VHELIO_AVG']= -9999.00
    #Fill in
    primarycat= cat[(cat['SPECPRIMARY'] ==  1)]
    nprimary= len(primarycat)
    for ii in range(nprimary):
        #Find all visits
        thesedataIndx= (cat['UNIQID'] == primarycat['SPECID'][ii])
        feh_avg= calc_avg(cat['FEH'][thesedataIndx])
        logg_avg= calc_avg(cat['LOGG'][thesedataIndx])
        teff_avg= calc_avg(cat['TEFF'][thesedataIndx])
        thesedata= cat[thesedataIndx]
        helioIndx= (thesedata['VRADERR'] != 0.)
        vhelio_avg= numpy.median(thesedata[helioIndx]['VHELIO'])
        cat['FEH_AVG'][thesedataIndx]= feh_avg
        cat['LOGG_AVG'][thesedataIndx]= logg_avg
        cat['TEFF_AVG'][thesedataIndx]= teff_avg
        cat['VHELIO_AVG'][thesedataIndx]= vhelio_avg
    #Save
    fitsio.write(outfilename,cat,clobber=True)


def calc_avg(quant):
    indx= quant > -1000.
    return numpy.mean(quant[indx])

if __name__ == '__main__':
    add_avg_ab(sys.argv[1],sys.argv[2])
