import sys
import os, os.path
import math
import numpy
from optparse import OptionParser
from readVclosData import readVclosData
from plot_bestfit import get_options
def fieldsTable(parser):
    (options,args)= parser.parse_args()
    #Read data
    print "Reading the data ..."
    data= readVclosData(postshutdown=options.postshutdown,
                        fehcut=options.fehcut,
                        cohort=options.cohort,
                        lmin=options.lmin,
                        bmax=options.bmax,
                        ak=True,
                        cutmultiples=options.cutmultiples,
                        jkmax=options.jkmax)
    #Parse data
    locs= numpy.array(sorted(list(set(data['LOCATION']))))
    nlocs= len(locs)
    platel= numpy.zeros(nlocs)
    for ii in range(nlocs):
        indx= (data['LOCATION'] == locs[ii])
        platel[ii]= numpy.mean(data['GLON'][indx])
    sortindx= numpy.argsort(platel)
    locs= locs[sortindx]
    outfile= open(options.plotfilename,'w')
    for ii in range(nlocs):
        indx= (data['LOCATION'] == locs[ii])
        #l
        printline= '$%i^\circ$ ' % int(round(numpy.mean(data['GLON'][indx])))
        # # of data
        printline+= '& %i' % numpy.sum(indx)
        # # of data H < 12.2
        printline+= '& %i' % numpy.sum((data['LOCATION'] == locs[ii])*(data['HMAG'] < 12.2))
        # # of data 12.2 <= H < 12.8
        nn= numpy.sum((data['LOCATION'] == locs[ii])*(data['HMAG'] >= 12.2)*(data['HMAG'] < 12.8))
        if nn > 0:
            printline+= '& %i' % nn
        else:
            printline+= '& 0'
        # # of data 12.8 <= H < 13.8
        nn= numpy.sum((data['LOCATION'] == locs[ii])*(data['HMAG'] >= 12.8)*(data['HMAG'] < 13.8))
        if nn > 0: printline+= '& %i' % nn
        else: printline+= '& 0 '
        #median ak
        printline+= '& %.1f ' % numpy.median(data['AK'][indx])
        #median visits
        #printline+= '& %i ' % numpy.median(data['NVISITS'][indx])
        #Write the line
        if not ii == (nlocs-1): printline+= '\\\\\n'
        else: printline+= '\n'
        outfile.write(printline)
    outfile.close()

if __name__ == '__main__':
    fieldsTable(get_options())
