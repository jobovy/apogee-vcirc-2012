import sys
import os, os.path
import math
import numpy
import cPickle as pickle
from optparse import OptionParser
from readVclosData import readVclosData
from fitvc import _PMSGRA, _VRSUN, _REFV0, _REFR0
import logl
from plot_pdfs import set_options
from resultsTable import sixtyeigthinterval
_INDIVCHI2= False
_FLATFIT= '../fits/all_simpledrift_noro_dwarf_vpec_sratio_hs.sav'
_FLATSAMPLES= '../fits/all_simpledrift_noro_dwarf_vpec_sratio_hs_10000samples.sav'
_NMOCKS= 5
_MOCKFIT= ['../fake_indivfeh_correct/all_simpledrift-dehnen_noro_dwarf_vpec_sratio_hs_hr3.sav',
           '../fake_indivfeh_correct/all2_simpledrift-dehnen_noro_dwarf_vpec_sratio_hs_hr3.sav',
           '../fake_indivfeh_correct/all3_simpledrift-dehnen_noro_dwarf_vpec_sratio_hs_hr3.sav',
           '../fake_indivfeh_correct/all4_simpledrift-dehnen_noro_dwarf_vpec_sratio_hs_hr3.sav',
           '../fake_indivfeh_correct/all5_simpledrift-dehnen_noro_dwarf_vpec_sratio_hs_hr3.sav']
_MOCKSAMPLES= ['../fake_indivfeh_correct/all_simpledrift-dehnen_noro_dwarf_vpec_sratio_hs_hr3_10000.sav',
               '../fake_indivfeh_correct/all2_simpledrift-dehnen_noro_dwarf_vpec_sratio_hs_hr3_10000.sav',
               '../fake_indivfeh_correct/all3_simpledrift-dehnen_noro_dwarf_vpec_sratio_hs_hr3_10000.sav',
               '../fake_indivfeh_correct/all4_simpledrift-dehnen_noro_dwarf_vpec_sratio_hs_hr3_10000.sav',
               '../fake_indivfeh_correct/all5_simpledrift-dehnen_noro_dwarf_vpec_sratio_hs_hr3_10000.sav']
_MOCKDATA= ['../fake_indivfeh_correct/fake_simpledrift-dehnen_noro_dwarf_vpec_sratio_hs.fits',
            '../fake_indivfeh_correct/fake2_simpledrift-dehnen_noro_dwarf_vpec_sratio_hs.fits',
            '../fake_indivfeh_correct/fake3_simpledrift-dehnen_noro_dwarf_vpec_sratio_hs.fits',
            '../fake_indivfeh_correct/fake4_simpledrift-dehnen_noro_dwarf_vpec_sratio_hs.fits',
            '../fake_indivfeh_correct/fake5_simpledrift-dehnen_noro_dwarf_vpec_sratio_hs.fits']            
def mockTable(parser):
    (options,args)= parser.parse_args()
    cmdline= '%python mockTable.py '+args[0]
    if len(args) == 0:
        parser.print_help()
        return
    #Open savefile for flat fit and samples
    savefile= open(_FLATFIT,'rb')
    flatparams= pickle.load(savefile)
    savefile.close()
    savefile= open(_FLATSAMPLES,'rb')
    flatsamples= pickle.load(savefile)
    savefile.close()
    #Open savefile for all mocks' fit and samples
    mockparams= []
    mocksamples= []
    for ii in range(_NMOCKS):
        savefile= open(_MOCKFIT[ii],'rb')
        mockparams.append(pickle.load(savefile))
        savefile.close()
        savefile= open(_MOCKSAMPLES[ii],'rb')
        mocksamples.append(pickle.load(savefile))
        savefile.close()
    #Set up sections
    names= ['$V_c(R_0)\ [\mathrm{km\ s}^{-1}]$',
            '$R_0\ [\mathrm{kpc}]$',
            '$V_{R,\odot}\ [\mathrm{km\ s}^{-1}]$',
            '$\Omega_{\odot}\ [\mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}]$',
            '$\\sigma_R(R_0)\ [\mathrm{km\ s}^{-1}]$',
            '$R_0/h_\sigma$',
            '$X^2 \equiv \sigma_\phi^2 / \sigma_R^2$',
            '$\Delta\chi^2/\mathrm{dof}$']
    scale= [_REFV0,
            _REFR0,_VRSUN,1.,_REFV0,1.,1.,1.]
    flatbestfits= []
    flatxs= []
    mockbestfits= []
    mockxs= []
    #V_c
    flatbestfits.append(flatparams[0])
    flatxs.append(numpy.array([s[0] for s in flatsamples]))
    #R_0
    flatbestfits.append(flatparams[1])
    flatxs.append(numpy.array([s[1] for s in flatsamples]))
    #VRsun
    flatbestfits.append(flatparams[6])
    flatxs.append(numpy.array([s[6] for s in flatsamples]))
    #Omega_sun
    flatbestfits.append(flatparams[7]*_PMSGRA)
    flatxs.append(numpy.array([s[7]*_PMSGRA for s in flatsamples]))
    #\sigma_R(R_0)
    flatbestfits.append(numpy.exp(flatparams[2]))
    flatxs.append(numpy.exp(numpy.array([s[2] for s in flatsamples])))
    #R_0/h_sigma
    flatbestfits.append(flatparams[1]*flatparams[9])
    flatxs.append(numpy.array([s[1]*s[9] for s in flatsamples]))
    #X^2
    flatbestfits.append(flatparams[8])
    flatxs.append(numpy.array([s[8] for s in flatsamples]))
    #chi2
    options= set_options(None)
    data= readVclosData(validfeh=True)
    options.isofile='../fits/isos.dat'
    thislogl= logl.logl(init=flatparams,options=options,data=data)
    flatbestfits.append(-2.*numpy.sum(thislogl)/(len(data)-len(flatparams)))
    if _INDIVCHI2:
        locs= numpy.array(sorted(list(set(data['LOCATION']))))
        nlocs= len(locs)
        platel= numpy.zeros(nlocs)
        for ii in range(nlocs):
            indx= (data['LOCATION'] == locs[ii])
            platel[ii]= numpy.mean(data['GLON'][indx])
            sortindx= numpy.argsort(platel)
            locs= locs[sortindx]
            platel= platel[sortindx]

            for jj in range(len(locs)):
                names.append('$\Delta\chi^2/\mathrm{dof}\ \mathrm{w/o}\ l=%i^\circ$' % int(round(platel[jj])))
                scale.append(1.)
                indx= (data['LOCATION'] != locs[jj])
                flatbestfits.append(-2.*numpy.sum(thislogl[indx])/(len(data[indx])-len(flatparams)))
    else: locs= []
    ##Same for all mocks
    for ii in range(_NMOCKS):
        thisbestfits= []
        thisxs= []
        #V_c
        thisbestfits.append(mockparams[ii][0])
        thisxs.append(numpy.array([s[0] for s in mocksamples[ii]]))
         #R_0
        thisbestfits.append(mockparams[ii][1])
        thisxs.append(numpy.array([s[1] for s in mocksamples[ii]]))
        #VRsun
        thisbestfits.append(mockparams[ii][6])
        thisxs.append(numpy.array([s[6] for s in mocksamples[ii]]))
        #Omega_sun
        thisbestfits.append(mockparams[ii][7]*_PMSGRA)
        thisxs.append(numpy.array([s[7]*_PMSGRA for s in mocksamples[ii]]))
        #\sigma_R(R_0)
        thisbestfits.append(numpy.exp(mockparams[ii][2]))
        thisxs.append(numpy.exp(numpy.array([s[2] for s in mocksamples[ii]])))
        #R_0/h_sigma
        thisbestfits.append(mockparams[ii][1]*mockparams[ii][9])
        thisxs.append(numpy.array([s[1]*s[9] for s in mocksamples[ii]]))
        #X^2
        thisbestfits.append(mockparams[ii][8])
        thisxs.append(numpy.array([s[8] for s in mocksamples[ii]]))
        #chi2
        #Load mock data
        data= readVclosData(datafilename=_MOCKDATA[ii])       
        thislogl= logl.logl(init=mockparams[ii],
                            data=data,
                            options=options)
        thisbestfits.append(-2.*numpy.sum(thislogl)/(len(data)-len(mockparams[ii])))
        if _INDIVCHI2:
            for jj in range(len(locs)):
                indx= (data['LOCATION'] != locs[jj])
                thisbestfits.append(-2.*numpy.sum(thislogl[indx])/(len(data[indx])-len(mockparams)))
        #Append all
        mockbestfits.append(thisbestfits)
        mockxs.append(thisxs)
    #Make table
    quantile= 0.68
    outfile= open(args[0],'w')
    for ii in range(len(names)):
        #Set up line
        printline= names[ii]
        #Flat
        printline+= ' & '
        if ii >= len(names)-(len(locs)+1): #chi2
            printline+= ' \ldots & '
        else:
            bestfit= flatbestfits[ii]*scale[ii]
            xs= flatxs[ii]*scale[ii]
            lowval, highval= sixtyeigthinterval(bestfit,xs,quantile=quantile)
            dlow, dhigh= bestfit-lowval, highval-bestfit
            #Prepare
            maxerr= numpy.amin(numpy.fabs([dlow,dhigh]))
            print names[ii], maxerr
            if math.log10(maxerr) >= 0.:
                value= '$%.0f$' % (bestfit)
                if numpy.fabs((dlow-dhigh)) < 1.:
                    err= '$\pm$%.0f' % ((dlow+dhigh)/2.)           
                else:
                    err= '$^{+%.0f}_{-%.0f}$' % (dhigh,dlow)           
            elif math.log10(maxerr) >= -1.:
                value= '$%.1f$' % (bestfit)
                if numpy.fabs((dlow-dhigh)) < 0.1:
                    err= '$\pm$%.1f' % ((dlow+dhigh)/2.)
                else:
                    err= '$^{+%.1f}_{-%.1f}$' % (dhigh,dlow)
            elif math.log10(maxerr) >= -2.:
                value= '$%.2f$' % (bestfit)
                if numpy.fabs((dlow-dhigh)) < 0.01:
                    err= '$\pm$%.2f' % ((dlow+dhigh)/2.)           
                else:
                    err= '$^{+%.2f}_{-%.2f}$' % (dhigh,dlow)
            elif math.log10(maxerr) >= -3.:
                value= '$%.3f$' % (bestfit)
                if numpy.fabs((dlow-dhigh)) < 0.001:
                    err= '$\pm$%.3f' % ((dlow+dhigh)/2.)           
                else:
                    err= '$^{+%.3f}_{-%.3f}$' % (dhigh,dlow)
            else:
                value= '$%.4f$' % (bestfit)
                if numpy.fabs((dlow-dhigh)) < 0.0001:
                    err= '$\pm$%.4f' % ((dlow+dhigh)/2.)           
                else:
                    err= '$^{+%.4f}_{-%.4f}$' % (dhigh,dlow)
            printline+= value+'&'+err
        #Mocks
        for jj in range(_NMOCKS):
            printline+= ' & '
            if ii >= len(names)-(len(locs)+1): #chi2
                bestfit= mockbestfits[jj][ii]*scale[ii]-flatbestfits[ii]*scale[ii]
                printline+= '%.2f & ' % bestfit
                if jj == _NMOCKS-1:
                    if not ii == (len(names)-1):
                        printline+= '\\\\'
                    outfile.write(printline+'\n')
                continue
            bestfit= mockbestfits[jj][ii]*scale[ii]
            xs= mockxs[jj][ii]*scale[ii]
            lowval, highval= sixtyeigthinterval(bestfit,xs,quantile=quantile)
            dlow, dhigh= bestfit-lowval, highval-bestfit
            #Prepare
            maxerr= numpy.amin(numpy.fabs([dlow,dhigh]))
            if math.log10(maxerr) >= 0.:
                value= '$%.0f$' % (bestfit)
                if numpy.fabs((dlow-dhigh)) < 1.:
                    err= '$\pm$%.0f' % ((dlow+dhigh)/2.)           
                else:
                    err= '$^{+%.0f}_{-%.0f}$' % (dhigh,dlow) 
            elif math.log10(maxerr) >= -1.:
                value= '$%.1f$' % (bestfit)
                if numpy.fabs((dlow-dhigh)) < 0.1:
                    err= '$\pm$%.1f' % ((dlow+dhigh)/2.)           
                else:
                    err= '$^{+%.1f}_{-%.1f}$' % (dhigh,dlow)
            elif math.log10(maxerr) >= -2.:
                value= '$%.2f$' % (bestfit)
                if numpy.fabs((dlow-dhigh)) < 0.01:
                    err= '$\pm$%.2f' % ((dlow+dhigh)/2.)           
                else:
                    err= '$^{+%.2f}_{-%.2f}$' % (dhigh,dlow)
            elif math.log10(maxerr) >= -3.:
                value= '$%.3f$' % (bestfit)
                if numpy.fabs((dlow-dhigh)) < 0.001:
                    err= '$\pm$%.3f' % ((dlow+dhigh)/2.)           
                else:
                    err= '$^{+%.3f}_{-%.3f}$' % (dhigh,dlow)
            else:
                value= '$%.4f$' % (bestfit)
                if numpy.fabs((dlow-dhigh)) < 0.0001:
                    err= '$\pm$%.4f' % ((dlow+dhigh)/2.)           
                else:
                    err= '$^{+%.4f}_{-%.4f}$' % (dhigh,dlow)
            printline+= value+'&'+err
        if not ii == (len(names)-1):
            printline+= '\\\\'
        if ii == (len(names)-(len(locs)+2)):
            printline+= '\\\\'
        if ii < (len(names)-(len(locs)+1)):      
            outfile.write(printline+'\n')
    outfile.write('\\enddata\n')
    outfile.write(cmdline+'\n')
    outfile.close()

def get_options():
    usage = "usage: %prog [options] <outputfilename>\n\noutputfilename= name of the file that the table will be saved to"
    parser = OptionParser(usage=usage)
    return parser

if __name__ == '__main__':
    mockTable(get_options())
