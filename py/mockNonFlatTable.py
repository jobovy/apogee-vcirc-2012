import sys
import os, os.path
import math
import numpy
import cPickle as pickle
from optparse import OptionParser
from fitvc import _PMSGRA, _VRSUN, _REFV0, _REFR0
from resultsTable import sixtyeigthinterval
_FLATFIT= '../fits/all_simpledrift_noro_dwarf_vpec_sratio_hs.sav'
_FLATSAMPLES= '../fits/all_simpledrift_noro_dwarf_vpec_sratio_hs_10000samples.sav'
_NSECTIONS= 3
_SECTIONS= ['powerlaw','linear','cubic']
_NMOCKS= 3
_MOCKFIT= ['../fake_indivfeh_correct/all_simpledrift-dehnen_noro_dwarf_FLAT_vpec_sratio_hs_hr3.sav',
           '../fake_indivfeh_correct/all2_simpledrift-dehnen_noro_dwarf_FLAT_vpec_sratio_hs_hr3.sav',
           '../fake_indivfeh_correct/all3_simpledrift-dehnen_noro_dwarf_FLAT_vpec_sratio_hs_hr3.sav']
_MOCKSAMPLES= ['../fake_indivfeh_correct/all_simpledrift-dehnen_noro_dwarf_FLAT_vpec_sratio_hs_hr3_10000.sav',
               '../fake_indivfeh_correct/all2_simpledrift-dehnen_noro_dwarf_FLAT_vpec_sratio_hs_hr3_10000.sav',
               '../fake_indivfeh_correct/all3_simpledrift-dehnen_noro_dwarf_FLAT_vpec_sratio_hs_hr3_10000.sav']
def mockNonFlatTable(parser):
    (options,args)= parser.parse_args()
    cmdline= '%python mockNonFlatTable.py '+args[0]
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
    #Set up sections
    names= ['$V_c(R_0)\ [\mathrm{km\ s}^{-1}]$',
            '$\\beta$',
            '$R_0\ [\mathrm{kpc}]$',
            '$V_{R,\odot}\ [\mathrm{km\ s}^{-1}]$',
            '$\Omega_{\odot}\ [\mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}]$',
            '$\\sigma_R(R_0)\ [\mathrm{km\ s}^{-1}]$',
            '$R_0/h_\sigma$',
            '$\\sigma_\phi^2 / \sigma_R^2$']
    missing= [0,1,0,0,0,0,0,0] #1 if missing for flat
    scale= [_REFV0,1.,
            _REFR0,_VRSUN,1.,_REFV0,1.,1.]
    flatbestfits= []
    flatxs= []
    #V_c
    flatbestfits.append(flatparams[0])
    flatxs.append(numpy.array([s[0] for s in flatsamples]))
    #beta
    flatbestfits.append(None)
    flatxs.append(None)
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
    #Make table
    quantile= 0.68
    outfile= open(args[0],'w')
    for kk in range(_NSECTIONS):
        #Open savefile for all mocks' fit and samples
        mockparams= []
        mocksamples= []
        for ii in range(_NMOCKS):
            savefile= open(_MOCKFIT[ii].replace('FLAT',_SECTIONS[kk]),'rb')
            mockparams.append(pickle.load(savefile))
            savefile.close()
            savefile= open(_MOCKSAMPLES[ii].replace('FLAT',_SECTIONS[kk]),'rb')
            mocksamples.append(pickle.load(savefile))
            savefile.close()
        mockbestfits= []
        mockxs= []
        ##Same for all mocks
        if _SECTIONS[kk] == 'powerlaw' or _SECTIONS[kk] == 'linear':
            skip= 1
        elif _SECTIONS[kk] == 'quadratic':
            skip= 2
        elif _SECTIONS[kk] == 'cubic':
            skip= 3
        if _SECTIONS[kk] == 'linear':
            names[1]= '$\\mathrm{d} V_c / \mathrm{d} R \left(R_0 \\right)\ [\mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}]$'
        elif _SECTIONS[kk] == 'quadratic':
            names.insert(2,'$\\mathrm{d}^2 V_c / \mathrm{d} R^2 \left(R_0 \\right)\ [\mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-2}]$')
            flatbestfits.insert(2,None)
            flatxs.insert(2,None)
            missing.insert(2,1)
            scale.insert(2,1.)
        elif _SECTIONS[kk] == 'cubic':
            if not 'quadratic' in _SECTIONS:
                #First insert quadratic
                names.insert(2,'$\\mathrm{d}^2 V_c / \mathrm{d} R^2 \left(R_0 \\right)\ [\mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-2}]$')
                flatbestfits.insert(2,None)
                flatxs.insert(2,None)
                missing.insert(2,1)
                scale.insert(2,1.)
            #Assumes that quadratic was already inserted
            names.insert(3,'$\\mathrm{d}^3 V_c / \mathrm{d} R^3 \left(R_0 \\right)\ [\mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-3}]$')
            flatbestfits.insert(3,None)
            flatxs.insert(3,None)
            missing.insert(3,1)
            scale.insert(3,1.)
        for ii in range(_NMOCKS):
            thisbestfits= []
            thisxs= []
            #V_c
            thisbestfits.append(mockparams[ii][0])
            thisxs.append(numpy.array([s[0] for s in mocksamples[ii]]))
            #Beta
            if _SECTIONS[kk] == 'powerlaw':
                thisbestfits.append(mockparams[ii][6])
                thisxs.append(numpy.array([s[6] for s in mocksamples[ii]]))
            elif _SECTIONS[kk] == 'linear' or _SECTIONS[kk] == 'quadratic' or _SECTIONS[kk] == 'cubic':
                thisbestfits.append(mockparams[ii][6]*_REFV0/_REFR0)
                thisxs.append(numpy.array([s[6]*_REFV0/_REFR0 for s in mocksamples[ii]]))
            if _SECTIONS[kk] == 'quadratic' or _SECTIONS[kk] == 'cubic':
                thisbestfits.append(mockparams[ii][7]*_REFV0/_REFR0**2.)
                thisxs.append(numpy.array([s[7]*_REFV0/_REFR0**2. for s in mocksamples[ii]]))
            if _SECTIONS[kk] == 'cubic':
                thisbestfits.append(mockparams[ii][8]*_REFV0/_REFR0**3.)
                thisxs.append(numpy.array([s[8]*_REFV0/_REFR0**3. for s in mocksamples[ii]]))
            #R_0
            thisbestfits.append(mockparams[ii][1])
            thisxs.append(numpy.array([s[1] for s in mocksamples[ii]]))
            #VRsun
            thisbestfits.append(mockparams[ii][6+skip])
            thisxs.append(numpy.array([s[6+skip] for s in mocksamples[ii]]))
            #Omega_sun
            thisbestfits.append(mockparams[ii][7+skip]*_PMSGRA)
            thisxs.append(numpy.array([s[7+skip]*_PMSGRA for s in mocksamples[ii]]))
            #\sigma_R(R_0)
            thisbestfits.append(numpy.exp(mockparams[ii][2]))
            thisxs.append(numpy.exp(numpy.array([s[2] for s in mocksamples[ii]])))
            #R_0/h_sigma
            thisbestfits.append(mockparams[ii][1]*mockparams[ii][9+skip])
            thisxs.append(numpy.array([s[1]*s[9+skip] for s in mocksamples[ii]]))
            #X^2
            thisbestfits.append(mockparams[ii][8+skip])
            thisxs.append(numpy.array([s[8+skip] for s in mocksamples[ii]]))
            #Append all
            mockbestfits.append(thisbestfits)
            mockxs.append(thisxs)
        for ii in range(len(names)):
            #Set up line
            if ii == 0:
                if _SECTIONS[kk] == 'powerlaw':
                    printline= 'power-law'
                else:
                    printline= _SECTIONS[kk]
            else:
                printline= " "
            printline+= ' & '
            printline+= names[ii]
            #Flat
            printline+= ' & '
            if missing[ii]:
                printline+= '\ldots & '
            else:
                bestfit= flatbestfits[ii]*scale[ii]
                xs= flatxs[ii]*scale[ii]
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
            #Mocks
            for jj in range(_NMOCKS):
                printline+= ' & '
                print _SECTIONS[kk], jj, ii, mockbestfits[jj][ii]
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
            if not (kk == (_NSECTIONS-1) and ii == (len(names)-1)):
                printline+= '\\\\'
            if not kk == (_NSECTIONS-1) and ii == (len(names)-1):
                printline+= '\\\\'
            outfile.write(printline+'\n')
    outfile.write('\\enddata\n')
    outfile.write(cmdline+'\n')
    outfile.close()

def get_options():
    usage = "usage: %prog [options] <outputfilename>\n\noutputfilename= name of the file that the table will be saved to"
    parser = OptionParser(usage=usage)
    return parser

if __name__ == '__main__':
    mockNonFlatTable(get_options())
