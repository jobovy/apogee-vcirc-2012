import sys
import os, os.path
import math
import numpy
import cPickle as pickle
from optparse import OptionParser
from fitvc import _PMSGRA, _VRSUN, _REFV0, _REFR0
_FLATFIT= '../fits/all_simpledrift_noro_dwarf_vpec_sratio_hs.sav'
_PLFIT= '../fits/all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav'
_FLATSAMPLES= '../fits/all_simpledrift_noro_dwarf_vpec_sratio_hs_10000samples.sav'
_PLSAMPLES= '../fits/all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_10000samples.sav'
def sixtyeigthinterval(bestfit,samples,quantile=.68):
    """68% interval around bestfit"""
    #Sort
    sortedsamples= sorted(samples)
    low= sortedsamples[int(numpy.floor((0.5-quantile/2.)*len(samples)))]
    high= sortedsamples[int(numpy.floor((0.5+quantile/2.)*len(samples)))]
    if low >= bestfit:
        low= bestfit-10.**(math.log10(high-bestfit)-1.) #HACK
    if high <= bestfit:
        high= bestfit+10.**(math.log10(bestfit-low)-1.) #HACK
    print low, bestfit, high
    if False:
       #Find bestfit in this
        bestfitindx= numpy.argmin(numpy.fabs(bestfit-samples))
        lowoffset= int(numpy.floor(quantile/2.*len(samples)))
        highoffset= int(numpy.floor(quantile/2.*len(samples)))
        print bestfitindx, lowoffset, highoffset
        if (bestfitindx-lowoffset) < 0:
            lowoffset= bestfitindx
            highoffset+= (highoffset-bestfitindx)
        elif (bestfitindx+highoffset) > (len(samples)-1):
            highoffset= len(samples)-1-bestfitindx
            lowoffset+= (lowoffset-bestfitindx)
        print lowoffset, highoffset
        low= sortedsamples[bestfitindx-lowoffset]
        high= sortedsamples[bestfitindx+highoffset]
    return (low,high)
def resultsTable(parser):
    (options,args)= parser.parse_args()
    cmdline= '%python resultsTable.py '+args[0]
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
    #Open savefile for power-law fit and samples
    savefile= open(_PLFIT,'rb')
    plparams= pickle.load(savefile)
    savefile.close()
    savefile= open(_PLSAMPLES,'rb')
    plsamples= pickle.load(savefile)
    savefile.close()
    #Set up sections
    names= ['$V_c(R_0)\ [\mathrm{km\ s}^{-1}]$',
            '$\\beta$',
            '$\\mathrm{d} V_c / \mathrm{d} R \left(R_0 \\right)\ [\mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}]$',
            '$A\ [\mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}]$',
            '$B\ [\mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}]$',
            '$(B^2-A^2)/(2\pi G)\ [\mathrm{M_{\odot}\ pc}^{-3}]$',
            '$\\Omega_0\ [\mathrm{km\ s}^{-1}\ \mathrm{kpc}^{-1}]$',
            '$R_0\ [\mathrm{kpc}]$',
            '$V_{R,\odot}\ [\mathrm{km\ s}^{-1}]$',
            '$V_{\phi,\odot}\ [\mathrm{km\ s}^{-1}]$',
            '$V_{\phi,\odot}-V_c\ [\mathrm{km\ s}^{-1}]$',
            '$\mu_{\mathrm{Sgr\ A}^{^*}}\ [\mathrm{mas\ yr}^{-1}]$',
            '$\\sigma_R(R_0)\ [\mathrm{km\ s}^{-1}]$',
            '$R_0/h_\sigma$',
            '$\\sigma_\phi^2 / \sigma_R^2$']
    missing= [0,1,1,0,0,1,0,0,0,0,0,0,0,0,0] #1 if missing for flat
    skip= [0,0,0,0,0,0,1,0,0,0,0,1,0,0,0] #1 if line skip after this parameter
    scale= [_REFV0,1.,1.,_REFV0/_REFR0,_REFV0/_REFR0,
            1./(10.**3.*4.302*2.*numpy.pi)*_REFV0/_REFR0*_REFV0/_REFR0,
            _REFV0/_REFR0, _REFR0,_VRSUN,_REFR0,_REFR0,
            1./4.74047,_REFV0,1.,1.]
    flatbestfits= []
    flatxs= []
    plbestfits= []
    plxs= []
    #V_c
    flatbestfits.append(flatparams[0])
    plbestfits.append(plparams[0])
    flatxs.append(numpy.array([s[0] for s in flatsamples]))
    plxs.append(numpy.array([s[0] for s in plsamples]))
    #beta
    flatbestfits.append(None)
    plbestfits.append(plparams[6])
    flatxs.append(None)
    plxs.append(numpy.array([s[6] for s in plsamples]))
    #dvcdr
    flatbestfits.append(None)
    plbestfits.append(plparams[6]*_REFV0/_REFR0)
    flatxs.append(None)
    plxs.append(numpy.array([s[6]*_REFV0/_REFR0 for s in plsamples]))
    #Oort A
    flatbestfits.append(0.5*flatparams[0]/flatparams[1])
    plbestfits.append(0.5*(plparams[0]/plparams[1]-plparams[6]))
    flatxs.append(numpy.array([0.5*s[0]/s[1] for s in flatsamples]))
    plxs.append(numpy.array([0.5*(s[0]/s[1]-s[6]) for s in plsamples]))
    #Oort B
    flatbestfits.append(-0.5*flatparams[0]/flatparams[1])
    plbestfits.append(-0.5*(plparams[0]/plparams[1]+plparams[6]))
    flatxs.append(numpy.array([-0.5*s[0]/s[1] for s in flatsamples]))
    plxs.append(numpy.array([-0.5*(s[0]/s[1]+s[6]) for s in plsamples]))
    #B^2-A^2
    flatbestfits.append(None)
    plbestfits.append((-0.5*(plparams[0]/plparams[1]+plparams[6]))**2.-
                      (0.5*(plparams[0]/plparams[1]-plparams[6]))**2.)
    flatxs.append(None)
    plxs.append(numpy.array([(-0.5*(s[0]/s[1]+s[6]))**2.-
                             (0.5*(s[0]/s[1]-s[6]))**2. for s in plsamples]))
    #Omega_0
    flatbestfits.append(flatparams[0]/flatparams[1])
    plbestfits.append(plparams[0]/plparams[1])
    flatxs.append(numpy.array([s[0]/s[1] for s in flatsamples]))
    plxs.append(numpy.array([s[0]/s[1] for s in plsamples]))
    #R_0
    flatbestfits.append(flatparams[1])
    plbestfits.append(plparams[1])
    flatxs.append(numpy.array([s[1] for s in flatsamples]))
    plxs.append(numpy.array([s[1] for s in plsamples]))
    #VRsun
    flatbestfits.append(flatparams[6])
    plbestfits.append(plparams[7])
    flatxs.append(numpy.array([s[6] for s in flatsamples]))
    plxs.append(numpy.array([s[7] for s in plsamples]))
    #Vphisun
    flatbestfits.append(flatparams[7]*flatparams[1]*_PMSGRA)
    plbestfits.append(plparams[8]*plparams[1]*_PMSGRA)
    flatxs.append(numpy.array([s[7]*s[1]*_PMSGRA for s in flatsamples]))
    plxs.append(numpy.array([s[8]*s[1]*_PMSGRA for s in plsamples]))
    #Vphisun-vc
    flatbestfits.append(flatparams[7]*flatparams[1]*_PMSGRA-flatparams[0]*_REFV0/_REFR0)
    plbestfits.append(plparams[8]*plparams[1]*_PMSGRA-plparams[0]*_REFV0/_REFR0)
    flatxs.append(numpy.array([s[7]*s[1]*_PMSGRA-s[0]*_REFV0/_REFR0 for s in flatsamples]))
    plxs.append(numpy.array([s[8]*s[1]*_PMSGRA-s[0]*_REFV0/_REFR0 for s in plsamples]))
    #pmsgra
    flatbestfits.append(flatparams[7]*_PMSGRA)
    plbestfits.append(plparams[8]*_PMSGRA)
    flatxs.append(numpy.array([s[7]*_PMSGRA for s in flatsamples]))
    plxs.append(numpy.array([s[8]*_PMSGRA for s in plsamples]))
    #\sigma_R(R_0)
    flatbestfits.append(numpy.exp(flatparams[2]))
    plbestfits.append(numpy.exp(plparams[2]))
    flatxs.append(numpy.exp(numpy.array([s[2] for s in flatsamples])))
    plxs.append(numpy.exp(numpy.array([s[2] for s in plsamples])))
    #R_0/h_sigma
    flatbestfits.append(flatparams[1]*flatparams[9])
    plbestfits.append(plparams[1]*plparams[10])
    flatxs.append(numpy.array([s[1]*s[9] for s in flatsamples]))
    plxs.append(numpy.array([s[1]*s[10] for s in plsamples]))
    #X^2
    flatbestfits.append(flatparams[8])
    plbestfits.append(plparams[9])
    flatxs.append(numpy.array([s[8] for s in flatsamples]))
    plxs.append(numpy.array([s[9] for s in plsamples]))
    #Make table
    quantile= 0.68
    outfile= open(args[0],'w')
    for ii in range(len(names)):
        #Set up line
        printline= names[ii]
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
            print names[ii], maxerr
            if math.log10(maxerr) >= 0.:
                value= '$%.0f$' % (bestfit)
                if numpy.fabs((dlow-dhigh)/dlow) < 0.1:
                    err= '$\pm$%.0f' % ((dlow+dhigh)/2.)           
                else:
                    err= '$^{+%.0f}_{-%.0f}$' % (dhigh,dlow)           
            elif math.log10(maxerr) >= -1.:
                value= '$%.1f$' % (bestfit)
                if numpy.fabs((dlow-dhigh)/dlow) < 0.1:
                    err= '$\pm$%.1f' % ((dlow+dhigh)/2.)
                else:
                    err= '$^{+%.1f}_{-%.1f}$' % (dhigh,dlow)
            elif math.log10(maxerr) >= -2.:
                value= '$%.2f$' % (bestfit)
                if numpy.fabs((dlow-dhigh)/dlow) < 0.1:
                    err= '$\pm$%.2f' % ((dlow+dhigh)/2.)           
                else:
                    err= '$^{+%.2f}_{-%.2f}$' % (dhigh,dlow)
            elif math.log10(maxerr) >= -3.:
                value= '$%.3f$' % (bestfit)
                if numpy.fabs((dlow-dhigh)/dlow) < 0.1:
                    err= '$\pm$%.3f' % ((dlow+dhigh)/2.)           
                else:
                    err= '$^{+%.3f}_{-%.3f}$' % (dhigh,dlow)
            else:
                value= '$%.4f$' % (bestfit)
                if numpy.fabs((dlow-dhigh)/dlow) < 0.1:
                    err= '$\pm$%.4f' % ((dlow+dhigh)/2.)           
                else:
                    err= '$^{+%.4f}_{-%.4f}$' % (dhigh,dlow)
            printline+= value+'&'+err
        #Power-law
        printline+= ' & '
        bestfit= plbestfits[ii]*scale[ii]
        xs= plxs[ii]*scale[ii]
        lowval, highval= sixtyeigthinterval(bestfit,xs,quantile=quantile)
        dlow, dhigh= bestfit-lowval, highval-bestfit
        #Prepare
        maxerr= numpy.amin(numpy.fabs([dlow,dhigh]))
        if math.log10(maxerr) >= 0.:
            value= '$%.0f$' % (bestfit)
            if numpy.fabs((dlow-dhigh)/dlow) < 0.1:
                err= '$\pm$%.0f' % ((dlow+dhigh)/2.)           
            else:
                err= '$^{+%.0f}_{-%.0f}$' % (dhigh,dlow) 
        elif math.log10(maxerr) >= -1.:
            value= '$%.1f$' % (bestfit)
            if numpy.fabs((dlow-dhigh)/dlow) < 0.1:
                err= '$\pm$%.1f' % ((dlow+dhigh)/2.)           
            else:
                err= '$^{+%.1f}_{-%.1f}$' % (dhigh,dlow)
        elif math.log10(maxerr) >= -2.:
            value= '$%.2f$' % (bestfit)
            if numpy.fabs((dlow-dhigh)/dlow) < 0.1:
                err= '$\pm$%.2f' % ((dlow+dhigh)/2.)           
            else:
                err= '$^{+%.2f}_{-%.2f}$' % (dhigh,dlow)
        elif math.log10(maxerr) >= -3.:
            value= '$%.3f$' % (bestfit)
            if numpy.fabs((dlow-dhigh)/dlow) < 0.1:
                err= '$\pm$%.3f' % ((dlow+dhigh)/2.)           
            else:
                err= '$^{+%.3f}_{-%.3f}$' % (dhigh,dlow)
        else:
            value= '$%.4f$' % (bestfit)
            if numpy.fabs((dlow-dhigh)/dlow) < 0.1:
                err= '$\pm$%.4f' % ((dlow+dhigh)/2.)           
            else:
                err= '$^{+%.4f}_{-%.4f}$' % (dhigh,dlow)
        #Print value+err
        printline+= value+'&'+err
        if not ii == (len(names)-1):
            printline+= '\\\\'
        if skip[ii]:
            printline+= '\\\\'
            #Write the line
        outfile.write(printline+'\n')
    outfile.write('\\enddata\n')
    outfile.write(cmdline+'\n')
    outfile.close()

def get_options():
    usage = "usage: %prog [options] <outputfilename>\n\noutputfilename= name of the file that the table will be saved to"
    parser = OptionParser(usage=usage)
    return parser

if __name__ == '__main__':
    resultsTable(get_options())
