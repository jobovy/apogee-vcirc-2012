#Calculate some theoretical stuff
import sys
import os, os.path
import cPickle as pickle
import math
import numpy
from galpy.df import dehnendf, shudf
from galpy.util import bovy_coords, save_pickles
_ERASESTR= "                                                                                "
#ignore annoying warnings
import warnings
warnings.filterwarnings("ignore")
#Overall setup
nls= 201
ls= numpy.linspace(0.,360.,nls)
nds= 51
_SAVEDIR= os.path.join(os.getenv('DATADIR'),
                       'bovy',
                       'vclos')
#Quick def
def safe_dl_to_rphi(d,l):
    R= math.sqrt(1.+d**2.-2.*d*math.cos(l))
    if R == 0.:
        R= 0.0001
        d+= 0.0001
    if 1./math.cos(l) < d and math.cos(l) > 0.:
        theta= math.pi-math.asin(d/R*math.sin(l))
    else:
        theta= math.asin(d/R*math.sin(l))
    return (R,theta,d,l)
def calc_pred(pred_file,dfc,nls,nds,ls,ds):
    if os.path.exists(pred_file):
        savefile= open(pred_file,'rb')
        ls= pickle.load(savefile)
        avg_pred= pickle.load(savefile)
        ii= pickle.load(savefile)
        savefile.close()
    else:
        avg_pred= numpy.zeros(nls)
        ii= 0
    while ii < len(ls):
        sys.stdout.write('\r'+"Working on %i / %i ...\r" %(ii+1,len(ls)))
        sys.stdout.flush()
        meanvlos= 0.
        norm= 0.
        for jj in range(nds):
            d,l= ds[jj], ls[ii]/180.*numpy.pi
            R,theta,d,l= safe_dl_to_rphi(d,l)
            surfmass= dfc.surfacemassLOS(d,l,deg=False)
            meanvlos+= surfmass*dfc.meanvT(R)*math.sin(theta+l)
            norm+= surfmass
        avg_pred[ii]= meanvlos/norm
        ii+= 1
        save_pickles(pred_file,ls,avg_pred,ii)
    sys.stdout.write('\r'+_ERASESTR+'\r')
    sys.stdout.flush()
    save_pickles(pred_file,ls,avg_pred,ii)

ds= numpy.linspace(0.,10./8.,nds)
#Start calculating
#Fiducial
pred_file= os.path.join(_SAVEDIR,'l_vhelio_fid.sav')
print "Working on fiducial ..."
dfc= dehnendf(beta=0.,correct=True,niter=20)
calc_pred(pred_file,dfc,nls,nds,ls,ds)
#betas
betas= [-0.2,-0.1,0.1,0.2]
for beta in betas:
    pred_file= os.path.join(_SAVEDIR,'l_vhelio_beta_%.2f.sav' % beta)
    print "Working on beta %.2f ..." % beta
    dfc= dehnendf(beta=beta,correct=True,niter=20)
    calc_pred(pred_file,dfc,nls,nds,ls,ds)
#hr
hrs= [1./4.,1./2.]
for hr in hrs:
    pred_file= os.path.join(_SAVEDIR,'l_vhelio_hr_%.6f.sav' % hr)
    print "Working on hr %.2f ..." % hr
    dfc= dehnendf(beta=0.,correct=True,niter=20,
                  profileParams=(hr,1.,0.2))
    calc_pred(pred_file,dfc,nls,nds,ls,ds)
#Alt dmax
altdmaxs= [5./8.,15./8.]
for altdmax in altdmaxs:
    thisds= numpy.linspace(0.,altdmax,nds)
    pred_file= os.path.join(_SAVEDIR,'l_vhelio_dmax_%.6f.sav' % altdmax)
    print "Working on dmax %.2f ..." % altdmax
    dfc= dehnendf(beta=0.,correct=True,niter=20)
    calc_pred(pred_file,dfc,nls,nds,ls,thisds)
#shu
pred_file= os.path.join(_SAVEDIR,'l_vhelio_shu.sav')
print "Working on Shu ..."
dfc= shudf(beta=0.,correct=True,niter=20)
calc_pred(pred_file,dfc,nls,nds,ls,ds)
#sig0.1
sig= 0.1
pred_file= os.path.join(_SAVEDIR,'l_vhelio_sig_%.6f.sav' % sig)
if not os.path.exists(pred_file):
    print "Working on sig %.2f ..." % sig
    dfc= dehnendf(beta=0.,correct=True,niter=20,
                  profileParams=(1./3.,1.,sig))
    calc_pred(pred_file,dfc,nls,nds,ls,ds)
#hs
hs= 2./3.
pred_file= os.path.join(_SAVEDIR,'l_vhelio_hs_%.6f.sav' % hs)
print "Working on hs %.2f ..." % hs
dfc= dehnendf(beta=0.,correct=True,niter=20,
              profileParams=(1./3.,hs,0.2))
calc_pred(pred_file,dfc,nls,nds,ls,ds)
