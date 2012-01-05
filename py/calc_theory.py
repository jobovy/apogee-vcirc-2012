#Calculate some theoretical stuff
import os, os.path
import cPickle as pickle
import math
import numpy
from galpy.df import dehnendf, shudf
from galpy.util import bovy_coords
#Overall setup
nls= 501
ls= numpy.linspace(0.,360.,nls)
nds= 101
ds= numpy.linspace(0.,10./8.,nds)
_SAVEDIR= os.path.join(os.getenv('DATADIR'),
                       'bovy',
                       'vclos')
#Quick def
def safe_dl_to_rphi(d,l):
    R= math.sqrt(1.+d**2.-2.*d*math.cos(l))
    if R == 0.:
        R+= 0.0001
        d+= 0.0001
    if 1./math.cos(l) < d and math.cos(l) > 0.:
        theta= math.pi-math.asin(d/R*math.sin(l))
    else:
        theta= math.asin(d/R*math.sin(l))
    return (R,theta)
def calc_pred(pred_file,dfc,nls,nds,ls):
    avg_pred= numpy.zeros(nls)
    for ii in range(nls):
        meanvlos= 0.
        norm= 0.
        for jj in range(nds):
            d,l= ds[jj], ls[ii]/180.*numpy.pi
            R,theta= safe_dl_to_rphi(d,l)
            surfmass= dfc.surfacemassLOS(d,l,deg=False)
            meanvlos+= surfmass*dfc.meanvT(R)*math.sin(theta+l)
            norm+= surfmass
        avg_pred[ii]= meanvlos/norm
    savefile= open(pred_file,'wb')
    pickle.dump(ls,savefile)
    pickle.dump(avg_pred,savefile)
    savefile.close()

#Start calculating
#Fiducial
pred_file= os.path.join(_SAVEDIR,'l_vhelio_fid.sav')
if not os.path.exists(pred_file):
    print "Working on fiducial ..."
    dfc= dehnendf(beta=0.,correct=True,niter=20)
    calc_pred(pred_file,dfc,nls,nds,ls)
#betas
betas= [-0.2,-0.1,0.1,0.2]
for beta in betas:
    pred_file= os.path.join(_SAVEDIR,'l_vhelio_beta_%.2f.sav' % beta)
    if not os.path.exists(pred_file):
        print "Working on beta %.2f ..." % beta
        dfc= dehnendf(beta=beta,correct=True,niter=20)
        calc_pred(pred_file,dfc,nls,nds,ls)
#hr
hrs= [1./4.,1./2.]
for hr in hrs:
    pred_file= os.path.join(_SAVEDIR,'l_vhelio_hr_%.6f.sav' % hr)
    if not os.path.exists(pred_file):
        print "Working on hr %.2f ..." % hr
        dfc= dehnendf(beta=0.,correct=True,niter=20,
                      profileParams=(hr,1.,0.2))
        calc_pred(pred_file,dfc,nls,nds,ls)
#hs
hs= 2./3.
pred_file= os.path.join(_SAVEDIR,'l_vhelio_hs_%.6f.sav' % hs)
if not os.path.exists(pred_file):
    print "Working on hs %.2f ..." % hs
    dfc= dehnendf(beta=0.,correct=True,niter=20,
                  profileParams=(1./3.,hs,0.2))
    calc_pred(pred_file,dfc,nls,nds,ls)
#shu
pred_file= os.path.join(_SAVEDIR,'l_vhelio_shu.sav')
if not os.path.exists(pred_file):
    print "Working on Shu ..."
    dfc= shudf(beta=0.,correct=True,niter=20)
    calc_pred(pred_file,dfc,nls,nds,ls)
altdmaxs= [5./8.,15./8.]
for altdmax in altdmaxs:
    thisds= numpy.linspace(0.,altdmax,nds)
    pred_file= os.path.join(_SAVEDIR,'l_vhelio_dmax_%.6f.sav' % altdmax)
    if not os.path.exists(pred_file):
        print "Working on dmax %.2f ..." % altdmax
        dfc= dehnendf(beta=0.,correct=True,niter=20)
        calc_pred(pred_file,dfc,nls,thisnds,ls)

