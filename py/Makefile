FITSDIR=../fits/
NSAMPLES=1000

.PHONY: all_simple all_simple_dwarf all_simple_noro all_simpledrift_noro \
	lgtr90_simpledrift_noro \
	all_simpledrift_noro_dwarf \
	all_simpledrift_noro_hs16 \
	all_simpledrift_noro_powerlaw all_simpledrift_noro_linear \
	all_simpledrift_noro_quadratic \
	all_simpledrift_noro_cubic \
	all_simpledrift_noro_vpec \
	all_simpledrift_noro_vpec_nods \
	all_simpledrift_noro_vpec_sratio \
	l180_simpledrift_noro \
	fourth_simpledrift_noro \
	nomult_simpledrift_noro \
	all_simpledrift_noro_gp \
	all_simpledrift_noro_gpemcee \
	all_simpledrift_noro_gp_fixtau0.5 \
	all_skeweddrift_noro \
	all_dehnen_noro

all_simple:
	$(PYTHON) fitvc.py $(FITSDIR)all_simple.sav
	$(PYTHON) fitvc.py $(FITSDIR)all_simple_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simple.sav --mcsample --nsamples=$(NSAMPLES)

all_simple_dwarf:
	$(PYTHON) fitvc.py $(FITSDIR)all_simple_dwarf.sav --dwarf
	$(PYTHON) fitvc.py $(FITSDIR)all_simple_dwarf_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simple_dwarf.sav --mcsample --nsamples=$(NSAMPLES) \
	--dwarf

all_simple_noro:
	$(PYTHON) fitvc.py $(FITSDIR)all_simple_noro.sav --noroprior
	$(PYTHON) fitvc.py $(FITSDIR)all_simple_noro_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simple_noro.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior

all_simpledrift_noro:
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro.sav --noroprior \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift
#	$(PYTHON) compareDataModel.py \
#	--init=$(FITSDIR)all_simpledrift_noro.sav \
#	--plottype=pvloslos -o ~/Desktop/test.png --nvlos=21 -m 4 \
#	 --dfmodel=simplegaussiandrift

all_simpledrift_noro_dwarf:
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_dwarf.sav --noroprior \
	--dfmodel=simplegaussiandrift --dwarf
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_dwarf_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift --dwarf

all_simpledrift_noro_hs16:
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_hs16.sav --noroprior \
	--dfmodel=simplegaussiandrift --hs=16.
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_hs16_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_hs16.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift --hs=16.

all_skeweddrift_noro:
	$(PYTHON) -W ignore fitvc.py $(FITSDIR)all_skeweddrift_noro.sav \
	--noroprior \
	--dfmodel=simpleskeweddrift --multi=16
	$(PYTHON) fitvc.py $(FITSDIR)all_skeweddrift_noro_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_skeweddrift_noro.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simpleskeweddrift

all_dehnen_noro:
	$(PYTHON) -W ignore fitvc.py $(FITSDIR)all_dehnen_noro.sav \
	--noroprior \
	--dfmodel=dehnen --multi=16
	$(PYTHON) fitvc.py $(FITSDIR)all_dehnen_noro_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_dehnen_noro.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=dehnen --multi=16

lgtr90_simpledrift_noro:
	$(PYTHON) fitvc.py $(FITSDIR)lgtr90_simpledrift_noro.sav --noroprior \
	--dfmodel=simplegaussiandrift --lmin=90.
	$(PYTHON) fitvc.py $(FITSDIR)lgtr90_simpledrift_noro_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)lgtr90_simpledrift_noro.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift --lmin=90.

l180_nomult_simpledrift_noro:
	$(PYTHON) fitvc.py $(FITSDIR)l180_nomult_simpledrift_noro.sav \
	--noroprior \
	--dfmodel=simplegaussiandrift --location=4154 --cutmultiples
	$(PYTHON) fitvc.py $(FITSDIR)l180_nomult_simpledrift_noro_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)l180_nomult_simpledrift_noro.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift --location=4154 --cutmultiples

fourth_simpledrift_noro:
	$(PYTHON) fitvc.py $(FITSDIR)fourth_simpledrift_noro.sav --noroprior \
	--dfmodel=simplegaussiandrift --downsample=4
	$(PYTHON) fitvc.py $(FITSDIR)fourth_simpledrift_noro_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)fourth_simpledrift_noro.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift --downsample=4

nomult_simpledrift_noro:
	$(PYTHON) fitvc.py $(FITSDIR)nomult_simpledrift_noro.sav --noroprior \
	--dfmodel=simplegaussiandrift --cutmultiples
	$(PYTHON) fitvc.py $(FITSDIR)nomult_simpledrift_noro_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)nomult_simpledrift_noro.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift --cutmultiples

all_simpledrift_noro_powerlaw:
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_powerlaw.sav --noroprior \
	--dfmodel=simplegaussiandrift --rotcurve=powerlaw
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_powerlaw_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_powerlaw.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift --rotcurve=powerlaw

all_simpledrift_noro_gp:
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_gp_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift --rotcurve=gp

all_simpledrift_noro_gp_fixtau0.5:
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_gp_fixtau0.5_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift --rotcurve=gp --gpfixtau=0.5

all_simpledrift_noro_gpemcee:
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_gpemcee_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift --rotcurve=gp --gpemcee \
	--gpnr=21

all_simpledrift_noro_linear:
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_linear.sav --noroprior \
	--dfmodel=simplegaussiandrift --rotcurve=linear
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_linear_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_linear.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift --rotcurve=linear

all_simpledrift_noro_quadratic:
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_quadratic.sav --noroprior \
	--dfmodel=simplegaussiandrift --rotcurve=quadratic
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_quadratic_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_quadratic.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift --rotcurve=quadratic

all_simpledrift_noro_cubic:
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_cubic.sav --noroprior \
	--dfmodel=simplegaussiandrift --rotcurve=cubic
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_cubic_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_cubic.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift --rotcurve=cubic

all_simpledrift_noro_vpec:
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_vpec.sav --noroprior \
	--dfmodel=simplegaussiandrift --fitvpec
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_vpec_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_vpec.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift --fitvpec

all_simpledrift_noro_vpec_nods:
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_vpec_nods.sav --noroprior \
	--dfmodel=simplegaussiandrift --fitvpec --nods
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_vpec_nods_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_vpec_nods.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift --fitvpec --nods

all_simpledrift_noro_vpec_sratio:
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_vpec_sratio.sav --noroprior \
	--dfmodel=simplegaussiandrift --fitvpec --fitsratio
	$(PYTHON) fitvc.py $(FITSDIR)all_simpledrift_noro_vpec_sratio_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_vpec_sratio.sav --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift --fitvpec --fitsratio

#simple analysis without errors
../fakedata/l30_500_d0.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 500 $@ --distuncertainty=0. --los=30.
../plots/l30_500_d0_vclike.ps: map_vc_like_simple.py ../fakedata/l30_500_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_500_d0.sav \
	-o $@ --nvcirc=201 --los=30. --distuncertainty=0.
../plots/l30_500_dinf_vclike.ps: map_vc_like_simple.py \
	../fakedata/l30_500_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_500_d0.sav \
	-o $@ --nvcirc=101 --los=30. --distuncertainty=10.
../plots/l30_500_d0_vclinear.ps: map_vc_like_simple.py \
	../fakedata/l30_500_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_500_d0.sav \
	-o $@ --los=30. --distuncertainty=0. --linearfit
../plots/l30_500_d0_vlosd.ps: plot_fake_data.py ../fakedata/l30_500_d0.sav
	$(PYTHON) plot_fake_data.py ../fakedata/l30_500_d0.sav \
	--scatterbins=10 --vlosdname=$@ --los=30.

#simple analysis without errors, 10k
../fakedata/l30_10k_d0.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 10000 $@ --distuncertainty=0. --los=30.
../plots/l30_10k_d0_vclike.ps: map_vc_like_simple.py ../fakedata/l30_10k_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_10k_d0.sav \
	-o $@ --nvcirc=201 --los=30. --distuncertainty=0.
../plots/l30_10k_d0_vclinear.ps: map_vc_like_simple.py \
	../fakedata/l30_10k_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_10k_d0.sav \
	-o $@ --los=30. --distuncertainty=0. --linearfit

#simple analysis without errors, sigma_R=0.1
../fakedata/l30_500_so.1_d0.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 500 $@ --distuncertainty=0. --los=30. \
	--so=0.1
../plots/l30_500_d0_so.1_vclike.ps: map_vc_like_simple.py \
	../fakedata/l30_500_so.1_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_500_so.1_d0.sav \
	-o $@ --nvcirc=201 --los=30. --distuncertainty=0. --so=0.1
../plots/l30_500_d0_so.1_fixdf_vclike.ps: map_vc_like_simple.py \
	../fakedata/l30_500_so.1_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_500_so.1_d0.sav \
	-o $@ --nvcirc=201 --los=30. --distuncertainty=0. --so=0.1 \
	--fixdfmoments
../plots/l30_500_d0_so.1_vclinear.ps: map_vc_like_simple.py \
	../fakedata/l30_500_so.1_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_500_so.1_d0.sav \
	-o $@ --so=0.1 --los=30. --linearfit

#simple analysis without errors, sigma_R=0.1, 10k
../fakedata/l30_10k_so.1_d0.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 10000 $@ \
	--distuncertainty=0. --los=30. \
	--so=0.1
../plots/l30_10k_d0_so.1_vclike.ps: map_vc_like_simple.py \
	../fakedata/l30_10k_so.1_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_10k_so.1_d0.sav \
	-o $@ --nvcirc=201 --los=30. --distuncertainty=0. --so=0.1
../plots/l30_10k_d0_so.1_fixdf_vclike.ps: map_vc_like_simple.py \
	../fakedata/l30_10k_so.1_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_10k_so.1_d0.sav \
	-o $@ --nvcirc=201 --los=30. --distuncertainty=0. --so=0.1 \
	--fixdfmoments
../plots/l30_10k_d0_so.1_vclinear.ps: map_vc_like_simple.py \
	../fakedata/l30_10k_so.1_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_10k_so.1_d0.sav \
	-o $@ --so=0.1 --los=30. --linearfit

#simple analysis without errors, sigma_R=0.1
../fakedata/l30_500_so.0125_d0.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 500 $@ --distuncertainty=0. --los=30. \
	--so=0.0125
../plots/l30_500_d0_so.0125_vclike.ps: map_vc_like_simple.py \
	../fakedata/l30_500_so.0125_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_500_so.0125_d0.sav \
	-o $@ --nvcirc=201 --los=30. --distuncertainty=0. --so=0.0125
../plots/l30_500_d0_so.0125_vclinear.ps: map_vc_like_simple.py \
	../fakedata/l30_500_so.0125_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_500_so.0125_d0.sav \
	-o $@ --so=0.0125 --los=30. --linearfit

#simple analysis without errors, l=120
../fakedata/l120_500_d0.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 500 $@ --distuncertainty=0. --los=120.
../plots/l120_500_d0_vclike.ps: map_vc_like_simple.py \
	../fakedata/l120_500_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l120_500_d0.sav \
	-o $@ --nvcirc=201 --los=120. --distuncertainty=0.
../plots/l120_500_dinf_vclike.ps: map_vc_like_simple.py \
	../fakedata/l120_500_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l120_500_d0.sav \
	-o $@ --nvcirc=201 --los=120. --distuncertainty=10.
../plots/l120_500_vlosd.ps: plot_fake_data.py ../fakedata/l120_500_d0.sav
	$(PYTHON) plot_fake_data.py ../fakedata/l120_500_d0.sav \
	--scatterbins=10 --vlosdname=$@ --los=120.
../plots/l120_500_d0_so.1_fixdf_vclike.ps: map_vc_like_simple.py \
	../fakedata/l120_500_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_500_d0.sav \
	-o $@ --nvcirc=201 --los=30. --distuncertainty=0. \
	--fixdfmoments

#simple analysis without errors, sigma_R=0.1, l=120
../fakedata/l120_500_so.1_d0.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 500 $@ --distuncertainty=0. \
	--los=120. \
	--so=0.1
../plots/l120_500_d0_so.1_vclike.ps: map_vc_like_simple.py \
	../fakedata/l120_500_so.1_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l120_500_so.1_d0.sav \
	-o $@ --nvcirc=201 --los=120. --distuncertainty=0. --so=0.1

#simple analysis without errors
../fakedata/l75_500_d0.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 500 $@ --distuncertainty=0. --los=75.
../fakedata/l75_500_d0_2.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 500 $@ --distuncertainty=0. --los=75. \
	--seed=2
../plots/l75_500_d0_vclike.ps: map_vc_like_simple.py ../fakedata/l75_500_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l75_500_d0.sav \
	-o $@ --nvcirc=201 --los=75. --distuncertainty=0.
../plots/l75_500_dinf_vclike.ps: map_vc_like_simple.py \
	../fakedata/l75_500_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l75_500_d0.sav \
	-o $@ --nvcirc=201 --los=75. --distuncertainty=10.
../plots/l75_500_dinf_2_vclike.ps: map_vc_like_simple.py \
	../fakedata/l75_500_d0_2.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l75_500_d0_2.sav \
	-o $@ --nvcirc=201 --los=75. --distuncertainty=10.
../plots/l75_500_d0_vclinear.ps: map_vc_like_simple.py \
	../fakedata/l75_500_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l75_500_d0.sav \
	-o $@ --los=75. --distuncertainty=0. --linearfit
../plots/l75_500_d0_vcbetalike.ps: map_vc_like_simple.py \
	../fakedata/l75_500_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l75_500_d0.sav \
	-o $@ --nvcirc=51 --los=75. --distuncertainty=0. \
	--nbeta=51
../plots/l75_500_vlosd.ps: plot_fake_data.py ../fakedata/l75_500_d0.sav
	$(PYTHON) plot_fake_data.py ../fakedata/l75_500_d0.sav \
	--scatterbins=10 --vlosdname=$@ --los=75.

#simple analysis without errors
../fakedata/l75_500_so.1_d0.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 500 $@ --distuncertainty=0. --los=75. \
	--so=0.1
../plots/l75_500_so.1_d0_vclike.ps: map_vc_like_simple.py \
	../fakedata/l75_500_so.1_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l75_500_so.1_d0.sav \
	-o $@ --nvcirc=201 --los=75. --distuncertainty=0. --so=0.1
../plots/l75_500_so.1_dinf_vclike.ps: map_vc_like_simple.py \
	../fakedata/l75_500_so.1_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l75_500_so.1_d0.sav \
	-o $@ --nvcirc=201 --los=75. --distuncertainty=10. --so=.1

#simple analysis without errors, fit for beta too
../plots/l30_500_so.1_d0_vcbetalike.ps: map_vc_like_simple.py \
	../fakedata/l30_500_so.1_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_500_so.1_d0.sav \
	-o $@ --nvcirc=51 --los=30. --distuncertainty=0. \
	--nbeta=51 --so=0.1
#simple analysis without errors, fit for beta too
../plots/l120_500_so.1_d0_vcbetalike.ps: map_vc_like_simple.py \
	../fakedata/l120_500_so.1_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l120_500_so.1_d0.sav \
	-o $@ --nvcirc=51 --los=120. --distuncertainty=0. \
	--nbeta=51 --so=0.1
#simple analysis without errors, fit for beta too
../plots/l30_500_d0_vcbetalike.ps: map_vc_like_simple.py \
	../fakedata/l30_500_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_500_d0.sav \
	-o $@ --nvcirc=51 --los=30. --distuncertainty=0. \
	--nbeta=51
#simple analysis without errors, fit for beta too
../plots/l120_500_d0_vcbetalike.ps: map_vc_like_simple.py \
	../fakedata/l120_500_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l120_500_d0.sav \
	-o $@ --nvcirc=51 --los=120. --distuncertainty=0. \
	--nbeta=51

#test data
../fakedata/l30_500.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 500 $@
../plots/l30_500_vlosd.ps: plot_fake_data.py ../fakedata/l30_500.sav
	$(PYTHON) plot_fake_data.py ../fakedata/l30_500.sav \
	--scatterbins=10 --vlosdname=$@

#local sample
../fakedata/local_1000.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 1000 $@ --local
../plots/local_1000_vt.ps: plot_fake_data.py ../fakedata/local_1000.sav
	$(PYTHON) plot_fake_data.py ../fakedata/local_1000.sav \
	--histbins=30 --vtname=$@
../plots/local_1000_vr.ps: plot_fake_data.py ../fakedata/local_1000.sav
	$(PYTHON) plot_fake_data.py ../fakedata/local_1000.sav \
	--histbins=30 --vrname=$@
../plots/local_1000_vrvt.ps: plot_fake_data.py ../fakedata/local_1000.sav
	$(PYTHON) plot_fake_data.py ../fakedata/local_1000.sav \
	--histbins=30 --vrvtname=$@

#proof that the sampling works
../fakedata/l30_10k.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 10000 $@
../plots/l30_10k_vlosd.ps: plot_fake_data.py ../fakedata/l30_10k.sav
	$(PYTHON) plot_fake_data.py ../fakedata/l30_10k.sav \
	--scatterbins=20 --vlosdname=$@