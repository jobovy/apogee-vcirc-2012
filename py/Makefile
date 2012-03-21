FITSDIR=../fits/
NSAMPLES=1000

.PHONY: all_simpledrift_noro_dwarf_vpec_sratio_hs \
	all_simpledrift_noro_dwarf_linear_vpec_sratio \
	all_simpledrift_noro_dwarf_quadratic_vpec_sratio \
	all_simpledrift_noro_dwarf_cubic_vpec_sratio \
	all_simpledrift_noro_dwarf_powerlaw_vpec_sratio \
	all_simpledrift_noro_dwarf_vpec_sratioinnerouter_hs \
	all_simpledrift_noro_dwarf_vpec_sratioinnerouter \
	all_simpledrift_noro_dwarf_vpec_sratio_hs_nooutliermean \
	all_simpledrift_noro_dwarf_vpec_sratio_m2

all_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fiths

all_simpledrift_noro_dwarf_vpec_sratio_m2:
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_m2.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fitm2 \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_m2_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_m2.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fitm2

all_simpledrift_noro_dwarf_vpec_sratio_hs_nooutliermean:
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_nooutliermean.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift --nooutliermean
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_nooutliermean$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_nooutliermean.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fiths --nooutliermean

all_simpledrift_noro_dwarf_vpec_sratioinnerouter_hs:
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratioinnerouter_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratioinnerouter --fiths \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratioinnerouter_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratioinnerouter_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratioinnerouter --fiths

all_simpledrift_noro_dwarf_vpec_sratioinnerouter:
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratioinnerouter.sav \
	--noroprior --dwarf --fitvpec --fitsratioinnerouter \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratioinnerouter_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratioinnerouter.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratioinnerouter

all_simpledrift_noro_dwarf_linear_vpec_sratio:
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_linear_vpec_sratio.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=linear
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_linear_vpec_sratio_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_linear_vpec_sratio.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=linear \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_simpledrift_noro_dwarf_quadratic_vpec_sratio:
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_quadratic_vpec_sratio.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=quadratic
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_quadratic_vpec_sratio_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_quadratic_vpec_sratio.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=quadratic \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_simpledrift_noro_dwarf_powerlaw_vpec_sratio:
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=powerlaw
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=powerlaw \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_simpledrift_noro_dwarf_cubic_vpec_sratio:
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_cubic_vpec_sratio.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=cubic
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_cubic_vpec_sratio_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_cubic_vpec_sratio.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=cubic \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

