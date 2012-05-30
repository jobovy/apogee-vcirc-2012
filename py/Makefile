FITSDIR=../fits/
NSAMPLES=1000

.PHONY: all_simpledrift_noro_dwarf_vpec_sratio_hs \
	all_simpledrift_noro_dwarf_vpec_sratio \
	all_simpledrift_noro_dwarf_linear_vpec_sratio \
	all_simpledrift_noro_dwarf_quadratic_vpec_sratio \
	all_simpledrift_noro_dwarf_cubic_vpec_sratio \
	all_simpledrift_noro_dwarf_powerlaw_vpec_sratio \
	all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs \
	all_simpledrift_noro_dwarf_vpec_sratioinnerouter_hs \
	all_simpledrift_noro_dwarf_vpec_sratioinnerouter \
	all_simpledrift_noro_dwarf_vpec_sratio_hs_nooutliermean \
	all_simpledrift_noro_dwarf_vpec_sratio_m2 \
	all_simpledrift_noro_dwarf_vpec_sratio_innerouter \
	all_simpledrift_noro_dwarf_vpec_sratio_ah_innerouterah \
	all_simpledrift_noro_dwarf_vpec_sratio_dm_innerouterdm \
	all_simpledrift_noro_dwarf_vpec_sratio_hs_innerouter \
	all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_innerouter \
	all_simpledrift_noro_dwarf_powerlaw_vpec_sratioinnerouter_innerouter \
	all_multiplepops_noro_dwarf_vpec_sratio \
	loc4240_simpledrift_noro_dwarf_vpec_sratio \
	allwoloc4240_simpledrift_noro_dwarf_vpec_sratio

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

all_simpledrift_noro_dwarf_vpec_sratio:
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

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

all_simpledrift_noro_dwarf_vpec_sratio_innerouter:
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_innerouter.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --fitsrinnerouter --dwarfinnerouter
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_innerouter_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_innerouter.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fitsrinnerouter --dwarfinnerouter

all_simpledrift_noro_dwarf_vpec_sratio_ah_innerouterah:
	$(PYTHON) fitvc.py --fitah --nooutliermean \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_ah_innerouterah.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fitahinnerouter \
	--dfmodel=simplegaussiandrift --fitsrinnerouter --dwarfinnerouter
	$(PYTHON) fitvc.py --fitah --nooutliermean \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_ah_innerouterah_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_ah_innerouterah.sav \
	 --mcsample --nsamples=$(NSAMPLES) --fitahinnerouter \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fitsrinnerouter --dwarfinnerouter

all_simpledrift_noro_dwarf_vpec_sratio_hs_innerouter:
	$(PYTHON) fitvc.py --nooutliermean \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_innerouter.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift --fitsrinnerouter --dwarfinnerouter
	$(PYTHON) fitvc.py --fiths --nooutliermean \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_innerouter_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_innerouter.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fitsrinnerouter --dwarfinnerouter

all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_innerouter:
	$(PYTHON) fitvc.py --nooutliermean \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_innerouter.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --fitsrinnerouter --dwarfinnerouter
	$(PYTHON) fitvc.py --nooutliermean \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_innerouter_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_innerouter.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=powerlaw \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fitsrinnerouter --dwarfinnerouter

all_simpledrift_noro_dwarf_powerlaw_vpec_sratioinnerouter_innerouter:
	$(PYTHON) fitvc.py --nooutliermean \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratioinnerouter_innerouter.sav \
	--noroprior --dwarf --fitvpec --fitsratioinnerouter --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --fitsrinnerouter --dwarfinnerouter
	$(PYTHON) fitvc.py --nooutliermean \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratioinnerouter_innerouter_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratioinnerouter_innerouter.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=powerlaw \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratioinnerouter --fitsrinnerouter --dwarfinnerouter

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
	$(PYTHON) fitvc.py --nooutliermean \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=powerlaw
	$(PYTHON) fitvc.py --nooutliermean \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=powerlaw \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs:
	$(PYTHON) fitvc.py --nooutliermean --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=powerlaw
	$(PYTHON) fitvc.py --nooutliermean --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
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

all_multiplepops_noro_dwarf_vpec_sratio:
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_multiplepops_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=multiplepops
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_multiplepops_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_multiplepops_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=multiplepops \
	--dwarf --fitvpec --fitsratio

loc4240_simpledrift_noro_dwarf_vpec_sratio:
	$(PYTHON) fitvc.py \
	$(FITSDIR)loc4240_simpledrift_noro_dwarf_vpec_sratio.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --location=4240
	$(PYTHON) fitvc.py \
	$(FITSDIR)loc4240_simpledrift_noro_dwarf_vpec_sratio_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)loc4240_simpledrift_noro_dwarf_vpec_sratio.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --location=4240

allwoloc4240_simpledrift_noro_dwarf_vpec_sratio:
	$(PYTHON) fitvc.py \
	$(FITSDIR)allwoloc4240_simpledrift_noro_dwarf_vpec_sratio.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --removelocation=4240
	$(PYTHON) fitvc.py \
	$(FITSDIR)allwoloc4240_simpledrift_noro_dwarf_vpec_sratio_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4240_simpledrift_noro_dwarf_vpec_sratio.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --removelocation=4240

allwoloc4240_simpledrift_noro_dwarf_powerlaw_vpec_sratio:
	$(PYTHON) fitvc.py --rotcurve=powerlaw \
	$(FITSDIR)allwoloc4240_simpledrift_noro_dwarf_powerlaw_vpec_sratio.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --removelocation=4240
	$(PYTHON) fitvc.py \
	$(FITSDIR)allwoloc4240_simpledrift_noro_dwarf_powerlaw_vpec_sratio_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4240_simpledrift_noro_dwarf_powerlaw_vpec_sratio.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=powerlaw \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --removelocation=4240

all_simpledrift_noro_dwarf_vpec_sratio_dm_innerouterdm:
	$(PYTHON) fitvc.py --fitdm --nooutliermean --Z=0.021 \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_dm_innerouterdm.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fitdminnerouter \
	--dfmodel=simplegaussiandrift --fitsrinnerouter --dwarfinnerouter
	$(PYTHON) fitvc.py --fitdm --nooutliermean \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_dm_innerouterdm_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_dm_innerouterdm.sav \
	 --mcsample --nsamples=$(NSAMPLES) --fitdminnerouter \
	--noroprior --dfmodel=simplegaussiandrift --Z=0.021 \
	--dwarf --fitvpec --fitsratio --fitsrinnerouter --dwarfinnerouter

