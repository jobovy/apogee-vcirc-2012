FITSDIR=../fits_indivfeh_correct/
NSAMPLES=10000

.PHONY: all_simpledrift_noro_dwarf_vpec_sratio_hs \
	all_simpledrift_noro_dwarf_vpec_sratio \
	all_simpledrift_noro_dwarf_linear_vpec_sratio_hs \
	all_simpledrift_noro_dwarf_linear_vpec_sratio_ah_hs_innerouterah \
	all_simpledrift_noro_dwarf_quadratic_vpec_sratio_hs \
	all_simpledrift_noro_dwarf_quadratic_vpec_sratio_ah_hs_innerouterah \
	all_simpledrift_noro_dwarf_cubic_vpec_sratio_hs \
	all_simpledrift_noro_dwarf_cubic_vpec_sratio_ah_hs_innerouterah \
	all_simpledrift_noro_dwarf_powerlaw_vpec_sratio \
	all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs \
	all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_ah_hs_innerouterah \
	all_simpledrift_noro_dwarf_vpec_sratioinnerouter_hs \
	all_simpledrift_noro_dwarf_vpec_sratioinnerouter \
	all_simpledrift_noro_dwarf_vpec_sratio_hs_nooutliermean \
	all_simpledrift_noro_dwarf_vpec_sratio_m2 \
	all_simpledrift_noro_dwarf_vpec_sratio_innerouter \
	all_simpledrift_noro_dwarf_vpec_sratio_ah_hs_innerouterah \
	all_simpledrift_noro_dwarf_vpec_sratio_dm_innerouterdm \
	all_simpledrift_noro_dwarf_vpec_sratio_ah_hs \
	all_simpledrift_noro_dwarf_vpec_sratio_dm_hs \
	all_simpledrift_noro_dwarf_vpec_sratio_hs_innerouter \
	all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_innerouter \
	all_simpledrift_noro_dwarf_powerlaw_vpec_sratioinnerouter_innerouter \
	all_multiplepops_noro_dwarf_vpec_sratio \
	all_twopops_noro_dwarf_vpec_sratio_hs \
	loc4240_simpledrift_noro_dwarf_vpec_sratio \
	allwoloc4240_simpledrift_noro_dwarf_vpec_sratio_hs \
	allwoloc4241_simpledrift_noro_dwarf_vpec_sratio_hs

all_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fiths

all_simpledrift_noro_dwarf_vpec_sratio_hs_dl:
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_dl.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift --fitdl
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_dl_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_dl.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fiths --fitdl

all_simpledrift_vo250_ro8_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_vo250_ro8_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift --fixvo=250. --fixro=8.
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_vo250_ro8_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_vo250_ro8_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fiths --fixvo=250. --fixro=8.

all_simpledrift_vo250_ro8.4_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_vo250_ro8.4_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift --fixvo=250. --fixro=8.4
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_vo250_ro8.4_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_vo250_ro8.4_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fiths --fixvo=250. --fixro=8.4

all_simpledrift_vo250_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_vo250_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift --fixvo=250.
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_vo250_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_vo250_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fiths --fixvo=250.

all_simpledrift_noro_dwarf_sbdvpec_sratio_hs:
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_sbdvpec_sratio_hs.sav \
	--noroprior --dwarf --sbdvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_sbdvpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_sbdvpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --sbdvpec --fitsratio --fiths

all_simpledrift_noro_dwarf_vpec_sratio_hs_hr2:
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat --hr=2. \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_hr2.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat --hr=2. \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_hr2_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_hr2.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fiths

all_simpledrift_noro_dwarf_vpec_sratio_hs_hr4:
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat --hr=4. \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_hr4.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat --hr=4. \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_hr4_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_hr4.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fiths

all_simpledrift_noro_dwarf_vpec_sratio_hs_justinner:
	$(PYTHON) fitvc.py --justinner \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_justinner.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py --justinner \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_justinner_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_justinner.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fiths

all_simpledrift_noro_dwarf_vpec_sratio_hs_justouter:
	$(PYTHON) fitvc.py --justouter \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_justouter.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py --justouter \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_justouter_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_justouter.sav \
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

all_simpledrift_noro_dwarf_vpec_sratio_ah_hs_innerouterah:
	$(PYTHON) fitvc.py --fitah --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_ah_hs_innerouterah.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fitahinnerouter \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py --fitah --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_ah_hs_innerouterah_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_ah_hs_innerouterah.sav \
	 --mcsample --nsamples=$(NSAMPLES) --fitahinnerouter \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_simpledrift_noro_dwarf_vpec_sratio_feh_hs_innerouterfeh:
	$(PYTHON) fitvc.py --fitfeh --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_feh_hs_innerouterfeh.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fitfehinnerouter \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py --fitfeh --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_feh_hs_innerouterfeh_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_feh_hs_innerouterfeh.sav \
	 --mcsample --nsamples=$(NSAMPLES) --fitfehinnerouter \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_simpledrift_noro_dwarf_vpec_sratio_dm_hs_innerouterdm:
	$(PYTHON) fitvc.py --fitdm --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_dm_hs_innerouterdm.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fitdminnerouter \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py --fitdm --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_dm_hs_innerouterdm_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_dm_hs_innerouterdm.sav \
	 --mcsample --nsamples=$(NSAMPLES) --fitdminnerouter \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_simpledrift_noro_dwarf_vpec_sratio_ah_hs:
	$(PYTHON) fitvc.py --fitah --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_ah_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py --fitah --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_ah_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_ah_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) --fiths \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio 

all_simpledrift_noro_dwarf_vpec_sratio_feh_hs:
	$(PYTHON) fitvc.py --fitfeh --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_feh_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py --fitfeh --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_feh_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_feh_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) --fiths \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio 

all_simpledrift_noro_dwarf_vpec_sratio_dm_hs:
	$(PYTHON) fitvc.py --fitdm --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_dm_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py --fitdm --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_dm_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_dm_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) --fiths \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio 

all_simpledrift_noro_dwarf_vpec_sratio_hs_innerouter:
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_innerouter.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift --fitsrinnerouter --dwarfinnerouter
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_innerouter_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_hs_innerouter.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fitsrinnerouter --dwarfinnerouter

all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_innerouter:
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_innerouter.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --fitsrinnerouter --dwarfinnerouter
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_innerouter_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_innerouter.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=powerlaw \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fitsrinnerouter --dwarfinnerouter

all_simpledrift_noro_dwarf_powerlaw_vpec_sratioinnerouter_innerouter:
	$(PYTHON) fitvc.py \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratioinnerouter_innerouter.sav \
	--noroprior --dwarf --fitvpec --fitsratioinnerouter --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --fitsrinnerouter --dwarfinnerouter
	$(PYTHON) fitvc.py \
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

all_simpledrift_noro_dwarf_linear_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=linear
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_linear_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=linear \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_simpledrift_noro_dwarf_linear_vpec_sratio_ah_hs_innerouterah:
	$(PYTHON) fitvc.py --fiths --fitah --fitahinnerouter \
	$(FITSDIR)all_simpledrift_noro_dwarf_linear_vpec_sratio_ah_hs_innerouterah.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=linear
	$(PYTHON) fitvc.py --fiths --fitah --fitahinnerouter \
	$(FITSDIR)all_simpledrift_noro_dwarf_linear_vpec_sratio_ah_hs_innerouterah_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_linear_vpec_sratio_ah_hs_innerouterah.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=linear \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_simpledrift_noro_dwarf_quadratic_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_quadratic_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=quadratic
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_quadratic_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_quadratic_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=quadratic \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_simpledrift_noro_dwarf_quadratic_vpec_sratio_ah_hs_innerouterah:
	$(PYTHON) fitvc.py --fiths --fitah --fitahinnerouter \
	$(FITSDIR)all_simpledrift_noro_dwarf_quadratic_vpec_sratio_ah_hs_innerouterah.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=quadratic
	$(PYTHON) fitvc.py --fiths --fitah --fitahinnerouter \
	$(FITSDIR)all_simpledrift_noro_dwarf_quadratic_vpec_sratio_ah_hs_innerouterah_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_quadratic_vpec_sratio_ah_hs_innerouterah.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=quadratic \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_simpledrift_noro_dwarf_powerlaw_vpec_sratio:
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=powerlaw
	$(PYTHON) fitvc.py --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=powerlaw \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=powerlaw
	$(PYTHON) fitvc.py --fiths --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=powerlaw \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_hr2:
	$(PYTHON) fitvc.py --fiths --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_hr2.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=powerlaw --hr=2.
	$(PYTHON) fitvc.py --fiths --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_hr2_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_hr2.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=powerlaw \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --hr=2.

all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_hr4:
	$(PYTHON) fitvc.py --fiths --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_hr4.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=powerlaw --hr=4.
	$(PYTHON) fitvc.py --fiths --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_hr4_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_hr4.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=powerlaw \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --hr=4.

all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_ah_hs_innerouterah:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_ah_hs_innerouterah.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=powerlaw --fitah --fitahinnerouter
	$(PYTHON) fitvc.py --fiths --fitah --fitahinnerouter \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_ah_hs_innerouterah_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_ah_hs_innerouterah.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=powerlaw \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_feh_hs_innerouterfeh:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_feh_hs_innerouterfeh.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=powerlaw --fitfeh --fitfehinnerouter
	$(PYTHON) fitvc.py --fiths --fitfeh --fitfehinnerouter \
	$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_feh_hs_innerouterfeh_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_powerlaw_vpec_sratio_feh_hs_innerouterfeh.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=powerlaw \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_simpledrift_noro_dwarf_cubic_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_cubic_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=cubic
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)all_simpledrift_noro_dwarf_cubic_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_cubic_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=cubic \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_simpledrift_noro_dwarf_cubic_vpec_sratio_ah_hs_innerouterah:
	$(PYTHON) fitvc.py --fiths --fitah --fitahinnerouter \
	$(FITSDIR)all_simpledrift_noro_dwarf_cubic_vpec_sratio_ah_hs_innerouterah.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --rotcurve=cubic
	$(PYTHON) fitvc.py --fiths --fitah --fitahinnerouter \
	$(FITSDIR)all_simpledrift_noro_dwarf_cubic_vpec_sratio_ah_hs_innerouterah_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_cubic_vpec_sratio_ah_hs_innerouterah.sav \
	 --mcsample --nsamples=$(NSAMPLES) --rotcurve=cubic \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio

all_multiplepops_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)all_multiplepops_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=multiplepops
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)all_multiplepops_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_multiplepops_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=multiplepops \
	--dwarf --fitvpec --fitsratio

all_twopops_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_twopops_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=twopops
	$(PYTHON) fitvc.py --fiths --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_twopops_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_twopops_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=twopops \
	--dwarf --fitvpec --fitsratio

all_twopops_noro_dwarf_vpec_sratio_hs_hr26:
	$(PYTHON) fitvc.py --fiths --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_twopops_noro_dwarf_vpec_sratio_hs_hr26.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=twopops --hr2=6.
	$(PYTHON) fitvc.py --fiths --isofile=$(FITSDIR)isos.dat \
	$(FITSDIR)all_twopops_noro_dwarf_vpec_sratio_hs_hr26_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_twopops_noro_dwarf_vpec_sratio_hs_hr26.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=twopops \
	--dwarf --fitvpec --fitsratio --hr2=6.

loc4240_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)loc4240_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --location=4240
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)loc4240_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)loc4240_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --location=4240

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

allwoloc4240_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4240_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --removelocation=4240
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4240_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4240_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --removelocation=4240

allwoloc4241_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4241_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --removelocation=4241
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4241_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4241_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --removelocation=4241

allwoloc4242_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4242_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --removelocation=4242
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4242_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4242_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --removelocation=4242

allwoloc4270_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4270_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --removelocation=4270
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4270_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4270_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --removelocation=4270

allwoloc4243_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4243_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --removelocation=4243
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4243_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4243_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --removelocation=4243

allwoloc4318_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4318_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --removelocation=4318
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4318_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4318_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --removelocation=4318

allwoloc4271_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4271_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --removelocation=4271
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4271_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4271_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --removelocation=4271

allwoloc4319_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4319_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --removelocation=4319
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4319_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4319_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --removelocation=4319

allwoloc4151_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4151_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --removelocation=4151
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4151_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4151_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --removelocation=4151

allwoloc4272_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4272_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --removelocation=4272
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4272_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4272_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --removelocation=4272

allwoloc4321_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4321_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --removelocation=4321
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4321_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4321_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --removelocation=4321

allwoloc4154_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4154_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --removelocation=4154
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4154_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4154_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --removelocation=4145

allwoloc4273_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4273_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --removelocation=4273
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4273_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4273_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --removelocation=4273

allwoloc4157_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4157_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio \
	--dfmodel=simplegaussiandrift --removelocation=4157
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4157_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4157_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --removelocation=4157

all_simpledrift_noro_dwarf_vpec_sratio_dm_z21:
	$(PYTHON) fitvc.py --fitdm --Z=0.021 \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_dm_hs_z21.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift 
	$(PYTHON) fitvc.py --fitdm \
	$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_dm_hs_z21_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)all_simpledrift_noro_dwarf_vpec_sratio_dm_hs_z21.sav \
	 --mcsample --nsamples=$(NSAMPLES) --fiths \
	--noroprior --dfmodel=simplegaussiandrift --Z=0.021 \
	--dwarf --fitvpec --fitsratio

short_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --cohort=short \
	$(FITSDIR)short_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py --cohort=short \
	$(FITSDIR)short_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)short_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fiths

medium_simpledrift_noro_dwarf_vpec_sratio_hs:
	$(PYTHON) fitvc.py --cohort=medium \
	$(FITSDIR)medium_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --fiths \
	--dfmodel=simplegaussiandrift
	$(PYTHON) fitvc.py --cohort=medium \
	$(FITSDIR)medium_simpledrift_noro_dwarf_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)medium_simpledrift_noro_dwarf_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --fiths

#powerlaw jackknife
allwoloc4240_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4240_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --removelocation=4240
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4240_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4240_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=powerlaw --removelocation=4240

allwoloc4241_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4241_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --removelocation=4241
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4241_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4241_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=powerlaw --removelocation=4241

allwoloc4242_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4242_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --removelocation=4242
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4242_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4242_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=powerlaw --removelocation=4242

allwoloc4270_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4270_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --removelocation=4270
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4270_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4270_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=powerlaw --removelocation=4270

allwoloc4243_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4243_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --removelocation=4243
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4243_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4243_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=powerlaw --removelocation=4243

allwoloc4318_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4318_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --removelocation=4318
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4318_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4318_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=powerlaw --removelocation=4318

allwoloc4271_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4271_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --removelocation=4271
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4271_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4271_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=powerlaw --removelocation=4271

allwoloc4319_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4319_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --removelocation=4319
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4319_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4319_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=powerlaw --removelocation=4319

allwoloc4151_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4151_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --removelocation=4151
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4151_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4151_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=powerlaw --removelocation=4151

allwoloc4272_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4272_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --removelocation=4272
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4272_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4272_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=powerlaw --removelocation=4272

allwoloc4321_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4321_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --removelocation=4321
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4321_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4321_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=powerlaw --removelocation=4321

allwoloc4154_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4154_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --removelocation=4154
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4154_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4154_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=powerlaw --removelocation=4145

allwoloc4273_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4273_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --removelocation=4273
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4273_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4273_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=powerlaw --removelocation=4273

allwoloc4157_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4157_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=powerlaw \
	--dfmodel=simplegaussiandrift --removelocation=4157
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4157_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4157_simpledrift_noro_dwarf_powerlaw_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=powerlaw --removelocation=4157


#Linear Jackknife
allwoloc4240_simpledrift_noro_dwarf_linear_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4240_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=linear \
	--dfmodel=simplegaussiandrift --removelocation=4240
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4240_simpledrift_noro_dwarf_linear_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4240_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=linear --removelocation=4240

allwoloc4241_simpledrift_noro_dwarf_linear_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4241_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=linear \
	--dfmodel=simplegaussiandrift --removelocation=4241
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4241_simpledrift_noro_dwarf_linear_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4241_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=linear --removelocation=4241

allwoloc4242_simpledrift_noro_dwarf_linear_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4242_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=linear \
	--dfmodel=simplegaussiandrift --removelocation=4242
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4242_simpledrift_noro_dwarf_linear_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4242_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=linear --removelocation=4242

allwoloc4270_simpledrift_noro_dwarf_linear_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4270_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=linear \
	--dfmodel=simplegaussiandrift --removelocation=4270
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4270_simpledrift_noro_dwarf_linear_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4270_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=linear --removelocation=4270

allwoloc4243_simpledrift_noro_dwarf_linear_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4243_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=linear \
	--dfmodel=simplegaussiandrift --removelocation=4243
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4243_simpledrift_noro_dwarf_linear_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4243_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=linear --removelocation=4243

allwoloc4318_simpledrift_noro_dwarf_linear_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4318_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=linear \
	--dfmodel=simplegaussiandrift --removelocation=4318
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4318_simpledrift_noro_dwarf_linear_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4318_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=linear --removelocation=4318

allwoloc4271_simpledrift_noro_dwarf_linear_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4271_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=linear \
	--dfmodel=simplegaussiandrift --removelocation=4271
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4271_simpledrift_noro_dwarf_linear_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4271_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=linear --removelocation=4271

allwoloc4319_simpledrift_noro_dwarf_linear_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4319_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=linear \
	--dfmodel=simplegaussiandrift --removelocation=4319
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4319_simpledrift_noro_dwarf_linear_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4319_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=linear --removelocation=4319

allwoloc4151_simpledrift_noro_dwarf_linear_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4151_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=linear \
	--dfmodel=simplegaussiandrift --removelocation=4151
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4151_simpledrift_noro_dwarf_linear_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4151_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=linear --removelocation=4151

allwoloc4272_simpledrift_noro_dwarf_linear_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4272_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=linear \
	--dfmodel=simplegaussiandrift --removelocation=4272
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4272_simpledrift_noro_dwarf_linear_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4272_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=linear --removelocation=4272

allwoloc4321_simpledrift_noro_dwarf_linear_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4321_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=linear \
	--dfmodel=simplegaussiandrift --removelocation=4321
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4321_simpledrift_noro_dwarf_linear_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4321_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=linear --removelocation=4321

allwoloc4154_simpledrift_noro_dwarf_linear_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4154_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=linear \
	--dfmodel=simplegaussiandrift --removelocation=4154
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4154_simpledrift_noro_dwarf_linear_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4154_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=linear --removelocation=4145

allwoloc4273_simpledrift_noro_dwarf_linear_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4273_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=linear \
	--dfmodel=simplegaussiandrift --removelocation=4273
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4273_simpledrift_noro_dwarf_linear_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4273_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=linear --removelocation=4273

allwoloc4157_simpledrift_noro_dwarf_linear_vpec_sratio_hs:
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4157_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	--noroprior --dwarf --fitvpec --fitsratio --rotcurve=linear \
	--dfmodel=simplegaussiandrift --removelocation=4157
	$(PYTHON) fitvc.py --fiths \
	$(FITSDIR)allwoloc4157_simpledrift_noro_dwarf_linear_vpec_sratio_hs_$(NSAMPLES)samples.sav \
	--init=$(FITSDIR)allwoloc4157_simpledrift_noro_dwarf_linear_vpec_sratio_hs.sav \
	 --mcsample --nsamples=$(NSAMPLES) \
	--noroprior --dfmodel=simplegaussiandrift \
	--dwarf --fitvpec --fitsratio --rotcurve=linear --removelocation=4157



