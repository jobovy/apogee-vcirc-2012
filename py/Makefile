FITSDIR=../fits/
NSAMPLES=1000

.PHONY: all_simpledrift_noro_dwarf_vpec_sratio_hs

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

