
#simple analysis without errors
../fakedata/l30_500_d0.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 500 $@ --distuncertainty=0. --los=30.
../plots/l30_500_d0_vclike.ps: map_vc_like_simple.py ../fakedata/l30_500_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_500_d0.sav \
	-o $@ --nvcirc=201 --los=30. --distuncertainty=0.

#simple analysis without errors, sigma_R=0.1
../fakedata/l30_500_so.1_d0.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 500 $@ --distuncertainty=0. --los=30. \
	--so=0.1
../plots/l30_500_d0_so.1_vclike.ps: map_vc_like_simple.py \
	../fakedata/l30_500_so.1_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l30_500_so.1_d0.sav \
	-o $@ --nvcirc=201 --los=30. --distuncertainty=0. --so=0.1

#simple analysis without errors, l=120
../fakedata/l120_500_d0.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 500 $@ --distuncertainty=0. --los=120.
../plots/l120_500_d0_vclike.ps: map_vc_like_simple.py \
	../fakedata/l120_500_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l120_500_d0.sav \
	-o $@ --nvcirc=201 --los=120. --distuncertainty=0.

#simple analysis without errors, sigma_R=0.1, l=120
../fakedata/l120_500_so.1_d0.sav: create_fake_los.py
	$(PYTHON) create_fake_los.py -n 500 $@ --distuncertainty=0. \
	--los=120. \
	--so=0.1
../plots/l120_500_d0_so.1_vclike.ps: map_vc_like_simple.py \
	../fakedata/l120_500_so.1_d0.sav
	$(PYTHON) map_vc_like_simple.py ../fakedata/l120_500_so.1_d0.sav \
	-o $@ --nvcirc=201 --los=120. --distuncertainty=0. --so=0.1

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