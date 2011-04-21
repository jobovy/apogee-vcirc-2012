
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