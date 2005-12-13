

.PHONY: all


all:
ifdef JILIO
	make -C JILIO
endif
	#make -C include
	make -C BCAL
	make -C CDC
	make -C CHERENKOV
	make -C DANA
	make -C FCAL
	make -C FDC
	make -C TAGGER
	make -C TOF
	make -C TRIGGER
	make -C UPV
	make -C TRACKING
	make -C HDDM

clean:
ifdef JILIO
	make -C JILIO clean
endif
	make -C include clean
	make -C BCAL clean
	make -C CDC clean
	make -C CHERENKOV clean
	make -C DANA clean
	make -C FCAL clean
	make -C FDC clean
	make -C HDDM clean
	make -C TAGGER clean
	make -C TOF clean
	make -C TRIGGER clean
	make -C UPV clean
	make -C TRACKING clean
	
