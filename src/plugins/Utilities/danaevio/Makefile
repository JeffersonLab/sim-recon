#  makes DANAEVIO factory and event processor

#  E.Wolin, JLab, 1-Apr-2010



PACKAGES = DANA:ROOT:EVIO



include $(HALLD_HOME)/src/BMS/Makefile.shlib


# to ensure header file made from xml before dependencies generated
.depends/$(BMS_OSNAME)/DDANAEVIO_factory.d: dana_evio_dict.h




# Rule to make header file containing string from dana evio dictionary
%.h: %.xml
	$(HALLD_HOME)/src/BMS/TextFile2String.py $<
