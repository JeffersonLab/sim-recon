
.PHONY: all clean

JLAB_WWW_DIR = /group/halld/WWW/software

all: clean
	doxygen Doxyfile
	cp images/* html

install:
	rm -rf $(JLAB_WWW_DIR)/junk/HDSoftware_Documentation
	touch $(JLAB_WWW_DIR)/HDSoftware_Documentation
	mv $(JLAB_WWW_DIR)/HDSoftware_Documentation $(JLAB_WWW_DIR)/junk
	cp -rf html $(JLAB_WWW_DIR)/HDSoftware_Documentation

clean:
	rm -rf html latex man
