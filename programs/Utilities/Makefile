

DIRS += root2email
# DIRS += bfield2root file2et hddm2cMsg patfind

.PHONY: all install clean pristine relink

all:
	$(foreach dir,$(DIRS),make -C $(dir);)

install:
	$(foreach dir,$(DIRS),make -C $(dir) install;)

clean:
	$(foreach dir,$(DIRS),make -C $(dir) clean;)

pristine:
	$(foreach dir,$(DIRS),make -C $(dir) pristine;)

relink:
	$(foreach dir,$(DIRS),make -C $(dir) relink;)

