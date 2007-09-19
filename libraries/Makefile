

DIRS += include BCAL CDC CHERENKOV DANA FCAL PID FDC START_COUNTER
DIRS += HDGEOMETRY TAGGER TOF TRIGGER UPV TRACKING HDDM

# N. B: the two blank lines between the define and endef are crucial
define newline


endef

.PHONY: all install clean pristine relink

all:
	$(foreach dir,$(DIRS),$(MAKE) -C $(dir)$(newline))

install:
	$(foreach dir,$(DIRS),$(MAKE) -C $(dir) install$(newline))

clean:
	$(foreach dir,$(DIRS),$(MAKE) -C $(dir) clean$(newline))

pristine:
	$(foreach dir,$(DIRS),$(MAKE) -C $(dir) pristine$(newline))

relink:
	$(foreach dir,$(DIRS),$(MAKE) -C $(dir) relink$(newline))

