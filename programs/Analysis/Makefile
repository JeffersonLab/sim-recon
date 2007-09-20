

DIRS += hddm hd_ana hd_dump hd_root mctrk_ana plugins hdview

# N. B: the two blank lines between the define and endef are crucial
define newline


endef

.PHONY: all install clean pristine relink

all:
	$(foreach dir,$(DIRS),make -C $(dir)$(newline))

install:
	$(foreach dir,$(DIRS),make -C $(dir) install$(newline))

clean:
	$(foreach dir,$(DIRS),make -C $(dir) clean$(newline))

pristine:
	$(foreach dir,$(DIRS),make -C $(dir) pristine$(newline))

relink:
	$(foreach dir,$(DIRS),make -C $(dir) relink$(newline))

