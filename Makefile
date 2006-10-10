
DIRS += external libraries programs

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

