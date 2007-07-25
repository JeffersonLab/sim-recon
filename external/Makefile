

DIRS += xstream/src

.PHONY: all install clean pristine relink

install:
	$(foreach dir,$(DIRS),make -C $(dir) install;)

all:
	$(foreach dir,$(DIRS),make -C $(dir);)

clean:
	$(foreach dir,$(DIRS),make -C $(dir) clean;)

pristine:
	$(foreach dir,$(DIRS),make -C $(dir) pristine;)

relink:
	$(foreach dir,$(DIRS),make -C $(dir) relink;)

