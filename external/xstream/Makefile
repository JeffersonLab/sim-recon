#
# Makefile.0 example
#
.DELETE_ON_ERROR :
.PHONY : all clean 

SHELL = /bin/sh

#
# get the list of subdirectories for the current directory only
#
subdirs := $(shell find . -name CVS -prune -o -name include -prune -o \
			  -maxdepth 1 -type d -print \
                          | grep -v '\.$$' | sed 's/\.\///')


submakefile  := $(sublevel)/Makefile

#
# rebuild in each subdirectory
#
all: $(subdirs) 

.PHONY: $(subdirs)
${subdirs}:
	if [ -r $@/Makefile ]; then \
		$(MAKE) -C $@; \
	fi

install:
	for a in $(subdirs); do \
		$(MAKE) -C $$a $@; \
	done

clean:
	rm -f *~ .*~
	for a in $(subdirs); do \
		$(MAKE) -C $$a $@; \
	done
