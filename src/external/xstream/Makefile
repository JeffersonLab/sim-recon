#
# Makefile.0 example
#
.DELETE_ON_ERROR :
.PHONY : all clean pristine

SHELL = /bin/sh

#
# get the list of subdirectories for the current directory only
#
#subdirs := $(shell find . -name CVS -prune -o -name include -prune -o \
			  -maxdepth 1 -type d -print \
                          | grep -v '\.$$' | sed 's/\.\///')

# This was modified to only build in the "src" directory. I'm
# not sure what the business of the the above is all about
# since only 2 of the subdirectories has a makefile!
# Feb 15, 2008 D.L.
subdirs := src

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

pristine:clean
	

clean:
	rm -f *~ .*~
	for a in $(subdirs); do \
		$(MAKE) -C $$a $@; \
	done
