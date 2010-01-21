ifndef OSNAME
OSNAME = $(shell uname)-$(shell uname -p)
endif

AR = ar
RANLIB = ranlib

top_srcdir = $(HALLD_HOME)/src
INCLUDES = -I$(top_srcdir)/include
EXTRA_DIST = debug.h md5_t.pl md5_t.h
libxstream_a_SOURCES = debug.cpp common.cpp xdr.cpp base64.cpp \
	digest.cpp tee.cpp md5.cpp bz.cpp dater.cpp fd.cpp \
	posix.cpp z_digest.cpp
libxstream_a_OBJECTS = debug.o common.o xdr.o base64.o \
	digest.o tee.o md5.o bz.o dater.o fd.o \
	posix.o z_digest.o
lib_LIBRARIES = libxstream.a

all: libxstream.a

.SUFFIXES:
.SUFFIXES: .cpp .o .obj

clean-libLIBRARIES:
	-test -z "$(lib_LIBRARIES)" || rm -f $(lib_LIBRARIES)
libxstream.a: $(libxstream_a_OBJECTS) $(libxstream_a_DEPENDENCIES) 
	-rm -f libxstream.a
	$(AR) cru libxstream.a $(libxstream_a_OBJECTS) $(libxstream_a_LIBADD)
	$(RANLIB) libxstream.a

install: $(lib_LIBRARIES)
	mkdir -p $(HALLD_HOME)/lib/$(OSNAME)
	cp $^ $(HALLD_HOME)/lib/$(OSNAME)
	cp -r ../include ./tmpinc
	rm -rf `find ./tmpinc -name ".svn"`
	mkdir -p $(HALLD_HOME)/include
	cp -r tmpinc/* $(HALLD_HOME)/include
	rm -rf tmpinc

mostlyclean-compile:
	-rm -f *.$(OBJEXT)

distclean-compile:
	-rm -f *.tab.c

pristine: clean

depclean:
	@echo "Nothing to do for make depclean in xstream"

clean:
	rm -f *.o *.a $(HALLD_HOME)/lib/$(OSNAME)/$(lib_LIBRARIES)

%.o: %.cpp
	g++ $(INCLUDES) -c -o $@ $^
