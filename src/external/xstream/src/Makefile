AR = ar
RANLIB = ranlib

OBJDIR = .obj/$(BMS_OSNAME)
LIBDIR = .lib/$(BMS_OSNAME)
INCLUDES = -I../include
EXTRA_DIST = debug.h md5_t.pl md5_t.h
libxstream_a_SOURCES = debug.cpp common.cpp xdr.cpp base64.cpp \
	digest.cpp tee.cpp md5.cpp bz.cpp dater.cpp fd.cpp \
	posix.cpp z_digest.cpp z.cpp
libxstream_a_OBJECTS_0 = debug.o common.o xdr.o base64.o \
	digest.o tee.o md5.o bz.o dater.o fd.o \
	posix.o z_digest.o z.o
libxstream_a_OBJECTS = $(addprefix $(OBJDIR)/,$(libxstream_a_OBJECTS_0))
lib_LIBRARIES = libxstream.a

ifdef HALLD_MY
INSTALL_DIR = $(HALLD_MY)
else
INSTALL_DIR = $(HALLD_HOME)
endif

all: mkdirs $(LIBDIR)/$(lib_LIBRARIES)

.SUFFIXES:
.SUFFIXES: .cpp .o .obj

mkdirs:
	mkdir -p $(OBJDIR) $(LIBDIR)

clean-libLIBRARIES:
	-test -z "$(lib_LIBRARIES)" || rm -f $(lib_LIBRARIES)
$(LIBDIR)/libxstream.a: $(libxstream_a_OBJECTS) $(libxstream_a_DEPENDENCIES) 
	-rm -f $@
	$(AR) cru $@ $(libxstream_a_OBJECTS) $(libxstream_a_LIBADD)
	$(RANLIB) $@

install: all
	mkdir -p $(INSTALL_DIR)/lib/$(BMS_OSNAME)
	cp -p $(LIBDIR)/$(lib_LIBRARIES) $(INSTALL_DIR)/lib/$(BMS_OSNAME)
	cp -pr ../include ./tmpinc
	rm -rf `find ./tmpinc -name ".svn"`
	mkdir -p $(INSTALL_DIR)/include
	cp -pr tmpinc/* $(INSTALL_DIR)/include
	rm -rf tmpinc

mostlyclean-compile:
	-rm -f *.$(OBJEXT)

distclean-compile:
	-rm -f *.tab.c

pristine: clean

depclean:
	@echo "Nothing to do for make depclean in xstream"

clean:
	rm -rf $(OBJDIR) $(LIBDIR)

$(OBJDIR)/%.o: %.cpp
	g++ $(CXXFLAGS) $(INCLUDES) -c -o $@ $^
