OStype = $(shell uname)
ARCHtype = $(shell uname -m)
BINDIR = ../bin.$(OStype)

#--------------------------------------------------------------------
# The symbol SPECIAL is used to active any of the following switches
# to enable custom versions of the simulation program.
#NTUPLE_OPTIONS = -DBACKGROUND_STUDIES
#NTUPLE_OPTIONS = -DBACKGROUND_PROFILING
#NTUPLE_OPTIONS = -DWERNERS_VTX_NTUPLE
NTUPLE_OPTIONS = -DCERENKOV_PID_NTUPLE
# Enable the following line to turn off writing out hddm data 
#IO_OPTIONS = -DDISABLE_OUTPUT
#--------------------------------------------------------------------

SWITCHES = $(IO_OPTIONS) $(NTUPLE_OPTIONS)

ifeq ($(OStype),Linux)
	ifeq ($(ARCHtype),alpha)
		CC	:= gcc
		CPP	:= g++
		F77	:= g77
		NETLIB := -lnsl
		XLIBS  := -L/usr/X11R6/lib -lXpm -lSM -lXm -lXt -lICE -lXext -lX11 -lXp -ldl
		COPTS = -g $(SWITCHES)
		FOPTS = -g -Wno-globals $(SWITCHES)
		GLIBS 	:= -L/usr/lib/gcc-lib/alpha-redhat-linux/egcs-2.91.66/ -lg2c 
	else
		CC	:= gcc
		CPP	:= g++
		F77	:= g77
		NETLIB := -lnsl
		XLIBS  := -L/usr/X11R6/lib -lXpm -lSM -lXm -lXt -lICE -lXext -lX11 -lXp -ldl -lpthread
		COPTS = -g $(SWITCHES)
		FOPTS = -g -Wno-globals $(SWITCHES)
		GLIBS	:= 
		LINUXEXPORTDYNAMIC = -Wl,-export-dynamic
		LINUXSTATIC = -Wl,-Bstatic
	endif
endif
ifeq ($(OStype),OSF1)
	CC	:= cc
	CPP	:= g++
	F77	:= f77
	NETLIB := 
	STATIC	:=
	XLIBS  := -L/usr/lib -lXm -lSM -lICE -lXt -lX11 -lm -lPW -ldnet_stub
	COPTS = -g -D_Tru64 $(SWITCHES)
	FOPTS = -g -fpe4 $(SWITCHES)
	LOPTS = -g -non_shared -fpe4
	GLIBS 	:= -L/r5d1/applications/gcc3/lib/gcc-lib/alphaev5-dec-osf4.0f/ -lg2c -L/r5d1/applications/gcc3/lib/gcc-lib/alphaev5-dec-osf4.0f/3.3/ -lgcc
endif

HDDS = ../hdds
HDDM = ../hddm

ifeq ($(OStype),OSF1)
AUXF = gufld.o guhadr.o gukine.o guout.o guphad.o gustep.o guxcs.o \
       gelhad/gtgama.o gelhad/caspim.o gelhad/caspip.o gpairg.o \
       hddsGeant3.o goptimize.o cobrems.o vunit.o gsrotm.o fint.o gthion.o
else
AUXF = gufld.o guhadr.o gukine.o guout.o guphad.o gustep.o guxcs.o \
       gelhad/gtgama.o gelhad/caspim.o gelhad/caspip.o gpairg.o \
       hddsGeant3.o goptimize.o cobrems.o
endif

AUXC = hddmInput.o hddmOutput.o hddm_s.o \
       hitStart.o hitCDC.o hitFDC.o hitBCal.o hitCerenkov.o hitFTOF.o \
       hitLGD.o bintree.o hitutil.o memcheck.o timel.o trapfpe.o

GXFIXED = gxcs.o gxphys.o
GELHAD = gelhad/libgelhad.a

CERNLIB = $(CERN_ROOT)/bin/cernlib

hdgeant++: gxint.o uginit.o uglast.o $(GELHAD) $(AUXF) $(AUXC) $(GXFIXED)
	$(CERN_ROOT)/bin/gxint -v $(CERN_LEVEL) -g 321 -m -d Motif \
	-o $@ gxint.o uginit.o uglast.o \
	$(AUXF) $(AUXC) $(GXFIXED) $(STATIC) \
	-Lgelhad -lgelhad $(XLIBS) $(GLIBS) \
	-- $(LINUXEXPORTDYNAMIC)

hdgeant: hdgeant.o uginit.o uglast.o $(GELHAD) $(AUXF) $(AUXC)
	$(F77) $(LOPTS) -o $@ -I$(CERN_ROOT)/include $(LINUXSTATIC) \
	 hdgeant.o uginit.o uglast.o $(AUXF) $(AUXC) -Lgelhad -lgelhad \
	`$(CERNLIB) geant321 mathlib graflib grafX11` \
	$(NETLIB) $(XLIBS)

install: hdgeant hdgeant++
	cp $^ $(HALLD_HOME)/bin/$(OStype)

.F.o:
	$(F77) $(FOPTS) -I$(CERN_ROOT)/include -DCERNLIB_MOTIF -c -o $@ $<

.f.o:
	$(F77) $(FOPTS) -I$(CERN_ROOT)/include -c -o $@ $<

.c.o:
	$(CC) $(COPTS) -I. -I../include -I$(CERN_ROOT)/include -c -o $@ $<

$(AUXC): hddm_s.h bintree.h geant3.h hddmOutput.h 

.cpp.o:
	$(CPP) $(COPTS) -I. -I../include -c -o $@ $<

hddsGeant3.f:
	$(HDDS)/hdds-geant $(HDDS)/main_HDDS.xml >hddsGeant3.f

hddm_s.h:
	$(HDDM)/hddm-c $(HDDM)/event.xml

hddm_s.c:
	$(HDDM)/hddm-c $(HDDM)/event.xml

hddm_s.o: 	hddm_s.c
	$(CC) $(COPTS) -I. -I../include -c -o $@ $<

$(GELHAD):
	$(MAKE) -C gelhad

clean:
	rm -f *.o core last.kumac* paw.metafile
