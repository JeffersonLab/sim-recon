OStype = $(shell uname)
ARCHtype = $(shell uname -m)
BINDIR = ../bin.$(OStype)

ifeq ($(OStype),Linux)
	ifeq ($(ARCHtype),alpha)
		CC	:= gcc
		CPP	:= g++
		F77	:= g77
		NETLIB := -lnsl
		XLIBS  := -L/usr/X11R6/lib -lXpm -lSM -lXm -lXt -lICE -lXext -lX11 -lXp
		COPTS = -g
		FOPTS = -g -Wno-globals
		GLIBS 	:= -L/usr/lib/gcc-lib/alpha-redhat-linux/egcs-2.91.66/ -lg2c 
	else
		CC	:= gcc
		CPP	:= g++
		F77	:= g77
		NETLIB := -lnsl
		XLIBS  := -L/usr/X11R6/lib -lXpm -lSM -lXm -lXt -lICE -lXext -lX11 -lXp
		COPTS = -g
		FOPTS = -g -Wno-globals
		GLIBS	:= 
	endif
endif
ifeq ($(OStype),OSF1)
	CC	:= cc
	CPP	:= g++
	F77	:= f77
	NETLIB := 
	STATIC	:=
	XLIBS  := -L/usr/lib -lXm -lSM -lICE -lXt -lX11 -lm -lPW -ldnet_stub
	COPTS = -g -D_Tru64
	FOPTS = -g -fpe4
	LOPTS = -g -non_shared -fpe4
	GLIBS 	:= -L/r5da/applications/gcc/lib/gcc-lib/alphaev5-dec-osf4.0f/2.95.3 -lg2c -lgcc
endif

HDDS = ../hdds
HDDM = ../hddm

AUXF = gufld.o guhadr.o gukine.o guout.o guphad.o gustep.o gdstep.o guxcs.o \
       hddsGeant3.o goptimize.o cobrems.o

AUXC = hddmInput.o hddmOutput.o hddm_s.o \
       hitStart.o hitCDC.o hitFDC.o hitBCal.o hitCerenkov.o hitFTOF.o \
       hitLGD.o bintree.o hitutil.o memcheck.o timel.o trapfpe.o

GXFIXED = gxcs.o gxphys.o

hdgeant++:	gxint.o uginit.o uglast.o $(AUXF) $(AUXC) $(GXFIXED)
	$(CERN_ROOT)/bin/gxint -g 321 -m -d Motif -o $@ $^ $(STATIC) $(XLIBS) $(GLIBS)

hdgeant:	hdgeant.o uginit.o uglast.o $(AUXF) $(AUXC)
	$(F77) $(LOPTS) -o $@ -I$(CERN_ROOT)/include $^ \
	    -L${CERN_ROOT}/lib -lgeant321\
	    -lpawlib -lgraflib -lgrafX11 -lpacklib -lmathlib \
	    $(NETLIB) $(XLIBS) 

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

clean:
	rm -f *.o core last.kumac* paw.metafile
