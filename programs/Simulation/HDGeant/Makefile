OStype = $(shell uname)
BINDIR = ../bin.$(OStype)

NETLIB := -lnsl
XLIBS  := -L/usr/X11R6/lib -lXpm -lSM -lXm -lXt -lICE -lXext -lX11 -lXp
 
COPTS = -g
FOPTS = -g -Wno-globals -fdebug-kludge

HDDS = ../hdds
HDDM = ../hddm

AUXF = gufld.o guhadr.o gukine.o guout.o guphad.o gustep.o gdstep.o guxcs.o \
       hddsGeant3.o goptimize.o cobrems.o

AUXC = hddmInput.o hddmOutput.o hddm_s.o \
       hitStart.o hitCDC.o hitFDC.o hitBCal.o hitCerenkov.o hitFTOF.o \
       hitLGD.o bintree.o hitutil.o memcheck.o

GXFIXED = gxcs.o gxphys.o

hdgeant++:	gxint.o uginit.o uglast.o $(AUXF) $(AUXC) $(GXFIXED)
	$(CERN_ROOT)/bin/gxint -g 321 -m -d Motif -o $@ $^

hdgeant:	hdgeant.o uginit.o uglast.o $(AUXF) $(AUXC)
	f77 $(FOPTS) -o $@ -I$(CERN_ROOT)/include $^ \
	    -L${CERN_ROOT}/lib -lgeant321\
	    -lpawlib -lgraflib -lgrafX11 -lpacklib -lmathlib \
	    $(NETLIB) $(XLIBS)

.F.o:
	f77 $(FOPTS) -I$(CERN_ROOT)/include -DCERNLIB_MOTIF -c -o $@ $<

.f.o:
	f77 $(FOPTS) -I$(CERN_ROOT)/include -c -o $@ $<

.c.o:
	gcc $(COPTS) -I. -I../include -c -o $@ $<

$(AUXC): hddm_s.h bintree.h geant3.h hddmOutput.h

.cpp.o:
	g++ $(COPTS) -I. -I../include -c -o $@ $<

hddsGeant3.f:
	$(HDDS)/hdds-geant $(HDDS)/main_HDDS.xml >hddsGeant3.f

hddm_s.h:
	$(HDDM)/hddm-c $(HDDM)/event.xml

clean:
	rm -f *.o core hddm_s.[ch] hddsGeant3.f last.kumac* paw.metafile
