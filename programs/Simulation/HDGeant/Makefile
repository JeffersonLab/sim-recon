OStype = $(shell uname)
BINDIR = ../bin.$(OStype)

NETLIB := -lnsl
XLIBS  := -L/usr/X11R6/lib -lXpm -lSM -lXm -lXt -lICE -lXext -lX11 -lXp
 
COPTS = -g
FOPTS = -g -Wno-globals

HDDS = ../hdds

AUXS = hddsGeant3.o goptimize.o

GXFIXED = gxcs.o gxphys.o

hdgeant++.x:	gxint.o uginit.o uglast.o $(AUXS) $(GXFIXED)
	gxint -g 321 -m -d Motif -o $@ $^

hdgeant.x:	hdgeant.F uginit.F uglast.F $(AUXS)
	f77 $(FOPTS) -o $@ -I/cern/pro/src/geant321 $^ \
	    -L${CERN_ROOT}/lib -lgeant321\
	    -lpawlib -lgraflib -lgrafX11 -lpacklib -lmathlib \
	    $(NETLIB) $(XLIBS)
.f.o:
	f77 $(FOPTS) -I/cern/pro/src/geant321 -c -o $@ $^

.F.o:
	f77 $(FOPTS) -I/cern/pro/src/geant321 -DCERNLIB_MOTIF -c -o $@ $^

hddsGeant3.f:
	$(HDDS)/hdds-geant $(HDDS)/main_HDDS.xml >hddsGeant3.f

clean:
	rm -f *.o core hddsGeant3.f last.kumac* paw.metafile
