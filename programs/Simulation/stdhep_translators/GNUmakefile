#
# makefile for stdhep translators
#
UNAME = $(shell uname)
ifeq "$(UNAME)" "AIX"
    CC=cc
    CPP = /lib/cpp
    AR=ar
    F77=xlf
    CFLAGS=  -I.
    FFLAGS=-qextname
    KPPFLAGS= -C -P -DAIX -DUNIX
endif
ifeq "$(UNAME)" "IRIX"
    IRIX6 = $(findstring 6,$(shell uname -r))
    CC=cc
    CPP = /lib/cpp
    AR=ar
    F77=f77
    KPPFLAGS= -C -P -DIRIX -DUNIX
    ifeq "$(IRIX6)" "6"
      CFLAGS= -n32 -I$(STDHEP_DIR)/src/inc
      FFLAGS= -n32
    else
      CFLAGS= 
      FFLAGS=
    endif
endif
ifeq "$(UNAME)" "IRIX64"
    CC=cc
    CPP = /lib/cpp
    AR=ar
    F77=f77
    CFLAGS= -n32 
    FFLAGS= -n32
    KPPFLAGS= -C -P -DIRIX -DUNIX
endif
ifeq "$(UNAME)" "Linux"
    CC=gcc
    CPP=gcc -E
    F77=g77
    AR=ar
    CFLAGS= 
    FFLAGS= -fno-second-underscore
    KPPFLAGS= -C -P -DLinux -DUNIX
endif
ifeq "$(UNAME)" "OSF1"
    CC=cc
    CPP = cpp
    F77=f77
    AR=ar
    CFLAGS=  -g
    FFLAGS=
    KPPFLAGS= -C -P -DOSF1 -DUNIX
endif
ifeq "$(UNAME)" "SunOS"
    CC=/opt/SUNWspro/bin/cc
    CPP = /usr/ccs/lib/cpp
    F77=f77
    AR=ar
    CFLAGS= 
    FFLAGS=
    KPPFLAGS= -C -P -DSunOS -DUNIX
    LOCAL= -lnsl -lsocket
endif
ifeq "$(UNAME)" "HP-UX"
    AR=ar
    CC=cc
    CPP = /lib/cpp
    F77=f77
    CFLAGS=  -Aa -D_HPUX_SOURCE
    FFLAGS=+ppu
    KPPFLAGS= -C -P -DUNIX
endif

INCS = -I$(STDHEP_DIR)/src/inc -I$(HALLD_HOME)/source/include

#------------------------------------------

all:  ascii2stdhep  stdhep2ascii 

stdhep2itape: stdhep2itape.o  $(STDHEP_DIR)/lib/libstdhepC.a
	$(CC) $(CFLAGS) $(INCS)  -o stdhep2itape stdhep2itape.o   \
		$(HALLDLIBS)/libdata.a \
		$(STDHEP_DIR)/lib/libstdhepC.a \
		$(STDHEP_DIR)/lib/libFmcfio.a $(LOCAL) -lm

itape2stdhep: itape2stdhep.o  $(STDHEP_DIR)/lib/libstdhepC.a
	$(CC) $(CFLAGS) $(INCS)  -o itape2stdhep itape2stdhep.o   \
		$(HALLDLIBS)/libdata.a \
		$(STDHEP_DIR)/lib/libstdhepC.a \
		$(STDHEP_DIR)/lib/libFmcfio.a $(LOCAL) -lm


ascii2stdhep: ascii2stdhep.o  $(STDHEP_DIR)/lib/libstdhepC.a
	$(CC) $(CFLAGS)  -o ascii2stdhep ascii2stdhep.o  \
		$(STDHEP_DIR)/lib/libstdhepC.a \
		$(STDHEP_DIR)/lib/libFmcfio.a $(LOCAL) -lm

stdhep2ascii: stdhep2ascii.o  $(STDHEP_DIR)/lib/libstdhepC.a
	$(CC) $(CFLAGS)  -o stdhep2ascii stdhep2ascii.o  \
		$(STDHEP_DIR)/lib/libstdhepC.a \
		$(STDHEP_DIR)/lib/libFmcfio.a $(LOCAL) -lm

install: all
	mv ascii2stdhep   $(HALLD_HOME)/bin/$(UNAME)/.  
	mv stdhep2ascii    $(HALLD_HOME)/bin/$(UNAME)/.  

clean:
	rm -f *.o

realclean:
	rm -f stdtst stdtsth stdtsti stdtstj stdtstp stdtstq stdlistj
	rm -f stdtstio stdtstiom stdtstiox stdtstxin stdtstxout
	rm -f testpdg stdtst_dpm
	rm -f *.bak *.o

#------------------------------------------

.SUFFIXES : 
.SUFFIXES : .o .c .F .inc .h


.F.o:
	$(CPP) $(KPPFLAGS) $(INCS) $< > $*.f 
	$(F77) $(FFLAGS) -c $*.f
	rm -f $*.f

.c.o:
	$(CC) $(CFLAGS) $(INCS) -c $*.c