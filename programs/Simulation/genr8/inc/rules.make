#
# Rules for Makefiles
#


#
# List of supported targets-
#
# irix4		(obsolete and untested for ages)
# irix5		(either 5.3 or 6.2)
# irix5o
# irix5mips1o   (only 5.3: 6.2 does not have the -mips1 libraries)
# irix5mips2o   (either 5.3 or 6.2)
# irix6		(irix6 targets are really irix6.1 on R8000 machines)
# irix6o
# irix6n32
# irix6n32o
# irix6n32mips3o (SGI IRIX 6.2, in the -mips3 -n32 mode, to get the most bang out of the R4000)
# irix6n32mips4r10ko (SGI IRIX 6.4 on Origin 200, optimized for the R10k)
# sunos5	(SUN SPARC running Solaris 5.4)
# sunos5o
# sunos5ss20o   SPARC-20 fully optimized, Sun compilers version "cc -V: SC4.0 18 Oct 1995 C 4.0"
# sunos5ultrao  same, but for the Ultra.
# sunos5ultra2o target is SUN Ultra2 (specifically, sun3.bnl.gov).
# sunos5p6o	Intel P6 (Pentium Pro) running Solaris
# rs6000	(IBM RS6000 running AIX3.2.5)
# rs6000o
# rs6000ppc601o	(same but built for the PowerPC 601 processor)
# rs6000p2o	(same but built for the POWER2 processor)
# aix4		(IBM RS6000 running AIX4.1)
# aix4o
# aix4ppc601o	(same but built for the PowerPC 601 processor)
# aix4p2o	(same but built for the POWER2 processor)
# hpux		(HP running some or other HPUX)
# hpuxo
# decunix	(DEC Unix V3.2)
# decunixo	(DEC Unix V3.2)
# convexp	Convex Exemplar running SPP-UX 1.0/HPUX A.09.03, instrumented for profiling
# linux		linux with gcc/g77
# linuxo	linux with gcc/g77, optimized
#

ifndef E852_MAKETARGET
E852_MAKETARGET=irix5
endif

PWD := $(shell pwd)

#
# XLIBS defines the location of the X11 libraries
#

XLIBS_irix4 =		-lFm_c -lXm_s -lXt_s -lX11_s -lsun -lmalloc -lPW
XLIBS_rs6000 =		-lFm_c -lXm -lXt -lX11 -lPW
XLIBS_rs6000o		:= $(XLIBS_rs6000)
XLIBS_rs6000ppc601o	:= $(XLIBS_rs6000)
XLIBS_rs6000p2o		:= $(XLIBS_rs6000)
XLIBS_aix4		:= $(XLIBS_rs6000)
XLIBS_aix4o		:= $(XLIBS_rs6000)
XLIBS_aix4ppc601o	:= $(XLIBS_rs6000)
XLIBS_aix4p2o		:= $(XLIBS_rs6000)
XLIBS_irix5		:= -lFm_c libXm.so libXt.so libXext.so libX11.so libPW.so
XLIBS_irix5o		:= $(XLIBS_irix5)
XLIBS_irix5mips1o	:= $(XLIBS_irix5)
XLIBS_irix5mips2o	:= $(XLIBS_irix5)
XLIBS_irix6 =		-lFm_c -lXm_s -lXt_s -lXext -lX11_s -lmalloc -lPW
XLIBS_irix6o =		-lFm_c -lXm_s -lXt_s -lXext -lX11_s -lmalloc -lPW
XLIBS_irix6n32  =	$(XLIBS_irix6)
XLIBS_irix6n32o =	$(XLIBS_irix6o)
XLIBS_sunos5 =		-lFm_c /usr/dt/lib/libXm.so libXt.so libX11.so libw.so
XLIBS_sunos5o		:= $(XLIBS_sunos5)
XLIBS_sunos5ss20o	:= $(XLIBS_sunos5)
XLIBS_sunos5ultrao	:= $(XLIBS_sunos5)
XLIBS_sunos5ultra2o	:= $(XLIBS_sunos5)
XLIBS_sunos5p6o		:= $(XLIBS_sunos5)
XLIBS_decunix=		-lFm_c -lXm -lXt -lX11 -lSM -lICE -lPW
XLIBS_linux=		-lFm_c libXm.so.1 libXt.so libX11.so libICE.so libSM.so libXext.so libXpm.so
XLIBS_linuxo=		$(XLIBS_linux)

XLIBS:= $(XLIBS_$(E852_MAKETARGET))

PEXLIBS= -lXmu -lPEX5

#
# INCLUDE defines the location of the C and C++ include files
#

INCLUDE_rs6000 =	-I. -I../include -I$(E852_HOME)/include -I$(E852_HOME)/include/CC -DFUNCPROTO=15  -DDONOT_USEDAQ
INCLUDE_rs6000o		:= $(INCLUDE_rs6000)
INCLUDE_rs6000ppc601o	:= $(INCLUDE_rs6000)
INCLUDE_rs6000p2o	:= $(INCLUDE_rs6000)
INCLUDE_aix4		:= $(INCLUDE_rs6000)
INCLUDE_aix4o		:= $(INCLUDE_rs6000)
INCLUDE_aix4ppc601o	:= $(INCLUDE_rs6000)
INCLUDE_aix4p2o		:= $(INCLUDE_rs6000)
INCLUDE_irix4 =		-I. -I../include -I$(E852_HOME)/include -I$(E852_HOME)/include/CC -DFUNCPROTO=15 -DDONOT_USEDAQ
INCLUDE_irix5		:= $(INCLUDE_irix4)
INCLUDE_irix5o		:= $(INCLUDE_irix4)
INCLUDE_irix5mips1o	:= $(INCLUDE_irix4)
INCLUDE_irix5mips2o	:= $(INCLUDE_irix4)
INCLUDE_irix6		:= $(INCLUDE_irix4)
INCLUDE_irix6o		:= $(INCLUDE_irix4)
INCLUDE_irix6n32	:= $(INCLUDE_irix4)
INCLUDE_irix6n32o	:= $(INCLUDE_irix4)
INCLUDE_irix6n32mips3o	:= $(INCLUDE_irix4)
INCLUDE_irix6n32mips4r10ko	:= $(INCLUDE_irix4)
INCLUDE_cray		:= $(INCLUDE_irix5)
#INCLUDE_hpux		:= $(INCLUDE_irix5) -DHPUX -D_INCLUDE_XOPEN_SOURCE -D_INCLUDE_POSIX_SOURCE  -DDONOT_USEDAQ
INCLUDE_hpux		:= $(INCLUDE_irix5) -DHPUX -DDONOT_USEDAQ
INCLUDE_hpuxo		:= $(INCLUDE_hpux)
INCLUDE_convexp		:= $(INCLUDE_hpux)
INCLUDE_sunos5		:= -I/usr/openwin/include -I/usr/dt/include $(INCLUDE_irix5)  -DDONOT_USEDAQ
INCLUDE_sunos5o		:= $(INCLUDE_sunos5)
INCLUDE_sunos5ss20o	:= $(INCLUDE_sunos5)
INCLUDE_sunos5ultrao	:= $(INCLUDE_sunos5)
INCLUDE_sunos5ultra2o	:= $(INCLUDE_sunos5)
INCLUDE_sunos5p6o	:= $(INCLUDE_sunos5)
INCLUDE_decunix		:= $(INCLUDE_rs6000) -D_XOPEN_SOURCE -D_OSF_SOURCE -DDECUNIX
INCLUDE_decunixo	:= $(INCLUDE_decunix)
INCLUDE_linux		:= $(INCLUDE_irix4)  -I/usr/X11R6/include
INCLUDE_linuxo		:= $(INCLUDE_irix4)  -I/usr/X11R6/include


INCLUDE:= $(INCLUDE_$(E852_MAKETARGET))

#
# FINCLUDE defines the location of F77 libraries
#    *and* any special f77 switches (such as -n32, -mips2, etc...)
#

FINCLUDE_rs6000= 	-qextname -I. -I../include -I$(E852_HOME)/include
FINCLUDE_rs6000o	:= $(FINCLUDE_rs6000)
FINCLUDE_rs6000ppc601o	:= $(FINCLUDE_rs6000)
FINCLUDE_rs6000p2o	:= $(FINCLUDE_rs6000)
FINCLUDE_aix4		:= $(FINCLUDE_rs6000)
FINCLUDE_aix4o		:= $(FINCLUDE_rs6000)
FINCLUDE_aix4ppc601o	:= $(FINCLUDE_rs6000)
FINCLUDE_aix4p2o	:= $(FINCLUDE_rs6000)
FINCLUDE_irix4=		-static -Wf,-I. -Wf,-I../include -Wf,-I$(E852_HOME)/include $(INCLUDE)
FINCLUDE_irix5=		$(FINCLUDE_irix4)
FINCLUDE_irix5o=	$(FINCLUDE_irix4)
FINCLUDE_irix5mips1o=	$(FINCLUDE_irix4)
FINCLUDE_irix5mips2o=	$(FINCLUDE_irix4)
FINCLUDE_irix6=		$(FINCLUDE_irix4)
FINCLUDE_irix6o=	$(FINCLUDE_irix4)
FINCLUDE_irix6n32=	$(FINCLUDE_irix4)
FINCLUDE_irix6n32o=	$(FINCLUDE_irix4)
FINCLUDE_irix6n32mips3o=	$(FINCLUDE_irix4)
FINCLUDE_irix6n32mips4r10ko=	$(FINCLUDE_irix4)
FINCLUDE_cray=		-I. -I../include -I$(E852_HOME)/include
FINCLUDE_hpux=		+ppu -C -K -I. -I../include -I$(E852_HOME)/include
FINCLUDE_hpuxo:=	$(FINCLUDE_hpux)
FINCLUDE_sunos5=	-e -I. -I../include -I$(E852_HOME)/include
FINCLUDE_sunos5o=	$(FINCLUDE_sunos5)
FINCLUDE_sunos5ss20o=	$(FINCLUDE_sunos5)
FINCLUDE_sunos5ultrao=	$(FINCLUDE_sunos5)
FINCLUDE_sunos5ultra2o=	$(FINCLUDE_sunos5)
FINCLUDE_sunos5p6o=	$(FINCLUDE_sunos5)
FINCLUDE_decunix:=	-fpe1 $(FINCLUDE_irix5)
FINCLUDE_decunixo:=	-fpe1 $(FINCLUDE_irix5)
FINCLUDE_linux:=	-fno-automatic -fno-second-underscore $(INCLUDE)
FINCLUDE_linuxo:=	-fno-automatic -fno-second-underscore $(INCLUDE)

FINCLUDE:= $(FINCLUDE_$(E852_MAKETARGET))

#
# Here: COPT- defines the flags used to build the %.oo objects (both C and f77)
#       CDBG- defines the flags used to build the %.o objects,
#			  note that for the optimised targets (such as irix5o)
#			  CDBG is explicitely set equal to COPT
#

COPT_default=		-O
COPT_irix5=		-g3 -O2 -Olimit 20000
COPT_irix5o=		-g3 -O2 -Olimit 20000
COPT_irix5mips1o=	-g3 -O2 -Olimit 20000
COPT_irix5mips2o=	-g3 -O2 -Olimit 20000
COPT_irix6=		-O3
COPT_irix6o=		-O3
COPT_irix6n32=		-O3
COPT_irix6n32o=		-O3
COPT_irix6n32mips3o=	-g3 -O3
COPT_irix6n32mips4r10ko=	-g3 -O3
COPT_rs6000=		-O3
COPT_rs6000o=		-O3
COPT_rs6000ppc601o=	-O3 -qarch=ppc  -qtune=601
COPT_rs6000p2o=		-O3 -qarch=pwr2 -qtune=pwr2
COPT_aix4=		-O3
COPT_aix4o=		-O3
COPT_aix4ppc601o=	-O3 -qarch=ppc  -qtune=601
COPT_aix4p2o=		-O3 -qarch=pwr2 -qtune=pwr2
COPT_sunos5=		-O
COPT_sunos5o=		-O
COPT_sunos5ss20o=	-fast -fsimple=2 -xlibmieee -xlibmil -xO5
COPT_sunos5ultrao=	-fast -fsimple=2 -xlibmieee -xlibmil -xO5
COPT_sunos5ultra2o=	-fast -fsimple=2 -xlibmieee -xlibmil -xO5
COPT_sunos5p6o=		-fast -fsimple=2 -xlibmieee -xlibmil -xO5
COPT_decunix=		-g3 -O2
COPT_decunixo=		-g3 -O2
COPT_hpuxo=		-O
COPT_linux=		-g -O2 -m486
COPT_linuxo=		-g -O2 -m486

CDBG_default=		-g
CDBG_irix5o=		$(COPT_irix5o)
CDBG_irix5mips1o=	$(COPT_irix5mips1o)
CDBG_irix5mips2o=	$(COPT_irix5mips2o)
CDBG_irix6o=		$(COPT_irix6o)
CDBG_irix6n32o=		$(COPT_irix6n32o)
CDBG_irix6n32mips3o=	$(COPT_irix6n32mips3o)
CDBG_irix6n32mips4r10ko=	$(COPT_irix6n32mips3o)
CDBG_rs6000o=		$(COPT_rs6000o)
CDBG_rs6000ppc601o=	$(COPT_rs6000ppc601o)
CDBG_rs6000p2o=		$(COPT_rs6000p2o)
CDBG_aix4o=		$(COPT_aix4o)
CDBG_aix4ppc601o=	$(COPT_aix4ppc601o)
CDBG_aix4p2o=		$(COPT_aix4p2o)
CDBG_sunos5o=		$(COPT_sunos5o)
CDBG_sunos5ss20o=	$(COPT_sunos5ss20o)
CDBG_sunos5ultrao=	$(COPT_sunos5ultrao)
CDBG_sunos5ultra2o=	$(COPT_sunos5ultrao)
CDBG_sunos5p6o=		$(COPT_sunos5p6o)
CDBG_hpuxo=		$(COPT_hpuxo)
CDBG_decunixo=		$(COPT_decunix)
CDBG_linuxo=		$(COPT_linux)

COPT:=$(COPT_$(E852_MAKETARGET))
CDBG:=$(CDBG_$(E852_MAKETARGET))

ifndef CDBG
CDBG:=$(CDBG_default)
endif
ifndef COPT
COPT:=$(COPT_default)
endif
ifdef Optimizer
CDBG:=$(COPT)
endif

#
# various XFM stuff
#

XFM_irix4=	xfm
#XFM_irix5=	/usr/bin/X11/xfm
XFM_irix5=	/home/lemond/e852/bin/xfm1
XFM_irix5mips1o=	/home/lemond/e852/bin/xfm1
XFM_irix5mips2o=	/home/lemond/e852/bin/xfm1
XFM_rs6000=	/home/lemond/e852/bin/xfm1-fake.rs6000
XFM_rs6000o	:= $(XFM_rs6000)

USEXFM_irix4=		yes
USEXFM_irix5=		yes
USEXFM_irix5o=		yes
USEXFM_irix5mips2o=	yes
USEXFM_irix5mips1o=	yes
USEXFM_irix6=		no
USEXFM_irix6o=		no
USEXFM_irix6n32=	no
USEXFM_irix6n32o=	no
USEXFM_rs6000=		no
USEXFM_rs6000o=		no
USEXFM_rs6000ppc601o=	no
USEXFM_rs6000p2o=	no
USEXFM_aix4=		no
USEXFM_aix4o=		no
USEXFM_aix4ppc601o=	no
USEXFM_aix4p2o=		no
USEXFM_sunos5=		no
USEXFM_sunos5o=		no

USEXFM:= $(USEXFM_$(E852_MAKETARGET))

ifndef USEXFM
USEXFM:=no
endif

ifndef NSLSERVER
USEXFM:=no
endif

XFM:= $(XFM_$(E852_MAKETARGET))

#
# Various architecture-dependant flags, common between C, f77
#

ARCHFLAGS_rs6000ppc601o=	-qtune=601  -qarch=ppc
ARCHFLAGS_rs6000p2o=		-qtune=pwr2 -qarch=pwr2
ARCHFLAGS_aix4ppc601o=		-qtune=601  -qarch=ppc
ARCHFLAGS_aix4p2o=		-qtune=pwr2 -qarch=pwr2
ARCHFLAGS_irix5mips1o=		-mips1
ARCHFLAGS_irix5mips2o=		-mips2
ARCHFLAGS_irix6o=		-WK,-inline -OPT:IEEE_arithmetic=3:roundoff=3
ARCHFLAGS_irix6n32=		-n32 -mips4
ARCHFLAGS_irix6n32o=		-n32 -mips4 -WK,-inline -OPT:IEEE_arithmetic=3:roundoff=3 
ARCHFLAGS_irix6n32mips3o=	-n32 -mips3
ARCHFLAGS_irix6n32mips4r10ko=	-n32 -mips4 -TARG:processor=r10000
ARCHFLAGS_sunos5ss20o=		-xtarget=ss20/61
ARCHFLAGS_sunos5ultrao=		-xtarget=ultra
ARCHFLAGS_sunos5ultra2o=	-xtarget=ultra2 -xcache=16/32/1:1024/64/1
ARCHFLAGS_sunos5p6o=		-xtarget=pentium_pro

ARCHFLAGS:= $(ARCHFLAGS_$(E852_MAKETARGET))

#
# Various architecture-dependant flags, used by GCC
#

GARCHFLAGS_irix5mips1o=		-mips1
GARCHFLAGS_irix5mips2o=		-mips2
GARCHFLAGS_irix6o=		-mips4
GARCHFLAGS_irix6n32=		-mips4
GARCHFLAGS_irix6n32o=		-mips4

GARCHFLAGS:= $(GARCHFLAGS_$(E852_MAKETARGET))

#
# Various target-dependant C compiler flags (such as -n32, -mips2, etc...)
#

CFLAGS_irix4=		-fullwarn -prototypes -woff 269
CFLAGS_irix5=		-fullwarn -prototypes $(E852_ANSI)
CFLAGS_irix5o=		-fullwarn -prototypes $(E852_ANSI)
CFLAGS_irix6=		-fullwarn -woff 1174,1209,1233 $(E852_ANSI)
CFLAGS_irix6o=		-fullwarn -woff 1174,1209,1233 $(E852_ANSI)
CFLAGS_irix6n32=	-fullwarn -Wl,-woff,105 -woff 1174,1209,1233 $(E852_ANSI)
CFLAGS_irix6n32o=	-fullwarn -Wl,-woff,105 -woff 1174,1209,1233 $(E852_ANSI)
CFLAGS_irix6n32mips3o=	-fullwarn -woff 1174,1209 $(E852_ANSI)
CFLAGS_irix6n32mips4r10ko=	-fullwarn -woff 1174,1209 $(E852_ANSI)
CFLAGS_decunix=		-ieee_with_inexact
CFLAGS_decunixo:=	$(CFLAGS_decunix)
CFLAGS_hpux:=		-Ae
CFLAGS_hpuxo:=		-Ae
CFLAGS_linux:=		-Wall
CFLAGS_linuxo:=		-Wall

CFLAGS:= $(CFLAGS_$(E852_MAKETARGET))

#
# C++ stuff
#

CCFLAGS_irix4=		+w
CCFLAGS_irix5=		+w
CCFLAGS_irix5mips1o=	+w -O2 -g3
CCFLAGS_irix5mips2o=	+w -O2 -g3
CCFLAGS_irix6=		
CCFLAGS_irix6o=		-O3
CCFLAGS_irix6n32=	
CCFLAGS_irix6n32o=	-O3
CCFLAGS_rs6000=		
CCFLAGS_rs6000o=	-O3
CCFLAGS_rs6000ppc601o=	-O3
CCFLAGS_rs6000p2o=	-O3
CCFLAGS_aix4=		
CCFLAGS_aix4o=		-O3
CCFLAGS_aix4ppc601o=	-O3
CCFLAGS_aix4p2o=	-O3

CCFLAGS:= $(CCFLAGS_$(E852_MAKETARGET))

#
# The irix6 C++ preprocessor had to be kludged
#

CCPPFLAGS_irix6 =	-32
CCPPFLAGS_irix6o =	-32
CCPPFLAGS_irix6n32 =	-32
CCPPFLAGS_irix6n32o =	-32

CCPPFLAGS = $(CCPPFLAGS_$(E852_MAKETARGET))

#
# Location of the Fortran run-time libraries
#

FLIB_irix4=		-lF77 -lI77 -lU77 -lisam -lm
FLIB_rs6000=		-lxlf90 -lxlf -lm
FLIB_rs6000o=		$(FLIB_rs6000)
FLIB_rs6000ppc601o=	$(FLIB_rs6000)
FLIB_rs6000p2o=		$(FLIB_rs6000)
FLIB_aix4=		-lxlf -lxlf90 -lm
FLIB_aix4o=		$(FLIB_aix4)
FLIB_aix4ppc601o=	$(FLIB_aix4)
FLIB_aix4p2o=		$(FLIB_aix4)
FLIB_irix5=		libftn.so libm.so
FLIB_irix5o=		libftn.so libm.so
FLIB_irix5mips1o=	
FLIB_irix5mips2o=	libftn.so libm.so
FLIB_irix6=		/usr/lib64/libftn.so /usr/lib64/libm.so
FLIB_irix6o=		/usr/lib64/libftn.so /usr/lib64/libm.so
FLIB_irix6n32=		/usr/lib32/libftn.so /usr/lib32/libm.so
FLIB_irix6n32o=		/usr/lib32/libftn.so /usr/lib32/libm.so
FLIB_irix6n32mips3o=	/usr/lib32/mips3/libftn.so
FLIB_irix6n32mips4r10ko=	/usr/lib32/mips4/libftn.so
FLIB_sunos5=		-lV77 -lM77 libF77.so.2
FLIB_sunos5o=		$(FLIB_sunos5)
FLIB_sunos5ss20o:=	libM77.so libF77.so
FLIB_sunos5ultrao:=	$(FLIB_sunos5ss20o)
FLIB_sunos5ultra2o:=	$(FLIB_sunos5ss20o)
FLIB_sunos5p6o:=	$(FLIB_sunos5ss20o)
FLIB_hpux=		-lf /usr/lib/libisamstub.sl
FLIB_hpuxo:=		$(FLIB_hpux)
FLIB_decunix:=		-lfor -lFutil -lUfor
FLIB_decunixo:=		$(FLIB_decunix)
FLIB_linux:=		/usr/lib/libf2c.a
FLIB_linuxo:=		/usr/lib/libf2c.a

FLIB:= $(FLIB_$(E852_MAKETARGET))

#
# Additional stuff put directly on the link command line
#	(not on the dependancies list)
#

SYSLIB_rs6000=		-Wl,-bI:/usr/lpp/xlf/lib/lowsys.exp
SYSLIB_rs6000o=		$(SYSLIB_rs6000)
SYSLIB_rs6000ppc601o=	$(SYSLIB_rs6000)
SYSLIB_rs6000p2o=	$(SYSLIB_rs6000)
SYSLIB_sunos5=		-xlibmieee -lnsl -ldl -lsocket -lresolv -lsunmath -lm -Wl,-R,/usr/lib:/usr/openwin/lib:/opt/SUNWspro/lib
SYSLIB_sunos5o=		$(SYSLIB_sunos5)
SYSLIB_sunos5ss20o=	$(SYSLIB_sunos5)
SYSLIB_sunos5ultrao=	$(SYSLIB_sunos5)
SYSLIB_sunos5ultra2o=	$(SYSLIB_sunos5)
SYSLIB_sunos5p6o=	$(SYSLIB_sunos5)
SYSLIB_irix5=		-lm
SYSLIB_irix5mips1o=	-nostdlib -L$(E852_HOME)/lib.$(E852_MAKETARGET) -mips1 -lftn -lm
SYSLIB_irix5mips2o=	-mips2 -lm
SYSLIB_irix6n32mips3o=	-n32 -mips3 -lm
SYSLIB_irix6n32mips4r10ko=	$(ARCHFLAGS) -lm
SYSLIB_hpux=		-lm
SYSLIB_hpuxo=		-lm
SYSLIB_decunix=		-ldnet_stub -lots -lm
SYSLIB_decunixo=	-ldnet_stub -lots -lm
SYSLIB_linux:=		-lm
SYSLIB_linuxo:=		-lm

SYSLIB:= $(SYSLIB_$(E852_MAKETARGET))

#
# Location of the CERN library is given by the CERNLIBDIR
# environment variable. If it is not set, the default
# is in $E852_HOME/lib.$E852_MAKETARGET
#

CERNLIB=		-lpacklib -lkernlib -lmathlib

ifndef CERNLIBDIR
CERNLIBDIR:=	$(E852_HOME)/lib.$(E852_MAKETARGET)
endif

#
# Additional C libraries
#

CLIB_irix4 =		-lsun -lm -lc_s
CLIB_irix5 =
CLIB_irix6o =		-lfastm
CLIB_rs6000 =		-lm
CLIB_rs6000o =		-lm
CLIB_rs6000ppc601o =	-lm
CLIB_rs6000p2o =	-lm

CLIB:= $(CLIB_$(E852_MAKETARGET))

#
# ACPP: the C preprocessor to create dependacies list
#

ACPP_default:=		acpp -M
ACPP_hpux:=		gcc -M
ACPP_hpuxo:=		gcc -M
ACPP_decunix:=		cc -M
ACPP_decunixo:=		cc -M
ACPP_sunos5:=		cc -xM1
ACPP_sunos5ss20o:=	cc -xM1
ACPP_sunos5ultrao:=	cc -xM1
ACPP_sunos5ultra2o:=	cc -xM1
ACPP_sunos5p6o:=	cc -xM1
ACPP_irix5:=		cc -M
ACPP_irix5o:=		cc -M
ACPP_irix5mips1o:=	cc -M
ACPP_irix5mips2o:=	cc -M
ACPP_irix6n32mips3o:=	cc -M
ACPP_irix6n32mips4r10ko:=	cc -M
ACPP_linux:=		gcc -M
ACPP_linuxo:=		gcc -M

ACPP:=$(ACPP_$(E852_MAKETARGET))

ifndef ACPP
ACPP=$(ACPP_default)
endif

CC_default:=		cc
CC_linux:=		gcc
CC_linuxo:=		gcc

CC:=$(CC_$(E852_MAKETARGET))

ifndef CC
CC=$(CC_default)
endif

FC_default:=		f77
FC_linux:=		g77
FC_linuxo:=		g77

FC:=$(FC_$(E852_MAKETARGET))

ifndef FC
FC=$(FC_default)
endif

#
# LIBDIR: all the library files are created/updated there
#		(instead of ~e852/lib)
#

ifndef LIBDIR
LIBDIR=$(E852_HOME)/lib.$(E852_MAKETARGET)
endif

#
# VPATH: library search path
#

VPATH_default=		/usr/lib:$(E852_HOME)/lib.$(E852_MAKETARGET)
VPATH_irix4=		/usr/irix4/usr/lib:$(E852_HOME)/lib
VPATH_irix5mips1o=	$(E852_HOME)/lib.$(E852_MAKETARGET)
VPATH_irix5mips2o=	/usr/lib/mips2:/usr/lib:$(E852_HOME)/lib.$(E852_MAKETARGET)
VPATH_irix6=		$(E852_HOME)/lib.irix6:/usr/lib64/mips4:/usr/lib64
VPATH_irix6o=		$(E852_HOME)/lib.irix6o:/usr/lib64/mips4:/usr/lib64
VPATH_irix6n32=		$(E852_HOME)/lib.irix6n32:/usr/lib32/mips4:/usr/lib32
VPATH_irix6n32o=	$(E852_HOME)/lib.irix6n32o:/usr/lib32/mips4:/usr/lib32
VPATH_irix6n32mips3o=	$(E852_HOME)/lib.irix6n32mips3o:/usr/lib32/mips3:/usr/lib32
VPATH_irix6n32mips4r10ko=	$(E852_HOME)/lib.irix6n32mips4r10ko:/usr/lib32/mips4:/usr/lib32
VPATH_sunos5=		/usr/lib:/usr/openwin/lib:/opt/SUNWspro/lib:$(E852_HOME)/lib.$(E852_MAKETARGET)
VPATH_sunos5o:=		$(VPATH_sunos5)
VPATH_sunos5ss20o:=	$(VPATH_sunos5)
VPATH_sunos5ultrao:=	$(VPATH_sunos5)
VPATH_sunos5ultra2o:=	$(VPATH_sunos5)
VPATH_sunos5p6o:=	$(VPATH_sunos5)
VPATH_linux:=		/usr/lib:/usr/X11R6/lib:$(E852_HOME)/lib.$(E852_MAKETARGET)
VPATH_linuxo:=		$(VPATH_linux)

VPATH:= $(VPATH_$(E852_MAKETARGET))

ifndef VPATH
VPATH:=$(VPATH_default)
endif

# add LIBDIR and the CERN library path

VPATH:=$(LIBDIR):$(CERNLIBDIR):$(VPATH)

#
# AR: the 'ar' command- to set group-writable permissions
#			on newly created libraries
#

AR= umask 2; ar

.PHONY : all

#
# Hairy rules to build XFM stuff
#

%_widget.o: %.fm xfm.startup
	-mv -f $*_temp.c $*_temp.c.sav
ifeq ($(USEXFM),yes)
	$(XFM) -startup xfm.startup -cfile $*_temp.c -cname $* -ansic -compilegroup $< -standc -cflags gva
	-mv -f SCCS/$*_temp.c SCCS/$*_temp.c.sav
	cp $*_temp.c SCCS/
	-mv -f SCCS/$*_tempP.h SCCS/$*_tempP.h.sav
	cp $*_tempP.h SCCS/
	rm -f $(E852_HOME)/include/$*.h
	cp $*_temp.h $(E852_HOME)/include/$*.h
	chmod a=r  $(E852_HOME)/include/$*.h
	rm -f SCCS/$*.h
	cp $*_temp.h SCCS/$*.h
else
	-rm -f $*_temp.c
	cp SCCS/$*_temp.c .
	-rm -f $*_tempP.h
	cp SCCS/$*_tempP.h .
	-rm -f $*_temp.h
	cp SCCS/$*.h $*_temp.h
endif
	perl -e 'while (<STDIN>) { s/\(int \*\)call_data/\(int \*\*\)call_data/go; print $$_; }; ' < $*_temp.c > $*_temp.c-tmp1
	mv $*_temp.c-tmp1 $*_temp.c
	$(CC) -g -c $*_temp.c -I. -DFUNCPROTO=15 $(INCLUDE)
	-rm -f $*_temp.c.sav
	mv $*_temp.h $*.h
	mv $*_temp.o $*_widget.o

%_cwidget.o: %.c %.h %P.h
	$(ACPP) $<  $(INCLUDE) > $@.depend
	$(CC) -c $*.c $(CDBG) $(ARCHFLAGS) $(INCLUDE) $(CFLAGS)
	rm -f $(E852_HOME)/include/$*.h
	cp $*.h $(E852_HOME)/include
	chmod a=r  $(E852_HOME)/include/$*.h
	mv $*.o $*_cwidget.o

%_fm.o: %.fm xfm.startup
	-mv -f $*_temp.c $*_temp.c.sav
ifeq ($(USEXFM),yes)
	$(XFM) -startup xfm.startup -ansic -compile $< -cfile $*_temp.c -cname $*_fm
	-mv -f SCCS/$*_temp.c SCCS/$*_temp.c.sav
	cp $*_temp.c SCCS/
else
	cp SCCS/$*_temp.c .
endif
	perl -e 'while (<STDIN>) { s/\(int \*\)call_data/\(int \*\*\)call_data/go; print $$_; }; ' < $*_temp.c > $*_temp.c-tmp1
	mv $*_temp.c-tmp1 $*_temp.c
	$(CC) -g -c $*_temp.c $(CDBG) $(ARCHFLAGS) $(CFLAGS) $(INCLUDE) -DFUNCPROTO=15
	-rm -f $*_temp.c.sav
	mv $*_temp.o $*_fm.o

%_fm1.o: %.fm
	rm -f $*_fm1.c
	$(XFM) -ansic -compile $< -cfile $*_fm1.c -cname $*
	$(CC) -g -c $*_fm1.c $(CDBG) $(ARCHFLAGS) $(CFLAGS) $(INCLUDE)

%: %.o
	-mv -f $@ $@.sav
	$(CC) -o $@ $^ $(SYSLIB)

%: %.c
	$(MAKE) $*.o
	$(MAKE) $*

%_p:
	-mv -f $@ $@.sav
	cc -p -o $@ $^

#
#
# These are the main rules to build object files-
#
# %.o --- build "normal" object files (usually with -g)
# %.oo -- build "optimised" object files (always with -O)
#
#

%.oo: %.c
	$(ACPP) $< $(INCLUDE) > $*.oo.depend
	$(CC) -c $< $(COPT) $(ARCHFLAGS) $(CFLAGS) $(INCLUDE)
	mv $*.o $*.oo

%.oo: %.f
	$(FC) -c $< $(COPT) $(ARCHFLAGS) $(FFLAGS) $(FINCLUDE)
	mv $*.o $*.oo

%.o: %.f
	$(FC) -c $< $(CDBG) $(ARCHFLAGS) $(FFLAGS) $(FINCLUDE)

%.o: %.c
	$(ACPP) $< $(INCLUDE) > $@.depend
	$(CC) -c $< $(CDBG) $(ARCHFLAGS) $(CFLAGS) $(INCLUDE)

#
# Rules to build object files that go into libaries
#

(%.o): %.c
	$(CC) -c $< $(CDBG) $(ARCHFLAGS) $(CFLAGS) $(INCLUDE)

(%.oo): %.c
	$(CC) -c $< $(COPT) $(ARCHFLAGS) $(CFLAGS) $(INCLUDE)
	mv $*.o $*.oo

(%.o): %.f
	$(FC) -c $< $(CDBG) $(ARCHFLAGS) $(FFLAGS) $(FINCLUDE)

(%.oo): %.f
	$(FC) -c $< $(COPT) $(ARCHFLAGS) $(FFLAGS) $(FINCLUDE)
	mv $*.o $*.oo

#
# Rules to build dependancy files
#
# Note: The perl ugly magic is needed to transform the output of
# "gcc -M" into what I want (remove the system (/usr/include/*.h)
# header files and add the "libfile.a(foo.o):" dependancy).

#
# Dependancies for normal .o and .oo files
#

#%.o.depend: %.c
#	gcc -M -MG $< $(INCLUDE) | perl -n -e 's#\S+\.o\ *:\ *(\S+)\.c#\1\.o:#; s#/usr/include/\S+h##g; print;' > $@
#
#%.oo.depend: %.c
#	gcc -M -MG $< $(INCLUDE) | perl -n -e 's#\S+\.o\ *:\ *(\S+)\.c#\1\.oo:#; s#/usr/include/\S+h##g; print;' > $@

#
# Dependancies for .o and .oo files that go into library file $(LIBFILE)

#%.o.a.depend: %.c
#	gcc -M -MG $< $(INCLUDE) | perl -n -e 's#\S+\.o\ *:\ *(\S+)\.c#$(LIBFILE)\(\1\.o\):#; s#/usr/include/\S+h##g; print;' > $@
#
#%.oo.a.depend: %.c
#	gcc -M -MG $< $(INCLUDE) | perl -n -e 's#\S+\.o\ *:\ *(\S+)\.c#$(LIBFILE)\(\1\.oo\):#; s#/usr/include/\S+h##g; print;' > $@

#
# Rules to stuff the compiled objects $(OBJS) into the library file $(LIBFILE)
#

YOBJS:=  $(wildcard *.o *.oo)
XOBJS:=  $(filter $(OBJS),$(YOBJS))

libcollect:
ifneq ($(words $(XOBJS)),0)
	ar rv $(LIBFILE) $(XOBJS)
	rm -f $(XOBJS)
else
	@echo \\c
endif


#
#
#    Here we define the C++ stuff- which compiler we use
#
#    and which libraries we link with
#
# 


ifndef GCC
GCC=1
endif

ifeq ($(GCC),0)
CXX= CC
CCOPT= $(COPT)
CCDBG= $(CDBG)
CCARCHFLAGS= $(ARCHFLAGS)

ifndef SGIcomplex
COMPLEX= -lcomplex
else
CCFLAGS:= -I$(E852_HOME)/include/CC/gcomplex $(CCFLAGS)
COMPLEX= -lgcomplex
endif
endif

ifeq ($(GCC),1)
CXX= g++
CCOPT= -O
CCDBG= $(CDBG)
CCARCHFLAGS= $(GARCHFLAGS)
CCFLAGS:= -Wall -I$(E852_HOME)/include/CC/ISO
COMPLEX=
endif

%.o: %.C
	$(CXX) $(CCPPFLAGS) -c -M $<  $(INCLUDE) > $@.depend
	$(CXX) -c $< $(CCDBG) $(CCARCHFLAGS) $(CCFLAGS) $(INCLUDE)

%.oo: %.C
	$(CXX) $(CCPPFLAGS) -c  -M $<  $(INCLUDE) > $@.depend
	$(CXX) -c $< $(CCOPT) $(CCARCHFLAGS) $(CCFLAGS) $(INCLUDE)
	mv $*.o $*.oo

%.o: %.c++
	$(CXX) $(CCPPFLAGS) -c -M $<  $(INCLUDE) > $@.depend
	$(CXX) -c $< $(CCDBG) $(CCARCHFLAGS) $(CCFLAGS) $(INCLUDE)

%.oo: %.c++
	$(CXX) $(CCPPFLAGS) -c -M $<  $(INCLUDE) > $@.depend
	$(CXX) -c $< $(CCOPT) $(CCARCHFLAGS) $(CCFLAGS) $(INCLUDE)
	mv $*.o $*.oo

%.o: %.cxx
	$(CXX) $(CCPPFLAGS) -c -M $<  $(INCLUDE) > $@.depend
	$(CXX) -c $< $(CCDBG) $(CCARCHFLAGS) $(CCFLAGS) $(INCLUDE)

%.oo: %.cxx
	$(CXX) $(CCPPFLAGS) -c -M $<  $(INCLUDE) > $@.depend
	$(CXX) -c $< $(CCOPT) $(CCARCHFLAGS) $(CCFLAGS) $(INCLUDE)
	mv $*.o $*.oo

%.o: %.cc
	$(CXX) $(CCPPFLAGS) -c -M $<  $(INCLUDE) > $@.depend
	$(CXX) -c $< $(CCDBG) $(CCARCHFLAGS) $(CCFLAGS) $(INCLUDE)

%.oo: %.cc
	$(CXX) $(CCPPFLAGS) -c -M $<  $(INCLUDE) > $@.depend
	$(CXX) -c $< $(CCOPT) $(CCARCHFLAGS) $(CCFLAGS) $(INCLUDE)
	mv $*.o $*.oo


#
# end of C++ stuff
# 

#
# Helper rules to install stuff in ~e852/bin
#

ifndef E852_ARCH
E852_ARCH=$(ARCH)
endif

%.perl.install: %.perl
	chmod a+x $^
	cp $^ $(E852_HOME)/bin/$^.new
	rm -f $(E852_HOME)/bin/$^
	mv $(E852_HOME)/bin/$^.new $(E852_HOME)/bin/$^

%.perl5.install: %.perl5
	chmod a+x $^
	cp $^ $(E852_HOME)/bin/$^.new
	rm -f $(E852_HOME)/bin/$^
	mv $(E852_HOME)/bin/$^.new $(E852_HOME)/bin/$^

%.sh.install: %.sh
	chmod a+x $^
	cp $^ $(E852_HOME)/bin/$^.new
	rm -f $(E852_HOME)/bin/$^
	mv $(E852_HOME)/bin/$^.new $(E852_HOME)/bin/$^

%.install2: %
	cp $^ $(E852_HOME)/bin.$(E852_ARCH)/$^.new
	rm -f $(E852_HOME)/bin.$(E852_ARCH)/$^
	mv $(E852_HOME)/bin.$(E852_ARCH)/$^.new $(E852_HOME)/bin.$(E852_ARCH)/$^

%.install2-save: %
	cp $^ $(E852_HOME)/bin.$(E852_ARCH)/$^.new
	-mv $(E852_HOME)/bin.$(E852_ARCH)/$^ $(E852_HOME)/bin.$(E852_ARCH)/$^.previous
	rm -f $(E852_HOME)/bin.$(E852_ARCH)/$^
	mv $(E852_HOME)/bin.$(E852_ARCH)/$^.new $(E852_HOME)/bin.$(E852_ARCH)/$^

%.install: %
	cp $^ $(E852_HOME)/bin.$(E852_MAKETARGET)/$^.new
	rm -f $(E852_HOME)/bin.$(E852_MAKETARGET)/$^
	mv $(E852_HOME)/bin.$(E852_MAKETARGET)/$^.new $(E852_HOME)/bin.$(E852_MAKETARGET)/$^

%.install-save: %
	cp $^ $(E852_HOME)/bin.$(E852_MAKETARGET)/$^.new
	-mv $(E852_HOME)/bin.$(E852_MAKETARGET)/$^ $(E852_HOME)/bin.$(E852_MAKETARGET)/$^.previous
	rm -f $(E852_HOME)/bin.$(E852_MAKETARGET)/$^
	mv $(E852_HOME)/bin.$(E852_MAKETARGET)/$^.new $(E852_HOME)/bin.$(E852_MAKETARGET)/$^


#
# Make everything precious- to avoid removal of source files after targets are built.
#

.PRECIOUS: %

debug:
	@echo E852_HOME\\t        $(E852_HOME)
	@echo LIBDIR\\t\\t         $(LIBDIR)
	@echo E852_MAKETARGET\\t $(E852_MAKETARGET)
	@echo VPATH\\t\\t        $(VPATH)
	@echo CC\\t\\t           $(CC)
	@echo CFLAGS\\t\\t       $(CFLAGS)
	@echo CCFLAGS\\t\\t       $(CCFLAGS)
	@echo INCLUDE\\t\\t      $(INCLUDE)


Debug:
	@echo E852_HOME\\t        $(E852_HOME)
	@echo LIBDIR\\t\\t         $(LIBDIR)
	@echo E852_MAKETARGET\\t $(E852_MAKETARGET)
	@echo VPATH\\t\\t        $(VPATH)
	@echo CC\\t\\t           $(CC)
	@echo CFLAGS\\t\\t       $(CFLAGS)
	@echo CCFLAGS\\t\\t       $(CCFLAGS)
	@echo INCLUDE\\t\\t      $(INCLUDE)


clean:
	rm -f *.o *.oo *.depend
	echo >dummy.depend


#
# end file
#
