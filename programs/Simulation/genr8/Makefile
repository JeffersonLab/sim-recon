#
#
#
# Makefile for mc "monetcarlo" subsystem
#

include inc/rules.make


#.c:
#        cc -o $* $(CFLAGS) -I$(INCLUDE) -I. $< -lm 

all: genr8 


SYSLIB=	-ldnet_stub -lots -lm
FLIB=	-lfor -lFutil -lUfor

CERNLIB=                -lpacklib -lmathlib -lkernlib 

INCLUDE= -I.  -Iinc  -I../include -DFUNCPROTO


#VPATH:=$(LIBDIR):$(CERNLIBDIR):$(VPATH)



kinematics.o:  kinematics.c
	 cc -c  kinematics.c   -g  $(INCLUDE)



genr8.o : genr8.c 
	 cc -c genr8.c   -g $(INCLUDE)

genr8 : genr8.o  kinematics.o
	cc -o $@ $@.o kinematics.o  libdata.a  -lm





# End of file

