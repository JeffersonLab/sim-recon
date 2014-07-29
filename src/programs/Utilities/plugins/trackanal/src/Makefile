#  makefile for check2.cc

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)


all: check2

check2.o: check2.cc 
	gcc -pg -c -g check2.cc -pthread -m64 -I/usr/share/cern/root5.28/include

check2: check2.o
	gcc -pg -o check2 check2.o -L/usr/share/cern/root5.28/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lGui -pthread -lm -ldl -rdynamic


clean:
	rm -f check2.o check2