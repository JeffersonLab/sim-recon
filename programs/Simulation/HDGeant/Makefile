OStype = $(shell uname)
BINDIR = ../bin.$(OStype)

COPTS = -g

gpp: gpp.C
	c++ -o $(BINDIR)/$@ $(COPTS) $^

clean:
	rm -f *.o core
