CXX := g++
CXX_FLAGS:= 

# use ANAGATE only
ANAGATE = 1
CXX_FLAGS += -DANAGATE
INC_DIR += -I. -IincludeAnaGate/ 
LIBS := libAnaGate/libCANDLL-1.8.so.0 libAnaGate/libAnaCommon-1.0.12.so.0

.PHONY: default 

default: beamtestHV

gluex_newCanlib.o:	gluex_newCanlib.cxx gluex_newCanlib.h
	$(CXX) $(CXX_FLAGS)  $(INC_DIR) -c gluex_newCanlib.cxx 
beamtestHV: beamtestHV.cxx gluex_newCanlib.o
	$(CXX) $(CXX_FLAGS)  $(INC_DIR) -c beamtestHV.cxx
	g++ gluex_newCanlib.o beamtestHV.o  $(LIBS) -o beamtestHV
clean:
	rm -f *.o *~ beamtestHV
