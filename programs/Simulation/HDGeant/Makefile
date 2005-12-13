
# Default makefile for HDGeant
# this just invokes make in the gelhad directory and
# then make with Makefile.bms in this directory.


all:
	make -C gelhad
	make -C hitutil
	make -f Makefile.bms

clean:
	make -C gelhad clean
	make -C hitutil clean
	make -f Makefile.bms clean

exec_clean:
	make -f Makefile.bms exec_clean
	
