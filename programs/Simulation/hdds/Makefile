DOM = $(BUILDS)/xerces-c

CC = /usr/bin/g++ -g

hdds-geant: hdds-geant.cpp
	$(CC) -I$(DOM)/include -o $@ $^ \
	-L$(DOM)/lib -lxerces-c1_5

hdds-mcfast: hdds-mcfast.cpp
	$(CC) -I$(DOM)/include -o $@ $^ \
	-L$(DOM)/lib -lxerces-c1_5

clean:
	/bin/rm -f *.o core *.depend
