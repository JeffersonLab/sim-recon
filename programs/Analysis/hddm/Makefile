DOM = $(BUILDS)/xerces-c

CC = /usr/bin/g++ -g

hddm-xml: hddm-xml.cpp
	$(CC) -I../include -I$(DOM)/include -o $@ $^ \
	-L$(DOM)/lib -lxerces-c1_5

hddm-c: hddm-c.cpp
	$(CC) -I../include -I$(DOM)/include -o $@ $^ \
	-L$(DOM)/lib -lxerces-c1_5

clean:
	/bin/rm -f *.o core *.depend
