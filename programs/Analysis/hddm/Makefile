DOM = /home/www/docs/xerces

CC = /usr/bin/g++ -g

hddm-xml: hddm-xml.cpp
	$(CC) -I$(DOM)/include -o $@ $^ \
	-L$(DOM)/lib -lxerces-c1_4

hddm-c: hddm-c.cpp
	$(CC) -I$(DOM)/include -o $@ $^ \
	-L$(DOM)/lib -lxerces-c1_4

clean:
	/bin/rm -f *.o core *.depend
