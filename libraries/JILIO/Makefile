
.PHONY: hd_serializers.xml

hd_serializers.xml:
	./make_serializers
	make -f Makefile.bms

clean:
	rm -f hd_serializers.xml hd_serializers.cc
	make -f Makefile.bms clean
