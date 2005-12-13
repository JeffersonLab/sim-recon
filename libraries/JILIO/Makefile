
.PHONY: hd_serializers.xml

hd_serializers.xml:
	./make_serializers
	make -f Makefile.bms

clean:
	rm -f hd_serializers.xml hd_serializers.cc ../include/hd_classes.h ../include/hd_serializers.h
	make -f Makefile.bms clean
