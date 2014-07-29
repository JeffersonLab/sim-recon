
all: doxygen/html/index.html

doxygen/html/index.html: doxygen.cfg
	/usr/bin/doxygen doxygen.cfg

devel: doxygen_devel/html/index.html

doxygen_devel/html/index.html: doxygen_devel.cfg
	/usr/bin/doxygen doxygen_devel.cfg

clean:
	rm -rf doxygen
	rm -rf doxygen_devel
