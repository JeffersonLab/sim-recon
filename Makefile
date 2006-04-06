

all:
	make -C external
	make -C libraries
	make -C programs

clean:
	make -C external clean
	make -C libraries clean
	make -C programs clean

