#! gnumake

#---------------------------------------------------
# The stripped down genr8 Makefile 
#---------------------------------------------------
#


INCLUDE := -I. -Iinc  
CC = cc

ifndef UNAME
UNAME = $(shell uname)
endif

PROGRAM= genr8

genr8: genr8.o genkin.o
	$(CC) -o $@  $^  -lm

genr8.o: genr8.c
	$(CC) -c -o genr8.o   genr8.c $(INCLUDE)

genkin.o: genkin.c
	$(CC) -c -o genkin.o   genkin.c $(INCLUDE)


clean:
	@rm  *.o 


install: genr8
	mv genr8 $(HALLD_HOME)/bin/$(UNAME)
















