
#include <stdlib.h>

#include <iostream>
#include <string>
using namespace std;


string INFILE;
string OUTFILE;
bool POSTSMEAR = false;
string MCSMEAROPTS;
bool DELETEUNSMEARED = false;
float BGGATE1=-200.0;
float BGGATE2= 200.0;

// Declare routines callable from FORTRAN
extern "C"{

	// The following allows one to specify that mcsmear should
	// automatically be run after hdgeant_ returns. This is
	// controlled by the POSTSMEAR and MCSMEAROPTS keyowrds
	// in the control.in file.
	void copytocplusplus_(char *infile, char *outfile, int *postsmear, char *mcsmearopts, int *deleteunsmeared){
		if(infile)INFILE = infile;
		if(outfile){
			OUTFILE= outfile;
			POSTSMEAR = *postsmear != 0;
			if(POSTSMEAR && mcsmearopts)MCSMEAROPTS = mcsmearopts;
			DELETEUNSMEARED = *deleteunsmeared != 0;
		}
	}
	
	// Copy the values of BGGATE card to global variables
	// visible from C++. This is to allow the bcalHit.cc
	// routines to use them to set the histogram limits
	void copygatetocplusplus_(float bggate1, float bggate2){
		BGGATE1 = bggate1;
		BGGATE2 = bggate2;
	}
}

