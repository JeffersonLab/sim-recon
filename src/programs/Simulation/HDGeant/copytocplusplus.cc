
#include <stdlib.h>

#include <iostream>
#include <string>
using namespace std;


string INFILE;
string OUTFILE;
bool POSTSMEAR = false;
string MCSMEAROPTS;
bool DELETEUNSMEARED = false;

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
}

