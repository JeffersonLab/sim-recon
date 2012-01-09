// $Id$
//

// Random number routines used in mcsmear. All random numbers
// should come from these routines.
//
#include <stdlib.h>
#include <math.h>

#include "DRandom2.h"

DRandom2 gDRandom; // declared extern in DRandom2.h

//--------------------------
// SampleGaussian
//--------------------------
double SampleGaussian(double sigma)
{
	return gDRandom.Gaus(0.0, sigma);
}

//--------------------------
// SamplePoisson
//--------------------------
double SamplePoisson(double lambda)
{	
	return gDRandom.Poisson(lambda);
}

//--------------------------
// SampleRange
//--------------------------
double SampleRange(double x1, double x2)
{
	double s, f;
	double xlo, xhi;
	
	if(x1<x2){
		xlo = x1;
		xhi = x2;
	}else{
		xlo = x2;
		xhi = x1;
	}

	s  = gDRandom.Rndm();
	f  = xlo + s*(xhi-xlo);
	
	return f;
}

