// $Id$
//
// Copied from PrimEx psim_digitize June 22, 2005  David Lawrence

// Randomly sample a gaussian distribution with a sigma of 1.0
// and mean of 0.0.
// This is done by creating a table of error function values
// on the first call to SampleGassian() and then inerpolating
// on that and all subsequent calls. 
// The sigma parameter is used to scale the result to the desired
// units.

#include <stdlib.h>
#include <math.h>

#include <TRandom3.h>

TRandom3 *rnd;

//--------------------------
// SampleGaussian
//--------------------------
double SampleGaussian(double sigma)
{
	if(!rnd)rnd = new TRandom3;
	
	return rnd->Gaus(0.0, sigma);
}

//--------------------------
// SamplePoisson
//--------------------------
double SamplePoisson(double lambda)
{
        if(!rnd)rnd = new TRandom3;
	
	return rnd->Poisson(lambda);
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

	if(!rnd)rnd = new TRandom3;
	s  = rnd->Rndm();
	f  = xlo + s*(xhi-xlo);
	
	return f;
}

