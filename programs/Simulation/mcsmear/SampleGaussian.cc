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

int INITIALIZED_ERF_TABLE = 0;
#define N_ERF_BINS 4096
static double ERF_TABLE[N_ERF_BINS];
#define RANDOM_MAX (float)(0x7FFFFFFF)

double FindERFFromIntegralFraction(double fraction);
void InitializeERF_Table(void);

//--------------------------
// SampleGaussian
//--------------------------
double SampleGaussian(double sigma)
{
	int r;
	double s, f;

	if(!INITIALIZED_ERF_TABLE)InitializeERF_Table();

	r  = random();
	s  = (float)r/RANDOM_MAX;
	f  = ERF_TABLE[(int)(s*(float)(N_ERF_BINS-1))];
	f *= sigma*M_SQRT2; // sqrt(2) comes from erf using sig=1/sqrt(2)
	f *= random()%2 ? -1.0:+1.0;

	return f;
}

//--------------------------
// SampleRange
//--------------------------
double SampleRange(double x1, double x2)
{
	int r;
	double s, f;
	double xlo, xhi;
	
	if(x1<x2){
		xlo = x1;
		xhi = x2;
	}else{
		xlo = x2;
		xhi = x1;
	}

	r  = random();
	s  = (float)r/RANDOM_MAX;
	f  = xlo + s*(xhi-xlo);
	
	return f;
}

//--------------------------
// FindERFFromIntegralFraction
//--------------------------
double FindERFFromIntegralFraction(double fraction)
{
	double x1=0.0, x2=1.0E6;
	double x;
	double if1,if2,ifm;

	do{
		x = (x1+x2)/2.0;
		if1 = erf(x1);
		if2 = erf(x2);
		ifm = erf(x);

		if(ifm>fraction){
			x2 = x;
		}else{
			x1 = x;
		}

	}while(fabs(1.0-(ifm/fraction))>1.0E-4);

	return x;
}

//--------------------------
// InitializeERF_Table
//--------------------------
void InitializeERF_Table(void)
{
	int i;
	double f;

	for(i=0;i<N_ERF_BINS;i++){
		f = (double)i/(double)(N_ERF_BINS-1);
		ERF_TABLE[i] = FindERFFromIntegralFraction(f);
	}
	
	INITIALIZED_ERF_TABLE = 1;
}
