// $Id$

#ifndef _FIT_UTILS_H_
#define _FIT_UTILS_H_

#define MAX_CHISQV_HITS 100

typedef struct{
	float x0;
	float y0;
	float theta;
	float phi;
	float q;
	float p;
	float p_trans;
	float chisq;
	float chisqv[MAX_CHISQV_HITS];
	int nhits;
}FitParms_t;


#endif //_FIT_UTILS_H_
