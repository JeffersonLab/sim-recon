#include <math.h>
#define JMAX 20

double rtnewt(void (*funcd)(double, double *, double *), double x1, double x2,
	double xacc)
{
	void nrerror(char error_text[]);
	int j;
	double df,dx,f,rtn;

	rtn=0.5*(x1+x2);
	for (j=1;j<=JMAX;j++) {
		(*funcd)(rtn,&f,&df);
		dx=f/df;
		rtn -= dx;
		if ((x1-rtn)*(rtn-x2) < 0.0)
			nrerror("Jumped out of brackets in rtnewt");
		if (fabs(dx) < xacc) return rtn;
	}
	nrerror("Maximum number of iterations exceeded in rtnewt");
	return 0.0;
}
#undef JMAX
