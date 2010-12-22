#if !defined(WIGNERD)
#define WIGNERD

#include <complex>

#include "GPUManager/GPUCustomTypes.h"

using std::complex;

GDouble wignerDSmall( GDouble aj, GDouble am, GDouble an, GDouble beta );
complex< GDouble > wignerD( int l, int m, int n, GDouble cosTheta, GDouble phi );
complex< GDouble > Y( int l, int m, GDouble cosTheta, GDouble phi );

#endif
