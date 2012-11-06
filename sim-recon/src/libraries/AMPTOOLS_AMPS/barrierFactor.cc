#include <cmath>
#include "AMPTOOLS_AMPS/breakupMomentum.h"
#include "AMPTOOLS_AMPS/barrierFactor.h"


// mass0 = mass of parent
// spin  = angular momentum of the decay
// mass1 = mass of first daughter
// mass2 = mass of second daughter

double barrierFactor( double mass0, int spin, double mass1, double mass2 ){
	
        double q;

        q = breakupMomentum(mass0, mass1, mass2);

        return barrierFactor( q, spin );

}


// q     = breakup momentum
// spin  = angular momentum of the decay

double barrierFactor ( double q, int spin ){

        double barrier;
        double z;

        z = ( (q*q) / (0.1973*0.1973) );

        switch (spin){

          case 0:
             barrier = 1.0;
             break;

          case 1:
             barrier = sqrt( (2.0*z) /
                             (z + 1.0) );
             break;
        
          case 2:
             barrier = sqrt( (13.0*z*z) /
                            ((z-3.0)*(z-3.0) + 9.0*z) );
             break;

          case 3:
             barrier = sqrt( (277.0*z*z*z) /
                            (z*(z-15.0)*(z-15.0) + 
                             9.0*(2.0*z-5.0)*(2.0*z-5.0)) );
             break;

          case 4:
             barrier = sqrt( (12746.0*z*z*z*z) /
                            ((z*z-45.0*z+105.0)*(z*z-45.0*z+105.0) +
                             25.0*z*(2.0*z-21.0)*(2.0*z-21.0)) );
             break;

          default:
             barrier = 0.0;
        }

        return barrier;
	
}
