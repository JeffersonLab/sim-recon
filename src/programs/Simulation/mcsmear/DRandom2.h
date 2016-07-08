// $Id$
//

// Random number generator used in mcsmear. All random numbers
// should come from the global "gDRandom" object declared here.
//
// Because we want to record the seeds used for every event,
// we use the TRandom2 class. This one has only 3 seed values
// (as opposed to 24 for TRandom1 and 624 for TRandom3) but 
// is fast with a large period (10^26).
//
// Because the seeds in TRandom2 are declared protected in
// the class with no methods to access/set them, we derive
// a new class, DRandom2 from TRandom2. This allows us access
// to the numbers for easy recording/retrieving. 

#ifndef _DRANDOM2_H_
#define _DRANDOM2_H_

#include <TRandom2.h>
#include <iostream>
using std::cerr;
using std::endl;

class DRandom2:public TRandom2{
	public:
		
		DRandom2(UInt_t seed=1):TRandom2(seed){}
	
		void GetSeeds(UInt_t &seed, UInt_t &seed1, UInt_t &seed2){
			seed = this->fSeed;
			seed1 = this->fSeed1;
			seed2 = this->fSeed2;
		}
		
		void SetSeeds(UInt_t &seed, UInt_t &seed1, UInt_t &seed2){
		
			// See the comments in TRandom2::SetSeed(int)
			if( (seed<2) | (seed1<8) | (seed2<16) ){
				cerr << endl;
				cerr << "*********************************************************" << endl;
				cerr << "WARNING: Random seeds passed to DRandom2::SetSeeds have" << endl;
				cerr << "forbidden values: " << endl;
				cerr << "  seed = " << seed << "  (must be at least 2)" <<endl;
				cerr << "  seed1 = " << seed1 << "  (must be at least 8)" <<endl;
				cerr << "  seed1 = " << seed2 << "  (must be at least 16)" <<endl;
				cerr << "See comments in source for TRandom2::SetSeed(int)" <<endl;
				cerr << "The seeds will all be adjusted to be in range." <<endl;
				cerr << "*********************************************************" << endl;
				cerr << endl;
				seed += 2;
				seed1 += 8;
				seed2 += 15;
			}
		
		
			this->fSeed = seed;		
			this->fSeed1 = seed1;		
			this->fSeed2 = seed2;		
		}
		
		// legacy mcsmear interface
		inline double SampleGaussian(double sigma) {
			return Gaus(0.0, sigma);
		}

		inline double SamplePoisson(double lambda) {	
			return Poisson(lambda);
		}

		inline double SampleRange(double x1, double x2) {
			double s, f;
			double xlo, xhi;
	
			if(x1<x2){
				xlo = x1;
				xhi = x2;
			}else{
				xlo = x2;
				xhi = x1;
			}

			s  = Rndm();
			f  = xlo + s*(xhi-xlo);
	
			return f;
		}


};

#endif  // _DRANDOM2_H_

extern DRandom2 gDRandom;


