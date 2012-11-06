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


#include <TRandom2.h>

class DRandom2:public TRandom2{
	public:
		
		DRandom2(UInt_t seed=1):TRandom2(seed){}
	
		void GetSeeds(UInt_t &seed, UInt_t &seed1, UInt_t &seed2){
			seed = this->fSeed;
			seed1 = this->fSeed1;
			seed2 = this->fSeed2;
		}
		
		void SetSeeds(UInt_t &seed, UInt_t &seed1, UInt_t &seed2){
			this->fSeed = seed;		
			this->fSeed1 = seed1;		
			this->fSeed2 = seed2;		
		}
};

extern DRandom2 gDRandom;


