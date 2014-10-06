#include <stdint.h>
#include <cmath>
using namespace std;

#include <TRandom2.h>

static TRandom2 randgen;
bool NO_PEDESTAL = true;              // flag to completely disable pedestal generation
bool NO_RANDOM_PEDESTAL = false;       // turn off random components of pedestals
float MEAN_PEDESTAL = 100.0;           // mean pedestal in single sample fADC counts
float SIGMA_COMMON_PEDESTAL = 4.0;     // common mode noise in single sample fADC counts
float SIGMA_INDIVIDUAL_PEDESTAL = 3.0; // stochastic noise in single sample fADC counts


extern "C" {
	void GetPedestals(uint32_t *peds, uint32_t Npeds);
}

void GetPedestals(uint32_t *peds, uint32_t Npeds)
{
	/// Generate a set of pedestal values based on the global variables.
	/// Ths will fill in the first Npeds-1 values in "peds" with a pedestal
	/// near zero and Npeds-th value will be near the mean pedestal setting.
	/// The idea is that the last value delivered will contain the mean
	/// pedestal for all channels in the module plus common mode noise.
	/// The first set of values will contain individual pedestal fluctuations
	/// for each channel. The f250 and f125 will use the one common pedestal
	/// to adjust the pulse integral and pulse peak values. For the measured
	/// pedestals output in Pulse Pedestal words, the stochastic part is added
	/// in as well to represent the inaccuracy of that measurement. Note that 
	/// this means the actual pedestal value used is not actually written
	/// to the output file anywhere.
	///
	/// If the global boolean NO_PEDESTAL is set to true, then all values
	/// will be set to zero.
	///
	/// If the global NO_RANDOM_PEDESTAL is set to true, then the first
	/// Npeds-1 values will be set to zero and the last value to the mean
	/// pedestal (=MEAN_PEDESTAL) with no random variation.
	///
	/// The values MEAN_PEDESTAL, SIGMA_COMMON_PEDESTAL, and
	/// SIGMA_INDIVIDUAL_PEDESTAL are used to calculate the values. Random
	/// sampling is done from Gaussian's with the specified sigmas.

	if(NO_PEDESTAL){
		for(uint32_t i=0; i<Npeds; i++) peds[i] = 0;
	}else if(NO_RANDOM_PEDESTAL){
		for(uint32_t i=0; i<Npeds-1; i++) peds[i] = 0;	
		peds[Npeds - 1] = (uint32_t)MEAN_PEDESTAL;
	}else{
		for(uint32_t i=0; i<Npeds-1; i++) peds[i] = round(randgen.Gaus(0.0, SIGMA_INDIVIDUAL_PEDESTAL));
		peds[Npeds - 1] = round(randgen.Gaus(MEAN_PEDESTAL, SIGMA_COMMON_PEDESTAL));
	}

}

