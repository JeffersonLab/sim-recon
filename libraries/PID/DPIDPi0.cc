// DPIDPhoton member functions

#include "DPIDPi0.h"


DPIDPi0::DPIDPi0()
{
   fMom4.SetXYZT(0., 0., 0., 0.) ;
}

DPIDPi0::~DPIDPi0()
{
}


// Set pi0 four momentum
void DPIDPi0::setMom4(const TLorentzVector gamma1, const TLorentzVector gamma2)
{
   fMom4 = gamma1 + gamma2;
}

// set Pi0 bits with respect to the photon detection 
void DPIDPi0::setOrig(const unsigned int tag1, const unsigned int tag2 )
{
   fTag1 = tag1;
   fTag2 = tag2;
}
