// DPhoton member functions

#include "DPi0.h"


DPi0::DPi0()
{
   fMom4.SetXYZT(0., 0., 0., 0.) ;
}

DPi0::~DPi0()
{
}


// Set pi0 four momentum
void DPi0::setMom4(const DLorentzVector gamma1, const DLorentzVector gamma2)
{
   fMom4 = gamma1 + gamma2;
}

// set Pi0 bits with respect to the photon detection 
void DPi0::setTags(const unsigned int tag1, const unsigned int tag2 )
{
   fTag1 = tag1;
   fTag2 = tag2;
}
