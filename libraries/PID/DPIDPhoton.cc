// DFCALPhoton member functions

#include "DPIDPhoton.h"


DPIDPhoton::DPIDPhoton()
{
   fMom4.SetXYZT(0., 0., 0., 0.) ;
}

DPIDPhoton::~DPIDPhoton()
{
}

// Set photon four momentum
void DPIDPhoton::setMom4(const TLorentzVector gamma)
{
//   fMom4.SetT(gamma.T());
//   fMom4.SetVect(gamma.Vect());
     fMom4 = gamma;
}


// Tag photon origin: 0/1 for FCAL/BCAL
void DPIDPhoton::setTag(const unsigned int tag)
{
   fTag = tag;
}


