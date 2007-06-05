// DFCALPhoton member functions

#include "DPhoton.h"


DPhoton::DPhoton()
{
   fMom4.SetXYZT(0., 0., 0., 0.) ;
}

DPhoton::~DPhoton()
{
}

// Set photon four momentum
void DPhoton::setMom4(const TLorentzVector gamma)
{
     fMom4 = gamma;
}

// Tag photon origin: 0/1 for FCAL/BCAL
void DPhoton::setTag(const unsigned int tag)
{
   fTag = tag;
}

// Distance to track's ReferenceTrajectory
void DPhoton::setDtRT(const double dtrt)
{
   fDtRT = dtrt;
}


