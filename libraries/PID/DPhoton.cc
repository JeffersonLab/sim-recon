// DFCALPhoton member functions

#include "DPhoton.h"


DPhoton::DPhoton()
{
   fPosition.SetXYZ(0., 0., 0.) ;
   fMomentum.SetXYZ(0., 0., 0.) ;
   fEnergy = 0 ;
}

DPhoton::~DPhoton()
{
}

// Set photon momentum
void DPhoton::setMomentum(const DVector3 aMom)
{
     fMomentum = aMom;
}

// Set photon position
void DPhoton::setPosition(const DVector3 aPosition)
{
     fPosition = aPosition;
}

// Set photon energy
void DPhoton::setEnergy(const double aEnergy)
{
     fEnergy = aEnergy;
}

// Tag photon origin: 0/1 for FCAL/BCAL
void DPhoton::setTag(const unsigned int aTag)
{
   fTag = aTag;
}

// Distance to track's ReferenceTrajectory
void DPhoton::setDtRT(const double aDtRT)
{
   fDtRT = aDtRT;
}

// set Error Matrix
// here is the place to calculate error matrix
void DPhoton::setErrorMatrix(const DVector3 aPosition, const DVector3 aVertex, const double aEnergy, DMatrixDSym* aSigmas )
{
   fErrorMatrix = aSigmas;
}


