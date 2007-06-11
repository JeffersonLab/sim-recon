// DFCALPhoton member functions

#include "DPhoton.h"
#include <DMatrix.h>


DPhoton::DPhoton()
{
   fTag = 0; // default is FCAL
   fDtRT = 10000; // in a galaxy far far away...
}

DPhoton::~DPhoton()
{
}

// obsolite functions
/* Set photon momentum
void DPhoton::setMomentum(const DVector3 aMom)
{
     fMomentum = aMom;
}

// Set photon position
void DPhoton::setPosition(const DVector3 aPosition)
{
     fPosition = aPosition;
}

// Set photon vertex
void DPhoton::setVertex(const DVector3& aVertex)
{
     fVertex = aVertex;
}

 Set photon energy
void DPhoton::setEnergy(const double aEnergy)
{
     fEnergy = aEnergy;
}
*/
// Tag photon origin: 0/1 for FCAL/BCAL
void DPhoton::setTag(unsigned int aTag)
{
   fTag = aTag;
}

// Distance to track's ReferenceTrajectory
void DPhoton::setDtRT(double aDtRT)
{
   fDtRT = aDtRT;
}

// A poton is described by momentum (p), position (r) and energy (E_g)
// (assuming that BCAL/FCAl cluster energy is calibrated E_g = E_c = E).
// Make photon errorMatrix (V) by rotating matrix of measured errors (V0) 
// in terms of vartex position (r_v), cluster position (r_c) and enrgy (E_c),
// Use 
//     V = A V0 A_transp
// where A is the matrix of first derivatives, with coeficients like 
//      pd_Px/pd_x_c ... 
// where pd_ standas for partial derivative.
// For example pd_Px/pd_x_c = E (r^2 - r_x^2) / |r|^3
// with r = r_c - r_v and p = E r /|r| being vectors of photon position and momentum. 
//
void DPhoton::makeErrorMatrix( const DMatrixDSym& aSigmas )
{
//   DVector3 r_c = getPosition();
   DVector3 r = position();
   double E = energy();
// init and do nothing ....
   DMatrix A(7,1);
   DMatrix At(A);
   At.T();
   

   setErrorMatrix( aSigmas );

}


