// DFCALPhoton member functions

#include "DPhoton.h"
#include <DMatrix.h>


DPhoton::DPhoton() : 
   fTag( kDefaultTag ),
   fTime( 0.),
   fDtRT( kDefaultDistance ), 
   fdThetaCharge( kDefaultDistance ), 
   fPositionCal( DVector3(0.0, 0.0, 0.0) )
{
   return;
}

DPhoton::DPhoton(const oid_t id) : DKinematicData(id),
   fTag( kDefaultTag ),
   fTime( 0. ),
   fDtRT( kDefaultDistance ),
   fdThetaCharge( kDefaultDistance ),
   fPositionCal( DVector3(0.0, 0.0, 0.0) )
{
   return;
}

DPhoton::~DPhoton()
{
}

// Set photon position
void DPhoton::setPositionCal( const DVector3& aPosition )
{
     fPositionCal = aPosition;
}

// Tag photon origin: 1/2/3 for FCAL/BCAL/Charge
void DPhoton::setTag( PhotonTag aTag )
{
   fTag = aTag;
}

void DPhoton::setTime( double aTime )
{
   // Time photon hit calorimeter
   fTime = aTime;
   
   // Time photon was at vertex. (n.b. this assumes fPosition
   // and fPositionCal are valid
   double d = (fPositionCal-position()).Mag();
   setT0(aTime - d/29.98, 2.0, SYS_NULL); /// FIXME!!
}

// Distance to track's ReferenceTrajectory
void DPhoton::setDtRT( double aDtRT )
{
   fDtRT = aDtRT;
}

// Polar angle distance to closest generated charge
void DPhoton::setdThetaCharge( double adTheta )
{
   fdThetaCharge = adTheta;
}

// A photon is described by momentum (p), position (r) and energy (E_g)
// (assuming that BCAL/FCAl cluster energy is calibrated E_g = E_c = E).
// Make photon errorMatrix (V) by rotating matrix of measured errors (V0) 
// in terms of vertex position (r_v), cluster position (r_c) and energy (E_c),
// Use 
//     V = A V0 A_transp
// where A is the matrix of first derivatives, with coefficients like 
//      pd_Px/pd_x_c ... 
// where pd_ stands for partial derivative.
// For example pd_Px/pd_x_c = E (r^2 - r_x^2) / R^3  = - pd_Px/pd_x_v
// with r = r_c - r_v and p = E r /R being vectors of photon position and momentum
// and R magnitude of r.
//
#define DELTA(i,j) ((i==j) ? 1 : 0)
void DPhoton::makeErrorMatrix( const DMatrixDSym& aSigmas )
{
   DVector3 r_c = getPositionCal();
   DVector3 r_v = position();
   DVector3 r = r_c - r_v;
   double R = r.Mag();
   double R2= r.Mag2();
   double R3 = R*R2;
   double E = energy();
   double f = E/R3;

// init and fill rotation matrix,
// first with momentum derivatives
   DMatrix A(7,7);
   for (int i = 0; i < 3; i++) {
	for ( int j = 0; j <3; j++) {
		A[i][j] = f*( R2*DELTA(i,j) - r(i)*r(j) );
                A[i][j+4] = - A[i][j];
 	}
  }


// fill energy part and remember: relation between energy and photon 
// position in calorimeter is neglected!
    A[3][3] = 1.;
    for (int j=0; j<3; j++) {
	A[j][3] = r(j)/R;
    } 

// fill spatial part where: dp_r_x/dp_x_c = - dp_r_x/dp_x_v ....
   for (int i = 0; i < 3; i++) {
	for ( int j = 0; j <3; j++) {
		int ik=i+4;
		int jk=j+4;
		A[ik][j] = DELTA(i,j);
                A[ik][jk] = - A[ik][j];
 	}
  }
  

   DMatrixDSym result = aSigmas; 
   
   result = result.Similarity(A); 

   setErrorMatrix( result );

}


