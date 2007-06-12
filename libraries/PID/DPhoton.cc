// DFCALPhoton member functions

#include "DPhoton.h"
#include <DMatrix.h>


DPhoton::DPhoton()
{
   fTag = 0; // default is FCAL
   fDtRT = 10000; // in a galaxy far far away...
   fPositionCal.SetXYZ(0.,0.,0.);
}

DPhoton::~DPhoton()
{
}

// Set photon position
void DPhoton::setPositionCal(const DVector3& aPosition)
{
     fPositionCal = aPosition;
}

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

// init to zeros
   DMatrix A(7,7);

// fill momentum derivatives
   for (int i = 0; i < 3; i++) {
	for ( int j = 0; j <3; j++) {
		
		A[i][j] = f*( R2*DELTA(i,j) - r(i)*r(j) );
                A[j][i] = A[i][j];
                A[i][j+4] = - A[i][j];

 	}
  }

// fill energy part and remember: relation between energy and photon 
// position in calorimeter is neglected!
    A[3][3] = 1.;
    for (int j=0; j<3; j++) {
	A[3][j] = r(j)/R;
    } 

// fill spatial part where: dp_r_x/dp_x_c = - dp_r_x/dp_x_v ....
   for (int i = 4; i < 7; i++) {
	for ( int j = 0; j <3; j++) {
		int k=j+4;
		A[i][j] = DELTA(i,k);
                A[j][i] = A[i][j];
                A[i][j+4] = - A[i][j];

 	}
  }

   DMatrixDSym result = aSigmas; 
//   result = result.Similarity(A); 
   
   setErrorMatrix( result.Similarity(A) );

}


