// DFCALPhoton member functions

#include "DFCALPhoton.h"

DFCALPhoton::DFCALPhoton()
{
   fEnergy = 0;
   fMom3.SetY(0) ;
   fMom3.SetX(0) ;
   fMom3.SetZ(0) ;
   fMom4.SetXYZT(0.,0.,0.,0.);
}

DFCALPhoton::~DFCALPhoton()
{
}

#define	FCAL_RADIATION_LENGTH  3.1
#define	FCAL_CRITICAL_ENERGY  0.01455
#define	FCAL_SHOWER_OFFSET   3.0 

// These change with the FCAL attenuation length (L=166cm)
#define	EGAMMA_NORM 0.6348
#define	EGAMMA_EPSILON   0.029
 
// Simple non-linear correction
void DFCALPhoton::setEnergy(const double energy) 
{

	double const A=1/EGAMMA_NORM; 				// Normalization factor 
        double const power = 1/(1+EGAMMA_EPSILON); 	        // Non-linear factor
 	fEnergy = pow(A*energy,power);
}

// Simple depth correction: parameters imported from Radphi. 
void DFCALPhoton::setMom3(const double energy, const DVector3 pos) 
{

	double x = pos.X();
	double y = pos.Y();

//  estimate shower depth based on shower maximum
        double zMax = (FCAL_RADIATION_LENGTH*( 
			FCAL_SHOWER_OFFSET + log(energy/FCAL_CRITICAL_ENERGY)));

	double z = FCAL_Zmin - Shower_Vertex_Z + zMax;

// normalization factor to momenta [GeV]
        double f = energy/sqrt(pow(x,2)+pow(y,2)+pow(z,2)); 
        fMom3.SetX(x*f) ;
    	fMom3.SetY(y*f) ;
	fMom3.SetZ(z*f) ;
}

// Set photon four momentum
void DFCALPhoton::setMom4()
{
   fMom4.SetT(fEnergy);
   fMom4.SetVect(fMom3);
}


