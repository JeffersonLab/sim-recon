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
#define	EGAMMA_NORM 0.627419
#define	EGAMMA_EPSILON   0.05259
 
// Simple non-linear correction
void DFCALPhoton::fixEnergy(const double energy) 
{

	double const A=1/EGAMMA_NORM; 				// Normalization factor 
        double const power = 1/(1+EGAMMA_EPSILON); 	        // Non-linear factor

}

// Simple depth correction: estimate shower depth based on shower maximum
// 			    parameters imported from Radphi. 
void DFCALPhoton::fixDepth(const double energy, const DVector3 pos) {

        double zMax = (FCAL_RADIATION_LENGTH*( 
			FCAL_SHOWER_OFFSET + log(energy/FCAL_CRITICAL_ENERGY)));

	double z = FCAL_Zmin - Shower_Vertex_Z + zMax;

        fPosition.SetXYZ(pos.X(), pos.Y(), z);

}

// Set photon momentum:
// make sure that energy and shower depth are already taken care off.
void DFCALPhoton::setMom3(const double energy, const DVector3 pos) 
{

	double x = pos.X();
	double y = pos.Y();
	double z = pos.Z();

// normalization factor to momenta [GeV]
        double f = energy/sqrt(pow(x,2)+pow(y,2)+pow(z,2)); 
        fMom3.SetXYZ(x*f, y*f, z*f);
}

// Set photon four momentum
void DFCALPhoton::setMom4()
{
   fMom4.SetT(fEnergy);
   fMom4.SetVect(fMom3);
}


