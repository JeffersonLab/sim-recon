// DFCALPhoton member functions

#include "DFCALPhoton.h"

DFCALPhoton::DFCALPhoton()
{
   fEnergy = 0;
   fTime = 0;
   fMom3.SetY(0) ;
   fMom3.SetX(0) ;
   fMom3.SetZ(0) ;
   fMom4.SetXYZT(0.,0.,0.,0.);
}

DFCALPhoton::~DFCALPhoton()
{
}

//
 
void DFCALPhoton::setEnergy(const double energy) 
{

        fEnergy = energy;

}

void DFCALPhoton::setTime(const double time) 
{

        fTime = time;

}

void DFCALPhoton::setPosition(const DVector3 aPosition ) 
{

        fPosition = aPosition;

}

// Set position errors
void DFCALPhoton::setPosError(const double aXerr, const double aYerr, const double aZerr) {

        fPositionError.SetX( aXerr);
        fPositionError.SetY( aYerr);
        fPositionError.SetZ( aZerr);

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


