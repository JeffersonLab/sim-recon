//
//    File: DPhoton_factory.cc
// Created: Tue Aprl 17 11:57:50 EST 2007
// Creator: kornicer (on Linux stan)
//

#include <math.h>
#include <TLorentzVector.h>

//#include "FCAL/DFCALPhoton.h"
#include "BCAL/DBCALShower.h"
//#include "TRACKING/DTrack.h"
//#include "TRACKING/DReferenceTrajectory.h"
#include "DPhoton.h"
#include "DPhoton_factory.h"
#include "JANA/JEvent.h"


//----------------
// Constructor
//----------------
DPhoton_factory::DPhoton_factory()
{
	// Set defaults
	
}


//------------------
// evnt
// PID Photons: FCAL photons are copied with 'tag' atribute set to zero.
//		BCAL photons are made from BCAl showers after applying
//              x-dependant energy corrections (right now based on 10 MeV
//		hit threshold in HDGent. (MK)
//------------------
jerror_t DPhoton_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{

        vector<const DTrack*> tracks;
	eventLoop->Get(tracks);

// loop over FCAL photons    
	vector<const DFCALPhoton*> fcalPhotons;
	eventLoop->Get(fcalPhotons);
	
        for ( unsigned int i=0; i < fcalPhotons.size(); i++ ) {

                const DFCALPhoton *gamma = fcalPhotons[i];
		DPhoton *photon =  makeFCalPhoton( gamma );
                
                double mdtrt = MinDistToRT(photon,tracks); 
                photon->setDtRT(mdtrt); 

		_data.push_back(photon);

        } 


// loop over BCAL photons and
// correct shower energy and position in makeBCalPhoton
	vector<const DBCALShower*> bcalShowers;
	eventLoop->Get(bcalShowers);
	
       for (unsigned int i=0; i< bcalShowers.size(); i++) {

		
/*#ifdef ECORR_Z
                TLorentzVector bcalShower((**shower).x,(**shower).y,(**shower).z,(**shower).E);
#else
                TLorentzVector bcalShower((**shower).x,(**shower).y,(**shower).z,(**shower).Ecorr);
#endif*/
                
		DPhoton *photon =  makeBCalPhoton(bcalShowers[i]);

                double mdtrt = MinDistToRT(photon,tracks); 
                photon->setDtRT(mdtrt); 

		_data.push_back(photon);

       } 

	return NOERROR;
}

//------------------
// toString
//------------------
const string DPhoton_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row:   E(GeV):   Px(GeV):	  Py(GeV):    Pz(GeV):    X(cm):    Y(cm):    Z(cm):    Tag: ");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DPhoton *phot = _data[i];
               
		printnewrow();
		printcol("%d",	i);
		printcol("%5.2f", phot->energy());
		printcol("%5.2f", phot->momentum().X());
		printcol("%5.2f", phot->momentum().Y());
		printcol("%5.2f", phot->momentum().Z());
		printcol("%7.2f", phot->position().X());
		printcol("%7.2f", phot->position().Y());
		printcol("%7.2f", phot->position().Z());
		printcol("%5i", phot->getTag());
		printrow();
	}

	return _table;
}


// at this time just copy data from DFCALPhoton and set 'tag' to zero
// final non-linear and depth corrections can be applied here
#define FCAL_BLOCK_WIDTH 4
#define TARGET_RADIUS 1.5
#define TARGET_LENGTH 30
DPhoton* DPhoton_factory::makeFCalPhoton(const DFCALPhoton* gamma) 
{

        DPhoton* photon = new DPhoton;
        
        double energy = gamma->getEnergy();
// default vertex is (0,0,65) and this has been taken into account in FCAL libraries
// during MC calibration, make sure FCALPhoton returns centroid position in 
// GlueX coordinates in the future..., in case we need it
//        DVector3 centroid = gamma->getPosition();  
//        centroid.SetZ(centroid.Z()+65.);

        DVector3 momentum = gamma->getMom3();
        DVector3 vertex(0., 0., 65.);

        photon->setPosition( vertex );
        photon->setMomentum( momentum );
        photon->setMass( 0. );
        photon->setTag(0);

// create the simplest error matrix:
// At this point, it is assumed that error matrix of measured quantites is diagonal,
// with elements like: sigma_Z_t = L/sqrt(12) sigma_X_t = sigma_Y_t = r0/2 
// L=target lenght, r0 = targer radius...
// This means that energy-depth-polar angle relation  is neglected.
// the order of sigmas is:  x_c, y_c, z_c, E, x_t, y_t, z_t
        DMatrixDSym sigmas(7,0);
        sigmas(1,1) = sigmas(2,2) = FCAL_BLOCK_WIDTH/sqrt(12); // x_c, y_c
        sigmas(3,3) = 2.54; //  z_c = rms of average depth for photons from 0-5 GeV

        sigmas(4,4) = 1.; // right now energy is 4.2%/sqrt(E)
        if (energy>0) sigmas(4,4) = 0.042/sqrt(energy);

        sigmas(5,5) = sigmas(6,6) = 0.5*TARGET_RADIUS; // x_t, y_t
        sigmas(7,7)  = TARGET_LENGTH/sqrt(12); // z_t

        photon->makeErrorMatrix( sigmas );

        return photon;
}

// for BCAL photons do energy correction here:
// The following definition is to chose normalization (A) and non-linear factor (B) extracted 
// from 10MeV hit-threshold simulation to implement z-dependent correction of shower energy.
// Disable to use shower Ecorr, which is shower energy multiplied 
// by just one callibration constant.
#define ECORR_Z
DPhoton* DPhoton_factory::makeBCalPhoton(const DBCALShower* shower) 
{

        DPhoton* photon = new DPhoton;

        DVector3 vertex(0.,0.,65.);
        DVector3 showerPosition(shower->x, shower->y, shower->z);        
        DVector3 position = showerPosition - vertex;
#ifdef ECORR_Z
        double Z =  position.Z();
           float A0 = 0.87971;
           float A1 = 134.52;
           float A2 = 1.822E-6;
           float B0 = 0.0585;
           float B1 = 199.9;
           float B2 = 7.784E-7;
           float A = A0 - A2*(Z-A1)*(Z-A1);
           float B = B0 + B2*(Z-B1)*(Z-B1);
        float E = pow((double)(shower->E/A),(double)(1/(1+B))); 
#else
        float E = shower->Ecorr; 
#endif
        float f = E/sqrt(pow(position.X(),2) + pow(position.Y(),2) + pow(position.Z(),2));        
        DVector3 momentum = f*position;        
           
        photon->setMomentum( f*position );
        photon->setPosition( vertex );
        photon->setMass( 0 );
        photon->setTag(1);
        return photon;
}

// Loop over tracks and look up its reference trajectory, which has a method to 
// calculate the distance from its point to any given 3-vector.
// Return the distance from the closest track.
double DPhoton_factory::MinDistToRT(const DPhoton* photon, vector<const DTrack*> tracks) 
{

   double dmin = 10000.; // cm
   DVector3 photonPosition( photon->position().X(), photon->position().Y(), photon->position().Z() );

   for (vector<const DTrack*>::const_iterator track  = tracks.begin(); 
    					      track != tracks.end(); 
					 			track++) {
	DReferenceTrajectory* rt = const_cast<DReferenceTrajectory*>((**track).rt);
        double dtrt = rt->DistToRT(photonPosition); 
        if (dtrt < dmin ) {
            dmin = dtrt;
        }
  }
  return dmin;
}

