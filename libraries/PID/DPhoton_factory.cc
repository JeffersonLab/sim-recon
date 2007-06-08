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
		printcol("%5.2f", phot->getEnergy());
		printcol("%5.2f", phot->getMomentum().X());
		printcol("%5.2f", phot->getMomentum().Y());
		printcol("%5.2f", phot->getMomentum().Z());
		printcol("%7.2f", phot->getPosition().X());
		printcol("%7.2f", phot->getPosition().Y());
		printcol("%7.2f", phot->getPosition().Z());
		printcol("%5i", phot->getTag());
		printrow();
	}

	return _table;
}


// at this time just copy data from DFCALPhoton and set 'tag' to zero
//DPhoton* DPhoton_factory::makeFCalPhoton(const TLorentzVector gamma) 
DPhoton* DPhoton_factory::makeFCalPhoton(const DFCALPhoton* gamma) 
{

        DPhoton* photon = new DPhoton;
        
        photon->setPosition(gamma->getPosition());
        photon->setMomentum(gamma->getMom3());
        photon->setEnergy(gamma->getEnergy());
// default vertex is (0,0,65) and this has been taken into account in FCAL libraries
//      photon->setVertex(0., 0. , 0); 
        photon->setTag(0);

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

        float x = shower->x;
        float y = shower->y;
        float z = shower->z - 65;
#ifdef ECORR_Z
           float A0 = 0.87971;
           float A1 = 134.52;
           float A2 = 1.822E-6;
           float B0 = 0.0585;
           float B1 = 199.9;
           float B2 = 7.784E-7;
           float A = A0 - A2*(z-A1)*(z-A1);
           float B = B0 + B2*(z-B1)*(z-B1);
        float E = pow((double)(shower->E/A),(double)(1/(1+B))); 
#else
        float E = shower->Ecorr; 
#endif
        float f = E/sqrt(pow(x,2) + pow(y,2) + pow(z,2));        

        DVector3 aMom(x*f, y*f, z*f);        
        DVector3 aPos(x, y, z);        
           
        photon->setMomentum( aMom );
        photon->setPosition( aPos );
        photon->setEnergy( E );
        photon->setTag(1);
        //photon->setVertex(0., 0. , 0);
        return photon;
}

// loop over tracks and look at the distance of its reference point from the photon
// return the distance from the closest track
double DPhoton_factory::MinDistToRT(const DPhoton* photon, vector<const DTrack*> tracks) 
{

   double dmin = 1000.; // cm
   DVector3 photonPosition( photon->getPosition().X(), photon->getPosition().Y(), photon->getPosition().Z() );

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

