//
//    File: DPhoton_factory.cc
// Created: Tue Aprl 17 11:57:50 EST 2007
// Creator: kornicer (on Linux stan)
//

#include <math.h>
#include <TLorentzVector.h>

#include "FCAL/DFCALPhoton.h"
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
	
       for (vector<const DFCALPhoton*>::const_iterator gamma  = fcalPhotons.begin(); 
                                                        gamma != fcalPhotons.end(); 
							gamma++) {

             
		DPhoton *photon =  makeFCalPhoton((**gamma).getMom4());
                
                double mdtrt = MinDistToRT(photon,tracks); 
                photon->setDtRT(mdtrt); 

		_data.push_back(photon);

       } 


// loop over BCAL photons and
// correct shower energy and position in makeBCalPhoton
	vector<const DBCALShower*> bcalShowers;
	eventLoop->Get(bcalShowers);
	
       for (vector<const DBCALShower*>::const_iterator shower  = bcalShowers.begin(); 
                                                        shower != bcalShowers.end(); 
							shower++) {

#ifdef ECORR_Z
                TLorentzVector bcalShower((**shower).x,(**shower).y,(**shower).z,(**shower).E);
#else
                TLorentzVector bcalShower((**shower).x,(**shower).y,(**shower).z,(**shower).Ecorr);
#endif
		DPhoton *photon =  makeBCalPhoton(bcalShower);

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

	printheader("row:   E(GeV):   Px(GeV):	  Py(GeV):    Pz(GeV):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DPhoton *phot = _data[i];
               
		printnewrow();
		printcol("%d",	i);
		printcol("%5.2f", phot->getMom4().T());
		printcol("%5.2f", phot->getMom4().X());
		printcol("%5.2f", phot->getMom4().Y());
		printcol("%5.2f", phot->getMom4().Z());
		printrow();
	}

	return _table;
}


// at this time just copy data from DFCALPhoton and set 'tag' to zero
DPhoton* DPhoton_factory::makeFCalPhoton(const TLorentzVector gamma) 
{

        DPhoton* photon = new DPhoton;
        
        photon->setMom4(gamma);
        photon->setTag(0);

        return photon;
}

// for BCAL photons do energy correction here:
DPhoton* DPhoton_factory::makeBCalPhoton(const TLorentzVector shower) 
{

        DPhoton* photon = new DPhoton;

        float z = shower.Z() - 65;
#ifdef ECORR_Z
// normalization (A) and non-linear factor (B) extracted from 10MeV 
// hit-threshold simulation
           float A0 = 0.87971;
           float A1 = 134.52;
           float A2 = 1.822E-6;
           float B0 = 0.0585;
           float B1 = 199.9;
           float B2 = 7.784E-7;
           float A = A0 - A2*(z-A1)*(z-A1);
           float B = B0 + B2*(z-B1)*(z-B1);
        float E = pow((double)(shower.E()/A),(double)(1/(1+B))); 
#else
        float E = shower.E(); 
#endif
        float f = E/sqrt(pow(shower.X(),2)+pow(shower.Y(),2)+pow(z,2));        
        TLorentzVector gamma(shower.X()*f,shower.Y()*f,z*f,E);        
        photon->setMom4(gamma);
        photon->setTag(1);

        return photon;
}

// loop over tracks and find the distance of its reference point from the photon
// return photon distance from the closest track
double DPhoton_factory::MinDistToRT(const DPhoton* photon, vector<const DTrack*> tracks) 
{

   double dmin = 1000.; // cm
// shange this to photnPostion, not momentum!!!!
   DVector3 photonPosition( photon->getMom4().X(), photon->getMom4().Y(), photon->getMom4().Z() );

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

