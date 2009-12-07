//
//    File: DPhoton_factory.cc
// Created: Tue Aprl 17 11:57:50 EST 2007
// Creator: kornicer (on Linux stan)
// Last modified: Blake Leverington Mon Nov 23

#include <math.h>
#include <TLorentzVector.h>

#include "DPhoton.h"
#include "DPhoton_factory.h"
#include "JANA/JEvent.h"


//----------------
// Constructor
//----------------
DPhoton_factory::DPhoton_factory()
{
	// Set defaults

        PHOTON_VERTEX_X =  0.0;
        PHOTON_VERTEX_Y =  0.0;
        PHOTON_VERTEX_Z = 65.0;  

	DELTA_PHI_SWUMCHARGE = 0.15;// Azimuthal angle separation between photon and swumcharged particle 
                                   // in radians
	DELTA_Z_SWUMCHARGE = 40;//Position separation between photon and swumcharged particle 
                                   // in cm
	DELTA_R_SWUMCHARGE = 25;//Position separation between photon and swumcharged particle 
                                   // in cm
                                  
	USE_BCAL_ONLY = 0;
	USE_FCAL_ONLY = 0;
	
	gPARMS->SetDefaultParameter( "PID:USE_BCAL_ONLY", USE_BCAL_ONLY );
	gPARMS->SetDefaultParameter( "PID:USE_FCAL_ONLY", USE_FCAL_ONLY );

	gPARMS->SetDefaultParameter( "PID:DELTA_PHI_SWUMCHARGE", DELTA_PHI_SWUMCHARGE );
	gPARMS->SetDefaultParameter( "PID:DELTA_Z_SWUMCHARGE", DELTA_Z_SWUMCHARGE );
	gPARMS->SetDefaultParameter( "PID:DELTA_R_SWUMCHARGE", DELTA_R_SWUMCHARGE );

        gPARMS->SetDefaultParameter( "PID:PHOTON_VERTEX_X", PHOTON_VERTEX_X );
        gPARMS->SetDefaultParameter( "PID:PHOTON_VERTEX_Y", PHOTON_VERTEX_Y );
        gPARMS->SetDefaultParameter( "PID:PHOTON_VERTEX_Z", PHOTON_VERTEX_Z );

}


//------------------
// evnt
// PID Photons: FCAL photons are copied with 'tag' attribute set to zero.
//		BCAL photons are made from BCAL showers after applying
//              x-dependent energy corrections (right now based on 10 MeV
//		hit threshold in HDGent. (MK)
//------------------
jerror_t DPhoton_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{

	vector<const DParticle*> chargedswum;
	eventLoop->Get(chargedswum);
// loop over FCAL photons    
	vector<const DFCALPhoton*> fcalPhotons;
	if( ! USE_BCAL_ONLY ) eventLoop->Get(fcalPhotons);
  
        JObject::oid_t nPhotons=0;
        for ( unsigned int i=0; i < fcalPhotons.size(); i++ ) {

		DPhoton *photon =  makeFCalPhoton(fcalPhotons[i], ++nPhotons);
                
		vector<double> dSwum;
                dSwum = dFromSwumChargeMC(photon,chargedswum);

		if (dSwum[0] <DELTA_PHI_SWUMCHARGE  && dSwum[1] <DELTA_Z_SWUMCHARGE && dSwum[2] <DELTA_R_SWUMCHARGE  ) photon->setTag( DPhoton::kCharge   );

		_data.push_back(photon);

        } 

// loop over BCAL photons and
// correct shower energy and position in makeBCalPhoton
	vector<const DBCALPhoton*> bcalPhotons;
	if( ! USE_FCAL_ONLY ) eventLoop->Get(bcalPhotons);
	
       for (unsigned int i=0; i< bcalPhotons.size(); i++) {

		DPhoton *photon =  makeBCalPhoton(bcalPhotons[i], ++nPhotons);

		vector<double> dSwum;
                dSwum = dFromSwumChargeMC(photon,chargedswum);

		if (dSwum[0] <DELTA_PHI_SWUMCHARGE  && dSwum[1] <DELTA_Z_SWUMCHARGE && dSwum[2] <DELTA_R_SWUMCHARGE  ) photon->setTag( DPhoton::kCharge   );

		_data.push_back(photon);

       } 

	return NOERROR;
}


// at this time just copy data from DFCALPhoton and set 'tag' to zero
// final non-linear and depth corrections can be applied here
#define FCAL_BLOCK_WIDTH 4
#define TARGET_RADIUS 1.5
#define TARGET_LENGTH 30.
DPhoton* DPhoton_factory::makeFCalPhoton(const DFCALPhoton* gamma, const JObject::oid_t id) 
{

        DPhoton* photon = new DPhoton( id );
        
        double energy = gamma->getEnergy();

// default vertex is (0,0,65) and this has been taken into account in FCAL libraries
// during MC calibration, make sure FCALPhoton returns centroid position in 
// GlueX coordinates in the future..., in case we need it
        DVector3 centroid = gamma->getPosition();  

        DVector3 vertex( PHOTON_VERTEX_X, PHOTON_VERTEX_Y, PHOTON_VERTEX_Z);
        DVector3 photonVector = centroid - vertex;

        double scale = energy
                     / sqrt( photonVector.X()*photonVector.X() 
                           + photonVector.Y()*photonVector.Y() 
                           + photonVector.Z()*photonVector.Z() );

        photon->setMomentum( DVector3( photonVector.X()*scale, 
                                       photonVector.Y()*scale, 
                                       photonVector.Z()*scale ) );

        photon->setPosition( vertex );
        photon->setPositionCal( centroid );
        photon->setMass( 0. );
        photon->setTag( DPhoton::kFcal );
        photon->setTime( gamma->getTime() );

// create the simplest error matrix:
// At this point, it is assumed that error matrix of measured quantities is diagonal,
// with elements like: sigma_Z_t = L/sqrt(12) sigma_X_t = sigma_Y_t = r0/2 
// L=target length, r0 = target radius...
// This means that energy-depth-polar angle relation  is neglected.
// the order of sigmas is:  x_c, y_c, z_c, E, x_t, y_t, z_t
        DMatrixDSym sigmas(7);

//        sigmas[0][0] = pow( FCAL_BLOCK_WIDTH/sqrt(12.0), 2.0 ); // x_c, y_c
//        sigmas[1][1] = pow( FCAL_BLOCK_WIDTH/sqrt(12.0), 2.0 ); // x_c, y_c
        sigmas[0][0] = pow( 0.7 , 2.0 ); // x_c, y_c
        sigmas[1][1] = pow( 0.7 , 2.0 ); // x_c, y_c

//        sigmas[2][2] = pow( 2.54, 2.0); //  z_c = rms of average depth for photons from 0-5 GeV
        sigmas[2][2] = pow( gamma->getPositionError().Z() , 2.0); 

        sigmas[3][3] =  (energy >= 0)  ? pow( 0.042*sqrt(energy) + 0.0001, 2.0 ) : 1e-6 ;

        sigmas[4][4] = pow( 0.5*TARGET_RADIUS, 2.0) ; // x_t, y_t
        sigmas[5][5] = pow( 0.5*TARGET_RADIUS, 2.0) ; // x_t, y_t
        sigmas[6][6] = pow( TARGET_LENGTH/sqrt(12.0), 2.0) ; // z_t

        photon->makeErrorMatrix( sigmas );
        return photon;
}

// Copy data from BCALPhoton
DPhoton* DPhoton_factory::makeBCalPhoton(const DBCALPhoton* gamma, const JObject::oid_t id) {

        DPhoton* photon = new DPhoton( id );

        double energy = gamma->lorentzMomentum().Energy();
        DVector3 centroid = gamma->showerPosition();
        DVector3 vertex(PHOTON_VERTEX_X, PHOTON_VERTEX_Y, PHOTON_VERTEX_Z);

        DVector3 photonVector = centroid - vertex;

        double scale = energy
                     / sqrt( photonVector.X()*photonVector.X() 
                           + photonVector.Y()*photonVector.Y() 
                           + photonVector.Z()*photonVector.Z() );

        photon->setMomentum( DVector3( photonVector.X()*scale, 
                                       photonVector.Y()*scale, 
                                       photonVector.Z()*scale) ) ;

        photon->setPosition( vertex );
        photon->setTime( gamma->showerTime() );

        DLorentzVector p = gamma->lorentzMomentum();
        photon->setMomentum( DVector3( p.X(), p.Y(), p.Z() ) );

        photon->setPositionCal( gamma->showerPosition() );
        photon->setMass( 0. );
        photon->setTag( DPhoton::kBcal );
        DMatrixDSym sigmas(7);

        sigmas(0,0) = pow( gamma->fitLayPointErr().X(), 2.0); // 
        sigmas(1,1) = pow( gamma->fitLayPointErr().Y(), 2.0); // 
        sigmas(2,2) = pow( gamma->fitLayPointErr().Z(), 2.0); // 

        sigmas[3][3] = 1.; // From Blake's simulation
        if ( energy>=0 ) sigmas[3][3] = pow( 0.0445*sqrt( energy ) + 0.009*energy, 2.0);

        sigmas[4][4] = pow( 0.5*TARGET_RADIUS, 2.0); // x_t, y_t
        sigmas[5][5] = pow( 0.5*TARGET_RADIUS, 2.0); // x_t, y_t
        sigmas[6][6] = pow( TARGET_LENGTH/sqrt(12.0), 2.0); // z_t

        photon->makeErrorMatrix( sigmas );
        return photon;
}

// Return the distance in azimuthal angle and position Z from the closest charged track.
vector<double>  
DPhoton_factory::dFromSwumChargeMC(const DPhoton* photon, vector<const DParticle*>  chargedswum) 
{

 DVector3 photonPoint =photon->getPositionCal();
 DVector3 diffVect,diffVectbcal,diffVectfcal;
 double dPhi = 10.0;
 double dPhiMin = 10.0;
 double dZ = 1000.0;
 double dZMin = 1000.0;
 double dR = 1000.0;
 double dRMin = 1000.0;
 vector<double> diffSwum(3);

 DApplication* dapp = dynamic_cast<DApplication*>(eventLoop->GetJApplication());
   if(!dapp){
     _DBG_<<"Cannot get DApplication from JEventLoop! (are you using a JApplication based program?)"<<endl;
     //   return 0;
   }
   DMagneticFieldMap *bfield = dapp->GetBfield();


 for( vector<const DParticle*>::const_iterator swum = chargedswum.begin();
	     swum != chargedswum.end(); ++swum ){
   if ( (**swum).charge() == 0 ) continue;

 

   bool hitbcal,hitfcal = false;
   double q = (**swum).charge(); 
   DVector3 pos = (**swum).position();
   DVector3 mom = (**swum).momentum();
   DMagneticFieldStepper *stepper1 = new DMagneticFieldStepper(bfield, q, &pos, &mom);
   DMagneticFieldStepper *stepper2 = new DMagneticFieldStepper(bfield, q, &pos, &mom);

   DVector3 pos_bcal = pos; 
   DVector3 mom_bcal = mom; 
   DVector3 pos_fcal = pos;
   DVector3 mom_fcal = mom;
   DVector3 origin(0.0, 0.0, 643.2);
   DVector3 norm(0.0, 0.0, 1.0);


   bool swimrad = stepper1->SwimToRadius(pos_bcal, mom_bcal, 65.0);

   if( swimrad || (pos_bcal.Z()>400 && !swimrad) ){
     bool swimplane = stepper2->SwimToPlane(pos_fcal, mom_fcal, origin, norm );//save a little time and do this here
     hitbcal = false;
   if( swimplane ){
     hitfcal = false;
   }else{
     hitfcal = true;
     diffVectfcal = photonPoint - pos_fcal;
   }
   }else{
     hitbcal = true;
     diffVectbcal = photonPoint - pos_bcal;
   }

   if( hitbcal && !hitfcal ){
     dPhi = photonPoint.Phi() - pos_bcal.Phi();
     dZ = photonPoint.Z() - pos_bcal.Z();
     dR = photonPoint.Perp() - pos_bcal.Perp();

    }else if(!hitbcal && hitfcal){
     dPhi = fabs(photonPoint.Phi() - pos_fcal.Phi());
     dZ = fabs(photonPoint.Z() - pos_fcal.Z());
     dR = photonPoint.Perp() - pos_fcal.Perp();

  }
   //assigns the minimum distance
 dPhiMin = dPhi < dPhiMin ? dPhi : dPhiMin; 
 dZMin = dZ < dZMin ? dZ : dZMin;  
 dRMin = dR < dRMin ? dR : dRMin; 
 }

 diffSwum[0] = dPhiMin;
 diffSwum[1] = dZMin; 
 diffSwum[2] = dRMin;

  return diffSwum;

}
