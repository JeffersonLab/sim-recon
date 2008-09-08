//
//    File: DPhoton_factory.cc
// Created: Tue Aprl 17 11:57:50 EST 2007
// Creator: kornicer (on Linux stan)
//

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
        DELTA_THETA_CHARGE = 0.05; // Polar angle separation between photon and charged particle 
                                   // in radians
	USE_BCAL_ONLY = 0;
	USE_FCAL_ONLY = 0;
	
	gPARMS->SetDefaultParameter( "PID:DELTA_THETA_CHARGE", DELTA_THETA_CHARGE);
	gPARMS->SetDefaultParameter( "PID:USE_BCAL_ONLY", USE_BCAL_ONLY );
	gPARMS->SetDefaultParameter( "PID:USE_FCAL_ONLY", USE_FCAL_ONLY );
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

// Disable this info until tracking is fixed
//        vector<const DTrack*> tracks;
//	eventLoop->Get(tracks);
// and use thrown info within DPi0 to identify photons from charged particles
        vector<const DMCThrown*> thrown;
	eventLoop->Get(thrown);

// loop over FCAL photons    
	vector<const DFCALPhoton*> fcalPhotons;
	if( ! USE_BCAL_ONLY ) eventLoop->Get(fcalPhotons);
  
        JObject::oid_t nPhotons=0;
        for ( unsigned int i=0; i < fcalPhotons.size(); i++ ) {

		DPhoton *photon =  makeFCalPhoton(fcalPhotons[i], ++nPhotons);
                
//                double mdtrt = MinDistToRT(photon,tracks); 
//                photon->setDtRT(mdtrt); 
                double dTheta = dThetaToChargeMC(photon,thrown);
                
                photon->setdThetaCharge( dTheta ); 
                if (dTheta < DELTA_THETA_CHARGE ) photon->setTag(3);

		_data.push_back(photon);

        } 


// loop over BCAL photons and
// correct shower energy and position in makeBCalPhoton
	vector<const DBCALPhoton*> bcalPhotons;
	if( ! USE_FCAL_ONLY ) eventLoop->Get(bcalPhotons);
	
       for (unsigned int i=0; i< bcalPhotons.size(); i++) {

		DPhoton *photon =  makeBCalPhoton(bcalPhotons[i], ++nPhotons);

//                double mdtrt = MinDistToRT(photon,tracks); 
//	        photon->setDtRT(mdtrt); 
                double dTheta = dThetaToChargeMC(photon,thrown);
                photon->setdThetaCharge( dTheta ); 

                if (dTheta < DELTA_THETA_CHARGE ) photon->setTag(3);

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
        centroid.SetZ(centroid.Z()+65.);

        DVector3 momentum = gamma->getMom3();
        DVector3 vertex(0., 0., 65.);

        photon->setPosition( vertex );
        photon->setMomentum( momentum );
        photon->setPositionCal( centroid );
        photon->setMass( 0. );
        photon->setTag(1);

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


        DVector3 vertex(0.,0.,65.);
        photon->setPosition( vertex );

        DLorentzVector p = gamma->lorentzMomentum();
        photon->setMomentum( DVector3( p.X(), p.Y(), p.Z() ) );

        photon->setPositionCal( gamma->showerPosition() );
        photon->setMass( p.M() );
        photon->setTag(2);
        DMatrixDSym sigmas(7);

        sigmas(0,0) = pow( gamma->fitLayPointErr().X(), 2.0); // 
        sigmas(1,1) = pow( gamma->fitLayPointErr().Y(), 2.0); // 
        sigmas(2,2) = pow( gamma->fitLayPointErr().Z(), 2.0); // 

        sigmas[3][3] = 1.; // From Blake's simulation
        if (p.T()>=0) sigmas[3][3] = pow( 0.0445*sqrt(p.T()) + 0.009*p.T(), 2.0);

        sigmas[4][4] = pow( 0.5*TARGET_RADIUS, 2.0); // x_t, y_t
        sigmas[5][5] = pow( 0.5*TARGET_RADIUS, 2.0); // x_t, y_t
        sigmas[6][6] = pow( TARGET_LENGTH/sqrt(12.0), 2.0); // z_t

        photon->makeErrorMatrix( sigmas );
        return photon;
}


// OLD CODE KEEP FOR THE MOMENT
// for BCAL photons do energy correction here:
// The following definition is to chose normalization (A) and non-linear factor (B) extracted 
// from 10MeV hit-threshold simulation to implement z-dependent correction of shower energy.
// Disable to use shower Ecorr, which is shower energy multiplied 
// by just one calibration constant.
// #define ECORR_Z
DPhoton* DPhoton_factory::makeBCalPhoton(const DBCALShower* shower) 
{

        DPhoton* photon = new DPhoton;

        DVector3 vertex(0.,0.,65.);
        DVector3 showerCentroid(shower->x, shower->y, shower->z);        
        DVector3 r = showerCentroid - vertex;
#ifdef ECORR_Z
        double Z =  r.Z();
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
        float E = shower->E; 
#endif
        float f = E/sqrt(pow(r.X(),2) + pow(r.Y(),2) + pow(r.Z(),2));        
           
        photon->setMomentum( f*r );
        photon->setPosition( vertex );
        photon->setPositionCal( showerCentroid );
        photon->setMass( 0 );
        photon->setTag(1);

// make errorMatrix here ....
// BCAL shower factory performs linear fits of cell position (x,y,z)
// across all layers within a shower (x_i = a_x + c_x * r_i, i=layer index)
// to obtain coefficients (Ax, Ay, Az) and (Cx, Cy, Cz) and their respective errors.
// Thus, A is the impact point of the photon at the inner BCAL layer and 
// C gives the photon direction (C_x = P_x/|P|). 
// At this time, only impact point errors will be used to init errors on shower position in the BCAL.
        DMatrixDSym sigmas(7);
        sigmas(0,0) = shower->error_Apx_x; // 
        sigmas(1,1) = shower->error_Apx_y; // 
        sigmas(2,2) = shower->error_Apx_z; // 

        sigmas[3][3] = 1.; // From Blake's simulation
        if (E>=0) sigmas[3][3] = 0.0445*sqrt(E) + 0.009*E;

        sigmas[4][4] = 0.5*TARGET_RADIUS; // x_t, y_t
        sigmas[5][5] = 0.5*TARGET_RADIUS; // x_t, y_t
        sigmas[6][6]  = TARGET_LENGTH/sqrt(12.0); // z_t

        photon->makeErrorMatrix( sigmas );
        return photon;
}

// Loop over tracks and look up its reference trajectory, which has a method to 
// calculate the distance from its point to any given 3-vector.
// Return the distance from the closest track.
double DPhoton_factory::MinDistToRT(const DPhoton* photon, vector<const DTrack*> tracks) 
{

   double dmin = 10000.; // cm
   DVector3 photonPoint( photon->getPositionCal().X(), photon->getPositionCal().Y(), photon->getPositionCal().Z() );

   for (vector<const DTrack*>::const_iterator track  = tracks.begin(); 
    					      track != tracks.end(); 
					 			track++) {
	DReferenceTrajectory* rt = const_cast<DReferenceTrajectory*>((**track).rt);
        double dtrt = rt->DistToRT(photonPoint); 
        if (dtrt < dmin ) {
            dmin = dtrt;
        }
  }
  return dmin;
}

// Return the distance in polar anlge from the closest charged track.
double DPhoton_factory::dThetaToChargeMC(const DPhoton* photon, vector<const DMCThrown*> thrown) 
{

   double dmin = 1000.; 

   double theta = atan2( sqrt( pow(photon->momentum().X(),2) + pow(photon->momentum().Y(),2) ), photon->momentum().Z() );

   for (vector<const DMCThrown*>::const_iterator mc  = thrown.begin(); 
					      mc != thrown.end(); 
								mc++) {

        if ( (**mc).charge() == 0 ) continue;
	double deltaT = fabs( (**mc).momentum().Theta() - theta); 
        dmin = deltaT < dmin ? deltaT : dmin; 

  }

  return dmin;

}

