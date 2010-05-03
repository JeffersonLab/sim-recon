//
//    File: DPhoton_factory.cc
// Created: Tue Aprl 17 11:57:50 EST 2007
// Creator: kornicer (on Linux stan)
// Last modified: Blake Leverington Mon Nov 23

#include <math.h>
#include <DLorentzVector.h>

#include "DPhoton.h"
#include "DPhoton_factory.h"
#include "JANA/JEvent.h"


//----------------
// Constructor
//----------------
DPhoton_factory::DPhoton_factory()
{
	// Initialize data members
	bfield = NULL;
	stepper = NULL;

	// Set defaults for configuration parameters
	USE_BCAL_ONLY = 0;
	USE_FCAL_ONLY = 0;
	PHOTON_VERTEX_X =  0.0;
	PHOTON_VERTEX_Y =  0.0;
	PHOTON_VERTEX_Z = 65.0;  

	gPARMS->SetDefaultParameter( "PID:USE_BCAL_ONLY", USE_BCAL_ONLY );
	gPARMS->SetDefaultParameter( "PID:USE_FCAL_ONLY", USE_FCAL_ONLY );
	gPARMS->SetDefaultParameter( "PID:PHOTON_VERTEX_X", PHOTON_VERTEX_X );
	gPARMS->SetDefaultParameter( "PID:PHOTON_VERTEX_Y", PHOTON_VERTEX_Y );
	gPARMS->SetDefaultParameter( "PID:PHOTON_VERTEX_Z", PHOTON_VERTEX_Z );
}

//----------------
// brun
//----------------
jerror_t DPhoton_factory::brun(JEventLoop *loop, int runnumber)
{
	// Get Magnetic field map so we can create magnetic field stepper
	DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	if(!dapp){
		_DBG_<<"Cannot get DApplication from JEventLoop! (are you using a JApplication based program?)"<<endl;
		return RESOURCE_UNAVAILABLE;
	}
	bfield = dapp->GetBfield();
	stepper = new DMagneticFieldStepper(bfield);
	
	// Get z-position of front face of FCAL from geometry
	vector<double> fcal_pos;
	loop->GetJGeometry()->Get("//posXYZ[@volume='ForwardEMcal']/@X_Y_Z", fcal_pos, " ");
	assert(fcal_pos.size()==3);
	fcal_origin.SetXYZ(0.0, 0.0, fcal_pos[2]);
	fcal_norm.SetXYZ(0.0, 0.0, 1.0);
	
	// Get calibration constants
	map<string, double> photon_track_matching;
	loop->GetCalib("PID/photon_track_matching", photon_track_matching);
	  DELTA_R_FCAL = photon_track_matching["DELTA_R_FCAL"];
	   MEAN_R_FCAL = photon_track_matching["MEAN_R_FCAL"];
	  DELTA_R_BCAL = photon_track_matching["DELTA_R_BCAL"];
	   MEAN_R_BCAL = photon_track_matching["MEAN_R_BCAL"];

	// Optionally notify user of values
	if(debug_level>0){
		cout<<"PID: Photon/Charged track matching parameters"<<endl;
		cout<<"     DELTA_R_FCAL = "<<photon_track_matching["DELTA_R_FCAL"]<<endl;
		cout<<"      MEAN_R_FCAL = "<<photon_track_matching["MEAN_R_FCAL"]<<endl;
		cout<<"     DELTA_R_BCAL = "<<photon_track_matching["DELTA_R_BCAL"]<<endl;
		cout<<"      MEAN_R_BCAL = "<<photon_track_matching["MEAN_R_BCAL"]<<endl;
		cout<<endl;
		cout<<"FCAL front face is at z="<<fcal_origin.Z()<<endl;
	}

	return NOERROR;
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
	// Get reconstructed photons from both BCAL and FCAL    
	vector<const DFCALPhoton*> fcalPhotons;
	vector<const DBCALPhoton*> bcalPhotons;
	if( ! USE_BCAL_ONLY ) eventLoop->Get(fcalPhotons);
	if( ! USE_FCAL_ONLY ) eventLoop->Get(bcalPhotons);
	
	// If no reconstructed photons exist, we can return now
	if(fcalPhotons.size()==0 && bcalPhotons.size()==0)return NOERROR;

	// Get charged tracks
	vector<const DTrackTimeBased*> tracks;
	eventLoop->Get(tracks);

#if 0	
	// Project charged tracks to FCAL and BCAL if needed
	vector<DVector3> track_projection_fcal;
	vector<DVector3> track_projection_bcal;
	if(fcalPhotons.size()!=0)ProjectToFCAL(tracks, track_projection_fcal);
	if(bcalPhotons.size()!=0)ProjectToBCAL(tracks, track_projection_bcal);
#endif

	// Loop over reconstructed FCAL photons
	for ( unsigned int i=0; i < fcalPhotons.size(); i++ ) {
		DPhoton *photon =  makeFCalPhoton(fcalPhotons[i], (JObject::oid_t)i);
		DVector3 photon_pos = photon->getPositionCal();

		// Loop over charged tracks to see if this matches any
		for(unsigned int j=0; j<tracks.size(); j++){
			const DReferenceTrajectory *rt = tracks[j]->rt;
			double dist = rt->DistToRT(photon_pos);

			if(fabs(dist - MEAN_R_FCAL) < DELTA_R_FCAL){
				photon->setTag( DPhoton::kCharge);
				break;
			}
		}

		// Add this FCAL photon to list of all reconstructed photons
		_data.push_back(photon);
	} 

	// Loop over reconstructed BCAL photons
	for ( unsigned int i=0; i < bcalPhotons.size(); i++ ) {
		DPhoton *photon =  makeBCalPhoton(bcalPhotons[i], (JObject::oid_t)i);
		DVector3 photon_pos = photon->getPositionCal();

		// Loop over charged tracks to see if this matches any
		for(unsigned int j=0; j<tracks.size(); j++){
			const DReferenceTrajectory *rt = tracks[j]->rt;
			double dist = rt->DistToRT(photon_pos);

			if(fabs(dist - MEAN_R_BCAL) < DELTA_R_BCAL){
				photon->setTag( DPhoton::kCharge);
				break;
			}
		}

		// Add this BCAL photon to list of all reconstructed photons
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

#if 0
//--------------------------
// ProjectToFCAL
//--------------------------
void DPhoton_factory::ProjectToFCAL(vector<const DTrackTimeBased*> &tracks, vector<DVector3> &track_projection)
{
	// Loop over charged tracks projecting each to the FCAL front surface
	for(unsigned int i=0; i<tracks.size(); i++){
		const DTrackTimeBased *track = tracks[i];
		if(track->charge()==0)continue;
		DVector3 pos = track->position();
		DVector3 mom = track->momentum();
		stepper->SetCharge(track->charge());
		if(!stepper->SwimToPlane(pos, mom, fcal_origin, fcal_norm )){
			// Stepper managed to swim to the FCAL! Record position at FCAL front face
			track_projection.push_back(pos);
		}
	}
}

//--------------------------
// ProjectToBCAL
//--------------------------
void DPhoton_factory::ProjectToBCAL(vector<const DTrackTimeBased*> &tracks, vector<DVector3> &track_projection)
{
	// Loop over charged tracks projecting each to the BCAL inner surface
	for(unsigned int i=0; i<tracks.size(); i++){
		const DTrackTimeBased *track = tracks[i];
		if(track->charge()==0)continue;
		DVector3 pos = track->position();
		DVector3 mom = track->momentum();
		stepper->SetCharge(track->charge());
		if(!stepper->SwimToRadius(pos, mom, 65.0)){
			// Stepper managed to swim to the BCAL! Record position at BCAL inner face
			track_projection.push_back(pos);
		}
	}
}
#endif

