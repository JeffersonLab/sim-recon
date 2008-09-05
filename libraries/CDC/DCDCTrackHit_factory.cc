// $Id$
//
//    File: DCDCTrackHit_factory.cc
// Created: Mon Oct 16 10:20:07 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.7.0 powerpc)
//

#include <cmath>
#include <pthread.h>
using namespace std;

#include "DCDCTrackHit_factory.h"
#include "DCDCHit.h"
#include "HDGEOMETRY/DGeometry.h"

// Static globals used by all instances of DCDCTrackHit_factory
static pthread_mutex_t wire_mutex = PTHREAD_MUTEX_INITIALIZER;
static bool wire_table_initialized = false;
static DCDCWire wire[CDC_MAX_RINGS][CDC_MAX_STRAWS];
static int Nstraws[CDC_MAX_RINGS];

//------------------
// init
//------------------
jerror_t DCDCTrackHit_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DCDCTrackHit_factory::brun(JEventLoop *loop, int runnumber)
{
	// Initialize the "wire" table that is a static global variable.
	// This table takes up about 200kB - 300kB so we don't want
	// to duplicate it for each thread which is why it is static.
	// It wastes some space since it allocates MAX_STRAWS for every ring even
	// even though only the outer one actually has that many. This is a
	// trade off to allow fast access to the array
	// during event time.
	
  // Get pointer to DGeometry object
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  dgeom  = dapp->GetDGeometry(runnumber);

  vector<double>cdc_origin;
  vector<double>cdc_half_length;
  dgeom->Get("//posXYZ[@volume='CentralDC']/@X_Y_Z",cdc_origin);
  dgeom->Get("//posXYZ[@volume='centralDC_option-1']/@X_Y_Z",cdc_half_length);


	Z_MIN = cdc_origin[2];
	Z_MAX = Z_MIN + 2.*cdc_half_length[2];

	gPARMS->SetDefaultParameter("CDC:Z_MIN",Z_MIN);
	gPARMS->SetDefaultParameter("CDC:Z_MAX",Z_MAX);
	
	// We use the mutex and wire_table_initialized flag to
	// make sure the table is initialized only once.
	pthread_mutex_lock(&wire_mutex);
	if(wire_table_initialized){
		pthread_mutex_unlock(&wire_mutex);
		return NOERROR;
	}
	
	// Note: This information should ideally come from HDDS, but there
	// is currently no trivial mechanism by which to access it from here.

	// We calculate the straw coordinates knowing:
	//
	// 1.) The number of straws in a ring given the ring number
	// 2.) The radius of the wire at the midplane of the specified ring
	// 3.) The first straw is always at phi=0 in the ring.
	//
	// Currently, the CDC design has 23 layers corresponding to ring=1-23.
	// The stereo angles are as follows:
	//
	// layers 5,6,14,15 are +6 degrees
	// layers 7,8,16,17 are -6 degrees
	// all other layers are zero degrees

	float degrees6 = 6.0*M_PI/180.0;
	float degrees0 = 0.0;
	double L = Z_MAX - Z_MIN; // full length of wire
	for(int ring=1; ring<=CDC_MAX_RINGS; ring++){
		int myNstraws=0;
		float radius = 0.0;
		float stereo=0.0;
		float phi_shift=0.0;

		switch(ring){
			// Phase shifted configuration
			// NOTE: These follow an XML file sent ot me by Beni Z. on July 18, 2008
			// but the sign of the stereo shift seems opposite to that used in GlueX-doc-990-v10.
			// We stick with Beni's convention for simplicity here on the assumption that
			// this choice is arbitrary.  7/21/2008 DL.
			case  1:	myNstraws=  43;	radius= 10.984;	stereo=  degrees0; phi_shift= 0.000;	break;
			case  2:	myNstraws=  50;	radius= 12.769;	stereo=  degrees0; phi_shift=+3.600;	break;
			case  3:	myNstraws=  57;	radius= 14.555;	stereo=  degrees0; phi_shift=-2.105;	break;
			case  4:	myNstraws=  63;	radius= 16.174;	stereo= -degrees6; phi_shift= 0.000;	break;
			case  5:	myNstraws=  70;	radius= 17.969;	stereo= -degrees6; phi_shift= 0.000;	break;
			case  6:	myNstraws=  77;	radius= 19.765;	stereo= +degrees6; phi_shift= 0.000;	break;
			case  7:	myNstraws=  84;	radius= 21.561;	stereo= +degrees6; phi_shift= 0.000;	break;
			case  8:	myNstraws=  98;	radius= 25.015;	stereo=  degrees0; phi_shift=+0.216;	break;
			case  9:	myNstraws= 105;	radius= 26.801;	stereo=  degrees0; phi_shift= 0.000;	break;
			case 10:	myNstraws= 112;	radius= 28.588;	stereo=  degrees0; phi_shift=+1.607;	break;
			case 11:	myNstraws= 119;	radius= 30.374;	stereo=  degrees0; phi_shift=-1.008;	break;
			case 12:	myNstraws= 126;	radius= 32.160;	stereo=  degrees0; phi_shift=+0.571;	break;
			case 13:	myNstraws= 132;	radius= 33.877;	stereo= -degrees6; phi_shift= 0.000;	break;
			case 14:	myNstraws= 139;	radius= 35.673;	stereo= -degrees6; phi_shift= 0.000;	break;
			case 15:	myNstraws= 146;	radius= 37.469;	stereo= +degrees6; phi_shift= 0.000;	break;
			case 16:	myNstraws= 153;	radius= 39.265;	stereo= +degrees6; phi_shift= 0.000;	break;
			case 17:	myNstraws= 165;	radius= 42.113;	stereo=  degrees0; phi_shift= 0.000;	break;
			case 18:	myNstraws= 172;	radius= 43.899;	stereo=  degrees0; phi_shift=+1.047;	break;
			case 19:	myNstraws= 179;	radius= 45.686;	stereo=  degrees0; phi_shift=-0.670;	break;
			case 20:	myNstraws= 186;	radius= 47.472;	stereo=  degrees0; phi_shift=+0.387;	break;
			case 21:	myNstraws= 193;	radius= 49.258;	stereo=  degrees0; phi_shift=-0.266;	break;
			case 22:	myNstraws= 200;	radius= 51.045;	stereo=  degrees0; phi_shift=+0.164;	break;
			case 23:	myNstraws= 207;	radius= 52.831;	stereo=  degrees0; phi_shift=-0.134;	break;
			case 24:	myNstraws= 214;	radius= 54.618;	stereo=  degrees0; phi_shift=+0.099;	break;
			default:
				cerr<<__FILE__<<":"<<__LINE__<<" Invalid value for CDC ring ("<<ring<<") should be 1-24 inclusive!"<<endl;
		}
		Nstraws[ring-1] = myNstraws;
		
		// I'm not sure why this is needed, but empirically, it just is.
		stereo = -stereo;
		
		// Convert phi_shift to radians
		phi_shift *= M_PI/180.0;
		
		float dphi = 2.0*M_PI/(float)myNstraws; // phi angle difference between straws
		for(int straw=1; straw<=myNstraws; straw++){
			DCDCWire *w = &wire[ring-1][straw-1];

			float phi = phi_shift + (float)(straw-1)*dphi;
			
			// Define position of midpoint of wire, phi of midpoint, and m_x, m_y
			w->ring = ring;
			w->straw = straw;
			w->stereo = stereo;
			w->origin.SetX(radius*cos(phi));
			w->origin.SetY(radius*sin(phi));
			w->origin.SetZ((Z_MAX + Z_MIN)/2.0);
			w->phi = phi;
			//w->L = L/cos(stereo);
			w->L = L + (stereo==degrees0 ? 0.0:1.5); // to make consistent with HDDS
			
			// Here, we need to define a coordinate system for the wire
			// in which the wire runs along one axis. We call the directions
			// of the axes in this coordinate system s,t, and u with
			// the wire running in the "u" direction. The "s" direction
			// will be defined by the direction pointing from the beamline
			// to the midpoint of the wire.
			w->sdir.SetXYZ(cos(w->phi), sin(w->phi), 0.0);
			w->sdir.SetMag(1.0); // This isn't really needed
			
			w->udir.SetXYZ(-sin(stereo)*sin(w->phi), sin(stereo)*cos(w->phi), cos(stereo));
			w->udir.SetMag(1.0); // This isn't really needed
			
			w->tdir = w->udir.Cross(w->sdir);
			w->tdir.SetMag(1.0); // This isn't really needed
		}
	}

	// Flag table as initialized and release the lock so other threads
	// can continue.
	wire_table_initialized=true;
	pthread_mutex_unlock(&wire_mutex);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DCDCTrackHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
	/// Convert from ring/straw indexing to x/y position
	/// of wire center and stereo angle.
	vector<const DCDCHit*> cdchits;
	loop->Get(cdchits);
	
	for(unsigned int i=0; i<cdchits.size(); i++){
		const DCDCHit* cdchit = cdchits[i];

		if(cdchit->ring>CDC_MAX_RINGS || cdchit->ring<1
			|| cdchit->straw>Nstraws[cdchit->ring-1] || cdchit->straw<1){
			cerr<<__FILE__<<":"<<__LINE__<<" Ring or straw out of range! ring="
				<<cdchit->ring<<" (should be 1-"<<CDC_MAX_RINGS<<")  straw="
				<<cdchit->straw<<" (should be 1-"<<Nstraws[cdchit->ring-1]<<")"<<endl;
			continue;
		}

		DCDCTrackHit *hit = new DCDCTrackHit;
		hit->wire = &wire[cdchit->ring-1][cdchit->straw-1];
		hit->tdrift = cdchit->t;
		hit->dist = hit->tdrift*55.0E-4; // Use number hardwired in simulation for now
		
		_data.push_back(hit);
	}


	return NOERROR;
}

