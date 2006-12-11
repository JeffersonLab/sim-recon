// $Id$
//
//    File: DCDCTrackHit_factory.cc
// Created: Mon Oct 16 10:20:07 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.7.0 powerpc)
//

#include <cmath>
#include <pthread.h>

#include "DCDCTrackHit_factory.h"
#include "DCDCHit.h"

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
	// Initialize the "wire" table that is a static global variable.
	// This table takes up about 200kB - 300kB so we don't want
	// to duplicate it for each thread which is why it is static.
	// It wastes some space since it allocates MAX_STRAWS for every ring even
	// even though only the outer one actually has that many. This is a
	// trade off to allow fast access to the array
	// during event time.
	
	Z_MIN = 17.0;
	Z_MAX = Z_MIN + 175.0;
	
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
		switch(ring){
			case  1:	myNstraws=  63;	radius= 16.049;	stereo=  degrees0; break;
			case  2:	myNstraws=  70;	radius= 17.831;	stereo=  degrees0; break;
			case  3:	myNstraws=  77;	radius= 19.613;	stereo=  degrees0; break;
			case  4:	myNstraws=  84;	radius= 21.395;	stereo=  degrees0; break;
			case  5:	myNstraws=  91;	radius= 23.178;	stereo= +degrees6; break;
			case  6:	myNstraws=  98;	radius= 24.960;	stereo= +degrees6; break;
			case  7:	myNstraws= 105;	radius= 26.742;	stereo= -degrees6; break;
			case  8:	myNstraws= 112;	radius= 28.524;	stereo= -degrees6; break;
			case  9:	myNstraws= 126;	radius= 32.089;	stereo=  degrees0; break;
			case 10:	myNstraws= 133;	radius= 33.871;	stereo=  degrees0; break;
			case 11:	myNstraws= 140;	radius= 35.654;	stereo=  degrees0; break;
			case 12:	myNstraws= 147;	radius= 37.435;	stereo=  degrees0; break;
			case 13:	myNstraws= 154;	radius= 39.218;	stereo=  degrees0; break;
			case 14:	myNstraws= 161;	radius= 41.001;	stereo= +degrees6; break;
			case 15:	myNstraws= 168;	radius= 42.783;	stereo= +degrees6; break;
			case 16:	myNstraws= 175;	radius= 44.566;	stereo= -degrees6; break;
			case 17:	myNstraws= 182;	radius= 46.348;	stereo= -degrees6; break;
			case 18:	myNstraws= 193;	radius= 49.149;	stereo=  degrees0; break;
			case 19:	myNstraws= 200;	radius= 50.932;	stereo=  degrees0; break;
			case 20:	myNstraws= 207;	radius= 52.714;	stereo=  degrees0; break;
			case 21:	myNstraws= 214;	radius= 54.497;	stereo=  degrees0; break;
			case 22:	myNstraws= 221;	radius= 56.279;	stereo=  degrees0; break;
			case 23:	myNstraws= 228;	radius= 58.062;	stereo=  degrees0; break;
			default:
				cerr<<__FILE__<<":"<<__LINE__<<" Invalid value for CDC ring ("<<ring<<") should be 1-23 inclusive!"<<endl;
		}
		Nstraws[ring-1] = myNstraws;
		
		float dphi = 2.0*M_PI/(float)myNstraws; // phi angle difference between straws
		for(int straw=1; straw<=myNstraws; straw++){
			DCDCWire *w = &wire[ring-1][straw-1];

			float phi = (float)(straw-1)*dphi;
			
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

//------------------
// toString
//------------------
const string DCDCTrackHit_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("ring: straw:    x(cm):     y(cm): stereo(rad):  tdrift(ns):  dist(cm):");

	for(unsigned int i=0; i<_data.size(); i++){
		DCDCTrackHit *hit = _data[i];
		const DCDCWire *w = hit->wire;

		printnewrow();
		printcol("%d",	w->ring);
		printcol("%d",	w->straw);
		printcol("%3.1f",	w->origin.x());
		printcol("%3.1f",	w->origin.y());
		printcol("%1.4f",	w->stereo);
		printcol("%3.1f",	hit->tdrift);
		printcol("%3.1f",	hit->dist);
		printrow();
	}

	return _table;

}
