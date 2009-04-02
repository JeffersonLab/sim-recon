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

	float deg2rad = M_PI/180.0;
	float degrees6 = 6.0*deg2rad;
	float degrees0 = 0.0;
	double L = Z_MAX - Z_MIN; // full length of wire
	for(int ring=1; ring<=CDC_MAX_RINGS; ring++){
		int myNstraws=0;
		float radius = 0.0;
		float stereo=0.0;
		float phi_shift=0.0;
		float deltaX=0.0;
		float deltaY=0.0;
		float rotX=0.0;
		float rotY=0.0;

		switch(ring){
			// Geometry "C" from http://www.jlab.org/Hall-D/software/wiki/index.php/CDC_Future_MC_Studies_Discussion%2C_September_29%2C_2008

			// axial
			case  1:	myNstraws=  43;	radius= 10.984;	stereo=  degrees0; phi_shift= 0.00000;	break;
			case  2:	myNstraws=  43;	radius= 12.341;	stereo=  degrees0; phi_shift= 4.18605;	break;
			case  3:	myNstraws=  55;	radius= 14.029;	stereo=  degrees0; phi_shift= 2.00000;	break;
			case  4:	myNstraws=  55;	radius= 15.410;	stereo=  degrees0; phi_shift= 5.27272727272;	break;

			// -stereo
			case  5:	myNstraws=  66;	radius= 17.085;	stereo= -degrees6; phi_shift= 0.33333;	break;
			case  6:	myNstraws=  66;	phi_shift= 0.33333;	deltaX= 18.450;	deltaY= 0.886027;	rotX=-6.500;	rotY=-0.30;	break;
			case  7:	myNstraws=  80;	radius= 20.581;	stereo= -degrees6; phi_shift= -0.5000;	break;
			case  8:	myNstraws=  80;	phi_shift= -0.5000;	deltaX= 21.950;	deltaY= 0.87224;	rotX=-6.475;	rotY=-0.27;	break;

			// +stereo
			case  9:	myNstraws=  93;	radius= 23.980;	stereo= +degrees6; phi_shift= 1.1000;	break;
			case 10:	myNstraws=  93;	phi_shift= 1.1000;	deltaX= 25.35;	deltaY= 0.857573;	rotX=+6.35;	rotY=+0.24;	break;
			case 11:	myNstraws= 106;	radius= 27.380;	stereo= +degrees6; phi_shift= -1.40;	break;
			case 12:	myNstraws= 106;	phi_shift= -1.400;	deltaX= 28.800;	deltaY= 0.835;	rotX=+6.3;	rotY=+0.21;	break;

			// axial
			case 13:	myNstraws= 124;	radius= 31.65;	stereo=  degrees0; phi_shift= 0.5000000;	break;
			case 14:	myNstraws= 124;	radius= 33.05;	stereo=  degrees0; phi_shift= 1.9516000;	break;
			case 15:	myNstraws= 133;	radius= 34.7;	stereo=  degrees0; phi_shift= 1.0000000;	break;
			case 16:	myNstraws= 133;	radius= 36.1;	stereo=  degrees0; phi_shift= 2.3533834;	break;

			// -stereo
			case 17:	myNstraws= 145;	radius= 37.720;	stereo= -degrees6; phi_shift= 0.2;	break;
			case 18:	myNstraws= 145;	phi_shift= 0.2;	deltaX= 39.12;	deltaY= 0.8321106;	rotX=-6.25;	rotY=-0.16;	break;
			case 19:	myNstraws= 158;	radius= 40.89;	stereo= -degrees6; phi_shift= 0.7;	break;
			case 20:	myNstraws= 158;	phi_shift= 0.7;	deltaX=42.29 ;	deltaY= 0.8391063;	rotX=-6.25;	rotY=-0.13;	break;

			// +stereo
			case 21:	myNstraws= 171;	radius= 44.05;	stereo= +degrees6; phi_shift= 1.1000;	break;
			case 22:	myNstraws= 171;	phi_shift= 1.1000;	deltaX=45.46 ;	deltaY= 0.8349124;	rotX=+6.20;	rotY=+0.11;	break;
			case 23:	myNstraws= 184;	radius= 47.20;	stereo= +degrees6; phi_shift= 1.40;	break;
			case 24:	myNstraws= 184;	phi_shift= 1.400;	deltaX= 48.63;	deltaY= 0.8303831;	rotX=+6.18;	rotY=+0.10;	break;

			// axial
			case 25:	myNstraws= 197;	radius= 51.0;	stereo=  degrees0; phi_shift= 0.200000000;	break;
			case 26:	myNstraws= 197;	radius= 52.4;	stereo=  degrees0; phi_shift= 1.113705000;	break;
			case 27:	myNstraws= 210;	radius= 54.05;	stereo=  degrees0; phi_shift= 0.800000000;	break;
			case 28:	myNstraws= 210;	radius= 55.45;	stereo=  degrees0; phi_shift= 1.657142857;	break;

			default:
				cerr<<__FILE__<<":"<<__LINE__<<" Invalid value for CDC ring ("<<ring<<") should be 1-28 inclusive!"<<endl;
		}
		Nstraws[ring-1] = myNstraws;
		
		// I'm not sure why this is needed, but empirically, it just is.
		stereo = -stereo;
		rotX = -rotX;
		rotY = -rotY;
		
		// For close-packed stereo layers, rotX is set. For normal stereo layers, stereo is set.
		// In the case of normal stereo wires, copy the rotation to the rotX field
		if(stereo!=0.0)rotX = stereo/deg2rad;
		
		// Convert phi_shift to radians
		phi_shift *= M_PI/180.0;
		
		// Close-packed stereos have and additional, initial phi shift in order to be close packed
		phi_shift += atan2(deltaY, deltaX);
		
		// If radius is 0 then it means this is a close-packed stereo whose position is set
		// by deltaX and deltaY. Calculate the initial radius from these
		if(radius==0.0)radius=sqrt(deltaX*deltaX + deltaY*deltaY);
		
		float dphi = 2.0*M_PI/(float)myNstraws; // phi angle difference between straws
		for(int straw=1; straw<=myNstraws; straw++){
			DCDCWire *w = &wire[ring-1][straw-1];

			float phi = phi_shift + (float)(straw-1)*dphi;
			
			// Define position of midpoint of wire, phi of midpoint, and m_x, m_y
			w->ring = ring;
			w->straw = straw;
			w->origin.SetX(radius*cos(phi));
			w->origin.SetY(radius*sin(phi));
			w->origin.SetZ((Z_MAX + Z_MIN)/2.0);
			w->phi = phi;
			//w->L = L/cos(stereo);
			//w->L = L + (stereo==degrees0 ? 0.0:1.5); // to make consistent with HDDS
			
			// Here, we need to define a coordinate system for the wire
			// in which the wire runs along one axis. We call the directions
			// of the axes in this coordinate system s,t, and u with
			// the wire running in the "u" direction. The "s" direction
			// will be defined by the direction pointing from the beamline
			// to the midpoint of the wire.
			
			w->udir.SetXYZ(0.0, 0.0, 1.0);
			w->udir.RotateX(-rotX*deg2rad);
			w->udir.RotateY(-rotY*deg2rad);
			w->udir.RotateZ(phi);
			
			// With the addition of close-packed stereo wires, the vector connecting
			// the center of the wire to the beamline ("s" direction) is not necessarily
			// perpendicular to the beamline. By definition, we want the "s" direction
			// to be perpendicular to the wire direction "u" and pointing at the beamline.
			// 
			// NOTE: This extensive comment is here because the result, when implmented
			// below caused a WORSE residual distribution in the close-packed stereo
			// layers. Owing to lack of time currently to track the issue down (most
			// likely in DReferenceTrajectory) I'm commenting out the "correct" calculation
			// of s, but leaving this comment so the issue can be revisited later. This
			// error leads to around 100 micron errors in the C.P.S. wires, but they
			// are completely washed out when the position smearing of 150 microns is applied
			// making the error unnoticable except when position smearing is not applied.
			//
			// April 2, 2009  D.L.
			//
			// Here is how this is calculated -- We define a vector equation with 2 unknowns
			// Z and S:
			//
			//    Zz + Ss = W
			//
			// where:  z = unit vector in z direction
			//         s = unit vector in "s" direction
			//         W = vector pointing to center of wire in lab coordinates
			//
			//  Rearranging, we get:
			//
			//     s = (W - Zz)/S
			//
			//  Since s must be perpendicular to u, we take a dot product of s and u
			// and set it equal to zero to determine Z:
			//
			//    u.s = 0 = u.(W - Zz)/S  =>  u.W = Zu.z
			//
			//   or
			//
			//     Z = u.W/u.z
			//
			//  Thus, the s direction is just given by (W - (u.W/u.z)z)
			//
			
			
			w->sdir=w->origin-TVector3(0,0,w->origin.Z());
			//w->sdir = w->origin - TVector3(0.0, 0.0, w->udir.Dot(w->origin)/w->udir.Z());  // see above comments
			w->sdir.SetMag(1.0);

			w->tdir = w->udir.Cross(w->sdir);
			w->tdir.SetMag(1.0); // This isn't really needed

			w->stereo = w->udir.Angle(TVector3(0,0,1));
			if(rotX>0.0)w->stereo = -w->stereo;
			w->L = L/cos(w->stereo);
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

