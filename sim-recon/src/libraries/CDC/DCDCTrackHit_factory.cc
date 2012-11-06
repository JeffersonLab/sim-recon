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
#include <HDGEOMETRY/DGeometry.h>
#include <TRACKING/DTrackHitSelectorTHROWN.h>
#include <TRACKING/DMCTrackHit.h>

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
  //vector<double>cdc_half_length;
  vector<double>cdc_length;
  vector<double>cdc_center;
  dgeom->Get("//posXYZ[@volume='CentralDC']/@X_Y_Z",cdc_origin);
  dgeom->Get("//tubs[@name='STRA']/@Rio_Z",cdc_length);

  Z_MIN = cdc_origin[2];
  Z_MAX = Z_MIN + cdc_length[2];
  
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
  // Currently, the CDC design has 28 layers corresponding to ring=1-28.

  // Layers 5-8 and 17-20 are at negative stereo angles (roughly -6 degrees)
  // Layers 9-12 and 21-24 are at positive stereo angles (roughly +6 degrees)
  // The remaining layers have no stereo angle
  
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
      // axial
    case  1:	myNstraws=  42;	radius= 10.7219; stereo=  degrees0; phi_shift= 0.00000;	break;
    case  2:	myNstraws=  42;	radius= 12.0970; stereo=  degrees0; phi_shift= 4.285714;	break;
    case  3:	myNstraws=  54;	radius= 13.7803; stereo=  degrees0; phi_shift= 2.00000;	break;
    case  4:	myNstraws=  54;	radius= 15.1621; stereo=  degrees0; phi_shift= 5.3333333;	break;
      
      // -stereo
    case  5:	myNstraws=  66;	radius= 16.9321; stereo= -degrees6; phi_shift= 0.33333;	break;
    case  6:	myNstraws=  66;	phi_shift= 0.33333;	deltaX= 18.2948;	deltaY= 0.871486;	rotX= -6.47674;	rotY= -0.302853;	break;
    case  7:	myNstraws=  80;	radius= 20.5213; stereo= -degrees6; phi_shift= -0.5000;	break;
    case  8:	myNstraws=  80;	phi_shift= -0.5000;	deltaX= 21.8912;	deltaY= 0.860106;	rotX=-6.39548;	rotY=-0.245615;	break;
      
      // +stereo
    case  9:	myNstraws=  93;	radius= 23.8544; stereo= +degrees6; phi_shift= 1.1000;	break;
    case 10:	myNstraws=  93;	phi_shift= 1.1000;	deltaX=25.229;	deltaY= 0.852573;	rotX=+6.34142;	rotY=+0.208647;	break;
    case 11:	myNstraws= 106;	radius= 27.1877; stereo= +degrees6; phi_shift= -1.40;	break;
    case 12:	myNstraws= 106;	phi_shift= -1.400;	deltaX= 28.5658;	deltaY= 0.846871;	rotX=+6.30035;	rotY=+0.181146;	break;
      
      // axial
    case 13:	myNstraws= 123;	radius= 31.3799; stereo=  degrees0; phi_shift= 0.5000000;	break;
    case 14:	myNstraws= 123;	radius= 32.7747; stereo=  degrees0; phi_shift= 1.9634146;	break;
    case 15:	myNstraws= 135;	radius= 34.4343; stereo=  degrees0; phi_shift= 1.0000000;	break;
    case 16:	myNstraws= 135;	radius= 35.8301; stereo=  degrees0; phi_shift= 2.3333333;	break;
      
      // -stereo
    case 17:	myNstraws= 146;	radius= 37.4446; stereo= -degrees6; phi_shift= 0.2;	break;
    case 18:	myNstraws= 146;	phi_shift= 0.2;	deltaX= 38.8295;	deltaY= 0.835653;	rotX=-6.21919;	rotY=-0.128247;	break;
    case 19:	myNstraws= 158;	radius= 40.5364; stereo= -degrees6; phi_shift= 0.7;	break;
    case 20:	myNstraws= 158;	phi_shift= 0.7;	deltaX=41.9225 ;	deltaY= 0.833676;	rotX=-6.20274;	rotY=-0.118271;	break;
      
      // +stereo
    case 21:	myNstraws= 170;	radius= 43.6152; stereo= +degrees6; phi_shift= 1.1000;	break;
    case 22:	myNstraws= 170;	phi_shift= 1.1000;	deltaX= 45.0025;	deltaY= 0.83173;	rotX=6.18859;	rotY=+0.109325;	break;
    case 23:	myNstraws= 182;	radius= 46.6849; stereo= +degrees6; phi_shift= 1.40;	break;
    case 24:	myNstraws= 182;	phi_shift= 1.400;	deltaX= 48.0733;	deltaY=  0.829899;	rotX=+6.1763;	rotY=+0.101315;	break;
      
      // axial
    case 25:	myNstraws= 197;	radius= 50.3747; stereo=  degrees0; phi_shift= 0.200000000;	break;
    case 26:	myNstraws= 197;	radius= 51.7722; stereo=  degrees0; phi_shift= 1.113705;	break;
    case 27:	myNstraws= 209;	radius= 53.3631; stereo=  degrees0; phi_shift= 0.800000000;	break;
    case 28:	myNstraws= 209;	radius= 54.7617; stereo=  degrees0; phi_shift= 1.661244;	break;
      
    default:
      cerr<<__FILE__<<":"<<__LINE__<<" Invalid value for CDC ring ("<<ring<<") should be 1-28 inclusive!"<<endl;
    }
    Nstraws[ring-1] = myNstraws;
    
    // I'm not sure why this is needed, but empirically, it just is.
    //stereo = -stereo;
    //rotX = -rotX;
    //rotY = -rotY;
    
    // Convert phi_shift to radians
    phi_shift *= M_PI/180.0;
    
    float dphi = 2.0*M_PI/(float)myNstraws; // phi angle difference between straws
    for(int straw=1; straw<=myNstraws; straw++){
      DCDCWire *w = &wire[ring-1][straw-1];
      
      float phi = phi_shift + (float)(straw-1)*dphi;
      
      // Define position of midpoint of wire, phi of midpoint, and m_x, m_y
      w->ring = ring;
      w->straw = straw;
      if (fabs(rotY) < 1.0E-9){
	w->origin.SetX(radius*cos(phi));
	w->origin.SetY(radius*sin(phi));
      }
      else{
	w->origin.SetX(deltaX);
	w->origin.SetY(deltaY);
	w->origin.RotateZ(phi);
      }
      w->origin.SetZ((Z_MAX + Z_MIN)/2.0);
      w->phi = phi;
      //w->L = L/cos(stereo);

      //w->L = L + ((fabs(stereo - degrees0) < 1.0E-9) ? 0.0:1.5); // to make consistent with HDDS
			
      // Here, we need to define a coordinate system for the wire
      // in which the wire runs along one axis. We call the directions
      // of the axes in this coordinate system s,t, and u with
      // the wire running in the "u" direction. The "s" direction
      // will be defined by the direction pointing from the beamline
      // to the midpoint of the wire.
      w->udir.SetXYZ(0.0, 0.0,1.0);
      if (fabs(rotY) < 1.0E-9){

	w->udir.RotateX(stereo);	
	w->udir.RotateZ(phi);
			}
      else{
	w->udir.SetXYZ(0.0, 0.0,1.0);
	//w->udir.RotateY(rotY*deg2rad);	
	w->udir.RotateX(rotX*deg2rad);
	w->udir.RotateY(rotY*deg2rad);	
	w->udir.RotateZ(phi);
      }
      
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
      
      //w->sdir=w->origin-DVector3(0,0,w->origin.Z());
      w->sdir = w->origin - DVector3(0.0, 0.0, w->udir.Dot(w->origin)/w->udir.Z());  // see above comments
      w->sdir.SetMag(1.0);
      
      w->tdir = w->udir.Cross(w->sdir);
      w->tdir.SetMag(1.0); // This isn't really needed
      
      w->stereo = w->udir.Angle(DVector3(0,0,1));
      //if(rotX>0.0)w->stereo = -w->stereo;
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
	
	// If this is simulated data then we want to match up the truth hit
	// with this "real" hit. Ideally, this would be done at the
	// DCDCHit object level, but the organization of the data in HDDM
	// makes that difficult. Here we have the full wire definition so
	// we make the connection here.
	vector<const DMCTrackHit*> mctrackhits;
	loop->Get(mctrackhits);
	
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
		hit->is_stereo=((cdchit->ring>4&&cdchit->ring<13)
				||(cdchit->ring>16&&cdchit->ring<25))
				?true:false;
		hit->tdrift = cdchit->t;
		double w_eff=29.5e-9;
		double gas_gain=1e5;
		double electron_charge=1.6022e-4; /* fC */
		hit->dE=cdchit->dE*w_eff/(gas_gain*electron_charge);
		hit->dist = hit->tdrift*55.0E-4; // Use number hardwired in simulation for now
		
		// Try matching truth hit with this "real" hit.
		const DMCTrackHit *mctrackhit = DTrackHitSelectorTHROWN::GetMCTrackHit(hit->wire, hit->dist, mctrackhits);
		
		hit->AddAssociatedObject(cdchit);
		if(mctrackhit)hit->AddAssociatedObject(mctrackhit);
		
		_data.push_back(hit);
	}
	/*
	int ring=0,straw=0;
	double t=0;
	for (unsigned int i=0;i<_data.size();i++){
	  if (_data[i]->wire->ring==ring && _data[i]->wire->straw==straw
	      && _data[i]->tdrift<t) printf("out of order\n");
	  printf("r %d s %d t %f\n",_data[i]->wire->ring,_data[i]->wire->straw,
		 _data[i]->tdrift);
	  ring=_data[i]->wire->ring;
	  straw=_data[i]->wire->straw;
	  t=_data[i]->tdrift;
	}
	*/
	
	return NOERROR;
}

//------------------
// GetCDCWire
//------------------
const DCDCWire* DCDCTrackHit_factory::GetCDCWire(int ring, int straw)
{
	if(ring<1 || ring>CDC_MAX_RINGS || straw<1 || straw>CDC_MAX_STRAWS)return NULL;
	
	return &wire[ring-1][straw-1];
}

