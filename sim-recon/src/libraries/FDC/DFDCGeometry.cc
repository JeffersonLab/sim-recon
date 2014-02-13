// $Id$
//
//    File: DFDCGeometry.h
// Created: Wed Nov 29 13:35 EST 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.7.0 powerpc)
//

#include "DFDCGeometry.h"

#if 0  // Old method for FDC geometry with hard coded positions!!!

// Static globals used by all instances of DCDCTrackHit_factory
static pthread_mutex_t wire_mutex = PTHREAD_MUTEX_INITIALIZER;
static bool wire_table_initialized = false;
static DFDCWire fdcwire[FDC_NUM_LAYERS][WIRES_PER_PLANE];

//-----------------
// DFDCGeometry   (Constructor)
//-----------------
DFDCGeometry::DFDCGeometry(void)
{
	// We use the mutex and wire_table_initialized flag to
	// make sure the table is initialized only once.
	pthread_mutex_lock(&wire_mutex);
	if(wire_table_initialized){
		pthread_mutex_unlock(&wire_mutex);
		return;
	}

	float degrees00  = 0.0;
	float degrees60  = 60.0*M_PI/180.0;
	for(int layer=1; layer<=FDC_NUM_LAYERS; layer++){
		
		float angle=0.0;
		float z_anode=188.5+85.0; 
		switch(layer){
			case  1: z_anode+= -76.0-4.5-3.0;	angle=  degrees00; break;
			case  2: z_anode+= -76.0-4.5-0.0;	angle= +degrees60; break;
			case  3: z_anode+= -76.0-4.5+3.0;	angle= -degrees60; break;
			case  4: z_anode+= -76.0+4.5-3.0;	angle=  degrees00; break;
			case  5: z_anode+= -76.0+4.5-0.0;	angle= +degrees60; break;
			case  6: z_anode+= -76.0+4.5+3.0;	angle= -degrees60; break;

			case  7: z_anode+= -25.33-4.5-3.0;	angle=  degrees00; break;
			case  8: z_anode+= -25.33-4.5-0.0;	angle= +degrees60; break;
			case  9: z_anode+= -25.33-4.5+3.0;	angle= -degrees60; break;
			case 10: z_anode+= -25.33+4.5-3.0;	angle=  degrees00; break;
			case 11: z_anode+= -25.33+4.5-0.0;	angle= +degrees60; break;
			case 12: z_anode+= -25.33+4.5+3.0;	angle= -degrees60; break;

			case 13: z_anode+= +25.33-4.5-3.0;	angle=  degrees00; break;
			case 14: z_anode+= +25.33-4.5-0.0;	angle= +degrees60; break;
			case 15: z_anode+= +25.33-4.5+3.0;	angle= -degrees60; break;
			case 16: z_anode+= +25.33+4.5-3.0;	angle=  degrees00; break;
			case 17: z_anode+= +25.33+4.5-0.0;	angle= +degrees60; break;
			case 18: z_anode+= +25.33+4.5+3.0;	angle= -degrees60; break;

			case 19: z_anode+= +76.0-4.5-3.0;	angle=  degrees00; break;
			case 20: z_anode+= +76.0-4.5-0.0;	angle= +degrees60; break;
			case 21: z_anode+= +76.0-4.5+3.0;	angle= -degrees60; break;
			case 22: z_anode+= +76.0+4.5-3.0;	angle=  degrees00; break;
			case 23: z_anode+= +76.0+4.5-0.0;	angle= +degrees60; break;
			case 24: z_anode+= +76.0+4.5+3.0;	angle= -degrees60; break;
		}
		
		// Somewhere between HDDS and here I'm missing a sign. I'm not
		// sure where it is, but empirically, it must be here for things
		// to be consistent.
		angle=-angle;

		for(int wire=1; wire<=WIRES_PER_PLANE; wire++){

			DFDCWire *w = &fdcwire[layer-1][wire-1];
			w->layer = layer;
			w->wire = wire;
			w->angle = angle;
			
			// find coordinates of center of wire in rotated system
			float u = U_OF_WIRE_ZERO + WIRE_SPACING*(float)(wire-1);
			
			// Rotate coordinates into lab system and set the wire's origin
			// Note that the FDC measures "angle" such that angle=0
			// corresponds to the anode wire in the vertical direction
			// (i.e. at phi=90 degrees).
			float x = u*sin(angle + M_PI/2.0);
			float y = u*cos(angle + M_PI/2.0);
			w->origin.SetXYZ(x,y,z_anode);
			
			// Length of wire is set by active radius
			w->L = 2.0*sqrt(pow(FDC_ACTIVE_RADIUS,2.0) - u*u);
			
			// Set directions of wire's coordinate system with "udir"
			// along wire.
			w->udir.SetXYZ(sin(angle),cos(angle),0.0);
			
			// "s" points in direction from beamline to midpoint of
			// wire. This happens to be the same direction as "origin"
			w->sdir = w->origin;
			w->sdir.SetMag(1.0);
			
			w->tdir = w->udir.Cross(w->sdir);
			w->tdir.SetMag(1.0); // This isn't really needed
		}
	}

	// Flag table as initialized and release the lock so other threads
	// can continue.
	wire_table_initialized=true;
	pthread_mutex_unlock(&wire_mutex);
}

//-----------------
// DFDCGeometry   (Constructor)
//-----------------
const DFDCWire* DFDCGeometry::GetDFDCWire(int layer, int wire)
{
	if(layer<1)return NULL;
	if(layer>FDC_NUM_LAYERS)return NULL;
	if(wire<1)return NULL;
	if(wire>WIRES_PER_PLANE)return NULL;

	return &fdcwire[layer-1][wire-1];
}
#endif

