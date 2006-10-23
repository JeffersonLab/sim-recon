// $Id$
//
//    File: DCDCWire.h
// Created: Wed Oct 18 11:39:35 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.7.0 powerpc)
//

#ifndef _DCDCWire_
#define _DCDCWire_

#include <TVector3.h>

// this is defined as a struct instead of a class because several
// thousand of these are statically allocated in DCDCTrackHit_factory
// and I want to make sure the compiler doesn't generate a
// constructor and invoke it for every one.

typedef struct{
	TVector3 wpos; // x,y,z-coordinate of wire at midplane in lab coordinates
	TVector3 sdir; // s direction of wire coord. system
	TVector3 tdir; // t direction of wire coord. system
	TVector3 udir; // u direction of wire coord. system (along the wire)

	float phi;		// phi angle of wire at midplane in lab coordinates

	int ring;
	int straw;
	float stereo;
}DCDCWire;

#endif // _DCDCWire_

