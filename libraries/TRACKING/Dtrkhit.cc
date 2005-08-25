// $Id$
//
//    File: Dtrkhit.cc
// Created: Wed Jul 20 13:43:55 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "Dtrkhit.h"

//---------------------------------
// Dtrkhit    (Constructor)
//---------------------------------
Dtrkhit::Dtrkhit(const DMCTrackHit* hit)
{
	r = hit->r;
	phi = hit->phi;
	z = hit->z;
	system = hit->system;
	id = hit->id;
	
	x = r*cos(phi);
	y = r*sin(phi);
	
	flags = 0x0;
}

//---------------------------------
// ~Dtrkhit    (Destructor)
//---------------------------------
Dtrkhit::~Dtrkhit()
{

}

