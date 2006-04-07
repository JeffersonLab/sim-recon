// $Id$
//
//    File: Dtrk_hit.cc
// Created: Wed Jul 20 13:43:55 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "Dtrk_hit.h"

//---------------------------------
// Dtrk_hit    (Constructor)
//---------------------------------
Dtrk_hit::Dtrk_hit(const DMCTrackHit* hit)
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
// ~Dtrk_hit    (Destructor)
//---------------------------------
Dtrk_hit::~Dtrk_hit()
{

}

