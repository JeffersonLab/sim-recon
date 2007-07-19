// $Id: Dmctrk_hit.cc 2673 2007-06-13 03:02:09Z davidl $
//
//    File: Dmctrk_hit.cc
// Created: Wed Jul 20 13:43:55 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "Dmctrk_hit.h"

//---------------------------------
// Dmctrk_hit    (Constructor)
//---------------------------------
Dmctrk_hit::Dmctrk_hit(const DMCTrackHit* hit)
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
// ~Dmctrk_hit    (Destructor)
//---------------------------------
Dmctrk_hit::~Dmctrk_hit()
{

}

