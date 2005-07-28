// $Id$
//
//    File: DTrkHit.cc
// Created: Wed Jul 20 13:43:55 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DTrkHit.h"

//---------------------------------
// DTrkHit    (Constructor)
//---------------------------------
DTrkHit::DTrkHit(float x, float y, float z, float r, float phi, int system, int track, int ihit)
{
	this->x = x;
	this->y = y;
	this->z = z;
	this->r = r;
	this->phi = phi;
	this->system = system;
	this->track = track;
	this->ihit = ihit;
	
	flags = 0x0;
}

//---------------------------------
// ~DTrkHit    (Destructor)
//---------------------------------
DTrkHit::~DTrkHit()
{

}

