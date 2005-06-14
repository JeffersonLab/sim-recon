// $Id$
//
//    File: DEvent.cc
// Created: Wed Jun  8 12:30:53 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DEvent.h"

//---------------------------------
// DEvent    (Constructor)
//---------------------------------
DEvent::DEvent()
{
	source = NULL;
	event_number = 0 ;
	run_number = 0;
	ref = NULL;
}

//---------------------------------
// ~DEvent    (Destructor)
//---------------------------------
DEvent::~DEvent()
{

}
