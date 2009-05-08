// $Id$
//
//    File: DTrackCandidate_factory_CDC_or_FDCpseudo.cc
// Created: Thu Apr 16 09:14:49 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DTrackCandidate_factory_CDC_or_FDCpseudo.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DTrackCandidate_factory_CDC_or_FDCpseudo::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackCandidate_factory_CDC_or_FDCpseudo::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackCandidate_factory_CDC_or_FDCpseudo::evnt(JEventLoop *loop, int eventnumber)
{

	/// This factory simply combines the list of candidates from the 
	/// DTrackCandidate:CDC and DTrackCandidate:FDCpseudo factories
	/// into a single lit. It simply copies the pointers and flags
	/// itself as not being the owner of any of these objects.
	vector<const DTrackCandidate*> cdc;
	vector<const DTrackCandidate*> fdc;

	loop->Get(cdc, "CDC");
	loop->Get(fdc, "FDCpseudo");
	
	for(unsigned int i=0; i<cdc.size(); i++)_data.push_back((DTrackCandidate*)cdc[i]);
	for(unsigned int i=0; i<fdc.size(); i++)_data.push_back((DTrackCandidate*)fdc[i]);
	
	SetFactoryFlag(NOT_OBJECT_OWNER); // make sure these aren't deleted twice

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DTrackCandidate_factory_CDC_or_FDCpseudo::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackCandidate_factory_CDC_or_FDCpseudo::fini(void)
{
	return NOERROR;
}

