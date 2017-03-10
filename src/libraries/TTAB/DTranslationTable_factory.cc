// $Id$
//
//    File: DTranslationTable_factory.cc
// Created: Thu Jun 27 15:33:38 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DTranslationTable_factory.h"
using namespace jana;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>

//------------------
// init
//------------------
jerror_t DTranslationTable_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTranslationTable_factory::brun(jana::JEventLoop *loop, int32_t runnumber)
{
	if(!_data.empty()){
		jout << "WARNING: Translation table for run " << runnumber << " requested but" << endl;
		jout << "         translation table already exists. Ignoring request and using" << endl;
		jout << "         table already loaded." << endl;
		tt->SetSystemsToParse(loop->GetJEvent().GetJEventSource());
		return NOERROR;
	}
	jout << "Creating DTranslationTable for run " << runnumber << endl;

	// Grab run-dependent translation table from CCDB
	tt = new DTranslationTable(loop);
	
	// Keep this translation table around and reuse it for
	// susequent events
	_data.push_back(tt);
	SetFactoryFlag(PERSISTANT);
	
	// If restricting parsing, make sure it is set for this source
	tt->SetSystemsToParse(loop->GetJEvent().GetJEventSource());

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTranslationTable_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DTranslationTable_factory::erun(void)
{
	// If we have a translation table already the delete it
	if(tt){
		delete tt;
	}
	_data.clear();

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTranslationTable_factory::fini(void)
{
	return NOERROR;
}

