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
#include "JFactoryGenerator_DTranslationTable.h"
using namespace jana;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddFactoryGenerator(new JFactoryGenerator_DTranslationTable());
}
} // "C"

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
jerror_t DTranslationTable_factory::brun(jana::JEventLoop *loop, int runnumber)
{
	// Grab run-dependent translation table from CCDB
	tt = new DTranslationTable(loop);
	
	// Keep this translation table around and reuse it for
	// susequent events
	_data.push_back(tt);
	SetFactoryFlag(PERSISTANT);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTranslationTable_factory::evnt(JEventLoop *loop, int eventnumber)
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

