// $Id$
//
//    File: JEventProcessor_dumpthrowns.cc
// Created: Tue Feb  4 09:29:35 EST 2014
// Creator: davidl (on Linux ifarm1102 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_dumpthrowns.h"
using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>

#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DMCThrown.h>

extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_dumpthrowns());
}
} // "C"


//------------------
// JEventProcessor_dumpthrowns (Constructor)
//------------------
JEventProcessor_dumpthrowns::JEventProcessor_dumpthrowns()
{
	events_written = 0;
	events_discarded = 0;
}

//------------------
// ~JEventProcessor_dumpthrowns (Destructor)
//------------------
JEventProcessor_dumpthrowns::~JEventProcessor_dumpthrowns()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_dumpthrowns::init(void)
{
	// 
	MAX_CANDIDATE_FILTER = 1000;
	gPARMS->SetDefaultParameter("MAX_CANDIDATE_FILTER", MAX_CANDIDATE_FILTER, "Maximum number of candidates allowed in event before any are written to file.");

	
	// Open output file
	ofs = new ofstream("gluex_throwns.txt");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_dumpthrowns::brun(JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_dumpthrowns::evnt(JEventLoop *loop, int eventnumber)
{
	// Get track candidates
	vector<const DTrackCandidate*> candidates;
	loop->Get(candidates);
	if(candidates.size()==0 || candidates.size()>MAX_CANDIDATE_FILTER){
		events_discarded++;
		return NOERROR;
	}
	
	// Get thrown particles
	vector<const DMCThrown*> throwns;
	loop->Get(throwns);

	// Write out thrown parameters	
	for(unsigned int i=0; i<throwns.size(); i++){
		const DMCThrown *thrown = throwns[i];
		
		// Write thrown parameters to string
		stringstream ss;
		ss << thrown->charge();
		ss << " " << thrown->x() << " " << thrown->y() << " " << thrown->z();
		ss << " " << thrown->px() << " " << thrown->py() << " " << thrown->pz();
		
		// Write thrown parameters string to file
		(*ofs) << ss.str() << endl;
		events_written++;
		
		// Sometimes, generated particles are added to the thrown
		// particles list. We want only the first MAX_CANDIDATE_FILTER
		// particles which *may* correspond to the first candidates.
		// (probably not, but at least it may work for the single particle
		// case which is what we are interested in at the moment).
		if((i+1) >= MAX_CANDIDATE_FILTER) break;
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_dumpthrowns::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_dumpthrowns::fini(void)
{
	if(ofs){
		ofs->close();
		delete ofs;
		ofs =NULL;
	}
	
	cout << endl;
	cout << "Wrote " << events_written << " thrown events to output file (discarded " << events_discarded << ")" << endl;
	cout << endl;

	return NOERROR;
}



