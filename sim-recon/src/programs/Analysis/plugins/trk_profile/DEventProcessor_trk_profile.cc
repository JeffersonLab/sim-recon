// $Id$
//
//    File: DEventProcessor_trk_profile.cc
// Created: Wed Jan 12 08:02:32 EST 2011
// Creator: davidl (on Darwin eleanor.jlab.org 10.6.0 i386)
//

#include "DEventProcessor_trk_profile.h"
#include <TRACKING/DTrackFitter.h>
#include <PID/DChargedTrack.h>
using namespace jana;


const DTrackFitter *fitter = NULL;

// Routine used to create our DEventProcessor
#include <JANA/JApplication.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_trk_profile());
}
} // "C"


//------------------
// DEventProcessor_trk_profile (Constructor)
//------------------
DEventProcessor_trk_profile::DEventProcessor_trk_profile()
{

}

//------------------
// ~DEventProcessor_trk_profile (Destructor)
//------------------
DEventProcessor_trk_profile::~DEventProcessor_trk_profile()
{

}

//------------------
// init
//------------------
jerror_t DEventProcessor_trk_profile::init(void)
{
	// Create histograms here
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_trk_profile::brun(JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_trk_profile::evnt(JEventLoop *loop, int eventnumber)
{
	// Get DChargedTrack objects to activate tracking
	vector<const DChargedTrack*> charged_tracks;
	loop->Get(charged_tracks);

	// Get the track fitter object for use later
	if(!fitter)loop->GetSingle(fitter);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_trk_profile::erun(void)
{
	if(fitter){
		map<string, prof_time::time_diffs> prof_times;
		fitter->GetProfilingTimes(prof_times);
		
		double Ntracks = prof_times["Ntracks"].real; // a special entry keeps count of tracks in "real" slot
		
		cout<<endl;
		cout<<"Printing profiling info for track fitter ---("<<Ntracks<<" tracks)---"<<endl;
		map<string, prof_time::time_diffs>::iterator iter = prof_times.begin();
		for(; iter!=prof_times.end(); iter++){
			if(iter->first == "Ntracks")continue; // skip Ntracks which is special
			cout<< " "<<iter->first<<endl
				<<"   real="<< iter->second.real/Ntracks*1000.0<<"ms"<<endl
				<<"   prof="<< iter->second.prof/Ntracks*1000.0<<"ms"<<endl
				<<"   virt="<< iter->second.virt/Ntracks*1000.0<<"ms"<<endl
				<< endl;
		}
		cout<<endl;
	}

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_trk_profile::fini(void)
{

	return NOERROR;
}

