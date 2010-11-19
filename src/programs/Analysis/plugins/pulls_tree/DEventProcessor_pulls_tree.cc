// $Id$
//
//    File: DEventProcessor_pulls_tree.cc
// Created: Fri Feb 19 13:04:29 EST 2010
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//

#include <PID/DChargedTrack.h>
#include <TRACKING/DTrackTimeBased.h>
#include <TRACKING/DTrackWireBased.h>
#include <TRACKING/DMCThrown.h>

#include "DEventProcessor_pulls_tree.h"
using namespace jana;


// Routine used to create our DEventProcessor
#include <JANA/JApplication.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_pulls_tree());
}
} // "C"


//------------------
// DEventProcessor_pulls_tree (Constructor)
//------------------
DEventProcessor_pulls_tree::DEventProcessor_pulls_tree()
{
	pthread_mutex_init(&mutex, NULL);
}

//------------------
// ~DEventProcessor_pulls_tree (Destructor)
//------------------
DEventProcessor_pulls_tree::~DEventProcessor_pulls_tree()
{

}

//------------------
// init
//------------------
jerror_t DEventProcessor_pulls_tree::init(void)
{
	pthread_mutex_lock(&mutex);
	
	pullWB_ptr = &pullWB;
	pullTB_ptr = &pullTB;

	pullsWB = new TTree("pullsWB","Wire-based hits");
	pullsWB->Branch("W","pull_t",&pullWB_ptr);

	pullsTB = new TTree("pullsTB","Time-based hits");
	pullsTB->Branch("W","pull_t",&pullTB_ptr);

	pthread_mutex_unlock(&mutex);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_pulls_tree::brun(JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_pulls_tree::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DChargedTrack*> tracks;
	const DMCThrown * thrown;

	try{
		loop->Get(tracks);
		loop->GetSingle(thrown);
	}catch(...){
		//return NOERROR;
	}
	pthread_mutex_lock(&mutex);
	loop->GetJApplication()->Lock();
	
	if(thrown){
		const DVector3 &x = thrown->momentum();
		TVector3 x_tvec(x.X(), x.Y(), x.Z());
		pullWB.pthrown = pullTB.pthrown = x_tvec;
	}else{
		pullWB.pthrown.SetXYZ(0,0,0);
		pullTB.pthrown.SetXYZ(0,0,0);
	}
	
	for(unsigned int i=0; i<tracks.size(); i++){
		if(tracks[i]->hypotheses.size()>0){
			const DTrackTimeBased *tbtrk = tracks[i]->hypotheses[0];
			if(tbtrk){
				pullTB.eventnumber = eventnumber;
				pullTB.trk_chisq = tbtrk->chisq;
				pullTB.trk_Ndof = tbtrk->Ndof;

				for(unsigned int j=0; j<tbtrk->pulls.size(); j++){
					pullTB.resi = tbtrk->pulls[j].resi;
					pullTB.err  = tbtrk->pulls[j].err;
					pullTB.s		= tbtrk->pulls[j].s;
					pullTB.pull = pullTB.resi/pullTB.err;
					pullsTB->Fill();
				}
				
				const DTrackWireBased* wbtrk;
				tbtrk->GetSingle(wbtrk);

				if(wbtrk){
					pullWB.eventnumber = eventnumber;
					pullWB.trk_chisq = wbtrk->chisq;
					pullWB.trk_Ndof = wbtrk->Ndof;

					for(unsigned int j=0; j<wbtrk->pulls.size(); j++){
						pullWB.resi = wbtrk->pulls[j].resi;
						pullWB.err  = wbtrk->pulls[j].err;
						pullWB.s		= wbtrk->pulls[j].s;
						pullWB.pull = pullWB.resi/pullWB.err;
						pullsWB->Fill();
					}
				}
			}
		}
	}
	
	loop->GetJApplication()->Unlock();
	pthread_mutex_unlock(&mutex);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_pulls_tree::erun(void)
{
	// Any final calculations on histograms (like dividing them)
	// should be done here. This may get called more than once.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_pulls_tree::fini(void)
{
	// Called at very end. This will be called only once

	return NOERROR;
}

