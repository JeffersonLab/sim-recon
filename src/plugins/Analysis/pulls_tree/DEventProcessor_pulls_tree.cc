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

	fitter = NULL;
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
	RECALCULATE_CHISQ = false;
	
	gPARMS->SetDefaultParameter("RECALCULATE_CHISQ", RECALCULATE_CHISQ, "Recalculate Chisq, Ndof, and pulls based on track parameters and hits rather than use what's recorded in track objects");

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_pulls_tree::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DChargedTrack*> tracks;
	const DMCThrown * thrown_single=NULL;
	fitter=NULL;

	try{
		loop->Get(tracks);
		loop->GetSingle(thrown_single);
		loop->GetSingle(fitter, "ALT1"); // use ALT1 as the stand-alone chisq calculator
	}catch(...){
		//return NOERROR;
	}
	pthread_mutex_lock(&mutex);
	loop->GetJApplication()->Lock();
	
	for(unsigned int i=0; i<tracks.size(); i++){

		// Initialize in case these don't get set below
		pullWB.pthrown.SetXYZ(0,0,0);
		pullTB.pthrown.SetXYZ(0,0,0);

		if(tracks[i]->hypotheses.size()>0){
			const DTrackTimeBased *tbtrk = tracks[i]->hypotheses[0];
			if(tbtrk){
				
				// Try filling in pthrown. First, by looking for an
				// associated object, then by checking if thrown_single
				// is not NULL. The latter case being for single
				// track events, the former for when DTrackTimeBased::THROWN
				// is used.
				const DMCThrown *thrown;
				tbtrk->GetSingle(thrown);
				if(!thrown)thrown = thrown_single;
				if(thrown){
					const DVector3 &x = thrown->momentum();
					pullTB.pthrown.SetXYZ(x.X(), x.Y(), x.Z());
				}
				
				// Copy chisq/Ndof/pull values to local variables in case we need to
				// replace them with ones we recalculate
				double chisq = tbtrk->chisq;
				int Ndof = tbtrk->Ndof;
				vector<DTrackFitter::pull_t> pulls = tbtrk->pulls;
				
				// Optionally replace chisq values
				if(RECALCULATE_CHISQ){
					RecalculateChisq(DTrackFitter::kTimeBased, tbtrk, chisq, Ndof, pulls);
				}
				
				pullTB.eventnumber = eventnumber;
				pullTB.trk_chisq = chisq;
				pullTB.trk_Ndof = Ndof;

				for(unsigned int j=0; j<pulls.size(); j++){
					pullTB.resi = pulls[j].resi;
					pullTB.err  = pulls[j].err;
					pullTB.s		= pulls[j].s;
					pullTB.pull = pullTB.resi/pullTB.err;
					pullsTB->Fill();
				}
				
				const DTrackWireBased* wbtrk;
				tbtrk->GetSingle(wbtrk);

				if(wbtrk){
					pullWB.eventnumber = eventnumber;
					pullWB.trk_chisq = wbtrk->chisq;
					pullWB.trk_Ndof = wbtrk->Ndof;

					// See comment above for TB tracks
					const DMCThrown *thrown;
					wbtrk->GetSingle(thrown);
					if(!thrown)thrown = thrown_single;
					if(thrown){
						const DVector3 &x = thrown->momentum();
						pullWB.pthrown.SetXYZ(x.X(), x.Y(), x.Z());
					}

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

//------------------
// RecalculateChisq
//------------------
void DEventProcessor_pulls_tree::RecalculateChisq(DTrackFitter::fit_type_t fit_type, const DKinematicData *kd, double &chisq, int &Ndof, vector<DTrackFitter::pull_t> &pulls)
{
	if(!fitter){
		_DBG_<<"DTrackFitter:ALT1 not available when RECALCULATE_CHISQ is set!"<<endl;
		return;
	}

	
}

