// $Id: $
//
//    File: DEventProcessor_cdc_hists.cc
//

#include <iostream>
#include <cmath>
using namespace std;

#include <TThread.h>

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DEventProcessor_cdc_hists.h"

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DMCTrackHit.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <CDC/DCDCHit.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCHit.h>
#include <DVector2.h>


// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_cdc_hists());
}
} // "C"


//------------------
// DEventProcessor_cdc_hists
//------------------
DEventProcessor_cdc_hists::DEventProcessor_cdc_hists()
{	
	cdc_ptr = &cdc;
	cdchit_ptr = &cdchit;
	
	pthread_mutex_init(&mutex, NULL);
}

//------------------
// ~DEventProcessor_cdc_hists
//------------------
DEventProcessor_cdc_hists::~DEventProcessor_cdc_hists()
{
}

//------------------
// init
//------------------
jerror_t DEventProcessor_cdc_hists::init(void)
{
	// open ROOT file (if needed)
	//if(ROOTfile != NULL) ROOTfile->cd();
	
	// Create THROWN directory
	TDirectory *dir = new TDirectoryFile("CDC","CDC");
	dir->cd();

	// Create Tree
	cdctree = new TTree("cdc","CDC Truth points");
	cdchittree = new TTree("cdchit","CDC Hits");
	cdcbranch = cdctree->Branch("T","CDC_branch",&cdc_ptr);
	cdchitbranch = cdchittree->Branch("H","CDChit_branch",&cdchit_ptr);
	
	idEdx = new TH1D("idEdx","Integrated dE/dx in CDC", 10000, 0.0, 1.0E-3);
	idEdx->SetXTitle("dE/dx (GeV/cm)");
	idEdx_vs_p = new TH2D("idEdx_vs_p","Integrated dE/dx vs. momentum in CDC", 100, 0.0, 1.0, 1000, 0.0, 1.0E1);
	idEdx->SetXTitle("momentum (GeV/c)");
	idEdx->SetYTitle("dE/dx (MeV/g^{-1}cm^{2})");

	// Go back up to the parent directory
	dir->cd("../");
	
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_cdc_hists::brun(JEventLoop *eventLoop, int runnumber)
{
	DApplication *dapp = dynamic_cast<DApplication*>(app);
	bfield = dapp->GetBfield();
	rt = new DReferenceTrajectory(bfield);
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_cdc_hists::erun(void)
{

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_cdc_hists::fini(void)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_cdc_hists::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DMCTrackHit*> mctrackhits;
	vector<const DCDCHit*> cdchits;
	vector<const DCDCTrackHit*> cdctrackhits;
	vector<const DFDCHit*> fdchits;
	vector<const DMCThrown*> mcthrowns;
	loop->Get(mctrackhits);
	loop->Get(cdchits);
	loop->Get(cdctrackhits);
	loop->Get(fdchits);
	loop->Get(mcthrowns);
	
	// Find number of wire hits in FDC
	int Nfdc_wire_hits = 0;
	for(unsigned int i=0; i<fdchits.size(); i++)if(fdchits[i]->type==0)Nfdc_wire_hits++;
	
	// Lock mutex
	pthread_mutex_lock(&mutex);
	
	// Swim reference trajectory for first thrown track
	const DMCThrown *mcthrown = mcthrowns.size()>0 ? mcthrowns[0]:NULL;
	if(mcthrown){
		rt->Swim(mcthrown->position(), mcthrown->momentum(), mcthrown->charge());
	}
	
	// Loop over all truth hits, ignoring all but CDC hits
	for(unsigned int i=0; i<mctrackhits.size(); i++){
		const DMCTrackHit *mctrackhit = mctrackhits[i];
		if(mctrackhit->system != SYS_CDC)continue;
		
		double r = mctrackhit->r;
		double phi = mctrackhit->phi;
		double x = r*cos(phi);
		double y = r*sin(phi);
		cdc.pos_truth.SetXYZ(x, y, mctrackhit->z);
		
		cdctree->Fill();
	}

	// Loop over all real hits
	double dEtot = 0.0;
	double dxtot = 0.0;
	for(unsigned int i=0; i<cdchits.size(); i++){
		const DCDCHit *hit = cdchits[i];
		const DCDCWire *wire = (cdchits.size() == cdctrackhits.size()) ? cdctrackhits[i]->wire:NULL;
		
		cdchit.ring		= hit->ring;
		cdchit.straw	= hit->straw;
		cdchit.dE		= hit->dE;
		cdchit.dx		= 0.0;
		cdchit.t			= hit->t;
		cdchit.pthrown = mcthrowns.size()>0 ? mcthrowns[0]->momentum().Mag():-1000.0;
		cdchit.ncdchits	= (int)cdchits.size();
		cdchit.ntothits	= (int)cdchits.size() + Nfdc_wire_hits;
		
		if(mcthrown && wire){
			cdchit.dx = rt->Straw_dx(wire, 0.8);
		}
		
		if(cdchit.dx!=0.0){
			dEtot += cdchit.dE;
			dxtot += cdchit.dx;
		}
		
		// Find residual of hit with "thrown" track (if present)
		if(mcthrown && wire){
			double s;
			double doca = rt->DistToRT(wire, &s);
			double mass = 0.13957; // assume pion
			double beta = 1.0/sqrt(1.0 + pow(mass/mcthrown->momentum().Mag(), 2.0))*2.998E10;
			double tof = s/beta/1.0E-9; // in ns
			double dist = (cdchit.t-tof)*55.0E-4;
			cdchit.resi_thrown = doca-dist;
		}else{
			cdchit.resi_thrown = 0.0;
		}

		cdchittree->Fill();
	}
	
	// Fill dE/dx histograms
	if(((Nfdc_wire_hits+cdchits.size()) >= 10) && (cdchits.size()>=5)){
		if(dxtot>0.0){
			idEdx->Fill(dEtot/dxtot);
			if(mcthrown){
				// The CDC gas is 85% Ar, 15% CO2 by mass.
				// density of  Ar: 1.977E-3 g/cm^3
				// density of CO2: 1.66E-3 g/cm^3
				double density = 0.85*1.66E-3 + 0.15*1.977E-3;
				idEdx_vs_p->Fill(mcthrown->momentum().Mag(), dEtot/dxtot*1000.0/density);
			}
		}
	}

	// Unlock mutex
	pthread_mutex_unlock(&mutex);
	
	return NOERROR;
}
