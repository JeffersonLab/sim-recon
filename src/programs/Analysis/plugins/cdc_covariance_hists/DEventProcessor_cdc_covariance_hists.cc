// $Id$
//
//    File: DEventProcessor_cdc_covariance_hists.cc
// Created: Thu Apr 23 08:30:24 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 Darwin Kernel Version 9.6.0)
//

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include "DEventProcessor_cdc_covariance_hists.h"

#include <TROOT.h>

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DMCTrackHit.h>
#include <TRACKING/DMCTrajectoryPoint.h>
#include <CDC/DCDCTrackHit_factory.h>
#include <HDGEOMETRY/DGeometry.h>
#include <DVector2.h>
#include <particleType.h>


// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_cdc_covariance_hists());
}
} // "C"


//------------------
// DEventProcessor_cdc_covariance_hists
//------------------
DEventProcessor_cdc_covariance_hists::DEventProcessor_cdc_covariance_hists()
{
	pthread_mutex_init(&mutex, NULL);

	Nevents = 0;
}

//------------------
// ~DEventProcessor_cdc_covariance_hists
//------------------
DEventProcessor_cdc_covariance_hists::~DEventProcessor_cdc_covariance_hists()
{

}

//------------------
// init
//------------------
jerror_t DEventProcessor_cdc_covariance_hists::init(void)
{
	// Create TRACKING directory
	TDirectory *dir = (TDirectory*)gROOT->FindObject("TRACKING");
	if(!dir)dir = new TDirectoryFile("TRACKING","TRACKING");
	dir->cd();

	cdc_cov = new TProfile2D("cdc_cov","CDC Covariance calculated from residuals", 27, 0.5, 27.5, 27, 0.5, 27.5);
	cdc_cov->SetStats(0);
	cdc_cov->SetXTitle("Ring number");
	cdc_cov->SetYTitle("Ring number");
	cdc_cov_calc = (TProfile2D*)cdc_cov->Clone("cdc_cov_calc");
	cdc_cov_calc->SetTitle("CDC Covariance calculated from materials");

	dir->cd("../");
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_cdc_covariance_hists::brun(JEventLoop *loop, int runnumber)
{
	// We need to make a DReferenceTrajectory which means we need the B-field
	DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	if(!dapp){
		_DBG_<<"Cannot get DApplication from JEventLoop! (are you using a JApplication based program perhaps?)"<<endl;
		return RESOURCE_UNAVAILABLE;
	}
	bfield=dapp->GetBfield();
	rt = new DReferenceTrajectory(bfield);
	rt->SetDRootGeom(dapp->GetRootGeom());
	
	// Get radius of innermost CDC layer
	//const DCDCWire *wire1=DCDCTrackHit_factory::GetCDCWire(1,1);
	//R_cdc1 = wire1->origin.Perp();
	R_cdc1 = 10.960;

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_cdc_covariance_hists::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_cdc_covariance_hists::fini(void)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_cdc_covariance_hists::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DMCThrown*> mcthrowns;
	vector<const DMCTrackHit*> mctrackhits;
	vector<const DMCTrajectoryPoint*> mctrajectorypoints;
	vector<const DCDCTrackHit*> cdctrackhits;	
	
	loop->Get(mcthrowns);
	loop->Get(mctrackhits);
	loop->Get(mctrajectorypoints);
	loop->Get(cdctrackhits);
	
	Nevents++;
	
	// Only look at events with one thrown and one reconstructed particle
	if(mcthrowns.size() !=1){
		_DBG_<<" mcthrowns.size()="<<mcthrowns.size()<<endl;
		return NOERROR;
	}
	
	// Look for truth point corresponding to hit of innermost
	// CDC layer. If we don't find one, skip this event.
	const DMCTrackHit *mctrackhit1 = NULL;
	for(unsigned int i=0; i<mctrackhits.size(); i++){
		const DMCTrackHit *hit = mctrackhits[i];
		if(hit->system == SYS_CDC){
			if((hit->r < (R_cdc1+0.8)) && (hit->r > (R_cdc1-0.8))){
				mctrackhit1 = hit;
				break;
			}
		}
	}
	if(!mctrackhit1)return NOERROR;
	DVector3 pos_mctrackhit1;
	pos_mctrackhit1.SetXYZ(mctrackhit1->r*cos(mctrackhit1->phi), mctrackhit1->r*sin(mctrackhit1->phi), mctrackhit1->z);
	
	// Look for trajectory point closest to first CDC layer hit.
	// Simultaneously, record distance to first layer.
	DVector3 pos_cdc1;
	DVector3 mom_cdc1;
	double s1 = 0.0;
	double diff_min = 1.0E6;
	double s1_min = 0.0;
	for(unsigned int i=0; i<mctrajectorypoints.size(); i++){
		const DMCTrajectoryPoint *traj = mctrajectorypoints[i];
		DVector3 pos_traj(traj->x, traj->y, traj->z);
		double diff = (pos_traj-pos_mctrackhit1).Mag();
		s1 += (double)traj->step;
		if(diff<diff_min){
			pos_cdc1.SetXYZ(traj->x, traj->y, traj->z);
			mom_cdc1.SetXYZ(traj->px, traj->py, traj->pz);
			s1_min = s1;
			diff_min = diff;
		}
	}
	s1 = s1_min;

	// Lock mutex
	pthread_mutex_lock(&mutex);

	rt->Swim(pos_cdc1, mom_cdc1);

	// Loop over CDC hits
	vector<double> resi_by_layer(27,-1000);
	vector<const DReferenceTrajectory::swim_step_t*> step_by_layer(27, (const DReferenceTrajectory::swim_step_t*)NULL);
	for(unsigned int i=0; i<cdctrackhits.size(); i++){
		const DCDCTrackHit *cdctrackhit = cdctrackhits[i];

		double s;
		DVector3 pos_doca, mom_doca;
		double doca = rt->DistToRT(cdctrackhit->wire, &s);		
		rt->GetLastDOCAPoint(pos_doca, mom_doca);

		// Estimate TOF assuming pion
		double mass = 0.13957;
		double beta = 1.0/sqrt(1.0 + pow(mass/mom_doca.Mag(), 2.0))*2.998E10;
		double tof = (s+s1)/beta/1.0E-9; // in ns
		double dist = (cdctrackhit->tdrift - tof)*55E-4;
		double resi = dist - doca;
		
		int layer = cdctrackhit->wire->ring;
		if(layer>=1 && layer<=27){
			resi_by_layer[layer-1] = resi;
			step_by_layer[layer-1] = rt->GetLastSwimStep();
		}
	}
	
	// Fill direct covariance matrix
	for(int layer1=1; layer1<=27; layer1++){
		double resi1 = resi_by_layer[layer1-1];
		if(!finite(resi1) || fabs(resi1)>100.0)continue;
		
		for(int layer2=layer1; layer2<=27; layer2++){
			double resi2 = resi_by_layer[layer2-1];
			if(!finite(resi2) || fabs(resi2)>100.0)continue;

			cdc_cov->Fill(layer1, layer2, resi1*resi2);
			cdc_cov->Fill(layer2, layer1, resi1*resi2);
		}
	}

	// Fill covariance matrix calculated from materials
	const DReferenceTrajectory::swim_step_t* step0 = step_by_layer[0];
	if(step0){
		for(int layerA=1; layerA<=27; layerA++){
			const DReferenceTrajectory::swim_step_t* stepA = step_by_layer[layerA-1];
			if(!stepA)continue;
			
			for(int layerB=layerA; layerB<=27; layerB++){
				const DReferenceTrajectory::swim_step_t* stepB = step_by_layer[layerB-1];
				if(!stepB)continue;

				// Correlations between A and B are due only to material between
				// the first detector and the most upstream of A or B.
				double sA = stepA->s;
				double sB = stepB->s;
				const DReferenceTrajectory::swim_step_t *step_end = sA<sB ? stepA:stepB;

				if(step0->s>step_end->s)continue; // Bullet proof
				
				double itheta02 = step_end->itheta02 - step0->itheta02;
				double itheta02s = step_end->itheta02s - step0->itheta02s;
				double itheta02s2 = step_end->itheta02s2 - step0->itheta02s2;
				double sigmaAB = sA*sB*itheta02 -(sA+sB)*itheta02s + itheta02s2;

				cdc_cov_calc->Fill(layerA, layerB, sigmaAB);
				cdc_cov_calc->Fill(layerB, layerA, sigmaAB);
			}
		}
	}

	// Unlock mutex
	pthread_mutex_unlock(&mutex);
	
	return NOERROR;
}
