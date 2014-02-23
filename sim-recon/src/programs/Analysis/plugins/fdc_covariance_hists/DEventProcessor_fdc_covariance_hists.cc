// $Id$
//
//    File: DEventProcessor_fdc_covariance_hists.cc
// Created: Mon Apr 20 10:18:30 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 Darwin Kernel Version 9.6.0)
//

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include "DEventProcessor_fdc_covariance_hists.h"

#include <TROOT.h>

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DMCTrackHit.h>
#include <TRACKING/DMCTrajectoryPoint.h>
#include <FDC/DFDCGeometry.h>
#include <FDC/DFDCPseudo.h>
#include <HDGEOMETRY/DGeometry.h>
#include <DVector2.h>
#include <particleType.h>


// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_fdc_covariance_hists());
}
} // "C"


//------------------
// DEventProcessor_fdc_covariance_hists
//------------------
DEventProcessor_fdc_covariance_hists::DEventProcessor_fdc_covariance_hists()
{
	pthread_mutex_init(&mutex, NULL);
	
	Nevents = 0;
}

//------------------
// ~DEventProcessor_fdc_covariance_hists
//------------------
DEventProcessor_fdc_covariance_hists::~DEventProcessor_fdc_covariance_hists()
{

}

//------------------
// init
//------------------
jerror_t DEventProcessor_fdc_covariance_hists::init(void)
{
	// Create TRACKING directory
	TDirectory *dir = (TDirectory*)gROOT->FindObject("TRACKING");
	if(!dir)dir = new TDirectoryFile("TRACKING","TRACKING");
	dir->cd();

	fdc_cov = new TProfile2D("fdc_cov","FDC Covariance calculated from residuals", 24, 0.5, 24.5, 24, 0.5, 24.5);	
	fdc_cov->SetStats(0);
	fdc_cov->SetXTitle("Layer number");
	fdc_cov->SetYTitle("Layer number");
	fdc_cov_calc = (TProfile2D*)fdc_cov->Clone("fdc_cov_calc");
	fdc_cov_calc->SetTitle("FDC Covariance calculated from materials");

	dir->cd("../");
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_fdc_covariance_hists::brun(JEventLoop *loop, int runnumber)
{	
	DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	if(!dapp){
		_DBG_<<"Cannot get DApplication from JEventLoop! (are you using a JApplication based program perhaps?)"<<endl;
		return RESOURCE_UNAVAILABLE;
	}
	bfield=dapp->GetBfield();
	rt = new DReferenceTrajectory(bfield);
	rt->SetDRootGeom(dapp->GetRootGeom());
	
	// Get z-position of most upstream FDC layer
	vector<double> z_wires;
	dapp->GetDGeometry(runnumber)->GetFDCZ(z_wires);
	Z_fdc1 = z_wires[0];

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_fdc_covariance_hists::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_fdc_covariance_hists::fini(void)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_fdc_covariance_hists::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DMCThrown*> mcthrowns;
	vector<const DMCTrackHit*> mctrackhits;
	vector<const DMCTrajectoryPoint*> mctrajectorypoints;
	vector<const DFDCPseudo*> fdcpseudohits;	
	
	loop->Get(mcthrowns);
	loop->Get(mctrackhits);
	loop->Get(mctrajectorypoints);
	loop->Get(fdcpseudohits);
	
	Nevents++;
	
	// Only look at events with one thrown and one reconstructed particle
	if(mcthrowns.size() !=1){
		_DBG_<<" mcthrowns.size()="<<mcthrowns.size()<<endl;
		return NOERROR;
	}
	
	// Look for truth point corresponding to hit of upstream most
	// FDC layer. If we don't find one, skip this event.
	const DMCTrackHit *mctrackhit1 = NULL;
	for(unsigned int i=0; i<mctrackhits.size(); i++){
		const DMCTrackHit *hit = mctrackhits[i];
		if(hit->system == SYS_FDC){
			if((hit->z < (Z_fdc1+0.5)) && (hit->z > (Z_fdc1-0.5))){
				mctrackhit1 = hit;
				break;
			}
		}
	}
	if(!mctrackhit1)return NOERROR;
	DVector3 pos_mctrackhit1;
	pos_mctrackhit1.SetXYZ(mctrackhit1->r*cos(mctrackhit1->phi), mctrackhit1->r*sin(mctrackhit1->phi), mctrackhit1->z);

	// Look for trajectory point closest to first FDC plane.
	// Simultaneously, record distance to first layer.
	DVector3 pos_fdc1;
	DVector3 mom_fdc1;
	double s1 = 0.0;
	double diff_min = 1.0E6;
	double s1_min = 0.0;
	for(unsigned int i=0; i<mctrajectorypoints.size(); i++){
		const DMCTrajectoryPoint *traj = mctrajectorypoints[i];
		DVector3 pos_traj(traj->x, traj->y, traj->z);
		double diff = (pos_traj-pos_mctrackhit1).Mag();
		s1 += (double)traj->step;
		if(diff<diff_min){
			pos_fdc1.SetXYZ(traj->x, traj->y, traj->z);
			mom_fdc1.SetXYZ(traj->px, traj->py, traj->pz);
			s1_min = s1;
			diff_min = diff;
		}
	}
	s1 = s1_min;

	// Lock mutex
	pthread_mutex_lock(&mutex);

	rt->Swim(pos_fdc1, mom_fdc1);

	// Loop over FDC hits
	vector<double> resi_by_layer(24,-1000);
	vector<const DReferenceTrajectory::swim_step_t*> step_by_layer(24, (const DReferenceTrajectory::swim_step_t*)NULL);
	for(unsigned int i=0; i<fdcpseudohits.size(); i++){
		const DFDCPseudo *fdcpseudohit = fdcpseudohits[i];

		double s;
		DVector3 pos_doca, mom_doca;
		double doca = rt->DistToRT(fdcpseudohit->wire, &s);		
		rt->GetLastDOCAPoint(pos_doca, mom_doca);

		// Estimate TOF assuming pion
		double mass = 0.13957;
		double beta = 1.0/sqrt(1.0 + pow(mass/mom_doca.Mag(), 2.0))*2.998E10;
		double tof = (s+s1)/beta/1.0E-9; // in ns
		double dist = (fdcpseudohit->time - tof)*55E-4;
		double resi = dist - doca;
		
		int layer = fdcpseudohit->wire->layer;
		if(layer>=1 && layer<=24){
			resi_by_layer[layer-1] = resi;
			step_by_layer[layer-1] = rt->GetLastSwimStep();
		}
	}

	// Fill direct covariance matrix
	for(int layer1=1; layer1<=24; layer1++){
		double resi1 = resi_by_layer[layer1-1];
		if(!finite(resi1) || fabs(resi1)>100.0)continue;
		
		for(int layer2=layer1; layer2<=24; layer2++){
			double resi2 = resi_by_layer[layer2-1];
			if(!finite(resi2) || fabs(resi2)>100.0)continue;

			fdc_cov->Fill(layer1, layer2, resi1*resi2);
			fdc_cov->Fill(layer2, layer1, resi1*resi2);
		}
	}

	// Fill covariance matrix calculated from materials
	const DReferenceTrajectory::swim_step_t* step0 = step_by_layer[0];
	if(step0){
		for(int layerA=1; layerA<=24; layerA++){
			const DReferenceTrajectory::swim_step_t* stepA = step_by_layer[layerA-1];
			if(!stepA)continue;
			
			for(int layerB=layerA; layerB<=24; layerB++){
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

				fdc_cov_calc->Fill(layerA, layerB, sigmaAB);
				fdc_cov_calc->Fill(layerB, layerA, sigmaAB);
			}
		}
	}

	// Unlock mutex
	pthread_mutex_unlock(&mutex);
	
	return NOERROR;
}

