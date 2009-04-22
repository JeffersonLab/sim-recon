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

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DMCTrackHit.h>
#include <TRACKING/DMCTrajectoryPoint.h>
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DTrackHitSelectorTHROWN.h>
#include <PID/DParticle.h>
#include <FDC/DFDCGeometry.h>
#include <CDC/DCDCTrackHit.h>
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
	fdchit_ptr = &fdchit;

	pthread_mutex_init(&mutex, NULL);
	
	NLRgood = 0;
	NLRbad = 0;
	NLRfit_unknown = 0;
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

	fdchits = new TTree("fdchit","FDC hits");
	fdchits->Branch("fdchit","dchit",&fdchit_ptr);
	
	fdc_cov = new TProfile2D("fdc_cov","", 24, 0.5, 24.5, 24, 0.5, 24.5);	
	fdc_cath_cov = new TProfile2D("fdc_cath_cov","", 24, 0.5, 24.5, 24, 0.5, 24.5);	
	//fdc_cov_numerator = new TH2D("fdc_cov_numerator","", 24, 0.5, 24.5, 24, 0.5, 24.5);
	//fdc_cov_denominator = (TH2D*)fdc_cov_numerator->Clone("fdc_cov_denominator");

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
	lorentz_def=dapp->GetLorentzDeflections();
	bfield=dapp->GetBfield();
	rt = new DReferenceTrajectory(bfield);
	
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
	char str[256];
	sprintf(str,"%3.4f%%", 100.0*(double)NLRbad/(double)(NLRbad+NLRgood));
	
	cout<<endl<<setprecision(4);
	cout<<"       NLRgood: "<<NLRgood<<endl;
	cout<<"        NLRbad: "<<NLRbad<<endl;
	cout<<"     NLRfit==0: "<<NLRfit_unknown<<endl;
	cout<<"Percentage bad: "<<str<<endl;
	cout<<"       Nevents: "<<Nevents<<endl;
	cout<<endl;

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_fdc_covariance_hists::evnt(JEventLoop *loop, int eventnumber)
{
	//vector<const DParticle*> particles;
	//vector<const DTrackCandidate*> candidates;
	vector<const DMCThrown*> mcthrowns;
	vector<const DMCTrackHit*> mctrackhits;
	vector<const DMCTrajectoryPoint*> mctrajectorypoints;
	vector<const DFDCPseudo*> fdcpseudohits;	
	
	//loop->Get(particles);
	//loop->Get(candidates);
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
	
	// Look for trajectory point closest to first FDC plane.
	DVector3 pos_fdc1;
	DVector3 mom_fdc1;
	for(unsigned int i=0; i<mctrajectorypoints.size(); i++){
		const DMCTrajectoryPoint *traj = mctrajectorypoints[i];
		if(traj->z > Z_fdc1){
			pos_fdc1.SetXYZ(traj->x, traj->y, traj->z);
			mom_fdc1.SetXYZ(traj->px, traj->py, traj->pz);
			break;
		}
	}

	// Lock mutex
	pthread_mutex_lock(&mutex);

	rt->Swim(pos_fdc1, mom_fdc1);

	// Loop over FDC hits
	vector<double> resi_by_layer(24,-1000);
	vector<double> resic_by_layer(24,-1000);
	for(unsigned int i=0; i<fdcpseudohits.size(); i++){
		const DFDCPseudo *fdcpseudohit = fdcpseudohits[i];

		// Lorentz corrected poisition along the wire is contained in x,y values.
		//DVector3 wpos(fdcpseudohit->x, fdcpseudohit->y, fdcpseudohit->wire->origin.Z());
		//DVector3 wdiff = wpos - fdcpseudohit->wire->origin;
		//double u_corr = fdcpseudohit->wire->udir.Dot(wdiff);

		// The hit info structure is used to pass info both in and out of FindLR()
		hit_info_t hit_info;
		hit_info.rt = rt;
		hit_info.wire = fdcpseudohit->wire;
		hit_info.tdrift = fdcpseudohit->time;
		hit_info.FindLR(mctrackhits, lorentz_def);
		
		fdchit.eventnumber = eventnumber;
		fdchit.wire = fdcpseudohit->wire->wire;
		fdchit.layer = fdcpseudohit->wire->layer;
		fdchit.t = fdcpseudohit->time;
		fdchit.tof = hit_info.tof;
		fdchit.doca = hit_info.doca;
		fdchit.resi = hit_info.dist - hit_info.doca;
		fdchit.resic = hit_info.u - (fdcpseudohit->s + hit_info.u_lorentz) ;
		fdchit.LRis_correct = hit_info.LRis_correct;
		fdchit.LRfit = hit_info.LRfit;
		fdchit.pos_doca = hit_info.pos_doca;
		fdchit.pos_wire = hit_info.pos_wire;
		
		fdchits->Fill();
		
		if(fdchit.layer>=1 && fdchit.layer<=24){
			resi_by_layer[fdchit.layer-1] = fdchit.resi;
			resic_by_layer[fdchit.layer-1] = fdchit.resic;
		}

		if(fdchit.LRfit!=0){
			if(fdchit.LRis_correct){
				NLRgood++;
			}else{
				NLRbad++;
			}
		}else{
			NLRfit_unknown++;
		}
	}
	
	// Fill covariance matrices
	for(int layer1=1; layer1<=24; layer1++){
		double resi1 = resi_by_layer[layer1-1];
		double resic1 = resic_by_layer[layer1-1];
		if(!finite(resi1) || fabs(resi1)>100.0)continue;
		if(!finite(resic1) || fabs(resic1)>100.0)continue;
		
//_DBG_<<"layer -- "<<layer1<<"  resi="<<resi1<<"  resic="<<resic1<<endl;
		for(int layer2=layer1; layer2<=24; layer2++){
			double resi2 = resi_by_layer[layer2-1];
			double resic2 = resic_by_layer[layer2-1];
			if(!finite(resi2) || fabs(resi2)>100.0)continue;
			if(!finite(resic2) || fabs(resic2)>100.0)continue;

			fdc_cov->Fill(layer1, layer2, resi1*resi2);
			fdc_cath_cov->Fill(layer1, layer2, resic1*resic2);
		}
	}

	// Unlock mutex
	pthread_mutex_unlock(&mutex);
	
	return NOERROR;
}

//------------------
// FindLR
//------------------
void DEventProcessor_fdc_covariance_hists::hit_info_t::FindLR(vector<const DMCTrackHit*> &mctrackhits, const DLorentzDeflections *lorentz_def)
{
	/// Decided on whether or not the reference trajectory is on the correct side of the wire
	/// based on truth information.
	///
	/// Upon entry, the rt, wire, and rdrift members should be set. The remaining fields are
	/// filled in on return from this method.
	///
	/// Essentially, this will look through the DMCTrackHit objects via the
	/// DTrackHitSelectorTHROWN::GetMCTrackHit method to find the truth point corresponding
	/// to this hit (if any). It then compares the vector pointing from the point of
	/// closest approach on the wire to the DOCA point so a similar vector pointing
	/// from the same place on the wire to the truth point. If the 2 vectors are within
	/// +/- 90 degrees, then the trajectory is said to be on the correct side of the wire.
	
	double s;
	doca = rt->DistToRT(wire, &s);
	rt->GetLastDOCAPoint(pos_doca, mom_doca);
	DVector3 shift = wire->udir.Cross(mom_doca);
	u = rt->GetLastDistAlongWire();
	pos_wire = wire->origin + u*wire->udir;

	// Estimate TOF assuming pion
	double mass = 0.13957;
	double beta = 1.0/sqrt(1.0 + pow(mass/mom_doca.Mag(), 2.0))*2.998E10;
	tof = s/beta/1.0E-9; // in ns
	dist = (tdrift - tof)*55E-4;
	
	// Find the Lorentz correction based on current track (if applicable)
	if(lorentz_def){
		DVector3 shift = wire->udir.Cross(mom_doca);
		shift.SetMag(1.0);
		double LRsign = shift.Dot(pos_doca-pos_wire)<0.0 ? +1.0:-1.0;
		double alpha = mom_doca.Angle(DVector3(0,0,1));
		u_lorentz = LRsign*lorentz_def->GetLorentzCorrection(pos_doca.X(), pos_doca.Y(), pos_doca.Z(), alpha, dist);
	}else{
		u_lorentz = 0.0;
	}
	
	// Look for a truth hit corresponding to this wire. If none is found, mark the hit as 
	// 0 (i.e. neither left nor right) and continue to the next hit.
	const DMCTrackHit *mctrackhit = DTrackHitSelectorTHROWN::GetMCTrackHit(wire, dist, mctrackhits);
	if(!mctrackhit){
		LRfit = 0;
		LRis_correct = false; // can't really tell what to set this to
	}else{
		DVector3 pos_truth(mctrackhit->r*cos(mctrackhit->phi), mctrackhit->r*sin(mctrackhit->phi), mctrackhit->z);
		DVector3 pos_diff_truth = pos_truth-pos_wire;
		DVector3 pos_diff_traj  = pos_doca-pos_wire;
		
		LRfit = shift.Dot(pos_diff_traj)<0.0 ? -1:1;
		LRis_correct = pos_diff_truth.Dot(pos_diff_traj)>0.0;
	}
}

