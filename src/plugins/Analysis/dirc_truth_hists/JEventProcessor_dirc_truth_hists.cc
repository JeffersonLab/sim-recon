// $Id$
//
//    File: JEventProcessor_dirc_truth_hists.cc
// Created: Thu Mar 19 10:25:58 EDT 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "JEventProcessor_dirc_truth_hists.h"
using namespace std;

#include "TLorentzVector.h"

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_dirc_truth_hists());
}
} // "C"


//------------------
// JEventProcessor_dirc_truth_hists (Constructor)
//------------------
JEventProcessor_dirc_truth_hists::JEventProcessor_dirc_truth_hists()
{

}

//------------------
// ~JEventProcessor_dirc_truth_hists (Destructor)
//------------------
JEventProcessor_dirc_truth_hists::~JEventProcessor_dirc_truth_hists()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_dirc_truth_hists::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	
	japp->RootWriteLock();
	{
		event = 0;

		// initialize tree
		dircTree = new TTree("dircTree", "dircTree");
		dircTree->Branch("event",&event);
		dircTree->Branch("nGenTrk",&nGenTrk);
		dircTree->Branch("genTrk_px",m_genTrk_px,Form("genTrk_px[%d]/F",kGenTrkMax));
		dircTree->Branch("genTrk_py",m_genTrk_py,Form("genTrk_py[%d]/F",kGenTrkMax));
		dircTree->Branch("genTrk_pz",m_genTrk_pz,Form("genTrk_pz[%d]/F",kGenTrkMax));
		dircTree->Branch("genTrk_E",m_genTrk_E,Form("genTrk_E[%d]/F",kGenTrkMax));
		dircTree->Branch("genTrk_x",m_genTrk_x,Form("genTrk_x[%d]/F",kGenTrkMax));
		dircTree->Branch("genTrk_y",m_genTrk_y,Form("genTrk_y[%d]/F",kGenTrkMax));
		dircTree->Branch("genTrk_z",m_genTrk_z,Form("genTrk_z[%d]/F",kGenTrkMax));
		dircTree->Branch("genTrk_t",m_genTrk_t,Form("genTrk_t[%d]/F",kGenTrkMax));

		dircTree->Branch("nTruthHit",&nTruthHit);
		dircTree->Branch("hit_x",m_hit_x,Form("hit_x[%d]/F",kTruthHitMax));
		dircTree->Branch("hit_y",m_hit_y,Form("hit_y[%d]/F",kTruthHitMax));
		dircTree->Branch("hit_z",m_hit_z,Form("hit_z[%d]/F",kTruthHitMax));
		dircTree->Branch("hit_t",m_hit_t,Form("hit_t[%d]/F",kTruthHitMax));
		dircTree->Branch("hitpixel_w",m_hitpixel_w,Form("hitpixel_w[%d]/F",kTruthHitMax));
		dircTree->Branch("hitpixel_y",m_hitpixel_y,Form("hitpixel_y[%d]/F",kTruthHitMax));

		hTruthPoint = new TH2F("hTruthPoint","hTruthPoint Y vx X; X; Y", 240, -120., 120., 200, -200., 200.);
		hTruthPointM = new TH1F("hTruthPointM","hTruthPoint Mass; M", 100, 0., 1.);
		hTruthHitT = new TH1F("hTruthHitT","hTruthHit Time; T", 200, 0, 200);
		hTruthHitDeltaT = new TH1F("hTruthHitDeltaT","hTruthHit #Delta Time; #DeltaT", 200, 0, 200);
		hTruthHitE = new TH1F("hTruthHitE","hTruthHit Energy; E", 200, 0, 10);
		hTruthHitLambda = new TH1F("hTruthHitLambda","hTruthHit Wavelength; #lambda", 200, 0, 1000);
		hTruthHitXY = new TH2F("hTruthHitXY","hTruthHit Y vx X; X; Y", 200, 218., 240., 240, -120., 120.);
		hTruthHitZX = new TH2F("hTruthHitZX","hTruthHit X vx Z; Z; X", 200, 541., 564., 200, 218., 240.);
		hTruthHitYLocW = new TH2F("hTruthHitYLocW","hTruthHit LocW vx Y; Y; LocW", 400, -120., 120., 60, -1.8, 34.2); // 6 x 6 mm pixel size
		hTruthHitMissingYLocW = new TH2F("hTruthHitMissingYLocW","hTruthHitMissing LocW vx Y; Y; LocW", 400, -120., 120., 60, -1.8, 34.2); // 6 x 6 mm pixel size

		hTruthHitIncidentAngleY = new TH2F("hTruthHitIncidentAngleY", "hTruthHitIncidentAngleY", 400, -120, 120, 180, 0., 45.);
		hTruthHitIncidentAngleLocW = new TH2F("hTruthHitIncidentAngleLocW", "hTruthHitIncidentAngleLocW", 60, -1.8, 34.2, 180, 0., 45.);
		
		hTruthHitYLocWT = new TH3F("hTruthHitYLocWT","hTruthHit LocW vs Y vs T; Y; LocW; T", 400, -120., 120., 60, -1.8, 34.2, 200, 0, 200);

		hNTruthHit = new TH2F("hNTruthHit", "hNTruthHit; # truth hits; # truth points", 100, 0, 100, 10, 0, 10);
	}
	japp->RootUnLock();

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_dirc_truth_hists::brun(JEventLoop *eventLoop, int runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_dirc_truth_hists::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	vector<const DDIRCTruthPoint*> locDIRCTruthPoint;
	vector<const DDIRCTruthHit*> locDIRCTruthHit;
	loop->Get(locDIRCTruthPoint);
	loop->Get(locDIRCTruthHit);

	nGenTrk = 0;
	nTruthHit = 0;

	for(unsigned int i=0; i<kGenTrkMax; i++){
                m_genTrk_px[i] = -999.;
                m_genTrk_py[i] = -999.;
		m_genTrk_pz[i] = -999.;
		m_genTrk_E[i] = -999.;
		m_genTrk_x[i] = -999.;
		m_genTrk_y[i] = -999.;
		m_genTrk_z[i] = -999.;
		m_genTrk_t[i] = -999.;
        }

	for(unsigned int i=0; i<kTruthHitMax; i++){
                m_hit_x[i] = -999.;
                m_hit_y[i] = -999.;
		m_hit_z[i] = -999.;
		m_hit_t[i] = -999.;

		m_hitpixel_w[i] = -999.;
                m_hitpixel_y[i] = -999.;
        }

	// truth hits are interactions of tracks with front of the DIRC bar
	double truthTime = 0.;
	for(size_t loc_i = 0; loc_i < locDIRCTruthPoint.size(); loc_i++) {
			
		TLorentzVector truthHit(locDIRCTruthPoint[loc_i]->px, locDIRCTruthPoint[loc_i]->py, locDIRCTruthPoint[loc_i]->pz, locDIRCTruthPoint[loc_i]->E);
		TLorentzVector truthHit_position(locDIRCTruthPoint[loc_i]->x, locDIRCTruthPoint[loc_i]->y, locDIRCTruthPoint[loc_i]->z, locDIRCTruthPoint[loc_i]->t);
		
		// only keep pi and K for truth hits
		if(truthHit.M() < 0.13 || truthHit.M() > 0.94 || locDIRCTruthPoint[loc_i]->itrack<1) 
			continue;
		
		japp->RootWriteLock();
		{
			hTruthPoint->Fill(locDIRCTruthPoint[loc_i]->x, locDIRCTruthPoint[loc_i]->y);
			if(loc_i == 0) truthTime = locDIRCTruthPoint[loc_i]->t;
			
			hTruthPointM->Fill(truthHit.M());

			nGenTrk++;
			m_genTrk_px[loc_i] = truthHit.X();
			m_genTrk_py[loc_i] = truthHit.Y();
			m_genTrk_pz[loc_i] = truthHit.Z();
			m_genTrk_E[loc_i] = truthHit.E();
			m_genTrk_x[loc_i] = truthHit_position.X();
			m_genTrk_y[loc_i] = truthHit_position.Y();
			m_genTrk_z[loc_i] = truthHit_position.Z();
			m_genTrk_t[loc_i] = truthHit_position.T();
		}
		japp->RootUnLock();
		
	}

	// if track didn't make it to DIRC plane, don't fill hits
	if(nGenTrk == 0) return NOERROR;
		
	// hits are photons incident on the photocathode
	for(size_t loc_i = 0; loc_i < locDIRCTruthHit.size(); loc_i++) {
		
		japp->RootWriteLock();
		{
			double time = locDIRCTruthHit[loc_i]->t;
			hTruthHitT->Fill(time);
			hTruthHitDeltaT->Fill(time - truthTime);
			hTruthHitE->Fill(locDIRCTruthHit[loc_i]->E*1e9);
			hTruthHitLambda->Fill(1.24e3/(locDIRCTruthHit[loc_i]->E*1e9));
			hTruthHitXY->Fill(locDIRCTruthHit[loc_i]->x, locDIRCTruthHit[loc_i]->y);
			hTruthHitZX->Fill(locDIRCTruthHit[loc_i]->z, locDIRCTruthHit[loc_i]->x);
			
			TVector3 hit(locDIRCTruthHit[loc_i]->x,locDIRCTruthHit[loc_i]->y,locDIRCTruthHit[loc_i]->z);
			TVector3 hit_mom(locDIRCTruthHit[loc_i]->px, locDIRCTruthHit[loc_i]->py, locDIRCTruthHit[loc_i]->pz);
			
			// approx reference for bottom edge of readout plane
			TVector3 ref(239.8,locDIRCTruthHit[loc_i]->y,563.9);
			
			// vector for normal to readout plane
			TVector3 ref_top(219.1,locDIRCTruthHit[loc_i]->y,541.0);
			TVector3 norm = ref_top - ref;
			norm.RotateY(-1.*TMath::Pi()/2.);
			
			// histogram the incident angle where photons where photons get reflected
			double incidentAngle = norm.Angle(hit_mom);
			double criticalAngle = 45. * TMath::Pi()/180.; // what is the real number?
			
			//cout<<endl;
			//cout<<norm.X()<<" "<<norm.Y()<<" "<<norm.Z()<<endl;
			//cout<<hit_mom.X()*1e9<<" "<<hit_mom.Y()*1e9<<" "<<hit_mom.Z()*1e9<<endl;
			//cout<<" "<<incidentAngle*180./TMath::Pi()<<endl;
			
			if(incidentAngle > criticalAngle)
				hTruthHitMissingYLocW->Fill(locDIRCTruthHit[loc_i]->y, (hit-ref).Mag());
			
			// plot position on local readout plane
			hTruthHitYLocW->Fill(locDIRCTruthHit[loc_i]->y, (hit-ref).Mag());
			hTruthHitYLocWT->Fill(locDIRCTruthHit[loc_i]->y, (hit-ref).Mag(), time);
			hTruthHitIncidentAngleY->Fill(locDIRCTruthHit[loc_i]->y, incidentAngle*180./TMath::Pi());
			hTruthHitIncidentAngleLocW->Fill((hit-ref).Mag(), incidentAngle*180./TMath::Pi());
			
			nTruthHit++;
			m_hit_x[loc_i] = hit.X();
			m_hit_y[loc_i] = hit.Y();
			m_hit_z[loc_i] = hit.Z();
			m_hit_t[loc_i] = time;
			
			int pixelBinY = hTruthHitYLocW->GetXaxis()->FindBin(locDIRCTruthHit[loc_i]->y);
			int pixelBinW = hTruthHitYLocW->GetYaxis()->FindBin((hit-ref).Mag());
			m_hitpixel_y[loc_i] = hTruthHitYLocW->GetXaxis()->GetBinCenter(pixelBinY);
			m_hitpixel_w[loc_i] = hTruthHitYLocW->GetYaxis()->GetBinCenter(pixelBinW);
		}
		japp->RootUnLock();
	}
	
	japp->RootWriteLock();
	{
		hNTruthHit->Fill(nTruthHit, nGenTrk);
		dircTree->Fill();
	}
	japp->RootUnLock();
	

	/*
	if(nGenTrk > 1) {
	cout<<"found new event with 2 or more pions"<<endl;

		vector<const DMCThrown*> locMCThrown;
		loop->Get(locMCThrown);
		
		for(size_t loc_i = 0; loc_i < locMCThrown.size(); loc_i++){
			cout<<locMCThrown[loc_i]->PID()<<endl;
		}
	}
	*/

	event++;
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_dirc_truth_hists::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_dirc_truth_hists::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

