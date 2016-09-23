#include "JEventProcessor_highlevel_online.h"
using namespace jana;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_highlevel_online());
  }
} // "C"

//------------------
// init
//------------------
jerror_t JEventProcessor_highlevel_online::init(void)
{
	// All histograms go in the "highlevel" directory
	TDirectory *main = gDirectory;

	gDirectory->mkdir("highlevel")->cd();

	/***************************************************************** RF *****************************************************************/

	int NRFbins = 1024;
	double ns_per_bin = 0.2;
	double RFlo = 0.0;
	double RFhi = RFlo + ns_per_bin*(double)NRFbins;
	double dft_min =  50.0; // min frequency in Mhz
	double dft_max = 900.0; // max frequency in Mhz
	int Ndft_bins = (1.0/dft_min - 1.0/dft_max)*1000.0/ns_per_bin;
	dHist_BeamBunchPeriod     = new TH1I("RFBeamBunchPeriod", "RF Beam Bunch Period;TAGH #Deltat (Within Same Counter) (ns)", NRFbins, RFlo, RFhi);
	dHist_BeamBunchPeriod_DFT = new TH1F("RFBeamBunchPeriod_DFT", "Fourier Transform of RF Beam Bunch Period;Beam bunch frequency (MHz)", Ndft_bins, dft_min, dft_max);

	/*************************************************************** TRIGGER **************************************************************/

	dHist_L1GTPRate = new TH2F("L1GTPRate","L1 GTP Rate by bit;Trigger Bit;L1 GTP Rate (kHz)", 8, 0.5, 8.5, 1000, 0.0, 100.0);

	dHist_BCALVsFCAL_TrigBit1 = new TH2I("BCALVsFCAL_TrigBit1","TRIG BIT 1;E (FCAL) (count);E (BCAL) (count)", 200, 0., 10000, 200, 0., 50000);
	dHist_BCALVsFCAL_TrigBit6 = new TH2I("BCALVsFCAL_TrigBit6","TRIG BIT 6;E (FCAL) (count);E (BCAL) (count)", 200, 0., 10000, 200, 0., 50000);

	/****************************************************** NUM RECONSTRUCTED OBJECTS *****************************************************/

	//2D Summary
	dHist_NumHighLevelObjects = new TH2I("NumHighLevelObjects", ";;# Objects / Event", 12, 0.5, 12.5, 101, -0.5, 100.5);
	dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(1, "DSCHit");
	dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(2, "DTOFPoint");
	dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(3, "DBCALShower");
	dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(4, "DFCALShower");
	dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(5, "DTimeBasedTrack");
	dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(6, "TrackSCMatches");
	dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(7, "TrackTOFMatches");
	dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(8, "TrackBCALMatches");
	dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(9, "TrackFCALMatches");
	dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(10, "DBeamPhoton");
	dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(11, "DChargedTrack");
	dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(12, "DNeutralShower");

	/************************************************************* KINEMATICS *************************************************************/

	// Beam Energy
	dHist_BeamEnergy = new TH1I("BeamEnergy", "Reconstructed Tagger Beam Energy;Beam Energy (GeV)", 240, 0.0, 12.0);

	// PVsTheta Time-Based Tracks
	dHist_PVsTheta_Tracks = new TH2I("PVsTheta_Tracks", "P vs. #theta for time-based tracks;#theta#circ;p (GeV/c)", 280, 0.0, 140.0, 150, 0.0, 12.0);

	// PhiVsTheta Time-Based Tracks
	dHist_PhiVsTheta_Tracks = new TH2I("PhiVsTheta_Tracks", "#phi vs. #theta for time-based tracks;#theta#circ;#phi#circ", 280, 0.0, 140.0, 60, -180.0, 180.0);

	/*************************************************************** VERTEX ***************************************************************/

	// Event Vertex-Z
	dEventVertexZ = new TH1I("EventVertexZ", "Reconstructed Event Vertex Z;Event Vertex-Z (cm)", 600, 0.0, 200.0);

	// Event Vertex-Y Vs Vertex-X
	dEventVertexYVsX = new TH2I("EventVertexYVsX", "Reconstructed Event Vertex X/Y;Event Vertex-X (cm);Event Vertex-Y (cm)", 400, -10.0, 10.0, 400, -10.0, 10.0);

	/*************************************************************** 2 gamma inv. mass ***************************************************************/
	// 2-gamma inv. mass
	d2gamma = new TH1I("TwoGammaMass", "2#gamma inv. mass;2#gamma inv. mass (GeV)", 400, 0.0, 1.2);

	// pi+ pi-
	dpip_pim = new TH1I("PiPlusPiMinus", "#pi^{+}#pi^{-} inv. mass w/ identified proton;#pi^{+}#pi^{-} inv. mass (GeV)", 400, 0.0, 2.0);

	// pi+ pi- pi0
	dpip_pim_pi0 = new TH1I("PiPlusPiMinusPiZero", "#pi^{+}#pi^{-}#pi^{o} inv. mass w/ identified proton;#pi^{+}#pi^{-}#pi^{o} inv. mass (GeV)", 200, 0.035, 2.0);
	dptrans = new TH1I("PiPlusPiMinusPiZeroProton_t", ";#pi^{+}#pi^{-}#pi^{o}p transverse momentum(GeV)", 500, 0.0, 1.0);
	
	dbeta_vs_p = new TH2I("BetaVsP", "#beta vs. p (best FOM all charged tracks);p (GeV);#beta", 200, 0.0, 2.0, 100, 0.0, 1.1);

	// back to main dir
	main->cd();
  
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_highlevel_online::brun(JEventLoop *locEventLoop, int32_t runnumber)
{
  // This is called whenever the run number changes

	fcal_cell_thr  =  64;
	bcal_cell_thr  =  20;

	fcal_row_mask_min = 26;
	fcal_row_mask_max = 32;
	fcal_col_mask_min = 26;
	fcal_col_mask_max = 32;

	if( runnumber < 11127 )
	{
		fcal_row_mask_min = 24;
		fcal_row_mask_max = 34;
		fcal_col_mask_min = 24;
		fcal_col_mask_max = 34;
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_highlevel_online::evnt(JEventLoop *locEventLoop, uint64_t eventnumber)
{
	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);

	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locFCALShowers);

	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers);

	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles, "PreSelect");

	vector<const DTOFPoint*> locTOFPoints;
	locEventLoop->Get(locTOFPoints);

	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locSCHits);

	vector<const DTAGHHit*> locTAGHHits;
	locEventLoop->Get(locTAGHHits);

    vector<const DBCALDigiHit*> locBCALDigiHits;
    locEventLoop->Get(locBCALDigiHits);

    vector<const DFCALDigiHit*> locFCALDigiHits;
    locEventLoop->Get(locFCALDigiHits);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	const DVertex* locVertex = NULL;
	locEventLoop->GetSingle(locVertex);

	vector<const DL1Trigger*> locL1Triggers;
	locEventLoop->Get(locL1Triggers);
	const DL1Trigger* locL1Trigger = locL1Triggers.empty() ? NULL : locL1Triggers[0];

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, "PreSelect");

	/********************************************************* PREPARE RF ********************************************************/
	
	// Do discrete Fourier transform for selected frequency range
	double *rf_dft = NULL;
	if(dHist_BeamBunchPeriod->GetEntries() > 100){
		
		int NRF_bins      = dHist_BeamBunchPeriod->GetNbinsX();
		int Ndft_bins     = dHist_BeamBunchPeriod_DFT->GetNbinsX();
		double ns_per_bin = dHist_BeamBunchPeriod->GetBinWidth(1);
		double rfdomain   = ns_per_bin*(double)NRF_bins;
		rf_dft = new double[Ndft_bins];
		for(int ibin=1; ibin<=Ndft_bins; ibin++){

			// n.b. the factor of 1000 below is to convert from MHz to GHz
			double f = dHist_BeamBunchPeriod_DFT->GetXaxis()->GetBinCenter(ibin);
			double k = rfdomain/(1000.0/f); // k is number of oscillations in whole x-range for this freq.

			double dphi = TMath::TwoPi()*k/(double)NRF_bins;

			double sumc = 0.0;
			double sums = 0.0;
			for(int jbin=1; jbin<=(int)NRF_bins; jbin++){

				double v = dHist_BeamBunchPeriod->GetBinContent(jbin);
				double theta = (double)jbin*dphi;
				sumc += v*cos(theta);
				sums += v*sin(theta);
			}

			rf_dft[ibin-1] = sqrt(sumc*sumc + sums*sums);
		}
	}

	/********************************************************* PREPARE TAGGER HITS ********************************************************/

	//organize TAGH hits
	map<int, set<double> > locTAGHHitMap; //double: TAGH tdc times (set allows time sorting)
	for(size_t loc_i = 0; loc_i < locTAGHHits.size(); ++loc_i)
	{
		if(locTAGHHits[loc_i]->has_TDC)
			locTAGHHitMap[locTAGHHits[loc_i]->counter_id].insert(locTAGHHits[loc_i]->time_tdc);
	}

	//collect TAGH delta-t's for RF beam bunch histogram
	map<int, set<double> >::iterator locTAGHIterator = locTAGHHitMap.begin();
	vector<double> locTAGHDeltaTs;
	for(; locTAGHIterator != locTAGHHitMap.end(); ++locTAGHIterator)
	{
		set<double>& locCounterHits = locTAGHIterator->second;
		if(locCounterHits.size() < 2)
			continue;
		set<double>::iterator locSetIterator = locCounterHits.begin();
		set<double>::iterator locPreviousSetIterator = locSetIterator;
		for(++locSetIterator; locSetIterator != locCounterHits.end(); ++locSetIterator, ++locPreviousSetIterator)
			locTAGHDeltaTs.push_back(*locSetIterator - *locPreviousSetIterator);
	}

	/********************************************************* PREPARE TRIGGER INFO *******************************************************/

	vector<unsigned int> locTrigBits(32, 0); //bit 1 -> 32 is index 0 -> 31

	if(locL1Trigger != NULL)
	{
		for(unsigned int bit = 0; bit < 32; ++bit)
			locTrigBits[bit] = (locL1Trigger->trig_mask & (1 << bit)) ? 1 : 0;
	}

	//Get total FCAL energy
	int fcal_tot_en = 0;
    for(unsigned int jj = 0; jj < locFCALDigiHits.size(); jj++)
    {
		const DFCALDigiHit *fcal_hit = locFCALDigiHits[jj];

		int row = fcal_hit->row;
		int col = fcal_hit->column;
		if((row >= fcal_row_mask_min && row <= fcal_row_mask_max) && (col >= fcal_col_mask_min && col <= fcal_col_mask_max))
			continue;

		uint32_t adc_time = (fcal_hit->pulse_time & 0x7FC0) >> 6;
		if(!((adc_time > 20) && (adc_time < 70)))
			continue;

		Int_t pulse_int = fcal_hit->pulse_integral - fcal_hit->nsamples_integral*100;
		if(pulse_int < 0)
			continue;

		const Df250PulsePedestal *pulsepedestal = NULL;
		fcal_hit->GetSingle(pulsepedestal);


		Int_t pulse_peak = -10;
		if(pulsepedestal)
			pulse_peak = pulsepedestal->pulse_peak - 100;

		if(pulse_peak <= fcal_cell_thr)
			continue;

		fcal_tot_en += pulse_int;
	}

	//Get total BCAL energy
	int bcal_tot_en = 0;
	for(unsigned int jj = 0; jj < locBCALDigiHits.size(); ++jj)
	{
		const DBCALDigiHit *bcal_hit = locBCALDigiHits[jj];
		Int_t pulse_int = bcal_hit->pulse_integral - bcal_hit->nsamples_integral*100;

		const Df250PulsePedestal *pulsepedestal;
		bcal_hit->GetSingle(pulsepedestal);

		Int_t pulse_peak = -10;
		if(pulsepedestal)
			pulse_peak = pulsepedestal->pulse_peak - 100;
		if(pulse_peak <= bcal_cell_thr)
			continue;

		bcal_tot_en += pulse_int;
	}
	

    /***************************************************************** RF *****************************************************************/

	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

	for(size_t loc_i = 0; loc_i < locTAGHDeltaTs.size(); ++loc_i)
		dHist_BeamBunchPeriod->Fill(locTAGHDeltaTs[loc_i]);
	
	if(rf_dft){
		for(int i=0; i<dHist_BeamBunchPeriod_DFT->GetNbinsX(); i++){
			dHist_BeamBunchPeriod_DFT->SetBinContent(i+1, rf_dft[i]);
		}
		delete[] rf_dft;
	}

	/*************************************************************** TRIGGER **************************************************************/

	if(locL1Trigger != NULL)
	{
		if(locTrigBits[0] == 1) //bit 1
			dHist_BCALVsFCAL_TrigBit1->Fill(Float_t(fcal_tot_en), Float_t(bcal_tot_en));
		if(locTrigBits[5] == 1) //bit 6
			dHist_BCALVsFCAL_TrigBit6->Fill(Float_t(fcal_tot_en), Float_t(bcal_tot_en));

		// Sync Events
		if(!locL1Trigger->gtp_rate.empty())
		{
			for(unsigned int ii = 0; ii < 8; ++ii)
				dHist_L1GTPRate->Fill(ii + 1, Float_t(locL1Trigger->gtp_rate[ii])/1000.0);
		}
	}

	/****************************************************** NUM RECONSTRUCTED OBJECTS *****************************************************/

	//# High-Level Objects
	dHist_NumHighLevelObjects->Fill(1, (Double_t)locSCHits.size());
	dHist_NumHighLevelObjects->Fill(2, (Double_t)locTOFPoints.size());
	dHist_NumHighLevelObjects->Fill(3, (Double_t)locBCALShowers.size());
	dHist_NumHighLevelObjects->Fill(4, (Double_t)locFCALShowers.size());
	dHist_NumHighLevelObjects->Fill(5, (Double_t)locTrackTimeBasedVector.size());
	dHist_NumHighLevelObjects->Fill(6, (Double_t)locDetectorMatches->Get_NumTrackSCMatches());
	dHist_NumHighLevelObjects->Fill(7, (Double_t)locDetectorMatches->Get_NumTrackTOFMatches());
	dHist_NumHighLevelObjects->Fill(8, (Double_t)locDetectorMatches->Get_NumTrackBCALMatches());
	dHist_NumHighLevelObjects->Fill(9, (Double_t)locDetectorMatches->Get_NumTrackFCALMatches());
	dHist_NumHighLevelObjects->Fill(10, (Double_t)locBeamPhotons.size());
	dHist_NumHighLevelObjects->Fill(11, (Double_t)locChargedTracks.size());
	dHist_NumHighLevelObjects->Fill(12, (Double_t)locNeutralShowers.size());

	/************************************************************* KINEMATICS *************************************************************/

	for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
		dHist_BeamEnergy->Fill(locBeamPhotons[loc_i]->energy());

	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		auto locChargedHypo = locChargedTracks[loc_i]->Get_BestTrackingFOM();
		double locP = locChargedHypo->momentum().Mag();
		double locTheta = locChargedHypo->momentum().Theta()*180.0/TMath::Pi();
		double locPhi = locChargedHypo->momentum().Phi()*180.0/TMath::Pi();
		double locBeta = locChargedHypo->measuredBeta();
		dHist_PVsTheta_Tracks->Fill(locTheta, locP);
		dHist_PhiVsTheta_Tracks->Fill(locTheta, locPhi);
		dbeta_vs_p->Fill(locP, locBeta);
	}

	/*************************************************************** VERTEX ***************************************************************/

	if(locChargedTracks.size() >= 2)
	{
		dEventVertexZ->Fill(locVertex->dSpacetimeVertex.Z());
		dEventVertexYVsX->Fill(locVertex->dSpacetimeVertex.X(), locVertex->dSpacetimeVertex.Y());
	}

	/*************************************************************** 2 gamma inv. mass ***************************************************************/
	vector<DLorentzVector> pi0s;
	for(uint32_t i=0; i<locNeutralParticles.size(); i++){

		auto hypoth1 = locNeutralParticles[i]->Get_Hypothesis(Gamma);
		if(!hypoth1) continue;
		const DLorentzVector &v1 = hypoth1->lorentzMomentum();

		for(uint32_t j=i+1; j<locNeutralParticles.size(); j++){

			auto hypoth2 = locNeutralParticles[j]->Get_Hypothesis(Gamma);
			if(!hypoth2) continue;
			const DLorentzVector &v2 = hypoth2->lorentzMomentum();

			DLorentzVector ptot(v1+v2);
			double mass = ptot.M();
			d2gamma->Fill(mass);
			if(mass>0.075 && mass<0.2) pi0s.push_back(ptot);
		}
	}

	/*************************************************************** pi+ pi- pi0 ***************************************************************/
	for(auto t1 : locChargedTracks){
		auto hypoth1 = t1->Get_BestTrackingFOM();
		if(!hypoth1) continue;
		if(hypoth1->PID() != PiPlus) continue;
		const DLorentzVector &pipmom = hypoth1->lorentzMomentum();

		for(auto t2 : locChargedTracks){
			if(t2 == t1) continue;
			auto hypoth2 = t2->Get_BestTrackingFOM();
			if(!hypoth2) continue;
			if(hypoth2->PID() != PiMinus) continue;
			const DLorentzVector &pimmom = hypoth2->lorentzMomentum();

			for(auto t3 : locChargedTracks){
				if(t3 == t1) continue;
				if(t3 == t2) continue;
				auto hypoth3 = t3->Get_BestTrackingFOM();
				if(!hypoth3) continue;
				if(hypoth3->PID() != Proton) continue;
				const DLorentzVector &protonmom = hypoth3->lorentzMomentum();

				DLorentzVector rhomom(pipmom + pimmom);
				dpip_pim->Fill(rhomom.M());

				for(uint32_t i=0; i<pi0s.size(); i++){
					DLorentzVector &pi0mom = pi0s[i];
	
					DLorentzVector omegamom(pi0mom + rhomom);
					DLorentzVector totmom(omegamom + protonmom);
					
					double ptrans = totmom.Perp();
					dptrans->Fill(ptrans);
					if(ptrans > 0.1) continue;
				
					dpip_pim_pi0->Fill(omegamom.M());
				}
			}
		}		
	}


	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_highlevel_online::erun(void)
{
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_highlevel_online::fini(void)
{
  // Called before program exit after event processing is finished.
  return NOERROR;
}

