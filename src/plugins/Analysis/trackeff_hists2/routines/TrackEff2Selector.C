#define TrackEff2Selector_cxx
// The class definition in TrackEff2Selector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("TrackEff2Selector.C")
// Root > T->Process("TrackEff2Selector.C","some options")
// Root > T->Process("TrackEff2Selector.C+")
//

#include "TrackEff2Selector.h"
#include <TH2.h>
#include <TStyle.h>


void TrackEff2Selector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

	dOutputFile = new TFile("dh_TrackEff2Hists.root", "RECREATE");

	unsigned int locNumPBins = 100;
	unsigned int locNumThetaBins = 140;
	unsigned int locNumDeltaPOverPBins = 1000;
	unsigned int locNumDeltaPtOverPtBins = 1000;
	unsigned int locNumDeltaThetaBins = 480;
	unsigned int locNumDeltaPhiBins = 800;
	unsigned int locNumConLevBins = 500;

	double locMinP = 0.0, locMaxP = 2.0;
	double locMinTheta = 0.0, locMaxTheta = 140.0;
	double locMinConLev = 0.0, locMaxConLev = 1.0;

	double locMinDeltaPOverP = -0.4, locMaxDeltaPOverP = 0.4; //more for wrong charge!
	double locMinDeltaPtOverPt = -0.4, locMaxDeltaPtOverPt = 0.4;
	double locMinDeltaTheta = -12.0, locMaxDeltaTheta = 12.0;
	double locMinDeltaPhi = -40.0, locMaxDeltaPhi = 40.0;

	double locMinDeltaPOverP_WrongCharge = -50.0, locMaxDeltaPOverP_WrongCharge = 50.0;
	double locMinDeltaPtOverPt_WrongCharge = -50.0, locMaxDeltaPtOverPt_WrongCharge = 50.0;
	double locMinDeltaTheta_WrongCharge = -360.0, locMaxDeltaTheta_WrongCharge = 360.0;
	double locMinDeltaPhi_WrongCharge = -360.0, locMaxDeltaPhi_WrongCharge = 360.0;

	//EFFICIENCIES
	dSelectorHist_TrackEff2_CandidateReconstructed = new TH2F("dSelectorHist_TrackEff2_CandidateReconstructed", "Reconstructed Track Candidates;#theta#circ;p (GeV/c)", locNumThetaBins, locMinTheta, locMaxTheta, locNumPBins, locMinP, locMaxP);
	dSelectorHist_TrackEff2_CandidateMissing = new TH2F("dSelectorHist_TrackEff2_CandidateMissing", "Missing Track Candidates;#theta#circ;p (GeV/c)", locNumThetaBins, locMinTheta, locMaxTheta, locNumPBins, locMinP, locMaxP);
	dSelectorHist_TrackEff2_WireBasedReconstructed = new TH2F("dSelectorHist_TrackEff2_WireBasedReconstructed", "Reconstructed Wire Based Tracks;#theta#circ;p (GeV/c)", locNumThetaBins, locMinTheta, locMaxTheta, locNumPBins, locMinP, locMaxP);
	dSelectorHist_TrackEff2_WireBasedMissing = new TH2F("dSelectorHist_TrackEff2_WireBasedMissing", "Missing Wire Based Tracks;#theta#circ;p (GeV/c)", locNumThetaBins, locMinTheta, locMaxTheta, locNumPBins, locMinP, locMaxP);
	dSelectorHist_TrackEff2_TimeBasedReconstructed = new TH2F("dSelectorHist_TrackEff2_TimeBasedReconstructed", "Reconstructed Time Based Tracks;#theta#circ;p (GeV/c)", locNumThetaBins, locMinTheta, locMaxTheta, locNumPBins, locMinP, locMaxP);
	dSelectorHist_TrackEff2_TimeBasedMissing = new TH2F("dSelectorHist_TrackEff2_TimeBasedMissing", "Missing Time Based Tracks;#theta#circ;p (GeV/c)", locNumThetaBins, locMinTheta, locMaxTheta, locNumPBins, locMinP, locMaxP);

	dSelectorHist_TrackEff2_CandidateChargeCorrect = new TH2F("dSelectorHist_TrackEff2_CandidateChargeCorrect", "Correct Charge Track Candidates;#theta#circ;p (GeV/c)", locNumThetaBins, locMinTheta, locMaxTheta, locNumPBins, locMinP, locMaxP);
	dSelectorHist_TrackEff2_CandidateChargeWrong = new TH2F("dSelectorHist_TrackEff2_CandidateChargeWrong", "Incorrect Charge Track Candidates;#theta#circ;p (GeV/c)", locNumThetaBins, locMinTheta, locMaxTheta, locNumPBins, locMinP, locMaxP);
	dSelectorHist_TrackEff2_WireBasedChargeCorrect = new TH2F("dSelectorHist_TrackEff2_WireBasedChargeCorrect", "Correct Charge Wire Based Tracks;#theta#circ;p (GeV/c)", locNumThetaBins, locMinTheta, locMaxTheta, locNumPBins, locMinP, locMaxP);
	dSelectorHist_TrackEff2_WireBasedChargeWrong = new TH2F("dSelectorHist_TrackEff2_WireBasedChargeWrong", "Incorrect Charge Wire Based Tracks;#theta#circ;p (GeV/c)", locNumThetaBins, locMinTheta, locMaxTheta, locNumPBins, locMinP, locMaxP);
	dSelectorHist_TrackEff2_TimeBasedChargeCorrect = new TH2F("dSelectorHist_TrackEff2_TimeBasedChargeCorrect", "Correct Charge Time Based Tracks;#theta#circ;p (GeV/c)", locNumThetaBins, locMinTheta, locMaxTheta, locNumPBins, locMinP, locMaxP);
	dSelectorHist_TrackEff2_TimeBasedChargeWrong = new TH2F("dSelectorHist_TrackEff2_TimeBasedChargeWrong", "Incorrect Charge Time Based Tracks;#theta#circ;p (GeV/c)", locNumThetaBins, locMinTheta, locMaxTheta, locNumPBins, locMinP, locMaxP);

	//CONFIDENCE LEVEL
	dSelectorHist_TrackEff2_WireBasedConLev = new TH1F("dSelectorHist_TrackEff2_WireBasedConLev", "Correct Charge Wire Based Tracks;Wire Based Track Reconstruction Confidence Level", locNumConLevBins, locMinConLev, locMaxConLev);
	dSelectorHist_TrackEff2_TimeBasedConLev = new TH1F("dSelectorHist_TrackEff2_TimeBasedConLev", "Correct Charge Time Based Tracks;Time Based Track Reconstruction Confidence Level", locNumConLevBins, locMinConLev, locMaxConLev);
	dSelectorHist_TrackEff2_WireBasedLowConLev = new TH2F("dSelectorHist_TrackEff2_WireBasedChargeCorrectLowConLev", "Correct Charge Wire Based Tracks, Confidence Level <= 0.01;#theta#circ;p (GeV/c)", locNumThetaBins, locMinTheta, locMaxTheta, locNumPBins, locMinP, locMaxP);
	dSelectorHist_TrackEff2_WireBasedNonLowConLev = new TH2F("dSelectorHist_TrackEff2_WireBasedChargeCorrectNonLowConLev", "Correct Charge Wire Based Tracks, Confidence Level > 0.01;#theta#circ;p (GeV/c)", locNumThetaBins, locMinTheta, locMaxTheta, locNumPBins, locMinP, locMaxP);
	dSelectorHist_TrackEff2_TimeBasedLowConLev = new TH2F("dSelectorHist_TrackEff2_TimeBasedChargeCorrectLowConLev", "Correct Charge Time Based Tracks, Confidence Level <= 0.01;#theta#circ;p (GeV/c)", locNumThetaBins, locMinTheta, locMaxTheta, locNumPBins, locMinP, locMaxP);
	dSelectorHist_TrackEff2_TimeBasedNonLowConLev = new TH2F("dSelectorHist_TrackEff2_TimeBasedChargeCorrectNonLowConLev", "Correct Charge Time Based Tracks, Confidence Level > 0.01;#theta#circ;p (GeV/c)", locNumThetaBins, locMinTheta, locMaxTheta, locNumPBins, locMinP, locMaxP);
	

	//CHARGE RECONSTRUCTION DETAIL
	dSelectorHist_TrackEff2_CandidateChargeCorrect_DeltaPOverPVsTheta = new TH2F("dSelectorHist_TrackEff2_CandidateChargeCorrect_DeltaPOverPVsTheta", "Correct Charge Track Candidates;#theta#circ;#Deltap/p", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPOverPBins, locMinDeltaPOverP, locMaxDeltaPOverP);
	dSelectorHist_TrackEff2_CandidateChargeWrong_DeltaPOverPVsTheta = new TH2F("dSelectorHist_TrackEff2_CandidateChargeWrong_DeltaPOverPVsTheta", "Incorrect Charge Track Candidates;#theta#circ;#Deltap/p", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPOverPBins, locMinDeltaPOverP_WrongCharge, locMaxDeltaPOverP_WrongCharge);
	dSelectorHist_TrackEff2_WireBasedChargeCorrect_DeltaPOverPVsTheta = new TH2F("dSelectorHist_TrackEff2_WireBasedChargeCorrect_DeltaPOverPVsTheta", "Correct Charge Wire Based Tracks;#theta#circ;#Deltap/p", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPOverPBins, locMinDeltaPOverP, locMaxDeltaPOverP);
	dSelectorHist_TrackEff2_WireBasedChargeWrong_DeltaPOverPVsTheta = new TH2F("dSelectorHist_TrackEff2_WireBasedChargeWrong_DeltaPOverPVsTheta", "Incorrect Charge Wire Based Tracks;#theta#circ;#Deltap/p", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPOverPBins, locMinDeltaPOverP_WrongCharge, locMaxDeltaPOverP_WrongCharge);
	dSelectorHist_TrackEff2_TimeBasedChargeCorrect_DeltaPOverPVsTheta = new TH2F("dSelectorHist_TrackEff2_TimeBasedChargeCorrect_DeltaPOverPVsTheta", "Correct Charge Time Based Tracks;#theta#circ;#Deltap/p", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPOverPBins, locMinDeltaPOverP, locMaxDeltaPOverP);
	dSelectorHist_TrackEff2_TimeBasedChargeWrong_DeltaPOverPVsTheta = new TH2F("dSelectorHist_TrackEff2_TimeBasedChargeWrong_DeltaPOverPVsTheta", "Incorrect Charge Time Based Tracks;#theta#circ;#Deltap/p", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPOverPBins, locMinDeltaPOverP_WrongCharge, locMaxDeltaPOverP_WrongCharge);

	dSelectorHist_TrackEff2_CandidateChargeCorrect_DeltaPtOverPtVsTheta = new TH2F("dSelectorHist_TrackEff2_CandidateChargeCorrect_DeltaPtOverPtVsTheta", "Correct Charge Track Candidates;#theta#circ;#Deltap_{T}/p_{T}", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPtOverPtBins, locMinDeltaPtOverPt, locMaxDeltaPtOverPt);
	dSelectorHist_TrackEff2_CandidateChargeWrong_DeltaPtOverPtVsTheta = new TH2F("dSelectorHist_TrackEff2_CandidateChargeWrong_DeltaPtOverPtVsTheta", "Incorrect Charge Track Candidates;#theta#circ;#Deltap_{T}/p_{T}", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPtOverPtBins, locMinDeltaPtOverPt_WrongCharge, locMaxDeltaPtOverPt_WrongCharge);
	dSelectorHist_TrackEff2_WireBasedChargeCorrect_DeltaPtOverPtVsTheta = new TH2F("dSelectorHist_TrackEff2_WireBasedChargeCorrect_DeltaPtOverPtVsTheta", "Correct Charge Wire Based Tracks;#theta#circ;#Deltap_{T}/p_{T}", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPtOverPtBins, locMinDeltaPtOverPt, locMaxDeltaPtOverPt);
	dSelectorHist_TrackEff2_WireBasedChargeWrong_DeltaPtOverPtVsTheta = new TH2F("dSelectorHist_TrackEff2_WireBasedChargeWrong_DeltaPtOverPtVsTheta", "Incorrect Charge Wire Based Tracks;#theta#circ;#Deltap_{T}/p_{T}", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPtOverPtBins, locMinDeltaPtOverPt_WrongCharge, locMaxDeltaPtOverPt_WrongCharge);
	dSelectorHist_TrackEff2_TimeBasedChargeCorrect_DeltaPtOverPtVsTheta = new TH2F("dSelectorHist_TrackEff2_TimeBasedChargeCorrect_DeltaPtOverPtVsTheta", "Correct Charge Time Based Tracks;#theta#circ;#Deltap_{T}/p_{T}", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPtOverPtBins, locMinDeltaPtOverPt, locMaxDeltaPtOverPt);
	dSelectorHist_TrackEff2_TimeBasedChargeWrong_DeltaPtOverPtVsTheta = new TH2F("dSelectorHist_TrackEff2_TimeBasedChargeWrong_DeltaPtOverPtVsTheta", "Incorrect Charge Time Based Tracks;#theta#circ;#Deltap_{T}/p_{T}", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPtOverPtBins, locMinDeltaPtOverPt_WrongCharge, locMaxDeltaPtOverPt_WrongCharge);

	dSelectorHist_TrackEff2_CandidateChargeCorrect_DeltaThetaVsTheta = new TH2F("dSelectorHist_TrackEff2_CandidateChargeCorrect_DeltaThetaVsTheta", "Correct Charge Track Candidates;#theta#circ;#Delta#theta#circ", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaThetaBins, locMinDeltaTheta, locMaxDeltaTheta);
	dSelectorHist_TrackEff2_CandidateChargeWrong_DeltaThetaVsTheta = new TH2F("dSelectorHist_TrackEff2_CandidateChargeWrong_DeltaThetaVsTheta", "Incorrect Charge Track Candidates;#theta#circ;#Delta#theta#circ", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaThetaBins, locMinDeltaTheta_WrongCharge, locMaxDeltaTheta_WrongCharge);
	dSelectorHist_TrackEff2_WireBasedChargeCorrect_DeltaThetaVsTheta = new TH2F("dSelectorHist_TrackEff2_WireBasedChargeCorrect_DeltaThetaVsTheta", "Correct Charge Wire Based Tracks;#theta#circ;#Delta#theta#circ", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaThetaBins, locMinDeltaTheta, locMaxDeltaTheta);
	dSelectorHist_TrackEff2_WireBasedChargeWrong_DeltaThetaVsTheta = new TH2F("dSelectorHist_TrackEff2_WireBasedChargeWrong_DeltaThetaVsTheta", "Incorrect Charge Wire Based Tracks;#theta#circ;#Delta#theta#circ", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaThetaBins, locMinDeltaTheta_WrongCharge, locMaxDeltaTheta_WrongCharge);
	dSelectorHist_TrackEff2_TimeBasedChargeCorrect_DeltaThetaVsTheta = new TH2F("dSelectorHist_TrackEff2_TimeBasedChargeCorrect_DeltaThetaVsTheta", "Correct Charge Time Based Tracks;#theta#circ;#Delta#theta#circ", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaThetaBins, locMinDeltaTheta, locMaxDeltaTheta);
	dSelectorHist_TrackEff2_TimeBasedChargeWrong_DeltaThetaVsTheta = new TH2F("dSelectorHist_TrackEff2_TimeBasedChargeWrong_DeltaThetaVsTheta", "Incorrect Charge Time Based Tracks;#theta#circ;#Delta#theta#circ", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaThetaBins, locMinDeltaTheta_WrongCharge, locMaxDeltaTheta_WrongCharge);

	dSelectorHist_TrackEff2_CandidateChargeCorrect_DeltaPhiVsTheta = new TH2F("dSelectorHist_TrackEff2_CandidateChargeCorrect_DeltaPhiVsTheta", "Correct Charge Track Candidates;#theta#circ;#Delta#phi#circ", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPhiBins, locMinDeltaPhi, locMaxDeltaPhi);
	dSelectorHist_TrackEff2_CandidateChargeWrong_DeltaPhiVsTheta = new TH2F("dSelectorHist_TrackEff2_CandidateChargeWrong_DeltaPhiVsTheta", "Incorrect Charge Track Candidates;#theta#circ;#Delta#phi#circ", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPhiBins, locMinDeltaPhi_WrongCharge, locMaxDeltaPhi_WrongCharge);
	dSelectorHist_TrackEff2_WireBasedChargeCorrect_DeltaPhiVsTheta = new TH2F("dSelectorHist_TrackEff2_WireBasedChargeCorrect_DeltaPhiVsTheta", "Correct Charge Wire Based Tracks;#theta#circ;#Delta#phi#circ", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPhiBins, locMinDeltaPhi, locMaxDeltaPhi);
	dSelectorHist_TrackEff2_WireBasedChargeWrong_DeltaPhiVsTheta = new TH2F("dSelectorHist_TrackEff2_WireBasedChargeWrong_DeltaPhiVsTheta", "Incorrect Charge Wire Based Tracks;#theta#circ;#Delta#phi#circ", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPhiBins, locMinDeltaPhi_WrongCharge, locMaxDeltaPhi_WrongCharge);
	dSelectorHist_TrackEff2_TimeBasedChargeCorrect_DeltaPhiVsTheta = new TH2F("dSelectorHist_TrackEff2_TimeBasedChargeCorrect_DeltaPhiVsTheta", "Correct Charge Time Based Tracks;#theta#circ;#Delta#phi#circ", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPhiBins, locMinDeltaPhi, locMaxDeltaPhi);
	dSelectorHist_TrackEff2_TimeBasedChargeWrong_DeltaPhiVsTheta = new TH2F("dSelectorHist_TrackEff2_TimeBasedChargeWrong_DeltaPhiVsTheta", "Incorrect Charge Time Based Tracks;#theta#circ;#Delta#phi#circ", locNumThetaBins, locMinTheta, locMaxTheta, locNumDeltaPhiBins, locMinDeltaPhi_WrongCharge, locMaxDeltaPhi_WrongCharge);
}

void TrackEff2Selector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t TrackEff2Selector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either TrackEff2Selector::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.


/*
	if(isreconstructable == false)
		return kTRUE; //not interested...
*/

	GetEntry(entry);
	if((entry % 1000) == 0)
		cout << "entry = " << entry << endl;

	TVector3 locP3_Thrown(pthrown_fX, pthrown_fY, pthrown_fZ);
	TVector3 locP3_TimeBased(pfit_fX, pfit_fY, pfit_fZ);
	TVector3 locP3_WireBased(pfit_wire_fX, pfit_wire_fY, pfit_wire_fZ);
	TVector3 locP3_Candidate(pcan_fX, pcan_fY, pcan_fZ);
	double locThrownTheta_Degrees = locP3_Thrown.Theta()*180.0/TMath::Pi();
	double locThrownMomentum = locP3_Thrown.Mag();
	double delta_p_over_p_can = (locP3_Candidate.Mag() - locThrownMomentum)/locThrownMomentum;
	double delta_p_over_p_wire = (locP3_WireBased.Mag() - locThrownMomentum)/locThrownMomentum;
	double delta_p_over_p = (locP3_TimeBased.Mag() - locThrownMomentum)/locThrownMomentum;
	double locFOM_WireBased = TMath::Prob(trk_chisq_wb, trk_Ndof_wb);
	double locFOM_TimeBased = TMath::Prob(trk_chisq, trk_Ndof);
	double locFOM_Cutoff = 0.01;
	double delta_phi_can_deg = delta_phi_can*180.0/(1000.0*TMath::Pi());
	double delta_phi_wire_deg = delta_phi_wire*180.0/(1000.0*TMath::Pi());
	double delta_phi_deg = delta_phi*180.0/(1000.0*TMath::Pi());
	double delta_theta_can_deg = delta_theta_can*180.0/(1000.0*TMath::Pi());
	double delta_theta_wire_deg = delta_theta_wire*180.0/(1000.0*TMath::Pi());
	double delta_theta_deg = delta_theta*180.0/(1000.0*TMath::Pi());

	//TRACK RECONSTRUCTION EFFICIENCIES
	if(dTrackReconstructedFlag_Candidate == true)
		dSelectorHist_TrackEff2_CandidateReconstructed->Fill(locThrownTheta_Degrees, locThrownMomentum);
	else
		dSelectorHist_TrackEff2_CandidateMissing->Fill(locThrownTheta_Degrees, locThrownMomentum);
	if(dTrackReconstructedFlag_Candidate == true){
		if(dTrackReconstructedFlag_WireBased == true)
			dSelectorHist_TrackEff2_WireBasedReconstructed->Fill(locThrownTheta_Degrees, locThrownMomentum);
		else
			dSelectorHist_TrackEff2_WireBasedMissing->Fill(locThrownTheta_Degrees, locThrownMomentum);
	}
	if(dTrackReconstructedFlag_WireBased == true){
		if(dTrackReconstructedFlag_TimeBased == true)
			dSelectorHist_TrackEff2_TimeBasedReconstructed->Fill(locThrownTheta_Degrees, locThrownMomentum);
		else
			dSelectorHist_TrackEff2_TimeBasedMissing->Fill(locThrownTheta_Degrees, locThrownMomentum);
	}

	//CHARGE RECONSTRUCTION EFFICIENCIES
	if(dTrackReconstructedFlag_Candidate == true){
		if(fabs(q_thrown - q_candidate) < 0.1)
			dSelectorHist_TrackEff2_CandidateChargeCorrect->Fill(locThrownTheta_Degrees, locThrownMomentum);
		else
			dSelectorHist_TrackEff2_CandidateChargeWrong->Fill(locThrownTheta_Degrees, locThrownMomentum);
	}
	if(dTrackReconstructedFlag_WireBased == true){
		if(fabs(q_thrown - q_wirebased) < 0.1)
			dSelectorHist_TrackEff2_WireBasedChargeCorrect->Fill(locThrownTheta_Degrees, locThrownMomentum);
		else
			dSelectorHist_TrackEff2_WireBasedChargeWrong->Fill(locThrownTheta_Degrees, locThrownMomentum);
	}
	if(dTrackReconstructedFlag_TimeBased == true){
		if(fabs(q_thrown - q_timebased) < 0.1)
			dSelectorHist_TrackEff2_TimeBasedChargeCorrect->Fill(locThrownTheta_Degrees, locThrownMomentum);
		else
			dSelectorHist_TrackEff2_TimeBasedChargeWrong->Fill(locThrownTheta_Degrees, locThrownMomentum);
	}

	//CONFIDENCE LEVEL
	if((dTrackReconstructedFlag_WireBased == true) && (fabs(q_thrown - q_wirebased) < 0.1)){
		dSelectorHist_TrackEff2_WireBasedConLev->Fill(locFOM_WireBased);
		if(locFOM_WireBased <= locFOM_Cutoff)
			dSelectorHist_TrackEff2_WireBasedLowConLev->Fill(locThrownTheta_Degrees, locThrownMomentum);
		else
			dSelectorHist_TrackEff2_WireBasedNonLowConLev->Fill(locThrownTheta_Degrees, locThrownMomentum);
	}
	if((dTrackReconstructedFlag_TimeBased == true) && (fabs(q_thrown - q_timebased) < 0.1)){
		dSelectorHist_TrackEff2_TimeBasedConLev->Fill(locFOM_TimeBased);
		if(locFOM_TimeBased <= locFOM_Cutoff)
			dSelectorHist_TrackEff2_TimeBasedLowConLev->Fill(locThrownTheta_Degrees, locThrownMomentum);
		else
			dSelectorHist_TrackEff2_TimeBasedNonLowConLev->Fill(locThrownTheta_Degrees, locThrownMomentum);
	}

	//CHARGE MIS-RECONSTRUCTION DELTAS
	if(dTrackReconstructedFlag_Candidate == true){
		if(fabs(q_thrown - q_candidate) < 0.1){
			dSelectorHist_TrackEff2_CandidateChargeCorrect_DeltaPOverPVsTheta->Fill(locThrownTheta_Degrees, delta_p_over_p_can);
			dSelectorHist_TrackEff2_CandidateChargeCorrect_DeltaPtOverPtVsTheta->Fill(locThrownTheta_Degrees, delta_pt_over_pt_can);
			dSelectorHist_TrackEff2_CandidateChargeCorrect_DeltaThetaVsTheta->Fill(locThrownTheta_Degrees, delta_theta_can_deg);
			dSelectorHist_TrackEff2_CandidateChargeCorrect_DeltaPhiVsTheta->Fill(locThrownTheta_Degrees, delta_phi_can_deg);
		}else{
			dSelectorHist_TrackEff2_CandidateChargeWrong_DeltaPOverPVsTheta->Fill(locThrownTheta_Degrees, delta_p_over_p_can);
			dSelectorHist_TrackEff2_CandidateChargeWrong_DeltaPtOverPtVsTheta->Fill(locThrownTheta_Degrees, delta_pt_over_pt_can);
			dSelectorHist_TrackEff2_CandidateChargeWrong_DeltaThetaVsTheta->Fill(locThrownTheta_Degrees, delta_theta_can_deg);
			dSelectorHist_TrackEff2_CandidateChargeWrong_DeltaPhiVsTheta->Fill(locThrownTheta_Degrees, delta_phi_can_deg);
		}
	}
	if(dTrackReconstructedFlag_WireBased == true){
		if(fabs(q_thrown - q_wirebased) < 0.1){
			dSelectorHist_TrackEff2_WireBasedChargeCorrect_DeltaPOverPVsTheta->Fill(locThrownTheta_Degrees, delta_p_over_p_wire);
			dSelectorHist_TrackEff2_WireBasedChargeCorrect_DeltaPtOverPtVsTheta->Fill(locThrownTheta_Degrees, delta_pt_over_pt_wire);
			dSelectorHist_TrackEff2_WireBasedChargeCorrect_DeltaThetaVsTheta->Fill(locThrownTheta_Degrees, delta_theta_wire_deg);
			dSelectorHist_TrackEff2_WireBasedChargeCorrect_DeltaPhiVsTheta->Fill(locThrownTheta_Degrees, delta_phi_wire_deg);
		}else{
			dSelectorHist_TrackEff2_WireBasedChargeWrong_DeltaPOverPVsTheta->Fill(locThrownTheta_Degrees, delta_p_over_p_wire);
			dSelectorHist_TrackEff2_WireBasedChargeWrong_DeltaPtOverPtVsTheta->Fill(locThrownTheta_Degrees, delta_pt_over_pt_wire);
			dSelectorHist_TrackEff2_WireBasedChargeWrong_DeltaThetaVsTheta->Fill(locThrownTheta_Degrees, delta_theta_wire_deg);
			dSelectorHist_TrackEff2_WireBasedChargeWrong_DeltaPhiVsTheta->Fill(locThrownTheta_Degrees, delta_phi_wire_deg);
		}
	}
	if(dTrackReconstructedFlag_TimeBased == true){
		if(fabs(q_thrown - q_timebased) < 0.1){
			dSelectorHist_TrackEff2_TimeBasedChargeCorrect_DeltaPOverPVsTheta->Fill(locThrownTheta_Degrees, delta_p_over_p);
			dSelectorHist_TrackEff2_TimeBasedChargeCorrect_DeltaPtOverPtVsTheta->Fill(locThrownTheta_Degrees, delta_pt_over_pt);
			dSelectorHist_TrackEff2_TimeBasedChargeCorrect_DeltaThetaVsTheta->Fill(locThrownTheta_Degrees, delta_theta_deg);
			dSelectorHist_TrackEff2_TimeBasedChargeCorrect_DeltaPhiVsTheta->Fill(locThrownTheta_Degrees, delta_phi_deg);
		}else{
			dSelectorHist_TrackEff2_TimeBasedChargeWrong_DeltaPOverPVsTheta->Fill(locThrownTheta_Degrees, delta_p_over_p);
			dSelectorHist_TrackEff2_TimeBasedChargeWrong_DeltaPtOverPtVsTheta->Fill(locThrownTheta_Degrees, delta_pt_over_pt);
			dSelectorHist_TrackEff2_TimeBasedChargeWrong_DeltaThetaVsTheta->Fill(locThrownTheta_Degrees, delta_theta_deg);
			dSelectorHist_TrackEff2_TimeBasedChargeWrong_DeltaPhiVsTheta->Fill(locThrownTheta_Degrees, delta_phi_deg);
		}
	}


/*
	//hist: deviation of measured values from true values at each level (in )

		bool dTrackReconstructedFlag_TimeBased;
		TVector3 pfit;
		double trk_chisq;
		int trk_Ndof;
		double delta_pt_over_pt;
		double delta_theta;	// mrad
		double delta_phi;		// mrad
		int num_timebased;


   Float_t         q_thrown;
   Int_t           PID_thrown;
   UInt_t          pthrown_fUniqueID;
   UInt_t          pthrown_fBits;
   Bool_t          isreconstructable;
   Int_t           Nstereo;
   Int_t           Ncdc;
   Int_t           Nfdc;
   Int_t           NLR_bad_stereo;
   Int_t           NLR_bad;
   Int_t           PID_hypothesized;
   Float_t         q_hypothesized;
   Float_t         FOM_hypothesized;
   Float_t         q_timebased;
   Int_t           PID_timebased;
   Bool_t          dTrackReconstructedFlag_TimeBased;
   Double_t        trk_chisq;
   Int_t           trk_Ndof;
   Double_t        delta_pt_over_pt;
   Double_t        delta_theta;
   Double_t        delta_phi;
   Int_t           num_timebased;
   Float_t         q_wirebased;
   Int_t           PID_wirebased;
   Bool_t          dTrackReconstructedFlag_WireBased;
   Double_t        trk_chisq_wb;
   Int_t           trk_Ndof_wb;
   Double_t        delta_pt_over_pt_wire;
   Double_t        delta_theta_wire;
   Double_t        delta_phi_wire;
   Int_t           num_wirebased;
   Float_t         q_candidate;
   Int_t           PID_candidate;
   Bool_t          dTrackReconstructedFlag_Candidate;
   Double_t        delta_pt_over_pt_can;
   Double_t        delta_theta_can;
   Double_t        delta_phi_can;
   Int_t           num_candidates;
*/
   return kTRUE;
}

void TrackEff2Selector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void TrackEff2Selector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

	dOutputFile->Write();
}
