//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec 20 14:16:04 2011 by ROOT version 5.30/00
// from TTree trkeff2/Tracking Efficiency
// found on file: cluster_run/hd_root_trackstudy_proton.root
//////////////////////////////////////////////////////////

#ifndef TrackEff2Selector_h
#define TrackEff2Selector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

class TrackEff2Selector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
 //track2          *E;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   ULong_t         event;
   Float_t         q_thrown;
   Int_t           PID_thrown;
   UInt_t          pthrown_fUniqueID;
   UInt_t          pthrown_fBits;
   Double_t        pthrown_fX;
   Double_t        pthrown_fY;
   Double_t        pthrown_fZ;
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
   UInt_t          pfit_fUniqueID;
   UInt_t          pfit_fBits;
   Double_t        pfit_fX;
   Double_t        pfit_fY;
   Double_t        pfit_fZ;
   Double_t        trk_chisq;
   Int_t           trk_Ndof;
   Double_t        delta_pt_over_pt;
   Double_t        delta_theta;
   Double_t        delta_phi;
   Int_t           num_timebased;
   Float_t         q_wirebased;
   Int_t           PID_wirebased;
   Bool_t          dTrackReconstructedFlag_WireBased;
   UInt_t          pfit_wire_fUniqueID;
   UInt_t          pfit_wire_fBits;
   Double_t        pfit_wire_fX;
   Double_t        pfit_wire_fY;
   Double_t        pfit_wire_fZ;
   Double_t        trk_chisq_wb;
   Int_t           trk_Ndof_wb;
   Double_t        delta_pt_over_pt_wire;
   Double_t        delta_theta_wire;
   Double_t        delta_phi_wire;
   Int_t           num_wirebased;
   Float_t         q_candidate;
   Int_t           PID_candidate;
   Bool_t          dTrackReconstructedFlag_Candidate;
   UInt_t          pcan_fUniqueID;
   UInt_t          pcan_fBits;
   Double_t        pcan_fX;
   Double_t        pcan_fY;
   Double_t        pcan_fZ;
   Double_t        delta_pt_over_pt_can;
   Double_t        delta_theta_can;
   Double_t        delta_phi_can;
   Int_t           num_candidates;

   // List of branches
   TBranch        *b_E_fUniqueID;   //!
   TBranch        *b_E_fBits;   //!
   TBranch        *b_E_event;   //!
   TBranch        *b_E_q_thrown;   //!
   TBranch        *b_E_PID_thrown;   //!
   TBranch        *b_E_pthrown_fUniqueID;   //!
   TBranch        *b_E_pthrown_fBits;   //!
   TBranch        *b_E_pthrown_fX;   //!
   TBranch        *b_E_pthrown_fY;   //!
   TBranch        *b_E_pthrown_fZ;   //!
   TBranch        *b_E_isreconstructable;   //!
   TBranch        *b_E_Nstereo;   //!
   TBranch        *b_E_Ncdc;   //!
   TBranch        *b_E_Nfdc;   //!
   TBranch        *b_E_NLR_bad_stereo;   //!
   TBranch        *b_E_NLR_bad;   //!
   TBranch        *b_E_PID_hypothesized;   //!
   TBranch        *b_E_q_hypothesized;   //!
   TBranch        *b_E_FOM_hypothesized;   //!
   TBranch        *b_E_q_timebased;   //!
   TBranch        *b_E_PID_timebased;   //!
   TBranch        *b_E_dTrackReconstructedFlag_TimeBased;   //!
   TBranch        *b_E_pfit_fUniqueID;   //!
   TBranch        *b_E_pfit_fBits;   //!
   TBranch        *b_E_pfit_fX;   //!
   TBranch        *b_E_pfit_fY;   //!
   TBranch        *b_E_pfit_fZ;   //!
   TBranch        *b_E_trk_chisq;   //!
   TBranch        *b_E_trk_Ndof;   //!
   TBranch        *b_E_delta_pt_over_pt;   //!
   TBranch        *b_E_delta_theta;   //!
   TBranch        *b_E_delta_phi;   //!
   TBranch        *b_E_num_timebased;   //!
   TBranch        *b_E_q_wirebased;   //!
   TBranch        *b_E_PID_wirebased;   //!
   TBranch        *b_E_dTrackReconstructedFlag_WireBased;   //!
   TBranch        *b_E_pfit_wire_fUniqueID;   //!
   TBranch        *b_E_pfit_wire_fBits;   //!
   TBranch        *b_E_pfit_wire_fX;   //!
   TBranch        *b_E_pfit_wire_fY;   //!
   TBranch        *b_E_pfit_wire_fZ;   //!
   TBranch        *b_E_trk_chisq_wb;   //!
   TBranch        *b_E_trk_Ndof_wb;   //!
   TBranch        *b_E_delta_pt_over_pt_wire;   //!
   TBranch        *b_E_delta_theta_wire;   //!
   TBranch        *b_E_delta_phi_wire;   //!
   TBranch        *b_E_num_wirebased;   //!
   TBranch        *b_E_q_candidate;   //!
   TBranch        *b_E_PID_candidate;   //!
   TBranch        *b_E_dTrackReconstructedFlag_Candidate;   //!
   TBranch        *b_E_pcan_fUniqueID;   //!
   TBranch        *b_E_pcan_fBits;   //!
   TBranch        *b_E_pcan_fX;   //!
   TBranch        *b_E_pcan_fY;   //!
   TBranch        *b_E_pcan_fZ;   //!
   TBranch        *b_E_delta_pt_over_pt_can;   //!
   TBranch        *b_E_delta_theta_can;   //!
   TBranch        *b_E_delta_phi_can;   //!
   TBranch        *b_E_num_candidates;   //!

	TFile* dOutputFile;

	TH2F *dSelectorHist_TrackEff2_CandidateReconstructed;
	TH2F *dSelectorHist_TrackEff2_CandidateMissing;
	TH2F *dSelectorHist_TrackEff2_WireBasedReconstructed;
	TH2F *dSelectorHist_TrackEff2_WireBasedMissing;
	TH2F *dSelectorHist_TrackEff2_TimeBasedReconstructed;
	TH2F *dSelectorHist_TrackEff2_TimeBasedMissing;

	TH2F *dSelectorHist_TrackEff2_CandidateChargeCorrect;
	TH2F *dSelectorHist_TrackEff2_CandidateChargeWrong;
	TH2F *dSelectorHist_TrackEff2_WireBasedChargeCorrect;
	TH2F *dSelectorHist_TrackEff2_WireBasedChargeWrong;
	TH2F *dSelectorHist_TrackEff2_TimeBasedChargeCorrect;
	TH2F *dSelectorHist_TrackEff2_TimeBasedChargeWrong;

	TH1F *dSelectorHist_TrackEff2_WireBasedConLev;
	TH1F *dSelectorHist_TrackEff2_TimeBasedConLev;
	TH2F *dSelectorHist_TrackEff2_WireBasedLowConLev;
	TH2F *dSelectorHist_TrackEff2_WireBasedNonLowConLev;
	TH2F *dSelectorHist_TrackEff2_TimeBasedLowConLev;
	TH2F *dSelectorHist_TrackEff2_TimeBasedNonLowConLev;

	TH2F *dSelectorHist_TrackEff2_CandidateChargeCorrect_DeltaPOverPVsTheta;
	TH2F *dSelectorHist_TrackEff2_CandidateChargeWrong_DeltaPOverPVsTheta;
	TH2F *dSelectorHist_TrackEff2_WireBasedChargeCorrect_DeltaPOverPVsTheta;
	TH2F *dSelectorHist_TrackEff2_WireBasedChargeWrong_DeltaPOverPVsTheta;
	TH2F *dSelectorHist_TrackEff2_TimeBasedChargeCorrect_DeltaPOverPVsTheta;
	TH2F *dSelectorHist_TrackEff2_TimeBasedChargeWrong_DeltaPOverPVsTheta;

	TH2F *dSelectorHist_TrackEff2_CandidateChargeCorrect_DeltaPtOverPtVsTheta;
	TH2F *dSelectorHist_TrackEff2_CandidateChargeWrong_DeltaPtOverPtVsTheta;
	TH2F *dSelectorHist_TrackEff2_WireBasedChargeCorrect_DeltaPtOverPtVsTheta;
	TH2F *dSelectorHist_TrackEff2_WireBasedChargeWrong_DeltaPtOverPtVsTheta;
	TH2F *dSelectorHist_TrackEff2_TimeBasedChargeCorrect_DeltaPtOverPtVsTheta;
	TH2F *dSelectorHist_TrackEff2_TimeBasedChargeWrong_DeltaPtOverPtVsTheta;

	TH2F *dSelectorHist_TrackEff2_CandidateChargeCorrect_DeltaThetaVsTheta;
	TH2F *dSelectorHist_TrackEff2_CandidateChargeWrong_DeltaThetaVsTheta;
	TH2F *dSelectorHist_TrackEff2_WireBasedChargeCorrect_DeltaThetaVsTheta;
	TH2F *dSelectorHist_TrackEff2_WireBasedChargeWrong_DeltaThetaVsTheta;
	TH2F *dSelectorHist_TrackEff2_TimeBasedChargeCorrect_DeltaThetaVsTheta;
	TH2F *dSelectorHist_TrackEff2_TimeBasedChargeWrong_DeltaThetaVsTheta;

	TH2F *dSelectorHist_TrackEff2_CandidateChargeCorrect_DeltaPhiVsTheta;
	TH2F *dSelectorHist_TrackEff2_CandidateChargeWrong_DeltaPhiVsTheta;
	TH2F *dSelectorHist_TrackEff2_WireBasedChargeCorrect_DeltaPhiVsTheta;
	TH2F *dSelectorHist_TrackEff2_WireBasedChargeWrong_DeltaPhiVsTheta;
	TH2F *dSelectorHist_TrackEff2_TimeBasedChargeCorrect_DeltaPhiVsTheta;
	TH2F *dSelectorHist_TrackEff2_TimeBasedChargeWrong_DeltaPhiVsTheta;

   TrackEff2Selector(TTree * /*tree*/ =0) { }
   virtual ~TrackEff2Selector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(TrackEff2Selector,0);
};

#endif

#ifdef TrackEff2Selector_cxx
void TrackEff2Selector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_E_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_E_fBits);
   fChain->SetBranchAddress("event", &event, &b_E_event);
   fChain->SetBranchAddress("q_thrown", &q_thrown, &b_E_q_thrown);
   fChain->SetBranchAddress("PID_thrown", &PID_thrown, &b_E_PID_thrown);
   fChain->SetBranchAddress("pthrown.fUniqueID", &pthrown_fUniqueID, &b_E_pthrown_fUniqueID);
   fChain->SetBranchAddress("pthrown.fBits", &pthrown_fBits, &b_E_pthrown_fBits);
   fChain->SetBranchAddress("pthrown.fX", &pthrown_fX, &b_E_pthrown_fX);
   fChain->SetBranchAddress("pthrown.fY", &pthrown_fY, &b_E_pthrown_fY);
   fChain->SetBranchAddress("pthrown.fZ", &pthrown_fZ, &b_E_pthrown_fZ);
   fChain->SetBranchAddress("isreconstructable", &isreconstructable, &b_E_isreconstructable);
   fChain->SetBranchAddress("Nstereo", &Nstereo, &b_E_Nstereo);
   fChain->SetBranchAddress("Ncdc", &Ncdc, &b_E_Ncdc);
   fChain->SetBranchAddress("Nfdc", &Nfdc, &b_E_Nfdc);
   fChain->SetBranchAddress("NLR_bad_stereo", &NLR_bad_stereo, &b_E_NLR_bad_stereo);
   fChain->SetBranchAddress("NLR_bad", &NLR_bad, &b_E_NLR_bad);
   fChain->SetBranchAddress("PID_hypothesized", &PID_hypothesized, &b_E_PID_hypothesized);
   fChain->SetBranchAddress("q_hypothesized", &q_hypothesized, &b_E_q_hypothesized);
   fChain->SetBranchAddress("FOM_hypothesized", &FOM_hypothesized, &b_E_FOM_hypothesized);
   fChain->SetBranchAddress("q_timebased", &q_timebased, &b_E_q_timebased);
   fChain->SetBranchAddress("PID_timebased", &PID_timebased, &b_E_PID_timebased);
   fChain->SetBranchAddress("dTrackReconstructedFlag_TimeBased", &dTrackReconstructedFlag_TimeBased, &b_E_dTrackReconstructedFlag_TimeBased);
   fChain->SetBranchAddress("pfit.fUniqueID", &pfit_fUniqueID, &b_E_pfit_fUniqueID);
   fChain->SetBranchAddress("pfit.fBits", &pfit_fBits, &b_E_pfit_fBits);
   fChain->SetBranchAddress("pfit.fX", &pfit_fX, &b_E_pfit_fX);
   fChain->SetBranchAddress("pfit.fY", &pfit_fY, &b_E_pfit_fY);
   fChain->SetBranchAddress("pfit.fZ", &pfit_fZ, &b_E_pfit_fZ);
   fChain->SetBranchAddress("trk_chisq", &trk_chisq, &b_E_trk_chisq);
   fChain->SetBranchAddress("trk_Ndof", &trk_Ndof, &b_E_trk_Ndof);
   fChain->SetBranchAddress("delta_pt_over_pt", &delta_pt_over_pt, &b_E_delta_pt_over_pt);
   fChain->SetBranchAddress("delta_theta", &delta_theta, &b_E_delta_theta);
   fChain->SetBranchAddress("delta_phi", &delta_phi, &b_E_delta_phi);
   fChain->SetBranchAddress("num_timebased", &num_timebased, &b_E_num_timebased);
   fChain->SetBranchAddress("q_wirebased", &q_wirebased, &b_E_q_wirebased);
   fChain->SetBranchAddress("PID_wirebased", &PID_wirebased, &b_E_PID_wirebased);
   fChain->SetBranchAddress("dTrackReconstructedFlag_WireBased", &dTrackReconstructedFlag_WireBased, &b_E_dTrackReconstructedFlag_WireBased);
   fChain->SetBranchAddress("pfit_wire.fUniqueID", &pfit_wire_fUniqueID, &b_E_pfit_wire_fUniqueID);
   fChain->SetBranchAddress("pfit_wire.fBits", &pfit_wire_fBits, &b_E_pfit_wire_fBits);
   fChain->SetBranchAddress("pfit_wire.fX", &pfit_wire_fX, &b_E_pfit_wire_fX);
   fChain->SetBranchAddress("pfit_wire.fY", &pfit_wire_fY, &b_E_pfit_wire_fY);
   fChain->SetBranchAddress("pfit_wire.fZ", &pfit_wire_fZ, &b_E_pfit_wire_fZ);
   fChain->SetBranchAddress("trk_chisq_wb", &trk_chisq_wb, &b_E_trk_chisq_wb);
   fChain->SetBranchAddress("trk_Ndof_wb", &trk_Ndof_wb, &b_E_trk_Ndof_wb);
   fChain->SetBranchAddress("delta_pt_over_pt_wire", &delta_pt_over_pt_wire, &b_E_delta_pt_over_pt_wire);
   fChain->SetBranchAddress("delta_theta_wire", &delta_theta_wire, &b_E_delta_theta_wire);
   fChain->SetBranchAddress("delta_phi_wire", &delta_phi_wire, &b_E_delta_phi_wire);
   fChain->SetBranchAddress("num_wirebased", &num_wirebased, &b_E_num_wirebased);
   fChain->SetBranchAddress("q_candidate", &q_candidate, &b_E_q_candidate);
   fChain->SetBranchAddress("PID_candidate", &PID_candidate, &b_E_PID_candidate);
   fChain->SetBranchAddress("dTrackReconstructedFlag_Candidate", &dTrackReconstructedFlag_Candidate, &b_E_dTrackReconstructedFlag_Candidate);
   fChain->SetBranchAddress("pcan.fUniqueID", &pcan_fUniqueID, &b_E_pcan_fUniqueID);
   fChain->SetBranchAddress("pcan.fBits", &pcan_fBits, &b_E_pcan_fBits);
   fChain->SetBranchAddress("pcan.fX", &pcan_fX, &b_E_pcan_fX);
   fChain->SetBranchAddress("pcan.fY", &pcan_fY, &b_E_pcan_fY);
   fChain->SetBranchAddress("pcan.fZ", &pcan_fZ, &b_E_pcan_fZ);
   fChain->SetBranchAddress("delta_pt_over_pt_can", &delta_pt_over_pt_can, &b_E_delta_pt_over_pt_can);
   fChain->SetBranchAddress("delta_theta_can", &delta_theta_can, &b_E_delta_theta_can);
   fChain->SetBranchAddress("delta_phi_can", &delta_phi_can, &b_E_delta_phi_can);
   fChain->SetBranchAddress("num_candidates", &num_candidates, &b_E_num_candidates);
}

Bool_t TrackEff2Selector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef TrackEff2Selector_cxx
