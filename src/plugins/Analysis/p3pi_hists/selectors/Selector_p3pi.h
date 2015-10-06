//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Oct  3 21:16:24 2015 by ROOT version 6.04/02
// from TTree p3pi_preco_any_Tree/p3pi_preco_any_Tree
// found on file: tree_p3pi_3185.root
//////////////////////////////////////////////////////////

#ifndef Selector_p3pi_h
#define Selector_p3pi_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"

// Added by me
#include "TFile.h"
#include "TH1I.h"
#include "TLorentzVector.h"

#include <set>

class Selector_p3pi : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

	//TARGET P4
	TLorentzVector dTargetP4;

	//DEFINE OUTPUT ROOT FILE
	TFile* dFile;

	//DEFINE HISTOGRAMS
	TH1I* dHist_Pi0Mass_Measured;
	TH1I* dHist_Pi0Mass_KinFit;
	TH1I* dHist_MissingMassSquared;
	TH1I* dHist_OmegaMass_Measured;
	TH1I* dHist_OmegaMass_KinFit;

	//FOR BELOW: NOTE THE FIXED ARRAY SIZES. BE CAREFUL WHEN USING THIS SELECTOR DIRECTLY!!!! 

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          RunNumber;
   ULong64_t       EventNumber;
   UInt_t          NumBeam;
   Int_t           Beam__PID[11];   //[NumBeam]
   TClonesArray    *Beam__X4_TargetCenter;
   TClonesArray    *Beam__P4_Measured;
   UInt_t          NumNeutralShowers;
   TClonesArray    *NeutralShower__X4_Shower;
   Float_t         NeutralShower__Energy_BCAL[8];   //[NumNeutralShowers]
   Float_t         NeutralShower__Energy_FCAL[8];   //[NumNeutralShowers]
   Float_t         NeutralShower__TrackBCAL_DeltaPhi[8];   //[NumNeutralShowers]
   Float_t         NeutralShower__TrackBCAL_DeltaZ[8];   //[NumNeutralShowers]
   Float_t         NeutralShower__TrackFCAL_DOCA[8];   //[NumNeutralShowers]
   Float_t         NeutralShower__PhotonRFDeltaTVar[8];   //[NumNeutralShowers]
   UInt_t          NumChargedHypos;
   Int_t           ChargedHypo__TrackID[17];   //[NumChargedHypos]
   Int_t           ChargedHypo__PID[17];   //[NumChargedHypos]
   TClonesArray    *ChargedHypo__X4_Measured;
   TClonesArray    *ChargedHypo__P4_Measured;
   UInt_t          ChargedHypo__NDF_Tracking[17];   //[NumChargedHypos]
   Float_t         ChargedHypo__ChiSq_Tracking[17];   //[NumChargedHypos]
   UInt_t          ChargedHypo__NDF_DCdEdx[17];   //[NumChargedHypos]
   Float_t         ChargedHypo__ChiSq_DCdEdx[17];   //[NumChargedHypos]
   Float_t         ChargedHypo__dEdx_CDC[17];   //[NumChargedHypos]
   Float_t         ChargedHypo__dEdx_FDC[17];   //[NumChargedHypos]
   Float_t         ChargedHypo__HitTime[17];   //[NumChargedHypos]
   Float_t         ChargedHypo__RFDeltaTVar[17];   //[NumChargedHypos]
   Float_t         ChargedHypo__dEdx_TOF[17];   //[NumChargedHypos]
   Float_t         ChargedHypo__dEdx_ST[17];   //[NumChargedHypos]
   Float_t         ChargedHypo__Energy_BCAL[17];   //[NumChargedHypos]
   Float_t         ChargedHypo__Energy_FCAL[17];   //[NumChargedHypos]
   Float_t         ChargedHypo__TrackBCAL_DeltaPhi[17];   //[NumChargedHypos]
   Float_t         ChargedHypo__TrackBCAL_DeltaZ[17];   //[NumChargedHypos]
   Float_t         ChargedHypo__TrackFCAL_DOCA[17];   //[NumChargedHypos]
   UInt_t          NumCombos;
   Bool_t          IsComboCut[23];   //[NumCombos]
   Float_t         RFTime_Measured[23];   //[NumCombos]
   Float_t         ChiSq_KinFit[23];   //[NumCombos]
   UInt_t          NDF_KinFit[23];   //[NumCombos]
   Int_t           ComboBeam__BeamIndex[23];   //[NumCombos]
   TClonesArray    *ComboBeam__X4_Measured;
   TClonesArray    *ComboBeam__P4_KinFit;
   Int_t           Proton__ChargedIndex[23];   //[NumCombos]
   Float_t         Proton__Beta_Timing_Measured[23];   //[NumCombos]
   Float_t         Proton__ChiSq_Timing_Measured[23];   //[NumCombos]
   UInt_t          Proton__NDF_Timing[23];   //[NumCombos]
   TClonesArray    *Proton__P4_KinFit;
   Int_t           PiPlus__ChargedIndex[23];   //[NumCombos]
   Float_t         PiPlus__Beta_Timing_Measured[23];   //[NumCombos]
   Float_t         PiPlus__ChiSq_Timing_Measured[23];   //[NumCombos]
   UInt_t          PiPlus__NDF_Timing[23];   //[NumCombos]
   TClonesArray    *PiPlus__P4_KinFit;
   Int_t           PiMinus__ChargedIndex[23];   //[NumCombos]
   Float_t         PiMinus__Beta_Timing_Measured[23];   //[NumCombos]
   Float_t         PiMinus__ChiSq_Timing_Measured[23];   //[NumCombos]
   UInt_t          PiMinus__NDF_Timing[23];   //[NumCombos]
   TClonesArray    *PiMinus__P4_KinFit;
   TClonesArray    *DecayingPi0__P4_KinFit;
   Int_t           Photon1__ShowerIndex[23];   //[NumCombos]
   TClonesArray    *Photon1__X4_Measured;
   TClonesArray    *Photon1__P4_Measured;
   Float_t         Photon1__Beta_Timing_Measured[23];   //[NumCombos]
   Float_t         Photon1__ChiSq_Timing_Measured[23];   //[NumCombos]
   UInt_t          Photon1__NDF_Timing[23];   //[NumCombos]
   TClonesArray    *Photon1__P4_KinFit;
   Int_t           Photon2__ShowerIndex[23];   //[NumCombos]
   TClonesArray    *Photon2__X4_Measured;
   TClonesArray    *Photon2__P4_Measured;
   Float_t         Photon2__Beta_Timing_Measured[23];   //[NumCombos]
   Float_t         Photon2__ChiSq_Timing_Measured[23];   //[NumCombos]
   UInt_t          Photon2__NDF_Timing[23];   //[NumCombos]
   TClonesArray    *Photon2__P4_KinFit;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_NumBeam;   //!
   TBranch        *b_Beam__PID;   //!
   TBranch        *b_Beam__X4_TargetCenter;   //!
   TBranch        *b_Beam__P4_Measured;   //!
   TBranch        *b_NumNeutralShowers;   //!
   TBranch        *b_NeutralShower__X4_Shower;   //!
   TBranch        *b_NeutralShower__Energy_BCAL;   //!
   TBranch        *b_NeutralShower__Energy_FCAL;   //!
   TBranch        *b_NeutralShower__TrackBCAL_DeltaPhi;   //!
   TBranch        *b_NeutralShower__TrackBCAL_DeltaZ;   //!
   TBranch        *b_NeutralShower__TrackFCAL_DOCA;   //!
   TBranch        *b_NeutralShower__PhotonRFDeltaTVar;   //!
   TBranch        *b_NumChargedHypos;   //!
   TBranch        *b_ChargedHypo__TrackID;   //!
   TBranch        *b_ChargedHypo__PID;   //!
   TBranch        *b_ChargedHypo__X4_Measured;   //!
   TBranch        *b_ChargedHypo__P4_Measured;   //!
   TBranch        *b_ChargedHypo__NDF_Tracking;   //!
   TBranch        *b_ChargedHypo__ChiSq_Tracking;   //!
   TBranch        *b_ChargedHypo__NDF_DCdEdx;   //!
   TBranch        *b_ChargedHypo__ChiSq_DCdEdx;   //!
   TBranch        *b_ChargedHypo__dEdx_CDC;   //!
   TBranch        *b_ChargedHypo__dEdx_FDC;   //!
   TBranch        *b_ChargedHypo__HitTime;   //!
   TBranch        *b_ChargedHypo__RFDeltaTVar;   //!
   TBranch        *b_ChargedHypo__dEdx_TOF;   //!
   TBranch        *b_ChargedHypo__dEdx_ST;   //!
   TBranch        *b_ChargedHypo__Energy_BCAL;   //!
   TBranch        *b_ChargedHypo__Energy_FCAL;   //!
   TBranch        *b_ChargedHypo__TrackBCAL_DeltaPhi;   //!
   TBranch        *b_ChargedHypo__TrackBCAL_DeltaZ;   //!
   TBranch        *b_ChargedHypo__TrackFCAL_DOCA;   //!
   TBranch        *b_NumCombos;   //!
   TBranch        *b_IsComboCut;   //!
   TBranch        *b_RFTime_Measured;   //!
   TBranch        *b_ChiSq_KinFit;   //!
   TBranch        *b_NDF_KinFit;   //!
   TBranch        *b_ComboBeam__BeamIndex;   //!
   TBranch        *b_ComboBeam__X4_Measured;   //!
   TBranch        *b_ComboBeam__P4_KinFit;   //!
   TBranch        *b_Proton__ChargedIndex;   //!
   TBranch        *b_Proton__Beta_Timing_Measured;   //!
   TBranch        *b_Proton__ChiSq_Timing_Measured;   //!
   TBranch        *b_Proton__NDF_Timing;   //!
   TBranch        *b_Proton__P4_KinFit;   //!
   TBranch        *b_PiPlus__ChargedIndex;   //!
   TBranch        *b_PiPlus__Beta_Timing_Measured;   //!
   TBranch        *b_PiPlus__ChiSq_Timing_Measured;   //!
   TBranch        *b_PiPlus__NDF_Timing;   //!
   TBranch        *b_PiPlus__P4_KinFit;   //!
   TBranch        *b_PiMinus__ChargedIndex;   //!
   TBranch        *b_PiMinus__Beta_Timing_Measured;   //!
   TBranch        *b_PiMinus__ChiSq_Timing_Measured;   //!
   TBranch        *b_PiMinus__NDF_Timing;   //!
   TBranch        *b_PiMinus__P4_KinFit;   //!
   TBranch        *b_DecayingPi0__P4_KinFit;   //!
   TBranch        *b_Photon1__ShowerIndex;   //!
   TBranch        *b_Photon1__X4_Measured;   //!
   TBranch        *b_Photon1__P4_Measured;   //!
   TBranch        *b_Photon1__Beta_Timing_Measured;   //!
   TBranch        *b_Photon1__ChiSq_Timing_Measured;   //!
   TBranch        *b_Photon1__NDF_Timing;   //!
   TBranch        *b_Photon1__P4_KinFit;   //!
   TBranch        *b_Photon2__ShowerIndex;   //!
   TBranch        *b_Photon2__X4_Measured;   //!
   TBranch        *b_Photon2__P4_Measured;   //!
   TBranch        *b_Photon2__Beta_Timing_Measured;   //!
   TBranch        *b_Photon2__ChiSq_Timing_Measured;   //!
   TBranch        *b_Photon2__NDF_Timing;   //!
   TBranch        *b_Photon2__P4_KinFit;   //!

   Selector_p3pi(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~Selector_p3pi() { }
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

   ClassDef(Selector_p3pi,0);
};

#endif

#ifdef Selector_p3pi_cxx
void Selector_p3pi::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Beam__X4_TargetCenter = 0;
   Beam__P4_Measured = 0;
   NeutralShower__X4_Shower = 0;
   ChargedHypo__X4_Measured = 0;
   ChargedHypo__P4_Measured = 0;
   ComboBeam__X4_Measured = 0;
   ComboBeam__P4_KinFit = 0;
   Proton__P4_KinFit = 0;
   PiPlus__P4_KinFit = 0;
   PiMinus__P4_KinFit = 0;
   DecayingPi0__P4_KinFit = 0;
   Photon1__X4_Measured = 0;
   Photon1__P4_Measured = 0;
   Photon1__P4_KinFit = 0;
   Photon2__X4_Measured = 0;
   Photon2__P4_Measured = 0;
   Photon2__P4_KinFit = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("NumBeam", &NumBeam, &b_NumBeam);
   fChain->SetBranchAddress("Beam__PID", Beam__PID, &b_Beam__PID);
   fChain->SetBranchAddress("Beam__X4_TargetCenter", &Beam__X4_TargetCenter, &b_Beam__X4_TargetCenter);
   fChain->SetBranchAddress("Beam__P4_Measured", &Beam__P4_Measured, &b_Beam__P4_Measured);
   fChain->SetBranchAddress("NumNeutralShowers", &NumNeutralShowers, &b_NumNeutralShowers);
   fChain->SetBranchAddress("NeutralShower__X4_Shower", &NeutralShower__X4_Shower, &b_NeutralShower__X4_Shower);
   fChain->SetBranchAddress("NeutralShower__Energy_BCAL", NeutralShower__Energy_BCAL, &b_NeutralShower__Energy_BCAL);
   fChain->SetBranchAddress("NeutralShower__Energy_FCAL", NeutralShower__Energy_FCAL, &b_NeutralShower__Energy_FCAL);
   fChain->SetBranchAddress("NeutralShower__TrackBCAL_DeltaPhi", NeutralShower__TrackBCAL_DeltaPhi, &b_NeutralShower__TrackBCAL_DeltaPhi);
   fChain->SetBranchAddress("NeutralShower__TrackBCAL_DeltaZ", NeutralShower__TrackBCAL_DeltaZ, &b_NeutralShower__TrackBCAL_DeltaZ);
   fChain->SetBranchAddress("NeutralShower__TrackFCAL_DOCA", NeutralShower__TrackFCAL_DOCA, &b_NeutralShower__TrackFCAL_DOCA);
   fChain->SetBranchAddress("NeutralShower__PhotonRFDeltaTVar", NeutralShower__PhotonRFDeltaTVar, &b_NeutralShower__PhotonRFDeltaTVar);
   fChain->SetBranchAddress("NumChargedHypos", &NumChargedHypos, &b_NumChargedHypos);
   fChain->SetBranchAddress("ChargedHypo__TrackID", ChargedHypo__TrackID, &b_ChargedHypo__TrackID);
   fChain->SetBranchAddress("ChargedHypo__PID", ChargedHypo__PID, &b_ChargedHypo__PID);
   fChain->SetBranchAddress("ChargedHypo__X4_Measured", &ChargedHypo__X4_Measured, &b_ChargedHypo__X4_Measured);
   fChain->SetBranchAddress("ChargedHypo__P4_Measured", &ChargedHypo__P4_Measured, &b_ChargedHypo__P4_Measured);
   fChain->SetBranchAddress("ChargedHypo__NDF_Tracking", ChargedHypo__NDF_Tracking, &b_ChargedHypo__NDF_Tracking);
   fChain->SetBranchAddress("ChargedHypo__ChiSq_Tracking", ChargedHypo__ChiSq_Tracking, &b_ChargedHypo__ChiSq_Tracking);
   fChain->SetBranchAddress("ChargedHypo__NDF_DCdEdx", ChargedHypo__NDF_DCdEdx, &b_ChargedHypo__NDF_DCdEdx);
   fChain->SetBranchAddress("ChargedHypo__ChiSq_DCdEdx", ChargedHypo__ChiSq_DCdEdx, &b_ChargedHypo__ChiSq_DCdEdx);
   fChain->SetBranchAddress("ChargedHypo__dEdx_CDC", ChargedHypo__dEdx_CDC, &b_ChargedHypo__dEdx_CDC);
   fChain->SetBranchAddress("ChargedHypo__dEdx_FDC", ChargedHypo__dEdx_FDC, &b_ChargedHypo__dEdx_FDC);
   fChain->SetBranchAddress("ChargedHypo__HitTime", ChargedHypo__HitTime, &b_ChargedHypo__HitTime);
   fChain->SetBranchAddress("ChargedHypo__RFDeltaTVar", ChargedHypo__RFDeltaTVar, &b_ChargedHypo__RFDeltaTVar);
   fChain->SetBranchAddress("ChargedHypo__dEdx_TOF", ChargedHypo__dEdx_TOF, &b_ChargedHypo__dEdx_TOF);
   fChain->SetBranchAddress("ChargedHypo__dEdx_ST", ChargedHypo__dEdx_ST, &b_ChargedHypo__dEdx_ST);
   fChain->SetBranchAddress("ChargedHypo__Energy_BCAL", ChargedHypo__Energy_BCAL, &b_ChargedHypo__Energy_BCAL);
   fChain->SetBranchAddress("ChargedHypo__Energy_FCAL", ChargedHypo__Energy_FCAL, &b_ChargedHypo__Energy_FCAL);
   fChain->SetBranchAddress("ChargedHypo__TrackBCAL_DeltaPhi", ChargedHypo__TrackBCAL_DeltaPhi, &b_ChargedHypo__TrackBCAL_DeltaPhi);
   fChain->SetBranchAddress("ChargedHypo__TrackBCAL_DeltaZ", ChargedHypo__TrackBCAL_DeltaZ, &b_ChargedHypo__TrackBCAL_DeltaZ);
   fChain->SetBranchAddress("ChargedHypo__TrackFCAL_DOCA", ChargedHypo__TrackFCAL_DOCA, &b_ChargedHypo__TrackFCAL_DOCA);
   fChain->SetBranchAddress("NumCombos", &NumCombos, &b_NumCombos);
   fChain->SetBranchAddress("IsComboCut", IsComboCut, &b_IsComboCut);
   fChain->SetBranchAddress("RFTime_Measured", RFTime_Measured, &b_RFTime_Measured);
   fChain->SetBranchAddress("ChiSq_KinFit", ChiSq_KinFit, &b_ChiSq_KinFit);
   fChain->SetBranchAddress("NDF_KinFit", NDF_KinFit, &b_NDF_KinFit);
   fChain->SetBranchAddress("ComboBeam__BeamIndex", ComboBeam__BeamIndex, &b_ComboBeam__BeamIndex);
   fChain->SetBranchAddress("ComboBeam__X4_Measured", &ComboBeam__X4_Measured, &b_ComboBeam__X4_Measured);
   fChain->SetBranchAddress("ComboBeam__P4_KinFit", &ComboBeam__P4_KinFit, &b_ComboBeam__P4_KinFit);
   fChain->SetBranchAddress("Proton__ChargedIndex", Proton__ChargedIndex, &b_Proton__ChargedIndex);
   fChain->SetBranchAddress("Proton__Beta_Timing_Measured", Proton__Beta_Timing_Measured, &b_Proton__Beta_Timing_Measured);
   fChain->SetBranchAddress("Proton__ChiSq_Timing_Measured", Proton__ChiSq_Timing_Measured, &b_Proton__ChiSq_Timing_Measured);
   fChain->SetBranchAddress("Proton__NDF_Timing", Proton__NDF_Timing, &b_Proton__NDF_Timing);
   fChain->SetBranchAddress("Proton__P4_KinFit", &Proton__P4_KinFit, &b_Proton__P4_KinFit);
   fChain->SetBranchAddress("PiPlus__ChargedIndex", PiPlus__ChargedIndex, &b_PiPlus__ChargedIndex);
   fChain->SetBranchAddress("PiPlus__Beta_Timing_Measured", PiPlus__Beta_Timing_Measured, &b_PiPlus__Beta_Timing_Measured);
   fChain->SetBranchAddress("PiPlus__ChiSq_Timing_Measured", PiPlus__ChiSq_Timing_Measured, &b_PiPlus__ChiSq_Timing_Measured);
   fChain->SetBranchAddress("PiPlus__NDF_Timing", PiPlus__NDF_Timing, &b_PiPlus__NDF_Timing);
   fChain->SetBranchAddress("PiPlus__P4_KinFit", &PiPlus__P4_KinFit, &b_PiPlus__P4_KinFit);
   fChain->SetBranchAddress("PiMinus__ChargedIndex", PiMinus__ChargedIndex, &b_PiMinus__ChargedIndex);
   fChain->SetBranchAddress("PiMinus__Beta_Timing_Measured", PiMinus__Beta_Timing_Measured, &b_PiMinus__Beta_Timing_Measured);
   fChain->SetBranchAddress("PiMinus__ChiSq_Timing_Measured", PiMinus__ChiSq_Timing_Measured, &b_PiMinus__ChiSq_Timing_Measured);
   fChain->SetBranchAddress("PiMinus__NDF_Timing", PiMinus__NDF_Timing, &b_PiMinus__NDF_Timing);
   fChain->SetBranchAddress("PiMinus__P4_KinFit", &PiMinus__P4_KinFit, &b_PiMinus__P4_KinFit);
   fChain->SetBranchAddress("DecayingPi0__P4_KinFit", &DecayingPi0__P4_KinFit, &b_DecayingPi0__P4_KinFit);
   fChain->SetBranchAddress("Photon1__ShowerIndex", Photon1__ShowerIndex, &b_Photon1__ShowerIndex);
   fChain->SetBranchAddress("Photon1__X4_Measured", &Photon1__X4_Measured, &b_Photon1__X4_Measured);
   fChain->SetBranchAddress("Photon1__P4_Measured", &Photon1__P4_Measured, &b_Photon1__P4_Measured);
   fChain->SetBranchAddress("Photon1__Beta_Timing_Measured", Photon1__Beta_Timing_Measured, &b_Photon1__Beta_Timing_Measured);
   fChain->SetBranchAddress("Photon1__ChiSq_Timing_Measured", Photon1__ChiSq_Timing_Measured, &b_Photon1__ChiSq_Timing_Measured);
   fChain->SetBranchAddress("Photon1__NDF_Timing", Photon1__NDF_Timing, &b_Photon1__NDF_Timing);
   fChain->SetBranchAddress("Photon1__P4_KinFit", &Photon1__P4_KinFit, &b_Photon1__P4_KinFit);
   fChain->SetBranchAddress("Photon2__ShowerIndex", Photon2__ShowerIndex, &b_Photon2__ShowerIndex);
   fChain->SetBranchAddress("Photon2__X4_Measured", &Photon2__X4_Measured, &b_Photon2__X4_Measured);
   fChain->SetBranchAddress("Photon2__P4_Measured", &Photon2__P4_Measured, &b_Photon2__P4_Measured);
   fChain->SetBranchAddress("Photon2__Beta_Timing_Measured", Photon2__Beta_Timing_Measured, &b_Photon2__Beta_Timing_Measured);
   fChain->SetBranchAddress("Photon2__ChiSq_Timing_Measured", Photon2__ChiSq_Timing_Measured, &b_Photon2__ChiSq_Timing_Measured);
   fChain->SetBranchAddress("Photon2__NDF_Timing", Photon2__NDF_Timing, &b_Photon2__NDF_Timing);
   fChain->SetBranchAddress("Photon2__P4_KinFit", &Photon2__P4_KinFit, &b_Photon2__P4_KinFit);
}

Bool_t Selector_p3pi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef Selector_p3pi_cxx
