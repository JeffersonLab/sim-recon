//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep  6 22:46:15 2011 by ROOT version 5.30/00
// from TTree dPluginTree_BCALMCComparison/BCAL MC Comparison
// found on file: /raid1/mattione/gluex/cluster_run/hd_root_bcal.root
//////////////////////////////////////////////////////////

#ifndef BCALSelector_h
#define BCALSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

class BCALSelector : public TSelector {
public :

	float Calc_BCALPathLengthCorrection(float locEnergy);
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
 //BCALMCComparison *dPluginBranch_BCALMCComparison;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Float_t         dTrueR;
   Float_t         dTrueZ;
   Float_t         dTruePhi;
   Float_t         dTrueE;
   Float_t         dTrueT;
   Float_t         dDeltaR;
   Float_t         dDeltaZ;
   Float_t         dDeltaPhi;
   Float_t         dDeltaE;
   Float_t         dDeltaT;
   Float_t         dShowerUncertaintyX;
   Float_t         dShowerUncertaintyY;
   Float_t         dShowerUncertaintyZ;
   Float_t         dShowerUncertaintyT;
   Float_t         dShowerUncertaintyE;
   Float_t         dPathLengthCorrection;

   // List of branches
   TBranch        *b_dPluginBranch_BCALMCComparison_fUniqueID;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_fBits;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_dTrueR;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_dTrueZ;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_dTruePhi;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_dTrueE;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_dTrueT;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_dDeltaR;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_dDeltaZ;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_dDeltaPhi;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_dDeltaE;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_dDeltaT;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_dShowerUncertaintyX;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_dShowerUncertaintyY;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_dShowerUncertaintyZ;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_dShowerUncertaintyT;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_dShowerUncertaintyE;   //!
   TBranch        *b_dPluginBranch_BCALMCComparison_dPathLengthCorrection;   //!

		TFile *dOutputFile;
		//Path Length Correction
		TH2F *dPluginHist_BCAL_PathLengthCorrection;
		TH2F *dPluginHist_BCAL_PathLengthCorrectionPostEVsZ;

		//DeltaR Dependence
		TH2F *dPluginHist_BCAL_DeltaRVsR;
		TH2F *dPluginHist_BCAL_DeltaRVsPhi;
		TH2F *dPluginHist_BCAL_DeltaRVsZ;
		TH2F *dPluginHist_BCAL_DeltaRVsE;
		TH2F *dPluginHist_BCAL_DeltaRVsT;

		//DeltaPhi Dependence
		TH2F *dPluginHist_BCAL_DeltaPhiVsR;
		TH2F *dPluginHist_BCAL_DeltaPhiVsPhi;
		TH2F *dPluginHist_BCAL_DeltaPhiVsZ;
		TH2F *dPluginHist_BCAL_DeltaPhiVsE;
		TH2F *dPluginHist_BCAL_DeltaPhiVsT;

		//DeltaZ Dependence
		TH2F *dPluginHist_BCAL_DeltaZVsR;
		TH2F *dPluginHist_BCAL_DeltaZVsPhi;
		TH2F *dPluginHist_BCAL_DeltaZVsZ;
		TH2F *dPluginHist_BCAL_DeltaZVsE;
		TH2F *dPluginHist_BCAL_DeltaZVsT;

		//DeltaE Dependence
		TH2F *dPluginHist_BCAL_DeltaEVsR;
		TH2F *dPluginHist_BCAL_DeltaEVsPhi;
		TH2F *dPluginHist_BCAL_DeltaEVsZ;
		TH2F *dPluginHist_BCAL_DeltaEVsE;
		TH2F *dPluginHist_BCAL_DeltaEVsT;

		//DeltaT Dependence
		TH2F *dPluginHist_BCAL_DeltaTVsR;
		TH2F *dPluginHist_BCAL_DeltaTVsPhi;
		TH2F *dPluginHist_BCAL_DeltaTVsZ;
		TH2F *dPluginHist_BCAL_DeltaTVsE;
		TH2F *dPluginHist_BCAL_DeltaTVsT;

		//Common Dependence
		TH2F *dPluginHist_BCAL_DeltaRVsDeltaPhi;
		TH2F *dPluginHist_BCAL_DeltaRVsDeltaZ;
		TH2F *dPluginHist_BCAL_DeltaRVsDeltaE;
		TH2F *dPluginHist_BCAL_DeltaRVsDeltaT;
		TH2F *dPluginHist_BCAL_DeltaPhiVsDeltaZ;
		TH2F *dPluginHist_BCAL_DeltaPhiVsDeltaE;
		TH2F *dPluginHist_BCAL_DeltaPhiVsDeltaT;
		TH2F *dPluginHist_BCAL_DeltaZVsDeltaE;
		TH2F *dPluginHist_BCAL_DeltaZVsDeltaT;
		TH2F *dPluginHist_BCAL_DeltaEVsDeltaT;

	//ShowerSigmaX Dependence
	TH2F *dPluginHist_BCAL_ShowerSigmaXVsX;
	TH2F *dPluginHist_BCAL_ShowerSigmaXVsY;
	TH2F *dPluginHist_BCAL_ShowerSigmaXVsZ;
	TH2F *dPluginHist_BCAL_ShowerSigmaXVsE;
	TH2F *dPluginHist_BCAL_ShowerSigmaXVsT;

	//ShowerSigmaY Dependence
	TH2F *dPluginHist_BCAL_ShowerSigmaYVsX;
	TH2F *dPluginHist_BCAL_ShowerSigmaYVsY;
	TH2F *dPluginHist_BCAL_ShowerSigmaYVsZ;
	TH2F *dPluginHist_BCAL_ShowerSigmaYVsE;
	TH2F *dPluginHist_BCAL_ShowerSigmaYVsT;

	//ShowerSigmaZ Dependence
	TH2F *dPluginHist_BCAL_ShowerSigmaZVsX;
	TH2F *dPluginHist_BCAL_ShowerSigmaZVsY;
	TH2F *dPluginHist_BCAL_ShowerSigmaZVsZ;
	TH2F *dPluginHist_BCAL_ShowerSigmaZVsE;
	TH2F *dPluginHist_BCAL_ShowerSigmaZVsT;

	//ShowerSigmaE Dependence
	TH2F *dPluginHist_BCAL_ShowerSigmaEVsX;
	TH2F *dPluginHist_BCAL_ShowerSigmaEVsY;
	TH2F *dPluginHist_BCAL_ShowerSigmaEVsZ;
	TH2F *dPluginHist_BCAL_ShowerSigmaEVsE;
	TH2F *dPluginHist_BCAL_ShowerSigmaEVsT;

	//ShowerSigmaT Dependence
	TH2F *dPluginHist_BCAL_ShowerSigmaTVsX;
	TH2F *dPluginHist_BCAL_ShowerSigmaTVsY;
	TH2F *dPluginHist_BCAL_ShowerSigmaTVsZ;
	TH2F *dPluginHist_BCAL_ShowerSigmaTVsE;
	TH2F *dPluginHist_BCAL_ShowerSigmaTVsT;


   BCALSelector(TTree * /*tree*/ =0) { }
   virtual ~BCALSelector() { }
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

   ClassDef(BCALSelector,0);
};

#endif

#ifdef BCALSelector_cxx
void BCALSelector::Init(TTree *tree)
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

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_dPluginBranch_BCALMCComparison_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_dPluginBranch_BCALMCComparison_fBits);
   fChain->SetBranchAddress("dTrueR", &dTrueR, &b_dPluginBranch_BCALMCComparison_dTrueR);
   fChain->SetBranchAddress("dTrueZ", &dTrueZ, &b_dPluginBranch_BCALMCComparison_dTrueZ);
   fChain->SetBranchAddress("dTruePhi", &dTruePhi, &b_dPluginBranch_BCALMCComparison_dTruePhi);
   fChain->SetBranchAddress("dTrueE", &dTrueE, &b_dPluginBranch_BCALMCComparison_dTrueE);
   fChain->SetBranchAddress("dTrueT", &dTrueT, &b_dPluginBranch_BCALMCComparison_dTrueT);
   fChain->SetBranchAddress("dDeltaR", &dDeltaR, &b_dPluginBranch_BCALMCComparison_dDeltaR);
   fChain->SetBranchAddress("dDeltaZ", &dDeltaZ, &b_dPluginBranch_BCALMCComparison_dDeltaZ);
   fChain->SetBranchAddress("dDeltaPhi", &dDeltaPhi, &b_dPluginBranch_BCALMCComparison_dDeltaPhi);
   fChain->SetBranchAddress("dDeltaE", &dDeltaE, &b_dPluginBranch_BCALMCComparison_dDeltaE);
   fChain->SetBranchAddress("dDeltaT", &dDeltaT, &b_dPluginBranch_BCALMCComparison_dDeltaT);
   fChain->SetBranchAddress("dShowerUncertaintyX", &dShowerUncertaintyX, &b_dPluginBranch_BCALMCComparison_dShowerUncertaintyX);
   fChain->SetBranchAddress("dShowerUncertaintyY", &dShowerUncertaintyY, &b_dPluginBranch_BCALMCComparison_dShowerUncertaintyY);
   fChain->SetBranchAddress("dShowerUncertaintyZ", &dShowerUncertaintyZ, &b_dPluginBranch_BCALMCComparison_dShowerUncertaintyZ);
   fChain->SetBranchAddress("dShowerUncertaintyT", &dShowerUncertaintyT, &b_dPluginBranch_BCALMCComparison_dShowerUncertaintyT);
   fChain->SetBranchAddress("dShowerUncertaintyE", &dShowerUncertaintyE, &b_dPluginBranch_BCALMCComparison_dShowerUncertaintyE);
   fChain->SetBranchAddress("dPathLengthCorrection", &dPathLengthCorrection, &b_dPluginBranch_BCALMCComparison_dPathLengthCorrection);
}

Bool_t BCALSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef BCALSelector_cxx
