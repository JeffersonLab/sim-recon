//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep  6 22:46:29 2011 by ROOT version 5.30/00
// from TTree dPluginTree_FCALMCComparison/FCAL MC Comparison
// found on file: /raid1/mattione/gluex/cluster_run/hd_root_fcal.root
//////////////////////////////////////////////////////////

#ifndef FCALSelector_h
#define FCALSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

class FCALSelector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
 //FCALMCComparison *dPluginBranch_FCALMCComparison;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Float_t         dTrueX;
   Float_t         dTrueY;
   Float_t         dTrueZ;
   Float_t         dTrueE;
   Float_t         dTrueT;
   Float_t         dDeltaX;
   Float_t         dDeltaY;
   Float_t         dDeltaZ;
   Float_t         dDeltaE;
   Float_t         dDeltaT;
   Float_t         dShowerUncertaintyX;
   Float_t         dShowerUncertaintyY;
   Float_t         dShowerUncertaintyZ;
   Float_t         dShowerUncertaintyT;
   Float_t         dShowerUncertaintyE;
   Float_t         dPathLengthCorrection;

   // List of branches
   TBranch        *b_dPluginBranch_FCALMCComparison_fUniqueID;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_fBits;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_dTrueX;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_dTrueY;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_dTrueZ;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_dTrueE;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_dTrueT;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_dDeltaX;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_dDeltaY;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_dDeltaZ;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_dDeltaE;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_dDeltaT;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_dShowerUncertaintyX;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_dShowerUncertaintyY;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_dShowerUncertaintyZ;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_dShowerUncertaintyT;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_dShowerUncertaintyE;   //!
   TBranch        *b_dPluginBranch_FCALMCComparison_dPathLengthCorrection;   //!

		TFile *dOutputFile;
		//Path Length Correction
		TH2F *dPluginHist_FCAL_PathLengthCorrection;

		//DeltaX Dependence
		TH2F *dPluginHist_FCAL_DeltaXVsX;
		TH2F *dPluginHist_FCAL_DeltaXVsY;
		TH2F *dPluginHist_FCAL_DeltaXVsZ;
		TH2F *dPluginHist_FCAL_DeltaXVsE;
		TH2F *dPluginHist_FCAL_DeltaXVsT;

		//DeltaY Dependence
		TH2F *dPluginHist_FCAL_DeltaYVsX;
		TH2F *dPluginHist_FCAL_DeltaYVsY;
		TH2F *dPluginHist_FCAL_DeltaYVsZ;
		TH2F *dPluginHist_FCAL_DeltaYVsE;
		TH2F *dPluginHist_FCAL_DeltaYVsT;

		//DeltaZ Dependence
		TH2F *dPluginHist_FCAL_DeltaZVsX;
		TH2F *dPluginHist_FCAL_DeltaZVsY;
		TH2F *dPluginHist_FCAL_DeltaZVsZ;
		TH2F *dPluginHist_FCAL_DeltaZVsE;
		TH2F *dPluginHist_FCAL_DeltaZVsT;

		//DeltaE Dependence
		TH2F *dPluginHist_FCAL_DeltaEVsX;
		TH2F *dPluginHist_FCAL_DeltaEVsY;
		TH2F *dPluginHist_FCAL_DeltaEVsZ;
		TH2F *dPluginHist_FCAL_DeltaEVsE;
		TH2F *dPluginHist_FCAL_DeltaEVsT;

		//DeltaT Dependence
		TH2F *dPluginHist_FCAL_DeltaTVsX;
		TH2F *dPluginHist_FCAL_DeltaTVsY;
		TH2F *dPluginHist_FCAL_DeltaTVsZ;
		TH2F *dPluginHist_FCAL_DeltaTVsE;
		TH2F *dPluginHist_FCAL_DeltaTVsT;

		//Common Dependence
		TH2F *dPluginHist_FCAL_DeltaXVsDeltaY;
		TH2F *dPluginHist_FCAL_DeltaXVsDeltaZ;
		TH2F *dPluginHist_FCAL_DeltaXVsDeltaE;
		TH2F *dPluginHist_FCAL_DeltaXVsDeltaT;
		TH2F *dPluginHist_FCAL_DeltaYVsDeltaZ;
		TH2F *dPluginHist_FCAL_DeltaYVsDeltaE;
		TH2F *dPluginHist_FCAL_DeltaYVsDeltaT;
		TH2F *dPluginHist_FCAL_DeltaZVsDeltaE;
		TH2F *dPluginHist_FCAL_DeltaZVsDeltaT;
		TH2F *dPluginHist_FCAL_DeltaEVsDeltaT;

	//ShowerSigmaX Dependence
	TH2F *dPluginHist_FCAL_ShowerSigmaXVsX;
	TH2F *dPluginHist_FCAL_ShowerSigmaXVsY;
	TH2F *dPluginHist_FCAL_ShowerSigmaXVsZ;
	TH2F *dPluginHist_FCAL_ShowerSigmaXVsE;
	TH2F *dPluginHist_FCAL_ShowerSigmaXVsT;

	//ShowerSigmaY Dependence
	TH2F *dPluginHist_FCAL_ShowerSigmaYVsX;
	TH2F *dPluginHist_FCAL_ShowerSigmaYVsY;
	TH2F *dPluginHist_FCAL_ShowerSigmaYVsZ;
	TH2F *dPluginHist_FCAL_ShowerSigmaYVsE;
	TH2F *dPluginHist_FCAL_ShowerSigmaYVsT;

	//ShowerSigmaZ Dependence
	TH2F *dPluginHist_FCAL_ShowerSigmaZVsX;
	TH2F *dPluginHist_FCAL_ShowerSigmaZVsY;
	TH2F *dPluginHist_FCAL_ShowerSigmaZVsZ;
	TH2F *dPluginHist_FCAL_ShowerSigmaZVsE;
	TH2F *dPluginHist_FCAL_ShowerSigmaZVsT;

	//ShowerSigmaE Dependence
	TH2F *dPluginHist_FCAL_ShowerSigmaEVsX;
	TH2F *dPluginHist_FCAL_ShowerSigmaEVsY;
	TH2F *dPluginHist_FCAL_ShowerSigmaEVsZ;
	TH2F *dPluginHist_FCAL_ShowerSigmaEVsE;
	TH2F *dPluginHist_FCAL_ShowerSigmaEVsT;

	//ShowerSigmaT Dependence
	TH2F *dPluginHist_FCAL_ShowerSigmaTVsX;
	TH2F *dPluginHist_FCAL_ShowerSigmaTVsY;
	TH2F *dPluginHist_FCAL_ShowerSigmaTVsZ;
	TH2F *dPluginHist_FCAL_ShowerSigmaTVsE;
	TH2F *dPluginHist_FCAL_ShowerSigmaTVsT;


   FCALSelector(TTree * /*tree*/ =0) { }
   virtual ~FCALSelector() { }
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

   ClassDef(FCALSelector,0);
};

#endif

#ifdef FCALSelector_cxx
void FCALSelector::Init(TTree *tree)
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

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_dPluginBranch_FCALMCComparison_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_dPluginBranch_FCALMCComparison_fBits);
   fChain->SetBranchAddress("dTrueX", &dTrueX, &b_dPluginBranch_FCALMCComparison_dTrueX);
   fChain->SetBranchAddress("dTrueY", &dTrueY, &b_dPluginBranch_FCALMCComparison_dTrueY);
   fChain->SetBranchAddress("dTrueZ", &dTrueZ, &b_dPluginBranch_FCALMCComparison_dTrueZ);
   fChain->SetBranchAddress("dTrueE", &dTrueE, &b_dPluginBranch_FCALMCComparison_dTrueE);
   fChain->SetBranchAddress("dTrueT", &dTrueT, &b_dPluginBranch_FCALMCComparison_dTrueT);
   fChain->SetBranchAddress("dDeltaX", &dDeltaX, &b_dPluginBranch_FCALMCComparison_dDeltaX);
   fChain->SetBranchAddress("dDeltaY", &dDeltaY, &b_dPluginBranch_FCALMCComparison_dDeltaY);
   fChain->SetBranchAddress("dDeltaZ", &dDeltaZ, &b_dPluginBranch_FCALMCComparison_dDeltaZ);
   fChain->SetBranchAddress("dDeltaE", &dDeltaE, &b_dPluginBranch_FCALMCComparison_dDeltaE);
   fChain->SetBranchAddress("dDeltaT", &dDeltaT, &b_dPluginBranch_FCALMCComparison_dDeltaT);
   fChain->SetBranchAddress("dShowerUncertaintyX", &dShowerUncertaintyX, &b_dPluginBranch_FCALMCComparison_dShowerUncertaintyX);
   fChain->SetBranchAddress("dShowerUncertaintyY", &dShowerUncertaintyY, &b_dPluginBranch_FCALMCComparison_dShowerUncertaintyY);
   fChain->SetBranchAddress("dShowerUncertaintyZ", &dShowerUncertaintyZ, &b_dPluginBranch_FCALMCComparison_dShowerUncertaintyZ);
   fChain->SetBranchAddress("dShowerUncertaintyT", &dShowerUncertaintyT, &b_dPluginBranch_FCALMCComparison_dShowerUncertaintyT);
   fChain->SetBranchAddress("dShowerUncertaintyE", &dShowerUncertaintyE, &b_dPluginBranch_FCALMCComparison_dShowerUncertaintyE);
   fChain->SetBranchAddress("dPathLengthCorrection", &dPathLengthCorrection, &b_dPluginBranch_FCALMCComparison_dPathLengthCorrection);
}

Bool_t FCALSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef FCALSelector_cxx
