//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep 23 12:57:29 2011 by ROOT version 5.30/00
// from TTree dPluginTree_TOFMCComparison/TOF MC Comparison
// found on file: cluster_run/test.root
//////////////////////////////////////////////////////////

#ifndef TOFSelector_h
#define TOFSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

class TOFSelector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
 //TOFMCComparison *dPluginBranch_TOFMCComparison;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Float_t         dTrueX;
   Float_t         dTrueY;
   Float_t         dTrueZ;
   Float_t         dTruedE;
   Float_t         dTrueT;
   Float_t         dTrueBetaGamma;
   Float_t         dDeltaX;
   Float_t         dDeltaY;
   Float_t         dDeltaZ;
   Float_t         dDeltadE;
   Float_t         dDeltaT;
   Float_t         dPathLengthCorrection;
   Bool_t          dHorizontalPlaneFlag;
   Bool_t          dVerticalPlaneFlag;

   // List of branches
   TBranch        *b_dPluginBranch_TOFMCComparison_fUniqueID;   //!
   TBranch        *b_dPluginBranch_TOFMCComparison_fBits;   //!
   TBranch        *b_dPluginBranch_TOFMCComparison_dTrueX;   //!
   TBranch        *b_dPluginBranch_TOFMCComparison_dTrueY;   //!
   TBranch        *b_dPluginBranch_TOFMCComparison_dTrueZ;   //!
   TBranch        *b_dPluginBranch_TOFMCComparison_dTruedE;   //!
   TBranch        *b_dPluginBranch_TOFMCComparison_dTrueT;   //!
   TBranch        *b_dPluginBranch_TOFMCComparison_dTrueBetaGamma;   //!
   TBranch        *b_dPluginBranch_TOFMCComparison_dDeltaX;   //!
   TBranch        *b_dPluginBranch_TOFMCComparison_dDeltaY;   //!
   TBranch        *b_dPluginBranch_TOFMCComparison_dDeltaZ;   //!
   TBranch        *b_dPluginBranch_TOFMCComparison_dDeltadE;   //!
   TBranch        *b_dPluginBranch_TOFMCComparison_dDeltaT;   //!
   TBranch        *b_dPluginBranch_TOFMCComparison_dPathLengthCorrection;   //!
   TBranch        *b_dPluginBranch_TOFMCComparison_dHorizontalPlaneFlag;   //!
   TBranch        *b_dPluginBranch_TOFMCComparison_dVerticalPlaneFlag;   //!

		TFile *dOutputFile;

		//dEVsBetaGamma
		TH2F *dPluginHist_TOF_dEVsBetaGamma;

		//DeltaX Dependence
		TH2F *dPluginHist_TOF_DeltaXVsX;
		TH2F *dPluginHist_TOF_DeltaXVsY;
		TH2F *dPluginHist_TOF_DeltaXVsZ;
		TH2F *dPluginHist_TOF_DeltaXVsdE;
		TH2F *dPluginHist_TOF_DeltaXVsT;

		//DeltaY Dependence
		TH2F *dPluginHist_TOF_DeltaYVsX;
		TH2F *dPluginHist_TOF_DeltaYVsY;
		TH2F *dPluginHist_TOF_DeltaYVsZ;
		TH2F *dPluginHist_TOF_DeltaYVsdE;
		TH2F *dPluginHist_TOF_DeltaYVsT;

		//DeltaZ Dependence
		TH2F *dPluginHist_TOF_DeltaZVsX;
		TH2F *dPluginHist_TOF_DeltaZVsY;
		TH2F *dPluginHist_TOF_DeltaZVsZ;
		TH2F *dPluginHist_TOF_DeltaZVsdE;
		TH2F *dPluginHist_TOF_DeltaZVsT;

		//DeltadE Dependence
		TH2F *dPluginHist_TOF_DeltadEVsX;
		TH2F *dPluginHist_TOF_DeltadEVsY;
		TH2F *dPluginHist_TOF_DeltadEVsZ;
		TH2F *dPluginHist_TOF_DeltadEVsdE;
		TH2F *dPluginHist_TOF_DeltadEVsT;

		//DeltaT Dependence
		TH2F *dPluginHist_TOF_DeltaTVsX;
		TH2F *dPluginHist_TOF_DeltaTVsY;
		TH2F *dPluginHist_TOF_DeltaTVsZ;
		TH2F *dPluginHist_TOF_DeltaTVsdE;
		TH2F *dPluginHist_TOF_DeltaTVsT;

		//Common Dependence
		TH2F *dPluginHist_TOF_DeltaXVsDeltaY;
		TH2F *dPluginHist_TOF_DeltaXVsDeltaZ;
		TH2F *dPluginHist_TOF_DeltaXVsDeltadE;
		TH2F *dPluginHist_TOF_DeltaXVsDeltaT;
		TH2F *dPluginHist_TOF_DeltaYVsDeltaZ;
		TH2F *dPluginHist_TOF_DeltaYVsDeltadE;
		TH2F *dPluginHist_TOF_DeltaYVsDeltaT;
		TH2F *dPluginHist_TOF_DeltaZVsDeltadE;
		TH2F *dPluginHist_TOF_DeltaZVsDeltaT;
		TH2F *dPluginHist_TOF_DeltadEVsDeltaT;


		//dEVsBetaGamma
		TH2F *dPluginHist_TOF_dEVsBetaGamma_HorizontalOnly;

		//DeltaX Dependence
		TH2F *dPluginHist_TOF_DeltaXVsX_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaXVsY_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaXVsZ_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaXVsdE_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaXVsT_HorizontalOnly;

		//DeltaY Dependence
		TH2F *dPluginHist_TOF_DeltaYVsX_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaYVsY_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaYVsZ_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaYVsdE_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaYVsT_HorizontalOnly;

		//DeltaZ Dependence
		TH2F *dPluginHist_TOF_DeltaZVsX_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaZVsY_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaZVsZ_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaZVsdE_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaZVsT_HorizontalOnly;

		//DeltadE Dependence
		TH2F *dPluginHist_TOF_DeltadEVsX_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltadEVsY_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltadEVsZ_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltadEVsdE_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltadEVsT_HorizontalOnly;

		//DeltaT Dependence
		TH2F *dPluginHist_TOF_DeltaTVsX_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaTVsY_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaTVsZ_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaTVsdE_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaTVsT_HorizontalOnly;

		//Common Dependence
		TH2F *dPluginHist_TOF_DeltaXVsDeltaY_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaXVsDeltaZ_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaXVsDeltadE_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaXVsDeltaT_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaYVsDeltaZ_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaYVsDeltadE_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaYVsDeltaT_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaZVsDeltadE_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltaZVsDeltaT_HorizontalOnly;
		TH2F *dPluginHist_TOF_DeltadEVsDeltaT_HorizontalOnly;


		//dEVsBetaGamma
		TH2F *dPluginHist_TOF_dEVsBetaGamma_VerticalOnly;

		//DeltaX Dependence
		TH2F *dPluginHist_TOF_DeltaXVsX_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaXVsY_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaXVsZ_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaXVsdE_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaXVsT_VerticalOnly;

		//DeltaY Dependence
		TH2F *dPluginHist_TOF_DeltaYVsX_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaYVsY_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaYVsZ_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaYVsdE_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaYVsT_VerticalOnly;

		//DeltaZ Dependence
		TH2F *dPluginHist_TOF_DeltaZVsX_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaZVsY_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaZVsZ_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaZVsdE_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaZVsT_VerticalOnly;

		//DeltadE Dependence
		TH2F *dPluginHist_TOF_DeltadEVsX_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltadEVsY_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltadEVsZ_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltadEVsdE_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltadEVsT_VerticalOnly;

		//DeltaT Dependence
		TH2F *dPluginHist_TOF_DeltaTVsX_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaTVsY_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaTVsZ_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaTVsdE_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaTVsT_VerticalOnly;

		//Common Dependence
		TH2F *dPluginHist_TOF_DeltaXVsDeltaY_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaXVsDeltaZ_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaXVsDeltadE_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaXVsDeltaT_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaYVsDeltaZ_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaYVsDeltadE_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaYVsDeltaT_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaZVsDeltadE_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltaZVsDeltaT_VerticalOnly;
		TH2F *dPluginHist_TOF_DeltadEVsDeltaT_VerticalOnly;

   TOFSelector(TTree * /*tree*/ =0) { }
   virtual ~TOFSelector() { }
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

   ClassDef(TOFSelector,0);
};

#endif

#ifdef TOFSelector_cxx
void TOFSelector::Init(TTree *tree)
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

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_dPluginBranch_TOFMCComparison_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_dPluginBranch_TOFMCComparison_fBits);
   fChain->SetBranchAddress("dTrueX", &dTrueX, &b_dPluginBranch_TOFMCComparison_dTrueX);
   fChain->SetBranchAddress("dTrueY", &dTrueY, &b_dPluginBranch_TOFMCComparison_dTrueY);
   fChain->SetBranchAddress("dTrueZ", &dTrueZ, &b_dPluginBranch_TOFMCComparison_dTrueZ);
   fChain->SetBranchAddress("dTruedE", &dTruedE, &b_dPluginBranch_TOFMCComparison_dTruedE);
   fChain->SetBranchAddress("dTrueT", &dTrueT, &b_dPluginBranch_TOFMCComparison_dTrueT);
   fChain->SetBranchAddress("dTrueBetaGamma", &dTrueBetaGamma, &b_dPluginBranch_TOFMCComparison_dTrueBetaGamma);
   fChain->SetBranchAddress("dDeltaX", &dDeltaX, &b_dPluginBranch_TOFMCComparison_dDeltaX);
   fChain->SetBranchAddress("dDeltaY", &dDeltaY, &b_dPluginBranch_TOFMCComparison_dDeltaY);
   fChain->SetBranchAddress("dDeltaZ", &dDeltaZ, &b_dPluginBranch_TOFMCComparison_dDeltaZ);
   fChain->SetBranchAddress("dDeltadE", &dDeltadE, &b_dPluginBranch_TOFMCComparison_dDeltadE);
   fChain->SetBranchAddress("dDeltaT", &dDeltaT, &b_dPluginBranch_TOFMCComparison_dDeltaT);
   fChain->SetBranchAddress("dPathLengthCorrection", &dPathLengthCorrection, &b_dPluginBranch_TOFMCComparison_dPathLengthCorrection);
   fChain->SetBranchAddress("dHorizontalPlaneFlag", &dHorizontalPlaneFlag, &b_dPluginBranch_TOFMCComparison_dHorizontalPlaneFlag);
   fChain->SetBranchAddress("dVerticalPlaneFlag", &dVerticalPlaneFlag, &b_dPluginBranch_TOFMCComparison_dVerticalPlaneFlag);
}

Bool_t TOFSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef TOFSelector_cxx
