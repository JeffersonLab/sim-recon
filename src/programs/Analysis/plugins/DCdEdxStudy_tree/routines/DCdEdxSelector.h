//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 30 12:10:49 2011 by ROOT version 5.30/00
// from TTree dPluginTree_DCdEdxInformation/DC dEdx Information
// found on file: cluster_run/hd_root_dcdedx_proton_fom.root
//////////////////////////////////////////////////////////

#ifndef DCdEdxSelector_h
#define DCdEdxSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

class DCdEdxSelector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
 //DCdEdxInformation *dPluginBranch_DCdEdxInformation;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Double_t        dBeta;
   Double_t        dMomentum;
   Double_t        dTheta;
   Double_t        dVertexZ;
   Double_t        ddEdx_FDC;
   Double_t        ddx_FDC;
   UInt_t          dNumHitsUsedFordEdx_FDC;
   Double_t        ddEdx_CDC;
   Double_t        ddx_CDC;
   UInt_t          dNumHitsUsedFordEdx_CDC;
   Double_t        dChiSq_DCdEdx;
   UInt_t          dNDF_DCdEdx;
   Double_t        dFOM;

	bool dCalcFOMManuallyFlag;
	TObjArray *dSigmaFuncArray_FDC_Proton;
	TObjArray *dSigmaFuncArray_CDC_Proton;
	TObjArray *dSigmaFuncArray_FDC_KPlus;
	TObjArray *dSigmaFuncArray_CDC_KPlus;
	TObjArray *dSigmaFuncArray_FDC_PiPlus;
	TObjArray *dSigmaFuncArray_CDC_PiPlus;
	vector<int> ddEdxSigmaNumHitsVector_FDC_Proton;
	vector<int> ddEdxSigmaNumHitsVector_CDC_Proton;
	vector<int> ddEdxSigmaNumHitsVector_FDC_PiPlus;
	vector<int> ddEdxSigmaNumHitsVector_CDC_PiPlus;
	vector<int> ddEdxSigmaNumHitsVector_FDC_KPlus;
	vector<int> ddEdxSigmaNumHitsVector_CDC_KPlus;

	bool GetdEdxMean_CDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locMeandEdx);
	bool GetdEdxMean_FDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locMeandEdx);
	bool GetdEdxSigma_FDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locSigmadEdx);
	bool GetdEdxSigma_CDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locSigmadEdx);
	bool Calc_FOM(double& locFOM);

	TF1* ddEdxMeanFunc_FDC_Proton;
	TF1* ddEdxMeanFunc_FDC_PiPlus;
	TF1* ddEdxMeanFunc_CDC_Proton;
	TF1* ddEdxMeanFunc_CDC_PiPlus;
	TF1* ddEdxMeanFunc_FDC_KPlus;
	TF1* ddEdxMeanFunc_CDC_KPlus;

	double dRhoZoverA_CDC;
	double dRhoZoverA_FDC;

	string dParticleName;

	TFile* dOutputFile;

	TH2F *dSelectorHist_dEdxVsBeta_HitsCutoff_FDC;
	TH2F *dSelectorHist_dEdxVsP_HitsCutoff_FDC;
	TH2F *dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC;
	TH2F *dSelectorHist_dEdxVsBeta_HitsCutoff_CDC;
	TH2F *dSelectorHist_dEdxVsP_HitsCutoff_CDC;
	TH2F *dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC;

	TH2F *dSelectorHist_NumHitsVsTheta_CDC;
	TH2F *dSelectorHist_dxVsNumHits_CDC;
	TH2F *dSelectorHist_dxVsTheta_CDC;
	TH2F *dSelectorHist_NumHitsVsTheta_FDC;
	TH2F *dSelectorHist_dxVsNumHits_FDC;
	TH2F *dSelectorHist_dxVsTheta_FDC;
	TH2F *dSelectorHist_NumHitsFDCVsNumHitsCDC;

	TH1F* dSelectorHist_NotEnoughHitsInTheta;
	TH1F* dSelectorHist_ThetaDistribution;
	TH1F* dSelectorHist_PercentageNotEnoughHitsInTheta;

	TH2F* dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_3Hits;
	TH2F* dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_6Hits;
	TH2F* dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_9Hits;
	TH2F* dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_12Hits;

	TH2F* dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_2Hits;
	TH2F* dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_4Hits;
	TH2F* dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_6Hits;
	TH2F* dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_8Hits;
	TH2F* dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_10Hits;
	TH2F* dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_12Hits;
	TH2F* dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_14Hits;

	TH1F* dSelectorHist_ConfidenceLevel_FDC;
	TH1F* dSelectorHist_ConfidenceLevel_CDC;
	TH1F* dSelectorHist_ConfidenceLevel_Both;

	TH1F* dSelectorHist_ConfidenceLevel_FDC_3Hits;
	TH1F* dSelectorHist_ConfidenceLevel_FDC_6Hits;
	TH1F* dSelectorHist_ConfidenceLevel_FDC_9Hits;
	TH1F* dSelectorHist_ConfidenceLevel_FDC_12Hits;

	TH1F* dSelectorHist_ConfidenceLevel_CDC_4Hits;
	TH1F* dSelectorHist_ConfidenceLevel_CDC_6Hits;
	TH1F* dSelectorHist_ConfidenceLevel_CDC_8Hits;
	TH1F* dSelectorHist_ConfidenceLevel_CDC_10Hits;
	TH1F* dSelectorHist_ConfidenceLevel_CDC_12Hits;
	TH1F* dSelectorHist_ConfidenceLevel_CDC_14Hits;

	TH1F* dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin1;
	TH1F* dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin2;
	TH1F* dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin3;
	TH1F* dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin4;
	TH1F* dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin5;
	TH1F* dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin6;
	TH1F* dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin7;
	TH1F* dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin8;
	TH1F* dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin9;
	TH1F* dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin10;
	TH1F* dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin11;
	TH1F* dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin12;

	TH1F* dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin1;
	TH1F* dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin2;
	TH1F* dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin3;
	TH1F* dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin4;
	TH1F* dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin5;
	TH1F* dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin6;
	TH1F* dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin7;
	TH1F* dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin8;
	TH1F* dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin9;
	TH1F* dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin10;
	TH1F* dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin11;
	TH1F* dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin12;

   // List of branches
   TBranch        *b_dPluginBranch_DCdEdxInformation_fUniqueID;   //!
   TBranch        *b_dPluginBranch_DCdEdxInformation_fBits;   //!
   TBranch        *b_dPluginBranch_DCdEdxInformation_dBeta;   //!
   TBranch        *b_dPluginBranch_DCdEdxInformation_dMomentum;   //!
   TBranch        *b_dPluginBranch_DCdEdxInformation_dTheta;   //!
   TBranch        *b_dPluginBranch_DCdEdxInformation_dVertexZ;   //!
   TBranch        *b_dPluginBranch_DCdEdxInformation_ddEdx_FDC;   //!
   TBranch        *b_dPluginBranch_DCdEdxInformation_ddx_FDC;   //!
   TBranch        *b_dPluginBranch_DCdEdxInformation_dNumHitsUsedFordEdx_FDC;   //!
   TBranch        *b_dPluginBranch_DCdEdxInformation_ddEdx_CDC;   //!
   TBranch        *b_dPluginBranch_DCdEdxInformation_ddx_CDC;   //!
   TBranch        *b_dPluginBranch_DCdEdxInformation_dNumHitsUsedFordEdx_CDC;   //!
   TBranch        *b_dPluginBranch_DCdEdxInformation_dChiSq_DCdEdx;   //!
   TBranch        *b_dPluginBranch_DCdEdxInformation_dNDF_DCdEdx;   //!
   TBranch        *b_dPluginBranch_DCdEdxInformation_dFOM;   //!

   DCdEdxSelector(TTree * /*tree*/ =0) { }
   virtual ~DCdEdxSelector() { }
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

   ClassDef(DCdEdxSelector,0);
};

#endif

#ifdef DCdEdxSelector_cxx
void DCdEdxSelector::Init(TTree *tree)
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

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_dPluginBranch_DCdEdxInformation_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_dPluginBranch_DCdEdxInformation_fBits);
   fChain->SetBranchAddress("dBeta", &dBeta, &b_dPluginBranch_DCdEdxInformation_dBeta);
   fChain->SetBranchAddress("dMomentum", &dMomentum, &b_dPluginBranch_DCdEdxInformation_dMomentum);
   fChain->SetBranchAddress("dTheta", &dTheta, &b_dPluginBranch_DCdEdxInformation_dTheta);
   fChain->SetBranchAddress("dVertexZ", &dVertexZ, &b_dPluginBranch_DCdEdxInformation_dVertexZ);
   fChain->SetBranchAddress("ddEdx_FDC", &ddEdx_FDC, &b_dPluginBranch_DCdEdxInformation_ddEdx_FDC);
   fChain->SetBranchAddress("ddx_FDC", &ddx_FDC, &b_dPluginBranch_DCdEdxInformation_ddx_FDC);
   fChain->SetBranchAddress("dNumHitsUsedFordEdx_FDC", &dNumHitsUsedFordEdx_FDC, &b_dPluginBranch_DCdEdxInformation_dNumHitsUsedFordEdx_FDC);
   fChain->SetBranchAddress("ddEdx_CDC", &ddEdx_CDC, &b_dPluginBranch_DCdEdxInformation_ddEdx_CDC);
   fChain->SetBranchAddress("ddx_CDC", &ddx_CDC, &b_dPluginBranch_DCdEdxInformation_ddx_CDC);
   fChain->SetBranchAddress("dNumHitsUsedFordEdx_CDC", &dNumHitsUsedFordEdx_CDC, &b_dPluginBranch_DCdEdxInformation_dNumHitsUsedFordEdx_CDC);
   fChain->SetBranchAddress("dChiSq_DCdEdx", &dChiSq_DCdEdx, &b_dPluginBranch_DCdEdxInformation_dChiSq_DCdEdx);
   fChain->SetBranchAddress("dNDF_DCdEdx", &dNDF_DCdEdx, &b_dPluginBranch_DCdEdxInformation_dNDF_DCdEdx);
   fChain->SetBranchAddress("dFOM", &dFOM, &b_dPluginBranch_DCdEdxInformation_dFOM);
}

Bool_t DCdEdxSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef DCdEdxSelector_cxx

