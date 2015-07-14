//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Oct  2 15:04:24 2011 by ROOT version 5.30/00
// from TTree dPluginTree_MCReconstructionStatuses/MC Reconstruction Statuses
// found on file: cluster_run/hd_root_pidstudies_proton.root
//////////////////////////////////////////////////////////

#ifndef PIDStudySelector_h
#define PIDStudySelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
//#include <../../../../../libraries/include/DLorentzVector.h>
#include <include/DLorentzVector.h>
#include <TH2F.h>
#include <TH1F.h>

class PIDStudySelector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
 //MCReconstructionStatuses *dPluginBranch_MCReconstructionStatuses;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   vector<MCReconstructionStatus*> dMCReconstructionStatusVector;

	TFile *dOutputFile;

   // List of branches
   TBranch        *b_dPluginBranch_MCReconstructionStatuses_fUniqueID;   //!
   TBranch        *b_dPluginBranch_MCReconstructionStatuses_fBits;   //!
   TBranch        *b_dPluginBranch_MCReconstructionStatuses_dMCReconstructionStatusVector;   //!

   PIDStudySelector(TTree * /*tree*/ =0) { }
   virtual ~PIDStudySelector() { }
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

   ClassDef(PIDStudySelector,0);
};

#endif

#ifdef PIDStudySelector_cxx
void PIDStudySelector::Init(TTree *tree)
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

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_dPluginBranch_MCReconstructionStatuses_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_dPluginBranch_MCReconstructionStatuses_fBits);
   fChain->SetBranchAddress("dMCReconstructionStatusVector", &dMCReconstructionStatusVector, &b_dPluginBranch_MCReconstructionStatuses_dMCReconstructionStatusVector);
}

Bool_t PIDStudySelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef PIDStudySelector_cxx
