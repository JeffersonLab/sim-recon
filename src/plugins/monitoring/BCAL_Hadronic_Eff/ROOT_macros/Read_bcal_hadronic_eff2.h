//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 15 08:54:12 2016 by ROOT version 5.34/34
// from TTree bcal_hadronic_eff/bcal_hadronic_eff
// found on file: /cache/halld/RunPeriod-2016-02/recon/ver99/tree_bcal_hadronic_eff/011529/tree_bcal_hadronic_eff_011529_000.root
//////////////////////////////////////////////////////////

#ifndef Read_bcal_hadronic_eff2_h
#define Read_bcal_hadronic_eff2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.

// include histogram declarations

    
    TH1F *h1_TrackDeltaPhiToShower;
    TH1F *h1_TrackDeltaZToShower;
    TH1F *h1_ProjectedBCALSectors_Layer;
    TH1F *h1_NearestBCALSectors_Layer_Downstream;
    TH1F *h1_NearestBCALSectors_Layer_Upstream;
	TH1F *h1_pid_pdg;

    TH1F *h1_layer_proj;
    TH1F *h1_layer_up;
    TH1F *h1_layer_down;
    TH1F *h1_layer_up_rand_left;
    TH1F *h1_layer_down_rand_left;
    TH1F *h1_layer_up_rand_right;
    TH1F *h1_layer_down_rand_right;
    TH1F *h1_layer_up_rand;
    TH1F *h1_layer_down_rand;
    
    TH1F *h1_eff_up;
    TH1F *h1_eff_down;
	TH1F *h1_eff_up_temp;
	TH1F *h1_eff_down_temp;

    TH1F *h1_diff_up;
    TH1F *h1_diff_down;
    TH1F *h1_diff_up_rand_left;
    TH1F *h1_diff_down_rand_left;
    TH1F *h1_diff_up_rand_right;
    TH1F *h1_diff_down_rand_right;

    TH1F *h1_mom_proj;
    TH1F *h1_mom_up;
    TH1F *h1_mom_down;
    TH1F *h1_mom_up_rand_left;
    TH1F *h1_mom_down_rand_left;
    TH1F *h1_mom_up_rand_right;
    TH1F *h1_mom_down_rand_right;
    TH1F *h1_mom_up_rand;
    TH1F *h1_mom_down_rand;

    TH1F *h1_eff_mom_up;
	TH1F *h1_eff_mom_down;
	TH1F *h1_eff_mom_up_temp;
	TH1F *h1_eff_mom_down_temp;

	TH1F *h1_z_proj;
	TH1F *h1_z_up;
	TH1F *h1_z_down;
	TH1F *h1_z_up_rand_left;
	TH1F *h1_z_down_rand_left;
	TH1F *h1_z_up_rand_right;
	TH1F *h1_z_down_rand_right;
	TH1F *h1_z_up_rand;
	TH1F *h1_z_down_rand;

	TH1F *h1_eff_z_up;
	TH1F *h1_eff_z_down;
	TH1F *h1_eff_z_up_temp;
	TH1F *h1_eff_z_down_temp;

	TH1F *h1_Evis_proj;
	TH1F *h1_Evis_up;
	TH1F *h1_Evis_down;
	TH1F *h1_Evis_up_rand_left;
	TH1F *h1_Evis_down_rand_left;
	TH1F *h1_Evis_up_rand_right;
	TH1F *h1_Evis_down_rand_right;
	TH1F *h1_Evis_up_rand;
	TH1F *h1_Evis_down_rand;

	TH1F *h1_eff_Evis_up;
	TH1F *h1_eff_Evis_down;
	TH1F *h1_eff_Evis_up_temp;
	TH1F *h1_eff_Evis_down_temp;

#include <TVector3.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

class Read_bcal_hadronic_eff2 : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           PID_PDG;
   Float_t         TrackVertexZ;
   TVector3        *TrackP3;
   UInt_t          TrackCDCRings;
   UInt_t          TrackFDCPlanes;
   Float_t         ShowerEnergy;
   Float_t         TrackDeltaPhiToShower;
   Float_t         TrackDeltaZToShower;
   UChar_t         TrackProjectedBCALSector;
   Float_t         ProjectedBCALHitPhi;
   Float_t         ProjectedBCALHitZ;
   UChar_t         BCALClusterLayers;
   UChar_t         IsHitInCluster;
   UChar_t         ProjectedBCALSectors_Layer1;
   UChar_t         NearestBCALSectors_Layer1_Downstream;
   UChar_t         NearestBCALSectors_Layer1_Upstream;
   Float_t         NearestBCALEnergy_Layer1_Downstream;
   Float_t         NearestBCALEnergy_Layer1_Upstream;
   UChar_t         ProjectedBCALSectors_Layer2;
   UChar_t         NearestBCALSectors_Layer2_Downstream;
   UChar_t         NearestBCALSectors_Layer2_Upstream;
   Float_t         NearestBCALEnergy_Layer2_Downstream;
   Float_t         NearestBCALEnergy_Layer2_Upstream;
   UChar_t         ProjectedBCALSectors_Layer3;
   UChar_t         NearestBCALSectors_Layer3_Downstream;
   UChar_t         NearestBCALSectors_Layer3_Upstream;
   Float_t         NearestBCALEnergy_Layer3_Downstream;
   Float_t         NearestBCALEnergy_Layer3_Upstream;
   UChar_t         ProjectedBCALSectors_Layer4;
   UChar_t         NearestBCALSectors_Layer4_Downstream;
   UChar_t         NearestBCALSectors_Layer4_Upstream;
   Float_t         NearestBCALEnergy_Layer4_Downstream;
   Float_t         NearestBCALEnergy_Layer4_Upstream;
   // List of branches
   TBranch        *b_PID_PDG;   //!
   TBranch        *b_TrackVertexZ;   //!
   TBranch        *b_TrackP3;   //!
   TBranch        *b_TrackCDCRings;   //!
   TBranch        *b_TrackFDCPlanes;   //!
   TBranch        *b_ShowerEnergy;   //!
   TBranch        *b_TrackDeltaPhiToShower;   //!
   TBranch        *b_TrackDeltaZToShower;   //!
   TBranch        *b_TrackProjectedBCALSector;   //!
   TBranch        *b_ProjectedBCALHitPhi;   //!
   TBranch        *b_ProjectedBCALHitZ;   //!
   TBranch        *b_BCALClusterLayers;   //!
   TBranch        *b_IsHitInCluster;   //!
   TBranch        *b_ProjectedBCALSectors_Layer1;   //!
   TBranch        *b_NearestBCALSectors_Layer1_Downstream;   //!
   TBranch        *b_NearestBCALSectors_Layer1_Upstream;   //!
   TBranch        *b_NearestBCALEnergy_Layer1_Downstream;   //!
   TBranch        *b_NearestBCALEnergy_Layer1_Upstream;   //!
   TBranch        *b_ProjectedBCALSectors_Layer2;   //!
   TBranch        *b_NearestBCALSectors_Layer2_Downstream;   //!
   TBranch        *b_NearestBCALSectors_Layer2_Upstream;   //!
   TBranch        *b_NearestBCALEnergy_Layer2_Downstream;   //!
   TBranch        *b_NearestBCALEnergy_Layer2_Upstream;   //!
   TBranch        *b_ProjectedBCALSectors_Layer3;   //!
   TBranch        *b_NearestBCALSectors_Layer3_Downstream;   //!
   TBranch        *b_NearestBCALSectors_Layer3_Upstream;   //!
   TBranch        *b_NearestBCALEnergy_Layer3_Downstream;   //!
   TBranch        *b_NearestBCALEnergy_Layer3_Upstream;   //!
   TBranch        *b_ProjectedBCALSectors_Layer4;   //!
   TBranch        *b_NearestBCALSectors_Layer4_Downstream;   //!
   TBranch        *b_NearestBCALSectors_Layer4_Upstream;   //!
   TBranch        *b_NearestBCALEnergy_Layer4_Downstream;   //!
   TBranch        *b_NearestBCALEnergy_Layer4_Upstream;   //!
    
   Read_bcal_hadronic_eff2(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~Read_bcal_hadronic_eff2() { }
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

   ClassDef(Read_bcal_hadronic_eff2,0);
};

#endif

#ifdef Read_bcal_hadronic_eff2_cxx
void Read_bcal_hadronic_eff2::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   TrackP3 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("PID_PDG", &PID_PDG, &b_PID_PDG);
   fChain->SetBranchAddress("TrackVertexZ", &TrackVertexZ, &b_TrackVertexZ);
   fChain->SetBranchAddress("TrackP3", &TrackP3, &b_TrackP3);
   fChain->SetBranchAddress("TrackCDCRings", &TrackCDCRings, &b_TrackCDCRings);
   fChain->SetBranchAddress("TrackFDCPlanes", &TrackFDCPlanes, &b_TrackFDCPlanes);
   fChain->SetBranchAddress("ShowerEnergy", &ShowerEnergy, &b_ShowerEnergy);
   fChain->SetBranchAddress("TrackDeltaPhiToShower", &TrackDeltaPhiToShower, &b_TrackDeltaPhiToShower);
   fChain->SetBranchAddress("TrackDeltaZToShower", &TrackDeltaZToShower, &b_TrackDeltaZToShower);
   fChain->SetBranchAddress("TrackProjectedBCALSector", &TrackProjectedBCALSector, &b_TrackProjectedBCALSector);
   fChain->SetBranchAddress("ProjectedBCALHitPhi", &ProjectedBCALHitPhi, &b_ProjectedBCALHitPhi);
   fChain->SetBranchAddress("ProjectedBCALHitZ", &ProjectedBCALHitZ, &b_ProjectedBCALHitZ);
   fChain->SetBranchAddress("BCALClusterLayers", &BCALClusterLayers, &b_BCALClusterLayers);
   fChain->SetBranchAddress("IsHitInCluster", &IsHitInCluster, &b_IsHitInCluster);
   fChain->SetBranchAddress("ProjectedBCALSectors_Layer1", &ProjectedBCALSectors_Layer1, &b_ProjectedBCALSectors_Layer1);
   fChain->SetBranchAddress("NearestBCALSectors_Layer1_Downstream", &NearestBCALSectors_Layer1_Downstream, &b_NearestBCALSectors_Layer1_Downstream);
   fChain->SetBranchAddress("NearestBCALSectors_Layer1_Upstream", &NearestBCALSectors_Layer1_Upstream, &b_NearestBCALSectors_Layer1_Upstream);
   fChain->SetBranchAddress("NearestBCALEnergy_Layer1_Downstream", &NearestBCALEnergy_Layer1_Downstream, &b_NearestBCALEnergy_Layer1_Downstream);
   fChain->SetBranchAddress("NearestBCALEnergy_Layer1_Upstream", &NearestBCALEnergy_Layer1_Upstream, &b_NearestBCALEnergy_Layer1_Upstream);
   fChain->SetBranchAddress("ProjectedBCALSectors_Layer2", &ProjectedBCALSectors_Layer2, &b_ProjectedBCALSectors_Layer2);
   fChain->SetBranchAddress("NearestBCALSectors_Layer2_Downstream", &NearestBCALSectors_Layer2_Downstream, &b_NearestBCALSectors_Layer2_Downstream);
   fChain->SetBranchAddress("NearestBCALSectors_Layer2_Upstream", &NearestBCALSectors_Layer2_Upstream, &b_NearestBCALSectors_Layer2_Upstream);
   fChain->SetBranchAddress("NearestBCALEnergy_Layer2_Downstream", &NearestBCALEnergy_Layer2_Downstream, &b_NearestBCALEnergy_Layer2_Downstream);
   fChain->SetBranchAddress("NearestBCALEnergy_Layer2_Upstream", &NearestBCALEnergy_Layer2_Upstream, &b_NearestBCALEnergy_Layer2_Upstream);
   fChain->SetBranchAddress("ProjectedBCALSectors_Layer3", &ProjectedBCALSectors_Layer3, &b_ProjectedBCALSectors_Layer3);
   fChain->SetBranchAddress("NearestBCALSectors_Layer3_Downstream", &NearestBCALSectors_Layer3_Downstream, &b_NearestBCALSectors_Layer3_Downstream);
   fChain->SetBranchAddress("NearestBCALSectors_Layer3_Upstream", &NearestBCALSectors_Layer3_Upstream, &b_NearestBCALSectors_Layer3_Upstream);
   fChain->SetBranchAddress("NearestBCALEnergy_Layer3_Downstream", &NearestBCALEnergy_Layer3_Downstream, &b_NearestBCALEnergy_Layer3_Downstream);
   fChain->SetBranchAddress("NearestBCALEnergy_Layer3_Upstream", &NearestBCALEnergy_Layer3_Upstream, &b_NearestBCALEnergy_Layer3_Upstream);
   fChain->SetBranchAddress("ProjectedBCALSectors_Layer4", &ProjectedBCALSectors_Layer4, &b_ProjectedBCALSectors_Layer4);
   fChain->SetBranchAddress("NearestBCALSectors_Layer4_Downstream", &NearestBCALSectors_Layer4_Downstream, &b_NearestBCALSectors_Layer4_Downstream);
   fChain->SetBranchAddress("NearestBCALSectors_Layer4_Upstream", &NearestBCALSectors_Layer4_Upstream, &b_NearestBCALSectors_Layer4_Upstream);
   fChain->SetBranchAddress("NearestBCALEnergy_Layer4_Downstream", &NearestBCALEnergy_Layer4_Downstream, &b_NearestBCALEnergy_Layer4_Downstream);
   fChain->SetBranchAddress("NearestBCALEnergy_Layer4_Upstream", &NearestBCALEnergy_Layer4_Upstream, &b_NearestBCALEnergy_Layer4_Upstream);
}

Bool_t Read_bcal_hadronic_eff2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef Read_bcal_hadronic_eff2_cxx
