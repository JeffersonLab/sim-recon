#define Read_bcal_hadronic_eff2_cxx
// The class definition in Read_bcal_hadronic_eff2.h has been generated automatically
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
// Root > T->Process("Read_bcal_hadronic_eff2.C")
// Root > T->Process("Read_bcal_hadronic_eff2.C","some options")
// Root > T->Process("Read_bcal_hadronic_eff2.C+")
//
// Example:
//
// Bcal_hadronic_eff>root -l tree_bcal_hadronic_eff_010492_apr4.root
// root [1] bcal_hadronic_eff->Process("Read_bcal_hadronic_eff2.C","010492_apr4 1")    Note: TString option = "010492_apr4 1", encodes "filerun layer". Note filerun should be consistent with input file name.
//
// Note: script depends on the following parameters, set at the start of the script
// layer - layer to be plotted  - now taken from option (optional argument)
// coinc_cut - difference between tracking position and sector allowed for a match
// filerun - run number or designation used for output: output filename = "R"+filerun+"_layer"+TString::Itoa(layer,10)+"_cut"+TString::Itoa(coinc_cut,10); - from option (optional argument)
// Output files created in dat, pdf and root subdirectories

#include "Read_bcal_hadronic_eff2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <vector>
#include <TVector3.h>
#include <iostream>
#include <fstream>


    char estring[256];
    
	// Int_t layer=1 ;
	Int_t layer;    // select layer for efficiencies - get from option
	Int_t coinc_cut = 3;   // nominal is 3
	Int_t nfiles=1;
	// TString filerun = "030890-mar23";
	TString filerun;	// get from option
	TString filename;

    
    Int_t ndx = 0;
    Int_t entry_max=1000000;
    Int_t rand_offset = -15;   // keep offset negative    
    Int_t count_proj = 0;
    Int_t count_down = 0;
    Int_t count_up = 0;



void Read_bcal_hadronic_eff2::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
   
  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
    
    TString option = GetOption();
    // TString filerun = option;
    
    TObjArray *tokens = option.Tokenize(" ");
    Int_t ntokens = tokens->GetEntries();
    if (ntokens != 2) {
    	cout << "*** Begin: Number of entries=" << ntokens << endl;
        exit(1);
    }
    // Float_t  = (((TObjString*)tokens->At(1))->GetString()).Atof();
    // cout << " token1=" << ((TObjString*)tokens->At(0))->GetString()  << " token2=" << ((TObjString*)tokens->At(1))->GetString() << endl;
    filerun = ((TObjString*)tokens->At(0))->GetString();
    layer = ((TObjString*)tokens->At(1))->GetString().Atoi();
    cout << "Begin: filerun=" << filerun << " layer=" << layer << endl;
    
    // Define histograms
    Int_t nbins=100;
    
    h1_TrackDeltaPhiToShower = new TH1F("TrackDeltaPhiToShower", "TrackDeltaPhiToShower",nbins, -10,10);
    h1_TrackDeltaZToShower = new TH1F("TrackDeltaZToShower", "TrackDeltaZToShower", nbins,-20,20);
    h1_ProjectedBCALSectors_Layer = new TH1F("TProjectedBCALSectors_Layer", "ProjectedBCALSectors_Layer", 2*nbins,0,200);
    h1_NearestBCALSectors_Layer_Downstream = new TH1F("NearestBCALSectors_Layer_Downstream", "NearestBCALSectors_Layer_Downstream", 2*nbins,0,200);
    h1_NearestBCALSectors_Layer_Upstream = new TH1F("NearestBCALSectors_Layer_Upstream", "NearestBCALSectors_Layer_Upstream", 2*nbins,0,200);
    h1_pid_pdg = new TH1F("PID_PDG", "PID_PDG", 2*nbins,-3000,3000);
    
    h1_layer_proj = new TH1F("h1_layer_proj","h1_layer_proj",2*nbins,0,200);
    h1_layer_up = new TH1F("h1_layer_up","h1_layer_up",2*nbins,0,200);
    h1_layer_down = new TH1F("h1_layer_down","h1_layer_down",2*nbins,0,200);
    h1_layer_up_rand_left = new TH1F("h1_layer_up_rand_left","h1_layer_up_rand_left",2*nbins,0,200);
    h1_layer_down_rand_left = new TH1F("h1_layer_down_rand_left","h1_layer_down_rand_left",2*nbins,0,200);
    h1_layer_up_rand_right = new TH1F("h1_layer_up_rand_right","h1_layer_up_rand_right",2*nbins,0,200);
    h1_layer_down_rand_right = new TH1F("h1_layer_down_rand_right","h1_layer_down_rand_right",2*nbins,0,200);
    h1_layer_up_rand = new TH1F("h1_layer_up_rand","h1_layer_up_rand",2*nbins,0,200);
    h1_layer_down_rand = new TH1F("h1_layer_down_rand","h1_layer_down_rand",2*nbins,0,200);
    
    h1_eff_up = new TH1F("h1_eff_up","h1_eff_up",2*nbins,0,200);
    h1_eff_down = new TH1F("h1_eff_down","h1_eff_down",2*nbins,0,200);
    h1_eff_up_temp = new TH1F("h1_eff_up_temp","h1_eff_up_temp",2*nbins,0,200);
    h1_eff_down_temp = new TH1F("h1_eff_down_temp","h1_eff_down_temp",2*nbins,0,200);
    
    h1_diff_up = new TH1F("h1_diff_up","h1_diff_up",100,-50,50);
    h1_diff_down = new TH1F("h1_diff_down","h1_diff_down",100,-50,50);
    h1_diff_up_rand_left = new TH1F("h1_diff_up_rand_left","h1_diff_up_rand_left",50,-40,10);
    h1_diff_down_rand_left = new TH1F("h1_diff_down_rand_left","h1_diff_down_rand_left",50,-40,10);
    h1_diff_up_rand_right = new TH1F("h1_diff_up_rand_right","h1_diff_up_rand_right",50,-10,40);
    h1_diff_down_rand_right = new TH1F("h1_diff_down_rand_right","h1_diff_down_rand_right",50,-10,40);

    h1_mom_proj = new TH1F("h1_mom_proj","h1_mom_proj",2*nbins,0,4);
    h1_mom_up = new TH1F("h1_mom_up","h1_mom_up",2*nbins,0,4);
    h1_mom_down = new TH1F("h1_mom_down","h1_mom_down",2*nbins,0,4);
    h1_mom_up_rand_left = new TH1F("h1_mom_up_rand_left","h1_mom_up_rand_left",2*nbins,0,4);
    h1_mom_down_rand_left = new TH1F("h1_mom_down_rand_left","h1_mom_down_rand_left",2*nbins,0,4);
    h1_mom_up_rand_right = new TH1F("h1_mom_up_rand_right","h1_mom_up_rand_right",2*nbins,0,4);
    h1_mom_down_rand_right = new TH1F("h1_mom_down_rand_right","h1_mom_down_rand_right",2*nbins,0,4);
    h1_mom_up_rand = new TH1F("h1_mom_up_rand","h1_mom_up_rand",2*nbins,0,4);
    h1_mom_down_rand = new TH1F("h1_mom_down_rand","h1_mom_down_rand",2*nbins,0,4);

    h1_eff_mom_up = new TH1F("h1_eff_mom_up","h1_eff_mom_up",2*nbins,0,4);
    h1_eff_mom_down = new TH1F("h1_eff_mom_down","h1_eff_mom_down",2*nbins,0,4);
    h1_eff_mom_up_temp = new TH1F("h1_eff_mom_up_temp","h1_eff_mom_up_temp",2*nbins,0,4);
    h1_eff_mom_down_temp = new TH1F("h1_eff_mom__down_temp","h1_eff_mom_down_temp",2*nbins,0,4);
    
    h1_z_proj = new TH1F("h1_z_proj","h1_z_proj",2*nbins,0,500);
    h1_z_up = new TH1F("h1_z_up","h1_z_up",2*nbins,0,500);
    h1_z_down = new TH1F("h1_z_down","h1_z_down",2*nbins,0,500);
    h1_z_up_rand_left = new TH1F("h1_z_up_rand_left","h1_z_up_rand_left",2*nbins,0,500);
    h1_z_down_rand_left = new TH1F("h1_z_down_rand_left","h1_z_down_rand_left",2*nbins,0,500);
    h1_z_up_rand_right = new TH1F("h1_z_up_rand_right","h1_z_up_rand_right",2*nbins,0,500);
    h1_z_down_rand_right = new TH1F("h1_z_down_rand_right","h1_z_down_rand_right",2*nbins,0,500);
    h1_z_up_rand = new TH1F("h1_z_up_rand","h1_z_up_rand",2*nbins,0,500);
    h1_z_down_rand = new TH1F("h1_z_down_rand","h1_z_down_rand",2*nbins,0,500);
    
    h1_eff_z_up = new TH1F("h1_eff_z_up","h1_eff_z_up",2*nbins,0,500);
    h1_eff_z_down = new TH1F("h1_eff_z_down","h1_eff_z_down",2*nbins,0,500);
    h1_eff_z_up_temp = new TH1F("h1_eff_z_up_temp","h1_eff_z_up_temp",2*nbins,0,500);
    h1_eff_z_down_temp = new TH1F("h1_eff_z__down_temp","h1_eff_z_down_temp",2*nbins,0,500);
    
    h1_Evis_proj = new TH1F("h1_Evis_proj","h1_Evis_proj",2*nbins,0,0.5);
    h1_Evis_up = new TH1F("h1_Evis_up","h1_Evis_up",2*nbins,0,0.5);
    h1_Evis_down = new TH1F("h1_Evis_down","h1_Evis_down",2*nbins,0,0.5);
    h1_Evis_up_rand_left = new TH1F("h1_Evis_up_rand_left","h1_Evis_up_rand_left",2*nbins,0,0.5);
    h1_Evis_down_rand_left = new TH1F("h1_Evis_down_rand_left","h1_Evis_down_rand_left",2*nbins,0,0.5);
    h1_Evis_up_rand_right = new TH1F("h1_Evis_up_rand_right","h1_Evis_up_rand_right",2*nbins,0,0.5);
    h1_Evis_down_rand_right = new TH1F("h1_Evis_down_rand_right","h1_Evis_down_rand_right",2*nbins,0,0.5);
    h1_Evis_up_rand = new TH1F("h1_Evis_up_rand","h1_Evis_up_rand",2*nbins,0,0.5);
    h1_Evis_down_rand = new TH1F("h1_Evis_down_rand","h1_Evis_down_rand",2*nbins,0,0.5);
    
    h1_eff_Evis_up = new TH1F("h1_eff_Evis_up","h1_eff_Evis_up",2*nbins,0,0.5);
    h1_eff_Evis_down = new TH1F("h1_eff_Evis_down","h1_eff_Evis_down",2*nbins,0,0.5);
    h1_eff_Evis_up_temp = new TH1F("h1_eff_Evis_up_temp","h1_eff_Evis_up_temp",2*nbins,0,0.5);
    h1_eff_Evis_down_temp = new TH1F("h1_eff_Evis__down_temp","h1_eff_Evis_down_temp",2*nbins,0,0.5);
    
    h1_eff_up->Sumw2();
    h1_layer_up->Sumw2();
    h1_layer_up_rand_left->Sumw2();
    h1_layer_up_rand_right->Sumw2();
    h1_layer_up_rand->Sumw2();
    h1_eff_down->Sumw2();
    h1_layer_down->Sumw2();
    h1_layer_down_rand_left->Sumw2();
    h1_layer_down_rand_right->Sumw2();
    h1_layer_down_rand->Sumw2();

    h1_eff_mom_up->Sumw2();
    h1_mom_up->Sumw2();
    h1_mom_up_rand_left->Sumw2();
    h1_mom_up_rand_right->Sumw2();
    h1_mom_up_rand->Sumw2();
    h1_eff_mom_down->Sumw2();
    h1_mom_down->Sumw2();
    h1_mom_down_rand_left->Sumw2();
    h1_mom_down_rand_right->Sumw2();
    h1_mom_down_rand->Sumw2();
    
    h1_eff_z_up->Sumw2();
    h1_z_up->Sumw2();
    h1_z_up_rand_left->Sumw2();
    h1_z_up_rand_right->Sumw2();
    h1_z_up_rand->Sumw2();
    h1_eff_z_down->Sumw2();
    h1_z_down->Sumw2();
    h1_z_down_rand_left->Sumw2();
    h1_z_down_rand_right->Sumw2();
    h1_z_down_rand->Sumw2();
    
    h1_eff_Evis_up->Sumw2();
    h1_Evis_up->Sumw2();
    h1_Evis_up_rand_left->Sumw2();
    h1_Evis_up_rand_right->Sumw2();
    h1_Evis_up_rand->Sumw2();
    h1_eff_Evis_down->Sumw2();
    h1_Evis_down->Sumw2();
    h1_Evis_down_rand_left->Sumw2();
    h1_Evis_down_rand_right->Sumw2();
    h1_Evis_down_rand->Sumw2();

    
    TString filebase = "/cache/halld/RunPeriod-2016-02/recon/ver01/tree_bcal_hadronic_eff/"+filerun+"/tree_bcal_hadronic_eff_";
    // TString filebase = "tree_bcal_hadronic_eff_";
    TFile *file = NULL;

    // file = new TFile(filebase+filerun+".root","read");
    // cout << "Opening root file: " << (filebase+filerun+".root").Data() << endl;
    filename = "R"+filerun+"_layer"+TString::Itoa(layer,10)+"_cut"+TString::Itoa(coinc_cut,10);

}

void Read_bcal_hadronic_eff2::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t Read_bcal_hadronic_eff2::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either Read_bcal_hadronic_eff2::GetEntry() or TBranch::GetEntry()
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


    Int_t ProjectedBCALSectors_Layer;
    Int_t NearestBCALSectors_Layer_Upstream;
    Int_t NearestBCALSectors_Layer_Downstream;

      
      Int_t proj4;
      Int_t up4;
      Int_t down4;
      Double_t updiff4;
      Double_t downdiff4;

      fChain->GetEntry(entry);

    if (entry >= entry_max) {
        // cout << " Event limit reached =" << entry << endl;
        return kFALSE;
    }
    
    Float_t Evis;
      // select layer of interest
      switch (layer) {
          case 1:
              ProjectedBCALSectors_Layer = (Int_t) ProjectedBCALSectors_Layer1;  
              NearestBCALSectors_Layer_Upstream = (Int_t) NearestBCALSectors_Layer1_Upstream;
              NearestBCALSectors_Layer_Downstream = (Int_t) NearestBCALSectors_Layer1_Downstream;
              Evis = NearestBCALEnergy_Layer1_Downstream;
              break;
          case 2:
              ProjectedBCALSectors_Layer = (Int_t) ProjectedBCALSectors_Layer2;  
              NearestBCALSectors_Layer_Upstream = (Int_t) NearestBCALSectors_Layer2_Upstream;
              NearestBCALSectors_Layer_Downstream = (Int_t) NearestBCALSectors_Layer2_Downstream;
              Evis = NearestBCALEnergy_Layer2_Downstream;
              break;
          case 3:
	       ProjectedBCALSectors_Layer = 0;
	      // check that layer 4 is hit, otherwise set ProjectedBCALSectors_Layer=0, so that this track will be ignored in the layer 3 efficiency
	      proj4 = (Int_t) ProjectedBCALSectors_Layer4; 
	      up4 = (Int_t) NearestBCALSectors_Layer4_Upstream;
	      down4 = (Int_t) NearestBCALSectors_Layer4_Downstream;
	      updiff4 = up4 - (proj4-1)%192 -1;
	      downdiff4 = down4 - (proj4-1)%192 -1;	      
	      if (proj4 &&  up4 && down4 && fabs(updiff4) <= coinc_cut && fabs(downdiff4) <= coinc_cut) { 
              ProjectedBCALSectors_Layer = (Int_t) ProjectedBCALSectors_Layer3;
              NearestBCALSectors_Layer_Upstream = (Int_t) NearestBCALSectors_Layer3_Upstream;
              NearestBCALSectors_Layer_Downstream = (Int_t) NearestBCALSectors_Layer3_Downstream;
              Evis = NearestBCALEnergy_Layer3_Downstream;
		  // cout << endl << " ProjectedBCALSectors_Layer=" << ProjectedBCALSectors_Layer << " NearestBCALSectors_Layer_Upstream=" << NearestBCALSectors_Layer_Upstream << " NearestBCALSectors_Layer_Downstream=" << NearestBCALSectors_Layer_Downstream << endl;
		  // cout << " proj4=" << proj4 << " up4=" << up4 << " down4=" << down4 << endl;
		  }
              break;
          case 4:
              ProjectedBCALSectors_Layer = (Int_t) ProjectedBCALSectors_Layer4; 
              NearestBCALSectors_Layer_Upstream = (Int_t) NearestBCALSectors_Layer4_Upstream;
              NearestBCALSectors_Layer_Downstream = (Int_t) NearestBCALSectors_Layer4_Downstream;
              Evis = NearestBCALEnergy_Layer4_Downstream;
              break;
          default:
              cout << "*** Read_bcal_hadronic_eff: illegal layer=" << layer << endl;
              return kFALSE;
      }
      
      // cout << " ProjectedBCALSectors_Layer=" << ProjectedBCALSectors_Layer << " NearestBCALSectors_Layer_Upstream=" <<  NearestBCALSectors_Layer_Upstream << " NearestBCALSectors_Layer_Downstream=" << NearestBCALSectors_Layer_Downstream << endl;


    Double_t mom = TrackP3->Mag();
    Double_t z = ProjectedBCALHitZ;
    Int_t pid_pdf = PID_PDG;
    
    // cout << "TrackP3->Mag()=" << mom << endl;

    
      // eff histograms
      h1_TrackDeltaPhiToShower->Fill(TrackDeltaPhiToShower);
      h1_TrackDeltaZToShower->Fill(TrackDeltaZToShower);
      h1_ProjectedBCALSectors_Layer->Fill(ProjectedBCALSectors_Layer);
      h1_NearestBCALSectors_Layer_Downstream->Fill(NearestBCALSectors_Layer_Downstream);
      h1_NearestBCALSectors_Layer_Upstream->Fill(NearestBCALSectors_Layer_Upstream);
      h1_pid_pdg->Fill(pid_pdf);
    
      Double_t diff_down = 1000.;
      Double_t diff_up = 1000.;
      Double_t diff_down_rand_left = 1000.;
      Double_t diff_up_rand_left = 1000.;
      Double_t diff_down_rand_right = 1000.;
      Double_t diff_up_rand_right = 1000.;
      Int_t subtr = (ProjectedBCALSectors_Layer-1)%192 +1;
      Int_t subtr_rand = (ProjectedBCALSectors_Layer+rand_offset-1)%192 +1;
      if (subtr_rand <=0) subtr_rand += 192;

      if (ProjectedBCALSectors_Layer && (subtr > 0)) {
          diff_down = NearestBCALSectors_Layer_Downstream - (ProjectedBCALSectors_Layer-1)%192 -1;
          diff_up = NearestBCALSectors_Layer_Upstream - (ProjectedBCALSectors_Layer-1)%192 -1;

	  h1_diff_down->Fill(diff_down);
	  h1_diff_up->Fill(diff_up);
      }

      if (ProjectedBCALSectors_Layer && (subtr_rand > 0)) {
	  diff_down_rand_left = NearestBCALSectors_Layer_Downstream - (ProjectedBCALSectors_Layer+rand_offset-1)%192 -1;
	  diff_up_rand_left = NearestBCALSectors_Layer_Upstream - (ProjectedBCALSectors_Layer+rand_offset-1)%192 -1 ;
	  diff_down_rand_right = NearestBCALSectors_Layer_Downstream - (ProjectedBCALSectors_Layer-rand_offset-1)%192 -1;
	  diff_up_rand_right = NearestBCALSectors_Layer_Upstream - (ProjectedBCALSectors_Layer-rand_offset-1)%192 -1;
      
	  h1_diff_down_rand_left->Fill(diff_down_rand_left);
	  h1_diff_up_rand_left->Fill(diff_up_rand_left);
	  h1_diff_down_rand_right->Fill(diff_down_rand_right);
	  h1_diff_up_rand_right->Fill(diff_up_rand_right);
      }

      
      if (ProjectedBCALSectors_Layer && subtr && subtr_rand) {
          h1_layer_proj->Fill(ProjectedBCALSectors_Layer);
          h1_mom_proj->Fill(mom);
          h1_z_proj->Fill(z);
          h1_Evis_proj->Fill(Evis);
          count_proj++;
          
            if (NearestBCALSectors_Layer_Downstream && ((fabs(diff_down)<=coinc_cut) || fabs(fabs(diff_down)-192)<=coinc_cut) ){
              // cout << " fabs(diff_down)=" << fabs(diff_down) << endl;
                h1_layer_down->Fill(ProjectedBCALSectors_Layer);
                h1_mom_down->Fill(mom);
                h1_z_down->Fill(z);
                h1_Evis_down->Fill(Evis);
                count_down++;
	    }
            if (NearestBCALSectors_Layer_Downstream && (fabs(diff_down_rand_left)<=coinc_cut || fabs(fabs(diff_down_rand_left)-192)<=coinc_cut) ){
	      // if (ProjectedBCALSectors_Layer < 20) cout << " sector=" <<ProjectedBCALSectors_Layer  << " ***fabs(diff_down_rand_left)=" << fabs(diff_down_rand_left) << " coinc_cut=" << coinc_cut << " nearest=" << NearestBCALSectors_Layer_Downstream << endl;
                h1_layer_down_rand_left->Fill(ProjectedBCALSectors_Layer);
                h1_mom_down_rand_left->Fill(mom);
                h1_z_down_rand_left->Fill(z);
                h1_Evis_down_rand_left->Fill(Evis);
            }
            if (NearestBCALSectors_Layer_Downstream && (fabs(diff_down_rand_right)<=coinc_cut || fabs(fabs(diff_down_rand_right)-192)<=coinc_cut) ){
	    // if (ProjectedBCALSectors_Layer < 20) cout  << " sector=" << ProjectedBCALSectors_Layer << " ***fabs(diff_down_rand_right)=" << fabs(diff_down_rand_right) << endl;
                h1_layer_down_rand_right->Fill(ProjectedBCALSectors_Layer);
                h1_mom_down_rand_right->Fill(mom);
                h1_z_down_rand_right->Fill(z);
                h1_Evis_down_rand_right->Fill(Evis);
            }
          
	    if (NearestBCALSectors_Layer_Upstream &&  (fabs(diff_up) <= coinc_cut || fabs(fabs(diff_up)-192)<=coinc_cut) ){
              // cout << " fabs(diff_up)=" << fabs(diff_up) << endl;
            h1_layer_up->Fill(ProjectedBCALSectors_Layer);
            h1_mom_up->Fill(mom);
            h1_z_up->Fill(z);
            h1_Evis_up->Fill(Evis);
            count_up++;
	    }
	    if (NearestBCALSectors_Layer_Upstream &&  (fabs(diff_up_rand_left) <= coinc_cut || fabs(fabs(diff_up_rand_left)-192)<=coinc_cut) ){
	    // if (ProjectedBCALSectors_Layer < 20)  cout  << " sector=" <<ProjectedBCALSectors_Layer  << " ***fabs(diff_up_rand_left)=" << fabs(diff_up_rand_left) << endl;
            h1_layer_up_rand_left->Fill(ProjectedBCALSectors_Layer);
            h1_mom_up_rand_left->Fill(mom);
            h1_z_up_rand_left->Fill(z);
            h1_Evis_up_rand_left->Fill(Evis);
        }
          if (NearestBCALSectors_Layer_Upstream &&  (fabs(diff_up_rand_right) <= coinc_cut || fabs(fabs(diff_up_rand_right)-192)<=coinc_cut) ){
	    // if (ProjectedBCALSectors_Layer < 20)  cout  << " sector=" << ProjectedBCALSectors_Layer  << " ***fabs(diff_up_rand_right)=" << fabs(diff_up_rand_right) << endl;
              h1_layer_up_rand_right->Fill(ProjectedBCALSectors_Layer);
              h1_mom_up_rand_right->Fill(mom);
              h1_z_up_rand_right->Fill(z);
              h1_Evis_up_rand_right->Fill(Evis);
          }
      
      }



   return kTRUE;
}

void Read_bcal_hadronic_eff2::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void Read_bcal_hadronic_eff2::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.


    TCanvas *c0 = new TCanvas("c0", "c0",200,10,1000,700);
    
    c0->Divide(3,2);
    c0->cd(1);
    // gPad->SetLogy();
    Double_t xmin = 0;
    Double_t xmax = 2;
    Double_t ymin = 100;
    Double_t ymax = 10000;
    
    h1_TrackDeltaPhiToShower->SetTitle(filename);
    // h1_TrackDeltaPhiToShower->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_TrackDeltaPhiToShower->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_TrackDeltaPhiToShower->GetXaxis()->SetTitleSize(0.05);
    h1_TrackDeltaPhiToShower->GetYaxis()->SetTitleSize(0.05);
    h1_TrackDeltaPhiToShower->GetXaxis()->SetTitle("Delta Phi to Shower");
    h1_TrackDeltaPhiToShower->SetMarkerColor(4);
    h1_TrackDeltaPhiToShower->Draw();
    
    c0->cd(2);
    
    h1_TrackDeltaZToShower->SetTitle(filename);
    // h1_TrackDeltaZToShower->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_TrackDeltaZToShower->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_TrackDeltaZToShower->GetXaxis()->SetTitleSize(0.05);
    h1_TrackDeltaZToShower->GetYaxis()->SetTitleSize(0.05);
    h1_TrackDeltaZToShower->GetXaxis()->SetTitle("Delta Z to Shower");
    h1_TrackDeltaZToShower->SetMarkerColor(4);
    h1_TrackDeltaZToShower->Draw();
    
    c0->cd(3);
    
    h1_ProjectedBCALSectors_Layer->SetTitle(filename);
    // h1_ProjectedBCALSectors_Layer->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_ProjectedBCALSectors_Layer->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_ProjectedBCALSectors_Layer->GetXaxis()->SetTitleSize(0.05);
    h1_ProjectedBCALSectors_Layer->GetYaxis()->SetTitleSize(0.05);
    h1_ProjectedBCALSectors_Layer->GetXaxis()->SetTitle("ProjectedBCALSectors_Layer");
    h1_ProjectedBCALSectors_Layer->SetMarkerColor(4);
    h1_ProjectedBCALSectors_Layer->Draw();
    
    c0->cd(4);
    
    h1_NearestBCALSectors_Layer_Downstream->SetTitle(filename);
    // h1_NearestBCALSectors_Layer->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_NearestBCALSectors_Layer->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_NearestBCALSectors_Layer_Downstream->GetXaxis()->SetTitleSize(0.05);
    h1_NearestBCALSectors_Layer_Downstream->GetYaxis()->SetTitleSize(0.05);
    h1_NearestBCALSectors_Layer_Downstream->GetXaxis()->SetTitle("NearestBCALSectors_Layer_Downstream");
    h1_NearestBCALSectors_Layer_Downstream->SetMarkerColor(4);
    h1_NearestBCALSectors_Layer_Downstream->Draw();
    
    c0->cd(5);
    
    h1_NearestBCALSectors_Layer_Upstream->SetTitle(filename);
    // h1_NearestBCALSectors_Layer_Upstream->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_NearestBCALSectors_Layer_Upstream->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_NearestBCALSectors_Layer_Upstream->GetXaxis()->SetTitleSize(0.05);
    h1_NearestBCALSectors_Layer_Upstream->GetYaxis()->SetTitleSize(0.05);
    h1_NearestBCALSectors_Layer_Upstream->GetXaxis()->SetTitle("NearestBCALSectors_Layer_Upstream");
    h1_NearestBCALSectors_Layer_Upstream->SetMarkerColor(4);
    h1_NearestBCALSectors_Layer_Upstream->Draw();
    
    c0->cd(6);
    
    h1_pid_pdg->SetTitle(filename);
    // h1_pid_pdg->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_pid_pdg->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_pid_pdg->GetXaxis()->SetTitleSize(0.05);
    h1_pid_pdg->GetYaxis()->SetTitleSize(0.05);
    h1_pid_pdg->GetXaxis()->SetTitle("PID_PDG");
    h1_pid_pdg->SetMarkerColor(4);
    h1_pid_pdg->Draw();
    
    TCanvas *c1 = new TCanvas("c1", "c1",200,10,1000,700);
    
    c1->Divide(3,2);
    c1->cd(1);
    gPad->SetLogy();
    
    h1_diff_up->SetTitle(filename);  
    // h1_diff_up->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_diff_up->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_diff_up->GetXaxis()->SetTitleSize(0.05);
    h1_diff_up->GetYaxis()->SetTitleSize(0.05);
    h1_diff_up->GetXaxis()->SetTitle("diff_up");
    h1_diff_up->SetMarkerColor(4);
    h1_diff_up->Draw();
    
    c1->cd(2);
    gPad->SetLogy();
    
    h1_diff_up_rand_left->SetTitle(filename);
    // h1_diff_up_rand_left->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_diff_up_rand_left->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_diff_up_rand_left->GetXaxis()->SetTitleSize(0.05);
    h1_diff_up_rand_left->GetYaxis()->SetTitleSize(0.05);
    h1_diff_up_rand_left->GetXaxis()->SetTitle("diff_up_rand_left");
    h1_diff_up_rand_left->SetMarkerColor(4);
    h1_diff_up_rand_left->Draw();
    
    c1->cd(3);
    gPad->SetLogy();
    
    h1_diff_up_rand_right->SetTitle(filename);
    // h1_diff_up_rand_right->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_diff_up_rand_right->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_diff_up_rand_right->GetXaxis()->SetTitleSize(0.05);
    h1_diff_up_rand_right->GetYaxis()->SetTitleSize(0.05);
    h1_diff_up_rand_right->GetXaxis()->SetTitle("diff_up_rand_right");
    h1_diff_up_rand_right->SetMarkerColor(4);
    h1_diff_up_rand_right->Draw();
    
    c1->cd(4);
    gPad->SetLogy();
    
    h1_diff_down->SetTitle(filename);
    // h1_diff_down->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_diff_down->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_diff_down->GetXaxis()->SetTitleSize(0.05);
    h1_diff_down->GetYaxis()->SetTitleSize(0.05);
    h1_diff_down->GetXaxis()->SetTitle("diff_down");
    h1_diff_down->SetMarkerColor(4);
    h1_diff_down->Draw();
    
    
    c1->cd(5);
    gPad->SetLogy();
    
    h1_diff_down_rand_left->SetTitle(filename);
    // h1_diff_down_rand_left->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_diff_down_rand_left->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_diff_down_rand_left->GetXaxis()->SetTitleSize(0.05);
    h1_diff_down_rand_left->GetYaxis()->SetTitleSize(0.05);
    h1_diff_down_rand_left->GetXaxis()->SetTitle("diff_down_rand_left");
    h1_diff_down_rand_left->SetMarkerColor(4);
    h1_diff_down_rand_left->Draw();
    
    c1->cd(6);
    gPad->SetLogy();
    
    h1_diff_down_rand_right->SetTitle(filename);
    // h1_diff_down_rand_right->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_diff_down_rand_right->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_diff_down_rand_right->GetXaxis()->SetTitleSize(0.05);
    h1_diff_down_rand_right->GetYaxis()->SetTitleSize(0.05);
    h1_diff_down_rand_right->GetXaxis()->SetTitle("diff_down_rand_right");
    h1_diff_down_rand_right->SetMarkerColor(4);
    h1_diff_down_rand_right->Draw();
    
    TCanvas *c2 = new TCanvas("c2", "c2",200,10,1000,700);
    
    c2->Divide(3,2);
    
    c2->cd(1);
    
    h1_layer_proj->SetTitle(filename);
    // h1_layer_proj->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_layer_proj->GetYaxis()->SetRangeUser(ymin,ymax);
    // h1_layer_proj->GetXaxis()->SetTitleSize(0.05);
    h1_layer_proj->GetYaxis()->SetTitleSize(0.05);
    h1_layer_proj->GetXaxis()->SetTitle("layer_proj");
    h1_layer_proj->SetMarkerColor(4);
    h1_layer_proj->Draw();

    c2->cd(2);
    // gPad->SetLogy();
    
    h1_layer_up_rand->Add(h1_layer_up_rand_left,h1_layer_up_rand_right,0.5,0.5);
    Double_t rand_sum_up = h1_layer_up_rand->Integral();
    Double_t rand_sum_up_left = h1_layer_up_rand_left->Integral();
    Double_t rand_sum_up_right = h1_layer_up_rand_right->Integral();
    cout << " ave=" << rand_sum_up << " left=" << rand_sum_up_left << " right=" << rand_sum_up_right << endl;
    
    ymin=0.1;
    ymin=100000;
    h1_layer_up->SetTitle(filename);
    // h1_layer_up->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_layer_up->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_layer_up->GetXaxis()->SetTitleSize(0.05);
    h1_layer_up->GetYaxis()->SetTitleSize(0.05);
    h1_layer_up->GetXaxis()->SetTitle("layer_up");
    h1_layer_up->SetMarkerColor(4);
    h1_layer_up->Draw();
    h1_layer_up_rand->SetLineColor(2);
    h1_layer_up_rand->Draw("same");
    
    c2->cd(3);
    // gPad->SetLogy();
    
    h1_layer_down_rand->Add(h1_layer_down_rand_left,h1_layer_down_rand_right,0.5,0.5);
    Double_t rand_sum_down = h1_layer_down_rand->Integral();
    Double_t rand_sum_down_left = h1_layer_down_rand_left->Integral();
    Double_t rand_sum_down_right = h1_layer_down_rand_right->Integral();
    cout << " ave=" << rand_sum_down << " left=" << rand_sum_down_left << " right=" << rand_sum_down_right << endl;
    
    h1_layer_down->SetTitle(filename);
    // h1_layer_down->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_layer_down->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_layer_down->GetXaxis()->SetTitleSize(0.05);
    h1_layer_down->GetYaxis()->SetTitleSize(0.05);
    h1_layer_down->GetXaxis()->SetTitle("layer_down");
    h1_layer_down->SetMarkerColor(4);
    h1_layer_down->Draw();
    h1_layer_down_rand->SetLineColor(2);
    h1_layer_down_rand->Draw("same");
    
    c2->cd(4);
    
    h1_eff_up_temp->Add(h1_layer_up,h1_layer_up_rand,1,-1);
    Double_t up_temp_diff =h1_eff_up_temp->Integral();
    h1_eff_up->Divide(h1_eff_up_temp,h1_layer_proj,1,1,"B");
    
    h1_eff_up->SetTitle(filename);
    // h1_eff_up->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_eff_up->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_eff_up->GetXaxis()->SetTitleSize(0.05);
    h1_eff_up->GetYaxis()->SetTitleSize(0.05);
    h1_eff_up->GetXaxis()->SetTitle("eff_up");
    h1_eff_up->SetMarkerColor(4);
    h1_eff_up->Draw();
    
    c2->cd(5);

    h1_eff_down_temp->Add(h1_layer_down,h1_layer_down_rand,1,-1);
    Double_t down_temp_diff =h1_eff_up_temp->Integral();
    h1_eff_down->Divide(h1_eff_down_temp,h1_layer_proj,1,1,"B");
    
    h1_eff_down->SetTitle(filename);
    // h1_eff_down->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_eff_down->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_eff_down->GetXaxis()->SetTitleSize(0.05);
    h1_eff_down->GetYaxis()->SetTitleSize(0.05);
    h1_eff_down->GetXaxis()->SetTitle("eff_down");
    h1_eff_down->SetMarkerColor(4);
    h1_eff_down->Draw();
    
    c2->cd(6);
    
    sprintf (estring,"Average Efficiency\n");
    printf("estring=%s",estring);
    TLatex *t1 = new TLatex(0.2,0.95,estring);    // t1->SetNDC();
    t1->SetTextSize(0.04);
    t1->Draw();
    
    Double_t count_up_diff = count_up - rand_sum_up;
    Double_t count_down_diff = count_down - rand_sum_down;
    cout << " ave=" << rand_sum_up << " left=" << rand_sum_up_left << " right=" << rand_sum_up_right << endl;
    cout << " ave2=" << rand_sum_down << " left2=" << rand_sum_down_left << " right2=" << rand_sum_down_right << endl;
    cout << " up_temp_diff=" << up_temp_diff << " down_temp_diff=" << down_temp_diff << endl;
    
    cout << endl << endl << endl << "count_proj=" << count_proj << " count_up=" << count_up << " count_down=" << count_down << " count_down_diff=" << count_down_diff << " count_up_diff=" << count_up_diff << " ratio_down=" <<
    count_down_diff/count_proj << " ratio_up=" << count_up_diff/count_proj << endl << endl;
    
    
    sprintf (estring,"Layer= %d\n",layer);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.90,estring);
    sprintf (estring,"Proj-Hit <= %d\n",coinc_cut);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.85,estring);
    sprintf (estring,"Projected counts= %d\n",count_proj);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.75,estring);
    sprintf (estring,"Upstream counts= %.0f\n",count_up_diff);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.70,estring);
    sprintf (estring,"Downstream counts= %.0f\n",count_down_diff);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.65,estring);
    sprintf (estring,"Upstream efficiency= %.3f\n",count_up_diff/count_proj);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.6,estring);
    sprintf (estring,"Downstream efficiency= %.3f\n",count_down_diff/count_proj);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.55,estring);
    
    
    TCanvas *c3 = new TCanvas("c3", "c3",200,10,1000,700);
    
    c3->Divide(3,2);
    
    c3->cd(1);
    
    h1_mom_proj->SetTitle(filename);
    // h1_mom_proj->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_mom_proj->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_mom_proj->GetXaxis()->SetTitleSize(0.05);
    h1_mom_proj->GetYaxis()->SetTitleSize(0.05);
    h1_mom_proj->GetXaxis()->SetTitle("mom_proj (GeV)");
    h1_mom_proj->SetMarkerColor(4);
    h1_mom_proj->Draw();
    
    c3->cd(2);
    gPad->SetLogy();
    
    h1_mom_up_rand->Add(h1_mom_up_rand_left,h1_mom_up_rand_right,0.5,0.5);
    Double_t mom_rand_sum_up = h1_mom_up_rand->Integral();
    Double_t mom_rand_sum_up_left = h1_mom_up_rand_left->Integral();
    Double_t mom_rand_sum_up_right = h1_mom_up_rand_right->Integral();
    cout << " ave=" << mom_rand_sum_up << " left=" << mom_rand_sum_up_left << " right=" << mom_rand_sum_up_right << endl;
    
    
    h1_mom_up->SetTitle(filename);
    // h1_mom_up->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_mom_up->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_mom_up->GetXaxis()->SetTitleSize(0.05);
    h1_mom_up->GetYaxis()->SetTitleSize(0.05);
    h1_mom_up->GetXaxis()->SetTitle("mom_up (GeV)");
    h1_mom_up->SetMarkerColor(4);
    h1_mom_up->Draw();
    h1_mom_up_rand->SetLineColor(2);
    h1_mom_up_rand->Draw("same");
    
    c3->cd(3);
    gPad->SetLogy();
    
    h1_mom_down_rand->Add(h1_mom_down_rand_left,h1_mom_down_rand_right,0.5,0.5);
    Double_t mom_rand_sum_down = h1_mom_down_rand->Integral();
    Double_t mom_rand_sum_down_left = h1_mom_down_rand_left->Integral();
    Double_t mom_rand_sum_down_right = h1_mom_down_rand_right->Integral();
    cout << " ave=" << mom_rand_sum_down << " left=" << mom_rand_sum_down_left << " right=" << mom_rand_sum_down_right << endl;
    
    h1_mom_down->SetTitle(filename);
    // h1_mom_down->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_mom_down->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_mom_down->GetXaxis()->SetTitleSize(0.05);
    h1_mom_down->GetYaxis()->SetTitleSize(0.05);
    h1_mom_down->GetXaxis()->SetTitle("mom_down (GeV)");
    h1_mom_down->SetMarkerColor(4);
    h1_mom_down->Draw();
    h1_mom_down_rand->SetLineColor(2);
    h1_mom_down_rand->Draw("same");
    
    c3->cd(4);
    
    h1_eff_mom_up_temp->Add(h1_mom_up,h1_mom_up_rand,1,-1);
    h1_eff_mom_up->Divide(h1_eff_mom_up_temp,h1_mom_proj,1,1,"B");
    
    h1_eff_mom_up->SetTitle(filename);
    // h1_eff_mom_up->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_eff_mom_up->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_eff_mom_up->GetXaxis()->SetTitleSize(0.05);
    h1_eff_mom_up->GetYaxis()->SetTitleSize(0.05);
    h1_eff_mom_up->GetXaxis()->SetTitle("eff_mom_up (GeV)");
    h1_eff_mom_up->SetMarkerColor(4);
    h1_eff_mom_up->Draw();
    
    c3->cd(5);
    
    h1_eff_mom_down_temp->Add(h1_mom_down,h1_mom_down_rand,1,-1);
    h1_eff_mom_down->Divide(h1_eff_mom_down_temp,h1_mom_proj,1,1,"B");
    
    h1_eff_mom_down->SetTitle(filename);
    // h1_eff_mom_down->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_eff_mom_down->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_eff_mom_down->GetXaxis()->SetTitleSize(0.05);
    h1_eff_mom_down->GetYaxis()->SetTitleSize(0.05);
    h1_eff_mom_down->GetXaxis()->SetTitle("eff_mom_down (GeV)");
    h1_eff_mom_down->SetMarkerColor(4);
    h1_eff_mom_down->Draw();
    
    c3->cd(6);
    
    sprintf (estring,"Average Efficiency\n");
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.95,estring);    // t1->SetNDC();
    t1->SetTextSize(0.04);
    t1->Draw();
    
    Double_t mom_count_up_diff = count_up - mom_rand_sum_up;
    Double_t mom_count_down_diff = count_down - mom_rand_sum_down;
    cout << endl << endl << endl << "count_proj=" << count_proj << " count_up=" << count_up << " count_down=" << count_down << " mom_count_down_diff=" << mom_count_down_diff << " mom_count_up_diff=" << mom_count_up_diff << " ratio_down=" <<
    mom_count_down_diff/count_proj << " ratio_up=" << mom_count_up_diff/count_proj << endl << endl;
    
    
    sprintf (estring,"Layer= %d\n",layer);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.90,estring);
    sprintf (estring,"Proj-Hit <= %d\n",coinc_cut);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.85,estring);
    sprintf (estring,"Projected counts= %d\n",count_proj);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.75,estring);
    sprintf (estring,"Upstream counts= %.0f\n",mom_count_up_diff);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.70,estring);
    sprintf (estring,"Downstream counts= %.0f\n",mom_count_down_diff);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.65,estring);
    sprintf (estring,"Upstream efficiency= %.3f\n",mom_count_up_diff/count_proj);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.6,estring);
    sprintf (estring,"Downstream efficiency= %.3f\n",mom_count_down_diff/count_proj);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.55,estring);
    
    
    TCanvas *c4 = new TCanvas("c4", "c4",200,10,1000,700);
    
    c4->Divide(3,2);
    
    c4->cd(1);
    
    h1_z_proj->SetTitle(filename);
    // h1_z_proj->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_z_proj->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_z_proj->GetXaxis()->SetTitleSize(0.05);
    h1_z_proj->GetYaxis()->SetTitleSize(0.05);
    h1_z_proj->GetXaxis()->SetTitle("z_proj (cm)");
    h1_z_proj->SetMarkerColor(4);
    h1_z_proj->Draw();
    
    c4->cd(2);
    gPad->SetLogy();
    
    h1_z_up_rand->Add(h1_z_up_rand_left,h1_z_up_rand_right,0.5,0.5);
    Double_t z_rand_sum_up = h1_z_up_rand->Integral();
    Double_t z_rand_sum_up_left = h1_z_up_rand_left->Integral();
    Double_t z_rand_sum_up_right = h1_z_up_rand_right->Integral();
    cout << " ave=" << z_rand_sum_up << " left=" << z_rand_sum_up_left << " right=" << z_rand_sum_up_right << endl;
    
    
    h1_z_up->SetTitle(filename);
    // h1_z_up->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_z_up->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_z_up->GetXaxis()->SetTitleSize(0.05);
    h1_z_up->GetYaxis()->SetTitleSize(0.05);
    h1_z_up->GetXaxis()->SetTitle("z_up (cm)");
    h1_z_up->SetMarkerColor(4);
    h1_z_up->Draw();
    h1_z_up_rand->SetLineColor(2);
    h1_z_up_rand->Draw("same");
    
    c4->cd(3);
    gPad->SetLogy();
    
    h1_z_down_rand->Add(h1_z_down_rand_left,h1_z_down_rand_right,0.5,0.5);
    Double_t z_rand_sum_down = h1_z_down_rand->Integral();
    Double_t z_rand_sum_down_left = h1_z_down_rand_left->Integral();
    Double_t z_rand_sum_down_right = h1_z_down_rand_right->Integral();
    cout << " ave=" << z_rand_sum_down << " left=" << z_rand_sum_down_left << " right=" << z_rand_sum_down_right << endl;
    
    h1_z_down->SetTitle(filename);
    // h1_z_down->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_z_down->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_z_down->GetXaxis()->SetTitleSize(0.05);
    h1_z_down->GetYaxis()->SetTitleSize(0.05);
    h1_z_down->GetXaxis()->SetTitle("z_down (cm)");
    h1_z_down->SetMarkerColor(4);
    h1_z_down->Draw();
    h1_z_down_rand->SetLineColor(2);
    h1_z_down_rand->Draw("same");
    
    c4->cd(4);
    
    h1_eff_z_up_temp->Add(h1_z_up,h1_z_up_rand,1,-1);
    h1_eff_z_up->Divide(h1_eff_z_up_temp,h1_z_proj,1,1,"B");
    
    h1_eff_z_up->SetTitle(filename);
    // h1_eff_z_up->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_eff_z_up->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_eff_z_up->GetXaxis()->SetTitleSize(0.05);
    h1_eff_z_up->GetYaxis()->SetTitleSize(0.05);
    h1_eff_z_up->GetXaxis()->SetTitle("eff_z_up (cm)");
    h1_eff_z_up->SetMarkerColor(4);
    h1_eff_z_up->Draw();
    
    c4->cd(5);
    
    h1_eff_z_down_temp->Add(h1_z_down,h1_z_down_rand,1,-1);
    h1_eff_z_down->Divide(h1_eff_z_down_temp,h1_z_proj,1,1,"B");
    
    h1_eff_z_down->SetTitle(filename);
    // h1_eff_z_down->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_eff_z_down->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_eff_z_down->GetXaxis()->SetTitleSize(0.05);
    h1_eff_z_down->GetYaxis()->SetTitleSize(0.05);
    h1_eff_z_down->GetXaxis()->SetTitle("eff_z_down (cm)");
    h1_eff_z_down->SetMarkerColor(4);
    h1_eff_z_down->Draw();
    
    c4->cd(6);
    
    sprintf (estring,"Average Efficiency\n");
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.95,estring);    // t1->SetNDC();
    t1->SetTextSize(0.04);
    t1->Draw();
    
    Double_t z_count_up_diff = count_up - z_rand_sum_up;
    Double_t z_count_down_diff = count_down - z_rand_sum_down;
    cout << endl << endl << endl << "count_proj=" << count_proj << " count_up=" << count_up << " count_down=" << count_down << " z_count_down_diff=" << z_count_down_diff << " z_count_up_diff=" << z_count_up_diff << " ratio_down=" <<
    z_count_down_diff/count_proj << " ratio_up=" << z_count_up_diff/count_proj << endl << endl;
    
    
    sprintf (estring,"Layer= %d\n",layer);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.90,estring);
    sprintf (estring,"Proj-Hit <= %d\n",coinc_cut);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.85,estring);
    sprintf (estring,"Projected counts= %d\n",count_proj);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.75,estring);
    sprintf (estring,"Upstream counts= %.0f\n",z_count_up_diff);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.70,estring);
    sprintf (estring,"Downstream counts= %.0f\n",z_count_down_diff);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.65,estring);
    sprintf (estring,"Upstream efficiency= %.3f\n",z_count_up_diff/count_proj);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.6,estring);
    sprintf (estring,"Downstream efficiency= %.3f\n",z_count_down_diff/count_proj);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.55,estring);
    
    
    TCanvas *c5 = new TCanvas("c5", "c5",200,10,1000,700);
    
    c5->Divide(3,2);
    
    c5->cd(1);
    
    h1_Evis_proj->SetTitle(filename);
    // h1_Evis_proj->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_Evis_proj->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_Evis_proj->GetXaxis()->SetTitleSize(0.05);
    h1_Evis_proj->GetYaxis()->SetTitleSize(0.05);
    h1_Evis_proj->GetXaxis()->SetTitle("Evis (GeV)");
    h1_Evis_proj->SetMarkerColor(4);
    h1_Evis_proj->Draw();
    
    c5->cd(2);
    gPad->SetLogy();
    
    h1_Evis_up_rand->Add(h1_Evis_up_rand_left,h1_Evis_up_rand_right,0.5,0.5);
    Double_t Evis_rand_sum_up = h1_Evis_up_rand->Integral();
    Double_t Evis_rand_sum_up_left = h1_Evis_up_rand_left->Integral();
    Double_t Evis_rand_sum_up_right = h1_Evis_up_rand_right->Integral();
    cout << " ave=" << Evis_rand_sum_up << " left=" << Evis_rand_sum_up_left << " right=" << Evis_rand_sum_up_right << endl;
    
    
    h1_Evis_up->SetTitle(filename);
    // h1_Evis_up->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_Evis_up->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_Evis_up->GetXaxis()->SetTitleSize(0.05);
    h1_Evis_up->GetYaxis()->SetTitleSize(0.05);
    h1_Evis_up->GetXaxis()->SetTitle("Evis_up (GeV)");
    h1_Evis_up->SetMarkerColor(4);
    h1_Evis_up->Draw();
    h1_Evis_up_rand->SetLineColor(2);
    h1_Evis_up_rand->Draw("same");
    
    c5->cd(3);
    gPad->SetLogy();
    
    h1_Evis_down_rand->Add(h1_Evis_down_rand_left,h1_Evis_down_rand_right,0.5,0.5);
    Double_t Evis_rand_sum_down = h1_Evis_down_rand->Integral();
    Double_t Evis_rand_sum_down_left = h1_Evis_down_rand_left->Integral();
    Double_t Evis_rand_sum_down_right = h1_Evis_down_rand_right->Integral();
    cout << " ave=" << Evis_rand_sum_down << " left=" << Evis_rand_sum_down_left << " right=" << Evis_rand_sum_down_right << endl;
    
    h1_Evis_down->SetTitle(filename);
    // h1_Evis_down->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_Evis_down->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_Evis_down->GetXaxis()->SetTitleSize(0.05);
    h1_Evis_down->GetYaxis()->SetTitleSize(0.05);
    h1_Evis_down->GetXaxis()->SetTitle("Evis_down (GeV)");
    h1_Evis_down->SetMarkerColor(4);
    h1_Evis_down->Draw();
    h1_Evis_down_rand->SetLineColor(2);
    h1_Evis_down_rand->Draw("same");
    
    c5->cd(4);
    
    h1_eff_Evis_up_temp->Add(h1_Evis_up,h1_Evis_up_rand,1,-1);
    h1_eff_Evis_up->Divide(h1_eff_Evis_up_temp,h1_Evis_proj,1,1,"B");
    
    h1_eff_Evis_up->SetTitle(filename);
    // h1_eff_Evis_up->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_eff_Evis_up->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_eff_Evis_up->GetXaxis()->SetTitleSize(0.05);
    h1_eff_Evis_up->GetYaxis()->SetTitleSize(0.05);
    h1_eff_Evis_up->GetXaxis()->SetTitle("eff_Evis_up (GeV)");
    h1_eff_Evis_up->SetMarkerColor(4);
    h1_eff_Evis_up->Draw();
    
    c5->cd(5);
    
    h1_eff_Evis_down_temp->Add(h1_Evis_down,h1_Evis_down_rand,1,-1);
    h1_eff_Evis_down->Divide(h1_eff_Evis_down_temp,h1_Evis_proj,1,1,"B");
    
    h1_eff_Evis_down->SetTitle(filename);
    // h1_eff_Evis_down->GetXaxis()->SetRangeUser(xmin,xmax);
    // h1_eff_Evis_down->GetYaxis()->SetRangeUser(ymin,ymax);
    h1_eff_Evis_down->GetXaxis()->SetTitleSize(0.05);
    h1_eff_Evis_down->GetYaxis()->SetTitleSize(0.05);
    h1_eff_Evis_down->GetXaxis()->SetTitle("eff_Evis_down (GeV)");
    h1_eff_Evis_down->SetMarkerColor(4);
    h1_eff_Evis_down->Draw();
    
    c5->cd(6);
    
    sprintf (estring,"Average Efficiency\n");
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.95,estring);    // t1->SetNDC();
    t1->SetTextSize(0.04);
    t1->Draw();
    
    Double_t Evis_count_up_diff = count_up - Evis_rand_sum_up;
    Double_t Evis_count_down_diff = count_down - Evis_rand_sum_down;
    cout << endl << endl << endl << "count_proj=" << count_proj << " count_up=" << count_up << " count_down=" << count_down << " Evis_count_down_diff=" << Evis_count_down_diff << " Evis_count_up_diff=" << Evis_count_up_diff << " ratio_down=" <<
    Evis_count_down_diff/count_proj << " ratio_up=" << Evis_count_up_diff/count_proj << endl << endl;
    
    
    sprintf (estring,"Layer= %d\n",layer);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.90,estring);
    sprintf (estring,"Proj-Hit <= %d\n",coinc_cut);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.85,estring);
    sprintf (estring,"Projected counts= %d\n",count_proj);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.75,estring);
    sprintf (estring,"Upstream counts= %.0f\n",Evis_count_up_diff);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.70,estring);
    sprintf (estring,"Downstream counts= %.0f\n",Evis_count_down_diff);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.65,estring);
    sprintf (estring,"Upstream efficiency= %.3f\n",Evis_count_up_diff/count_proj);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.6,estring);
    sprintf (estring,"Downstream efficiency= %.3f\n",Evis_count_down_diff/count_proj);
    printf("estring=%s",estring);
    t1->DrawLatex(0.2,0.55,estring);
    

    // open output histogram 

   TString outhist_name = "root/"+filename+"_out.root";
   TFile *outhist = new TFile(outhist_name,"recreate");

   h1_eff_up->Write();
    h1_eff_down->Write();
    h1_eff_mom_up->Write();
    h1_eff_mom_down->Write();
    h1_eff_z_up->Write();
    h1_eff_z_down->Write();
    h1_eff_Evis_up->Write();
    h1_eff_Evis_down->Write();
    outhist->Close();

    // open file for output
    
   TString outfile = "dat/"+filename+".dat";
   cout << "Opening file: " << outfile.Data() << endl;

   ofstream filedat;
   filedat.open (outfile.Data(), ios::out);

   filedat << " run= "  << filerun << " layer= "  << layer << " coinc_cut= " << coinc_cut <<  endl;
    filedat << " count_proj= " << count_proj << " count_up_diff= " << count_up_diff << " count_down_diff= " << count_down_diff << " eff_up= "
           << count_up_diff/count_proj << " eff_down= " << count_down_diff/count_proj << endl;
    filedat << " count_proj= " << count_proj << " mom_count_up_diff= " << mom_count_up_diff << " mom_count_down_diff= " << mom_count_down_diff << " mom_eff_up= "
    << mom_count_up_diff/count_proj << " eff_down= " << mom_count_down_diff/count_proj << endl;
    filedat << " count_proj= " << count_proj << " z_count_up_diff= " << z_count_up_diff << " z_count_down_diff= " << z_count_down_diff << " z_eff_up= "
    << z_count_up_diff/count_proj << " eff_down= " << z_count_down_diff/count_proj << endl;
    filedat << " count_proj= " << count_proj << " Evis_count_up_diff= " << Evis_count_up_diff << " Evis_count_down_diff= " << Evis_count_down_diff << " Evis_eff_up= "
    << Evis_count_up_diff/count_proj << " eff_down= " << Evis_count_down_diff/count_proj << endl;
    filedat.close();
    
    c0->SaveAs("pdf/"+filename+".pdf(");
    c1->SaveAs("pdf/"+filename+".pdf");
    c2->SaveAs("pdf/"+filename+".pdf");
    c3->SaveAs("pdf/"+filename+".pdf");
    c4->SaveAs("pdf/"+filename+".pdf");
    c5->SaveAs("pdf/"+filename+".pdf)");


}
