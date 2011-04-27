//
// fill histograms from trackanal.root tree
// the histograms are created in setupHits.C according to
// the number of Thrown tracks of the first event expecting
// all events to be the same Thrown topology.
//
// include files */

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

// root include files
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h" 
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"

#include <string.h>


TH1F *hist[99];   
TH2F *hist2d[10];
void setupHists(TH1F** , TH2F** , Int_t , Int_t *, Int_t *);
void moreHists(TH1F** , TH2F** , Int_t , Int_t *, Int_t *, Float_t *);

int main(int argc, char **argv) {

  Int_t DEBUG = 0;
  char InFile[128] = "trackanal.root";
  TFile *INF = NULL;
  INF = new TFile(InFile);
  INF->ls();
  
  TTree *TrackTree = (TTree*)INF->Get("TrackTree");
  
  Int_t EventNum;
  Int_t NTrThrown;
  Int_t ThrownPType[200]; // particle type of thrown tracks
  Float_t ThrownPp[200];  // particle momentum of thrown tracks
  Float_t ThrownQ[200];   // electric charge of thrown particle
  Int_t NTrCand;
  Float_t TrCandP[200];   // momentum of track candidate
  Float_t TrCandQ[200];   // charge of track candidate
  Float_t TrCandN[200];   // number of hits(Ndof) of track candidate
  Int_t NTrCandHits;
  Int_t NTrFit;
  Int_t NTrHits;
  Int_t MaxT, MaxC, MaxF;
  Int_t trlistPtype[250]; // particle type of track candidate with best FOM
  Int_t trlistcand[250]; // track candidate number
  Float_t trlistPp[250];  // particle momentum of track candidate with best FOM
  Float_t trlistFOM[250];  //track fit FOM
  Float_t trlistchisq[250];  //chisq
  Float_t nh[900] ;        // for each track candidate the chamber hits of each thrown particle track     
  Float_t ptypes[900];     // for each track candidate the chamer hits for each particle type
  
  TrackTree->SetBranchAddress("EventNum", &EventNum);
  TrackTree->SetBranchAddress("MaxT", &MaxT);
  TrackTree->SetBranchAddress("MaxC", &MaxC);
  TrackTree->SetBranchAddress("MaxF", &MaxF);

  TrackTree->SetBranchAddress("NTrThrown", &NTrThrown);
  TrackTree->SetBranchAddress("ThrownPType", ThrownPType);
  TrackTree->SetBranchAddress("ThrownPp", ThrownPp);
  TrackTree->SetBranchAddress("ThrownQ", ThrownQ);

  TrackTree->SetBranchAddress("NTrCand", &NTrCand);
  TrackTree->SetBranchAddress("TrCandQ", TrCandQ);
  TrackTree->SetBranchAddress("TrCandP", TrCandP);
  TrackTree->SetBranchAddress("TrCandN", TrCandN);

  TrackTree->SetBranchAddress("NTrFit", &NTrFit);
  TrackTree->SetBranchAddress("trlistPtype", trlistPtype);
  TrackTree->SetBranchAddress("trlistFOM", trlistFOM);
  TrackTree->SetBranchAddress("trlistchisq", trlistchisq);
  TrackTree->SetBranchAddress("trlistPp", trlistPp);
  TrackTree->SetBranchAddress("trlistcand", trlistcand);

  TrackTree->SetBranchAddress("NTrCandHits", &NTrCandHits);
  TrackTree->SetBranchAddress("nh", nh);
  TrackTree->SetBranchAddress("ptypes", ptypes);

  Int_t Offset[5] = {0, 0, 0, 0, 0};
  for (int j=0;j<MaxF;j++){
    for (int i=0;i<MaxF;i++){
      ptypes[j*MaxF+i] = 0.0 ;
      nh[j*MaxF+i] = 0.0;
    }
  }

  Long64_t nentries = TrackTree->GetEntries();
  cout<<"Number of Entries: "<<nentries<<endl;


  Float_t PMax[50];
  for (Int_t j = 0;j<50;j++){
    PMax[j] = 0;
  }
  
  
  for (Long64_t i=0; i<nentries;i++) {
    //for (Long64_t i=0; i<5;i++) {
    
    //    TrackTree->GetEntry(i++);
    TrackTree->GetEntry(i);
    
    // take the first event to know how many tracks are thrown
    // and generate the histograms accordingly.cd ../../
    if (i==0) { 
      setupHists(hist, hist2d, NTrThrown, Offset, ThrownPType);
    }

    for (Int_t k = 0;k<NTrThrown;k++){
      if (ThrownPp[k]>PMax[k])
	PMax[k] = ThrownPp[k];
    }
    

    hist[0]->Fill((Float_t)NTrThrown);
    hist[2]->Fill((Float_t)NTrCand);  
    
    Int_t SZ = MaxF;
    Int_t idx = 0;
    Int_t MostHits[MaxC];
    Int_t TotHits[MaxC];
    memset(MostHits,0,4*MaxC);
    memset(TotHits,0,4*MaxC);
    
    Int_t MultCount[NTrThrown];
    Float_t NQ = 0;
    for (Int_t j=0;j<NTrThrown;j++){
      MultCount[j] = 0;
      if (ThrownQ[j]!=0){
	NQ += 1.;
      }
    }
    hist[1]->Fill(NQ);


    for (Int_t j=0;j<NTrCand;j++){
      Int_t MaxHits = 0;
      Int_t idx = -1;
      Int_t SumHits = 0;
      Int_t MaxHits1 = 0;
      Int_t ptype = 0;

      for (Int_t k=0;k<SZ;k++){
	//cout<<j<<" and "<<k<<endl;
	SumHits += nh[j*SZ+k];
	if (((Int_t)nh[j*SZ+k])>MaxHits){ // original thrown track with most hits
	  MaxHits = nh[j*SZ+k];
	  idx = k;
	}
	if (((Int_t)ptypes[j*SZ+k])>MaxHits1){ // particle type with most hits
	  MaxHits1 = ptypes[j*SZ+k];
	  ptype = k;
	}
      }
      if (idx>-1){
	MostHits[j] = nh[j*SZ+idx];
      }

      TotHits[j]  = SumHits;      
      hist[3]->Fill((Float_t)SumHits);
      if (idx>-1){
	hist[4+idx]->Fill((Float_t)SumHits);
      }

      if (idx>0){
	Float_t rat = (Float_t)MostHits[j]/(Float_t) TotHits[j]*100.;
	hist[Offset[2]+idx-1]->Fill(rat);
	MultCount[idx-1]++;
      }
     
      // loop over track fits that match the candidate
      Float_t chi2=9e99;
      Int_t index  = -1;
      Int_t index1 = -1;
      Int_t index2 = -1;
      Float_t FOM = 0;
      for (Int_t k=0;k<NTrFit;k++){
	if (trlistcand[k]==j+1){
	  if (chi2>trlistchisq[k]){
	    chi2 = trlistchisq[k];
	    index1 = k;
	  }
	  if (FOM<trlistFOM[k]){
	    FOM = trlistFOM[k];
	    index2 = k;
	  }
	}
      }

      if (index2>-1){
	index = index2;
      } else {
	index = index1;
      }
      
      if (index<0) // all fits for this candidate failed.
	continue;
      
      // for the thrown tracks it is always track 1=pi+ track2=pi- and track3=proton
      Int_t test1 = 0;
      for (Int_t kk=0;kk<NTrThrown;kk++){
	test1 += ((idx==kk+1) && (ptype==ThrownPType[kk]));
      }

      if ( test1 && idx>0) {	
	Int_t trptype = trlistPtype[index];
	Float_t prat = ((Float_t)trptype)/((Float_t)ptype);
	Float_t ptrack = trlistPp[index];
	Float_t ptrue = ThrownPp[idx-1];
	hist[Offset[0]+idx-1]->Fill(ptrack/ptrue);
	hist[Offset[1]+idx-1]->Fill(prat);
	hist2d[0+idx-1]->Fill(ptrack/ptrue, (Float_t)nh[j*SZ+idx]);
	hist2d[Offset[5]+idx]->Fill(prat, (Float_t)nh[j*SZ+idx]);
      } else if (idx==0) {  //most hits for this track candidate come from secondary
	Int_t trptype = trlistPtype[index];
	Float_t prat = ((Float_t)trptype)/((Float_t)ptype);
 	hist[Offset[2]-1]->Fill(prat);
	hist2d[Offset[5]]->Fill(prat, (Float_t)nh[j*SZ+idx]);
      }
      
      if (DEBUG){
	cout<<j<<"  "<<idx<<"  ";
	if (idx>-1) {
	  cout<<nh[j*SZ+idx]/SumHits*100.<<"  "<<SumHits<<endl;
	} else {
	  cout<<"  ?/SumHits           "<<SumHits<<endl;
	}
      }
    }
    
    if (DEBUG==2){
      for (Int_t j=0;j<NTrFit;j++){    
	cout<< i<<"  "<< j << "  " << trlistcand[j]<< "  " <<trlistFOM[j] << "  " <<
	  trlistchisq[j]<<"   "<<trlistPp[j]<<"  "<<trlistPtype[j]<<endl;
      }
    }
    
    //if (i==20)
    //  break;
    
    for (Int_t k=0; k<NTrThrown; k++){
      hist[Offset[3]+k]->Fill((Float_t)MultCount[k]);
    }
    
  }
  
  //cout<<"End looping over data"<<endl;
  // now loop once more the make nice momentum distributions
  // but first create the histograms.
  
  moreHists(hist, hist2d, NTrThrown, 
	    Offset, ThrownPType, PMax);

  for (Long64_t i=0; i<nentries;i++) {
    //for (Long64_t i=0; i<5;i++) {
    
    //    TrackTree->GetEntry(i++);
    TrackTree->GetEntry(i);

    for (Int_t k=0;k<NTrThrown;k++){
      hist[Offset[4]+k]->Fill(ThrownPp[k]);
    }

  }

  

  char fnam[128];
  sprintf(fnam,"histograms_trackanal.root");
  TFile *fout = new TFile(fnam,"RECREATE");

  TList * hist_list= new TList();
  Int_t Nhist = Offset[4] + NTrThrown;
  cout<<"Number of 1D histos: "<<Nhist<<endl;
  for (Int_t j=0;j<Nhist;j++){
    hist_list->Add(hist[j]);
  }

  TList * hist2d_list= new TList();
  Int_t Nhist2d = Offset[5] + NTrThrown;
  cout<<"Number of 2D histos: "<<Nhist2d<<endl;
  for (Int_t j=0;j<Nhist2d;j++){
    hist2d_list->Add(hist2d[j]);
  }

  hist_list->Write("hist_list", TObject::kSingleKey);
  hist2d_list->Write("hist2d_list", TObject::kSingleKey);
  fout->Close();
 

  return 0;
}
// -------------- END OF MAIN -------------------------------

//
// ------------------------------------------------------------
//

void setupHists(TH1F** hist, TH2F** hist2d, Int_t NTrThrown, 
		Int_t *Offset, Int_t *ThrownPType){

  Int_t k = 0;

  hist[0] = new TH1F("hist0","All Tracks Thrown", 16, -0.5, 15.5);
  hist[1] = new TH1F("hist1","Charged Tracks Thrown", 16, -0.5, 15.5);
  hist[2] = new TH1F("hist2","Track Candidates", 16, -0.5, 15.5);
  hist[3] = new TH1F("hist3","Tracking Hits of Candidates ", 41, -0.5, 40.5);
  hist[0]->GetXaxis()->SetTitle("Number of track per event");
  hist[1]->GetXaxis()->SetTitle("Number of track per event");
  hist[2]->GetXaxis()->SetTitle("Number of track per event");
  hist[3]->GetXaxis()->SetTitle("Number of hits for track candidates");

  for (Int_t i=0;i<NTrThrown+1;i++) {
    char str1[128];
    char str2[128];
    sprintf(str1,"hist%d",i+4);
    sprintf(str2,"Thrown Track %d, Ptype=%d : Cand Hits ",i,ThrownPType[i]);
    if (i==0){
      sprintf(str2,"Track X Cand Hits ");
    }
    hist[i+4] = new TH1F(str1, str2, 41, -0.5, 40.5);
    hist[i+4]->GetXaxis()->SetTitle("Number of hits per track");
  }

  Offset[0] = NTrThrown + 4 + 1;
  for (Int_t i=0;i<NTrThrown;i++) {
    char str1[128];
    char str2[128];
    sprintf(str1,"hist%d",i+Offset[0]);
    sprintf(str2,"Track %d Candidate p_det/p_true",i+1);
    hist[i+Offset[0]] = new TH1F(str1, str2,160, -2., 2.);
    hist[i+Offset[0]]->GetXaxis()->SetTitle("p_det/p_true");
    k = i;
  }

  Offset[1] = Offset[0] + NTrThrown;
  
  for (Int_t i=0;i<NTrThrown;i++) {
    char str1[128];
    char str2[128];
    sprintf(str1,"hist%d",i+Offset[1]);
    sprintf(str2,"PID_Tracking/True_PID Thrown Track %d",i+1);
    hist[i+Offset[1]] = new TH1F(str1, str2, 50, -0.5, 2.5); 
    hist[i+Offset[1]]->GetXaxis()->SetTitle("pid_det/pid_gen");
    k = i;
  } 
  char str1[128];
  sprintf(str1,"hist%d",k+1+Offset[1]);
  hist[k+1+Offset[1]] = new TH1F(str1, "PID_Tracking/True_PID secondary Track X", 50, -0.5, 2.5); 
  hist[k+1+Offset[1]]->GetXaxis()->SetTitle("pid_det/pid_gen");
  
  Offset[2] = Offset[1] + NTrThrown + 1;

  for (Int_t i=0;i<NTrThrown;i++) {
    char str1[128];
    char str2[128];
    sprintf(str1,"hist%d",i+Offset[2]);
    sprintf(str2,"Thrown Track %d PType=%d : NHits/TotHits",i+1,ThrownPType[i]);
    hist[i+Offset[2]] = new TH1F(str1, str2, 20, 40., 100.1);
    hist[i+Offset[2]]->GetXaxis()->SetTitle("% of hits from thrown particle");
    k = i;
  }

  Offset[3] = Offset[2] + NTrThrown;

  for (Int_t i=0;i<NTrThrown;i++) {
    char str1[128];
    char str2[128];
    sprintf(str1,"hist%d",i+Offset[3]);
    sprintf(str2,"Thrown Track %d PType=%d : Candidates",i+1,ThrownPType[i] );
    hist[i+Offset[3]] = new TH1F(str1, str2, 6, -0.5, 5.5);
    hist[i+Offset[3]]->GetXaxis()->SetTitle("Candidate Multiplicity");
    k = i;
  }

  Offset[4] =  Offset[3] + NTrThrown;

  // now setup the 2 dimensional histograms
  for (Int_t i=0;i<NTrThrown;i++) {
    char str1[128];
    char str2[128];
    sprintf(str1,"hist2d%d",i);
    sprintf(str2,"Thrown track %d NHits vs. #Delta p/p", i+1);
    hist2d[i] = new TH2F(str1, str2, 160, -2.,2.,41, -0.5,40.5);
    hist2d[i]->GetXaxis()->SetTitle("(p_det/p_gen)");
    hist2d[i]->GetYaxis()->SetTitle("Number of Hits");
    k = i;
  }
  Offset[5] = NTrThrown;

  for (Int_t i=0;i<NTrThrown;i++) {
    char str1[128];
    char str2[128];
    sprintf(str1,"hist2d%d",i+Offset[5]+1);
    sprintf(str2,"Thrown track %d NHits vs.ID/TRUE_ID", i+1);
    hist2d[i+Offset[5]+1] = new TH2F(str1, str2, 50, -.5,2.5, 41, -0.5,40.5);
    hist2d[i+Offset[5]+1]->GetXaxis()->SetTitle("ID_det/ID_true");
    hist2d[i+Offset[5]+1]->GetYaxis()->SetTitle("Number of Hits");
    k = i;
  }
  sprintf(str1,"hist2d%d",Offset[5]);
  hist2d[Offset[5]] = new TH2F(str1,"TrX NHits vs. ID/TRUE_ID X",50, -.5,2.5, 41, -0.5,40.5);
  hist2d[Offset[5]]->GetXaxis()->SetTitle("ID_det/ID_true");
  hist2d[Offset[5]]->GetYaxis()->SetTitle("Number of Hits");

}
//
//
void moreHists(TH1F** hist, TH2F** hist2d, Int_t NTrThrown, 
	       Int_t *Offset, Int_t *ThrownPType, Float_t *PMax){


  for (Int_t i=0;i<NTrThrown;i++) {
    char str1[128];
    char str2[128];
    sprintf(str1,"hist%d",i+Offset[4]);
    sprintf(str2,"Thrown Track %d PType=%d",i+1,ThrownPType[i]);
    hist[i+Offset[4]] = new TH1F(str1, str2, 50, 0.0, PMax[i]*1.1);
    hist[i+Offset[4]]->GetXaxis()->SetTitle("Track Momentum [GeV/c]");
  }

}
