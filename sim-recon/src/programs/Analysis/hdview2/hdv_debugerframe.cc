
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <fstream>
using namespace std;

#include <pthread.h>

#include <TRACKING/DMCThrown.h>
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DTrackTimeBased_factory.h>
#include "hdview2.h"
#include "MyProcessor.h"
#include "hdv_debugerframe.h"
#include "FDC/DFDCGeometry.h"
#include "FCAL/DFCALGeometry.h"
#include "DVector2.h"
#include "HDGEOMETRY/DGeometry.h"
#include "PID/DNeutralParticle.h"
#include <PID/DParticleSet.h>
#include <PID/DPhysicsEvent.h>
#include <PID/DKinematicData.h>

#include <TPolyMarker.h>
#include <TLine.h>
#include <TMarker.h>
#include <TBox.h>
#include <TVector3.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include <TGLabel.h>
#include <TGComboBox.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGTextEntry.h>
#include <TArrow.h>
#include <TLatex.h>
#include <TColor.h>
#include <TG3DLine.h>
#include <TMath.h>

extern JApplication *japp;
//TGeoVolume *MOTHER = NULL;
//TGeoCombiTrans *MotherRotTrans = NULL;

extern int GO;

//-------------------
// Constructor
//-------------------
hdv_debugerframe::hdv_debugerframe(hdv_mainframe *hdvmf, const TGWindow *p, UInt_t w, UInt_t h):TGMainFrame(p,w,h)
{
  this->hdvmf = hdvmf;
  NTrCand = 0;
  NTrWireBased = 0;
  NTrTimeBased = 0;
  InitMid1Frame = 0;
  InitMid2Frame = 0;

  // First, define all of the of the graphics objects. Below that, make all
  // of the connections to the methods so these things will work!
  
  // The main GUI window is divided into three sections, top, middle, and bottom.
  // Create those frames here.
  TGLayoutHints *hints = new TGLayoutHints(kLHintsNormal|kLHintsExpandX|kLHintsExpandY|kLHintsLeft, 5,5,5,5);
  TGLayoutHints *lhints = new TGLayoutHints(kLHintsNormal, 2,2,2,2);
  TGLayoutHints *chints = new TGLayoutHints(kLHintsCenterY|kLHintsCenterX, 2,2,2,2);
  TGLayoutHints *xhints = new TGLayoutHints(kLHintsNormal|kLHintsExpandX, 2,2,2,2);

  topframe = new TGHorizontalFrame(this,  800, 500);
  mid1frame = new TGHorizontalFrame(this, 800, 500);
  mid2frame = new TGHorizontalFrame(this, 800, 500);
  botframe = new TGHorizontalFrame(this,  400, 50);

  AddFrame(topframe, lhints);
  AddFrame(mid1frame, lhints);
  AddFrame(mid2frame, lhints);
  AddFrame(botframe, chints);


  // top frame for track Candidates
  hitdrawopts = new TGGroupFrame(topframe, "Track Candidate Hits", kVerticalFrame);
  topframe->AddFrame(hitdrawopts, lhints);

  for (int k=0;k<10;k++){
    char str1[128];
    sprintf(str1,"Candidate%d",k+1);
    char str2[128];
    if (k<NTrCand){ 
      sprintf(str2,"Hits Track Candidate  %d",k+1);
      if (k>9)
	sprintf(str2,"Hits Track Candidate %d",k+1);
    } else {
        sprintf(str2,".......................");
    }
    checkbuttons[str1] = new TGCheckButton(hitdrawopts,	str2); 
    hitdrawopts->AddFrame(checkbuttons[str1],lhints);
  }
  
  // Track Candidate Info
  TGGroupFrame *trackinfo = new TGGroupFrame(topframe, "     trk:     type:     p:     theta:   phi:       z:  ", kHorizontalFrame);
  topframe->AddFrame(trackinfo, hints);
  // Column names
  vector<string> colnames;
  colnames.push_back("trk");
  colnames.push_back("type");
  colnames.push_back("p");
  colnames.push_back("theta");
  colnames.push_back("phi");
  colnames.push_back("z");
  
  for(unsigned int i=0; i<colnames.size(); i++){
    // create frames
    tf[colnames[i]] = new TGVerticalFrame(trackinfo);
    trackinfo->AddFrame( tf[colnames[i]], xhints);
    //string lab = colnames[i]+":";
    //TGLabel *tl = new TGLabel(tf[colnames[i]], lab.c_str());
    //tf[colnames[i]]->AddFrame(tl, chints);

    vector<TGLabel*> tv;
    //tv.push_back(tl);
    for (int k=0;k<10;k++){
      TGLabel *lab = new TGLabel(tf[colnames[i]],"-----"); 
      tf[colnames[i]]->AddFrame(lab, chints);
      tv.push_back(lab);
    }
    candlabs[colnames[i]] = tv;
  }

  for (int k=0;k<10;k++){
    if (k<(int)TrackCandidates.size()){ 
      const DKinematicData *trk = TrackCandidates[k];
      stringstream trkno, type, p, theta, phi, z;
      trkno<<setprecision(4)<<k+1;
      candlabs["trk"][k]->SetText(trkno.str().c_str());
 
      double mass = trk->mass();
      if(fabs(mass-0.13957)<1.0E-4)type<<"pi";
      else if(fabs(mass-0.93827)<1.0E-4)type<<"proton";
      else if(fabs(mass-0.493677)<1.0E-4)type<<"K";
      else if(fabs(mass-0.000511)<1.0E-4)type<<"e";
      else if (fabs(mass)<1.e-4 && fabs(trk->charge())<1.e-4) type << "gamma";
      else type<<"q=";
      if (fabs(trk->charge())>1.e-4){
	type<<(trk->charge()>0 ? "+":"-");
      }
      int row = k;
      candlabs["type"][row]->SetText(type.str().c_str());

      p<<setprecision(3)<<fixed<<trk->momentum().Mag();
      candlabs["p"][row]->SetText(p.str().c_str());
      
      theta<<setprecision(2)<<fixed<<trk->momentum().Theta()*TMath::RadToDeg();
      candlabs["theta"][row]->SetText(theta.str().c_str());
      
      double myphi = trk->momentum().Phi();
      if(myphi<0.0)myphi+=2.0*M_PI;
      phi<<setprecision(2)<<fixed<<myphi;
      candlabs["phi"][row]->SetText(phi.str().c_str());
      
      z<<setprecision(2)<<fixed<<trk->position().Z();
      candlabs["z"][row]->SetText(z.str().c_str());
    } else {
      int row = k;
      candlabs["trk"][row]->SetText("------");
      candlabs["type"][row]->SetText("------");
      candlabs["p"][row]->SetText("------");
      candlabs["theta"][row]->SetText("------");
      candlabs["phi"][row]->SetText("------");
      candlabs["z"][row]->SetText("------");
    }
  
  }

  // mid1 frame for Wire based Tracking
  hitdrawoptsWB = new TGGroupFrame(mid1frame, "Wire Based Hits", kVerticalFrame);
  mid1frame->AddFrame(hitdrawoptsWB, lhints);
  trackinfoWB = new TGGroupFrame(mid1frame, "    trk:       type:       p:       theta:     phi:       z:  chisq/Ndof: Ndof:    cand:", 
					       kHorizontalFrame);
  mid1frame->AddFrame(trackinfoWB, hints);
  SetUpMid1Frame();



  hitdrawoptsTB = new TGGroupFrame(mid2frame, "Time Based Hits", kVerticalFrame);
  mid2frame->AddFrame(hitdrawoptsTB, lhints);
  trackinfoTB = new TGGroupFrame(mid2frame, "    trk:       type:       p:       theta:     phi:       z:  chisq/Ndof: Ndof:    cand:", 
					       kHorizontalFrame);
  mid2frame->AddFrame(trackinfoTB, hints);
  SetUpMid2Frame();

       
  //========== Done Button ===========
  done = new TGTextButton(botframe,"Done");
  botframe->AddFrame(done, chints);
  
  //&&&&&&&&&&&&&&&& Connections
  map<string, TGCheckButton*>::iterator iter = checkbuttons.begin();
  for(; iter!=checkbuttons.end(); iter++){
    iter->second->Connect("Clicked()","hdv_mainframe", hdvmf, "DoMyRedraw()");
  }
  
  // Add out checkbuttons to the list kept in hdv_mainframe
  hdvmf->AddCheckButtons(checkbuttons);
  
  // Finish up and map the window
  SetWindowName("Hall-D Event Viewer Debuger");
  SetIconName("HDView");
  
  done->Connect("Clicked()","hdv_debugerframe", this, "DoDone()");
  
  SetWindowName("Hall-D Event View Debuger");
  SetIconName("HDDebugView");
  
  MapSubwindows();
  Resize(GetDefaultSize());

}

//-------------------
// DoDone
//-------------------
void hdv_debugerframe::DoDone(void)
{
	//LowerWindow();
	UnmapWindow();
}

//-------------------
// UpdateTrackLables
//-------------------

void hdv_debugerframe::UpdateTrackLabels()
{
  
  int size = NTrCand; 
  for (int k=0;k<10;k++){
    char str1[128];
    sprintf(str1,"Candidate%d",k+1);
    char str2[128];
    if (k<size){ 
      sprintf(str2,"Hits Track Candidate  %d",k+1);
      if (k>9)
	sprintf(str2,"Hits Track Candidate %d",k+1);
    } else {
        sprintf(str2,".......................");
    }
    checkbuttons[str1]->SetText(str2);

  }
  for (int k=0;k<10;k++){

    if (k<(int)TrackCandidates.size()){ 
      const DKinematicData *trk = TrackCandidates[k];
      stringstream trkno, type, p, theta, phi, z;
      trkno<<setprecision(4)<<k+1;
      int row = k;
      candlabs["trk"][row]->SetText(trkno.str().c_str());
 
      double mass = trk->mass();
      if(fabs(mass-0.13957)<1.0E-4)type<<"pi";
      else if(fabs(mass-0.93827)<1.0E-4)type<<"proton";
      else if(fabs(mass-0.493677)<1.0E-4)type<<"K";
      else if(fabs(mass-0.000511)<1.0E-4)type<<"e";
      else if (fabs(mass)<1.e-4 && fabs(trk->charge())<1.e-4) type << "gamma";
      else type<<"q=";
      if (fabs(trk->charge())>1.e-4){
	type<<(trk->charge()>0 ? "+":"-");
      }
      candlabs["type"][row]->SetText(type.str().c_str());

      p<<setprecision(3)<<fixed<<trk->momentum().Mag();
      candlabs["p"][row]->SetText(p.str().c_str());
      
      theta<<setprecision(2)<<fixed<<trk->momentum().Theta()*TMath::RadToDeg();
      candlabs["theta"][row]->SetText(theta.str().c_str());
      
      double myphi = trk->momentum().Phi();
      if(myphi<0.0)myphi+=2.0*M_PI;
      phi<<setprecision(2)<<fixed<<myphi;
      candlabs["phi"][row]->SetText(phi.str().c_str());
      
      z<<setprecision(2)<<fixed<<trk->position().Z();
      candlabs["z"][row]->SetText(z.str().c_str());
    } else {
      int row = k;  
      candlabs["trk"][row]->SetText("------");
      candlabs["type"][row]->SetText("------");
      candlabs["p"][row]->SetText("------");
      candlabs["theta"][row]->SetText("------");
      candlabs["phi"][row]->SetText("------");
      candlabs["z"][row]->SetText("------");
    }
  }

  SetUpMid1Frame();
  SetUpMid2Frame();

  map<string, TGCheckButton*>::iterator iter = checkbuttons.begin();
  for(; iter!=checkbuttons.end(); iter++){
    iter->second->Connect("Clicked()","hdv_mainframe", hdvmf, "DoMyRedraw()");
  }
  
  MapSubwindows();
  Resize(GetDefaultSize());

}

//************************
//*
//*    SetUpMid1Frame
//*
//************************
void hdv_debugerframe::SetUpMid1Frame(){

  TGLayoutHints *lhints = new TGLayoutHints(kLHintsNormal, 2,2,2,2);
  TGLayoutHints *xhints = new TGLayoutHints(kLHintsNormal|kLHintsExpandX, 2,2,2,2);
  TGLayoutHints *chints = new TGLayoutHints(kLHintsCenterY|kLHintsCenterX, 2,2,2,2);

  
  for (int k=0;k<10;k++){
    char str1[128];
    sprintf(str1,"WireBased%d",k+1);
    char str2[128];
    if (k<NTrWireBased){
      if (k<10){ 
	sprintf(str2,"Hits Wire Based Track  %d",k+1);
	if (k>9)
	  sprintf(str2,"Hits Wire Based Track %d",k+1);
      }
    } else{ 
      sprintf(str2,".......................");
    }
    if (!InitMid1Frame) {
      checkbuttons[str1] = new TGCheckButton(hitdrawoptsWB, str2); 
      hitdrawoptsWB->AddFrame(checkbuttons[str1], lhints);
    } else {
      checkbuttons[str1]->SetText(str2);  
    }
  } 

  if (!InitMid1Frame) {    
    vector<string> colnamesw;
    colnamesw.push_back("trk");
    colnamesw.push_back("type");
    colnamesw.push_back("p");
    colnamesw.push_back("theta");
    colnamesw.push_back("phi");
    colnamesw.push_back("z");
    colnamesw.push_back("chisq/Ndof");
    colnamesw.push_back("Ndof");
    colnamesw.push_back("cand");
    
    for(unsigned int i=0; i<colnamesw.size(); i++){
      // create frames
      tfWB[colnamesw[i]] = new TGVerticalFrame(trackinfoWB);
      trackinfoWB->AddFrame( tfWB[colnamesw[i]], xhints);
      //string lab = colnamesw[i]+":";
      //TGLabel *tl = new TGLabel(tfWB[colnamesw[i]], lab.c_str());
      //tfWB[colnamesw[i]]->AddFrame(tl, chints);
      
      vector<TGLabel*> tv;
      //tv.push_back(tl);
      for (int k=0;k<10;k++){
	TGLabel *lab = new TGLabel(tfWB[colnamesw[i]],"-----"); 
	tfWB[colnamesw[i]]->AddFrame(lab, chints);
	tv.push_back(lab);
      }
      wblabs[colnamesw[i]] = tv;
    }
    
    InitMid1Frame = 1;
  }
  
  for (int k=0;k<10;k++){
    
    if (k<NTrWireBased) {

      const DTrackWireBased *trk = subTrackWireBased[k];
      stringstream trkno, type, p, theta, phi, z, chisq_per_dof, Ndof, cand;
      trkno<<setprecision(4)<<k+1;
      
      int row = k;
      wblabs["trk"][row]->SetText(trkno.str().c_str());
      
      double mass = trk->mass();
      if(fabs(mass-0.13957)<1.0E-4)type<<"pi";
      else if(fabs(mass-0.93827)<1.0E-4)type<<"proton";
      else if(fabs(mass-0.493677)<1.0E-4)type<<"K";
      else if(fabs(mass-0.000511)<1.0E-4)type<<"e";
      else if (fabs(mass)<1.e-4 && fabs(trk->charge())<1.e-4) type << "gamma";
      else type<<"q=";
      if (fabs(trk->charge())>1.e-4){
	type<<(trk->charge()>0 ? "+":"-");
      }
      
      wblabs["type"][row]->SetText(type.str().c_str());
      
      p<<setprecision(3)<<fixed<<trk->momentum().Mag();
      wblabs["p"][row]->SetText(p.str().c_str());
      
      theta<<setprecision(2)<<fixed<<trk->momentum().Theta()*TMath::RadToDeg();
      wblabs["theta"][row]->SetText(theta.str().c_str());
      
      double myphi = trk->momentum().Phi();
      if(myphi<0.0)myphi+=2.0*M_PI;
      phi<<setprecision(2)<<fixed<<myphi;
      wblabs["phi"][row]->SetText(phi.str().c_str());
      
      z<<setprecision(2)<<fixed<<trk->position().Z();
      wblabs["z"][row]->SetText(z.str().c_str());
      
      chisq_per_dof<<setprecision(1)<<fixed<<trk->chisq/trk->Ndof;
      Ndof<<trk->Ndof;
      cand << trk->candidateid;
      
      wblabs["chisq/Ndof"][row]->SetText(chisq_per_dof.str().c_str());
      wblabs["Ndof"][row]->SetText(Ndof.str().c_str());
      wblabs["cand"][row]->SetText(cand.str().c_str());

    } else {
      wblabs["trk"][k]->SetText("-----");
      wblabs["type"][k]->SetText("-----");
      wblabs["p"][k]->SetText("-----");
      wblabs["theta"][k]->SetText("-----");
      wblabs["phi"][k]->SetText("-----");
      wblabs["z"][k]->SetText("-----");
      wblabs["chisq/Ndof"][k]->SetText("-----");
      wblabs["Ndof"][k]->SetText("-----");
      wblabs["cand"][k]->SetText("-----");

    }
      
  }  

}
//************************
//*
//*    SetUpMid2Frame
//*
//************************
void hdv_debugerframe::SetUpMid2Frame(){

  TGLayoutHints *lhints = new TGLayoutHints(kLHintsNormal, 2,2,2,2);
  TGLayoutHints *xhints = new TGLayoutHints(kLHintsNormal|kLHintsExpandX, 2,2,2,2);
  TGLayoutHints *chints = new TGLayoutHints(kLHintsCenterY|kLHintsCenterX, 2,2,2,2);

  
  for (int k=0;k<10;k++){
    char str1[128];
    sprintf(str1,"TimeBased%d",k+1);
    char str2[128];
    if (k<NTrTimeBased){
      if (k<10){ 
	sprintf(str2,"Hits Time Based Track  %d",k+1);
	if (k>9)
	  sprintf(str2,"Hits Time Based Track %d",k+1);
      }
    } else{ 
      sprintf(str2,".......................");
    }
    if (!InitMid2Frame) {
      checkbuttons[str1] = new TGCheckButton(hitdrawoptsTB, str2); 
      hitdrawoptsTB->AddFrame(checkbuttons[str1], lhints);
    } else {
      checkbuttons[str1]->SetText(str2);  
    }
  } 

  if (!InitMid2Frame) {    
    vector<string> colnamesw;
    colnamesw.push_back("trk");
    colnamesw.push_back("type");
    colnamesw.push_back("p");
    colnamesw.push_back("theta");
    colnamesw.push_back("phi");
    colnamesw.push_back("z");
    colnamesw.push_back("chisq/Ndof");
    colnamesw.push_back("Ndof");
    colnamesw.push_back("cand");
    
    for(unsigned int i=0; i<colnamesw.size(); i++){
      // create frames
      tfTB[colnamesw[i]] = new TGVerticalFrame(trackinfoTB);
      trackinfoTB->AddFrame( tfTB[colnamesw[i]], xhints);
      //string lab = colnamesw[i]+":";
      //TGLabel *tl = new TGLabel(tfWB[colnamesw[i]], lab.c_str());
      //tfWB[colnamesw[i]]->AddFrame(tl, chints);
      
      vector<TGLabel*> tv;
      //tv.push_back(tl);
      for (int k=0;k<10;k++){
	TGLabel *lab = new TGLabel(tfTB[colnamesw[i]],"-----"); 
	tfTB[colnamesw[i]]->AddFrame(lab, chints);
	tv.push_back(lab);
      }
      tblabs[colnamesw[i]] = tv;
    }
    
    InitMid2Frame = 1;
  }
  
  for (int k=0;k<10;k++){
    
    if (k<NTrTimeBased) {

      const DTrackTimeBased *trk = subTrackTimeBased[k];
      stringstream trkno, type, p, theta, phi, z, chisq_per_dof, Ndof, cand;
      trkno<<setprecision(4)<<k+1;
      
      int row = k;
      tblabs["trk"][row]->SetText(trkno.str().c_str());
      
      double mass = trk->mass();
      if(fabs(mass-0.13957)<1.0E-4)type<<"pi";
      else if(fabs(mass-0.93827)<1.0E-4)type<<"proton";
      else if(fabs(mass-0.493677)<1.0E-4)type<<"K";
      else if(fabs(mass-0.000511)<1.0E-4)type<<"e";
      else if (fabs(mass)<1.e-4 && fabs(trk->charge())<1.e-4) type << "gamma";
      else type<<"q=";
      if (fabs(trk->charge())>1.e-4){
	type<<(trk->charge()>0 ? "+":"-");
      }
      
      tblabs["type"][row]->SetText(type.str().c_str());
      
      p<<setprecision(3)<<fixed<<trk->momentum().Mag();
      tblabs["p"][row]->SetText(p.str().c_str());
      
      theta<<setprecision(2)<<fixed<<trk->momentum().Theta()*TMath::RadToDeg();
      tblabs["theta"][row]->SetText(theta.str().c_str());
      
      double myphi = trk->momentum().Phi();
      if(myphi<0.0)myphi+=2.0*M_PI;
      phi<<setprecision(2)<<fixed<<myphi;
      tblabs["phi"][row]->SetText(phi.str().c_str());
      
      z<<setprecision(2)<<fixed<<trk->position().Z();
      tblabs["z"][row]->SetText(z.str().c_str());
      
      chisq_per_dof<<setprecision(1)<<fixed<<trk->chisq/trk->Ndof;
      Ndof<<trk->Ndof;
      cand << trk->candidateid;
      
      tblabs["chisq/Ndof"][row]->SetText(chisq_per_dof.str().c_str());
      tblabs["Ndof"][row]->SetText(Ndof.str().c_str());
      tblabs["cand"][row]->SetText(cand.str().c_str());

    } else {
      tblabs["trk"][k]->SetText("-----");
      tblabs["type"][k]->SetText("-----");
      tblabs["p"][k]->SetText("-----");
      tblabs["theta"][k]->SetText("-----");
      tblabs["phi"][k]->SetText("-----");
      tblabs["z"][k]->SetText("-----");
      tblabs["chisq/Ndof"][k]->SetText("-----");
      tblabs["Ndof"][k]->SetText("-----");
      tblabs["cand"][k]->SetText("-----");

    }
      
  }  

}
