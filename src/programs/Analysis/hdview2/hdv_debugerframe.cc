
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
#include "PID/DNeutralTrack.h"
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

  // First, define all of the of the graphics objects. Below that, make all
  // of the connections to the methods so these things will work!
  
  // The main GUI window is divided into three sections, top, middle, and bottom.
  // Create those frames here.
  TGLayoutHints *lhints = new TGLayoutHints(kLHintsNormal, 2,2,2,2);
  TGLayoutHints *chints = new TGLayoutHints(kLHintsCenterY|kLHintsCenterX, 2,2,2,2);
  TGLayoutHints *xhints = new TGLayoutHints(kLHintsNormal|kLHintsExpandX, 2,2,2,2);
  topframe = new TGHorizontalFrame(this, 400, h);
  botframe = new TGHorizontalFrame(this, 400, h);
  AddFrame(topframe, lhints);
  AddFrame(botframe, chints);
  
  

  hitdrawopts = new TGGroupFrame(topframe, "Drawing Debuger Options", kVerticalFrame);
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
  TGGroupFrame *trackinfo = new TGGroupFrame(topframe, "Track Candidate Info", kHorizontalFrame);
  topframe->AddFrame(trackinfo, xhints);
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
    string lab = colnames[i]+":";
    TGLabel *tl = new TGLabel(tf[colnames[i]], lab.c_str());
    tf[colnames[i]]->AddFrame(tl, chints);

    vector<TGLabel*> tv;
    tv.push_back(tl);
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
      int row = k+1;
      candlabs["type"][k]->SetText(type.str().c_str());

      p<<setprecision(4)<<trk->momentum().Mag();
      candlabs["p"][row]->SetText(p.str().c_str());
      
      theta<<setprecision(4)<<trk->momentum().Theta()*TMath::RadToDeg();
      candlabs["theta"][row]->SetText(theta.str().c_str());
      
      double myphi = trk->momentum().Phi();
      if(myphi<0.0)myphi+=2.0*M_PI;
      phi<<setprecision(4)<<myphi;
      candlabs["phi"][row]->SetText(phi.str().c_str());
      
      z<<setprecision(4)<<trk->position().Z();
      candlabs["z"][row]->SetText(z.str().c_str());
    } else {
      candlabs["trk"][k+1]->SetText("------");
      candlabs["type"][k+1]->SetText("------");
      candlabs["p"][k+1]->SetText("------");
      candlabs["theta"][k+1]->SetText("------");
      candlabs["phi"][k+1]->SetText("------");
      candlabs["z"][k+1]->SetText("------");
    }
  
  }


      
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
  
  MapSubwindows();
  Resize(GetDefaultSize());
  //MapWindow();
  //LowerWindow();
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
  
  int size = GetNTrCand(); 
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
      int row = k+1;
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

      p<<setprecision(4)<<trk->momentum().Mag();
      candlabs["p"][row]->SetText(p.str().c_str());
      
      theta<<setprecision(4)<<trk->momentum().Theta()*TMath::RadToDeg();
      candlabs["theta"][row]->SetText(theta.str().c_str());
      
      double myphi = trk->momentum().Phi();
      if(myphi<0.0)myphi+=2.0*M_PI;
      phi<<setprecision(4)<<myphi;
      candlabs["phi"][row]->SetText(phi.str().c_str());
      
      z<<setprecision(4)<<trk->position().Z();
      candlabs["z"][row]->SetText(z.str().c_str());
    } else {
      candlabs["trk"][k+1]->SetText("------");
      candlabs["type"][k+1]->SetText("------");
      candlabs["p"][k+1]->SetText("------");
      candlabs["theta"][k+1]->SetText("------");
      candlabs["phi"][k+1]->SetText("------");
      candlabs["z"][k+1]->SetText("------");
    }
  }
  map<string, TGCheckButton*>::iterator iter = checkbuttons.begin();
  for(; iter!=checkbuttons.end(); iter++){
    iter->second->Connect("Clicked()","hdv_mainframe", hdvmf, "DoMyRedraw()");
  }
  
  MapSubwindows();
  //Resize(GetDefaultSize());
  Resize(600,300);

}
