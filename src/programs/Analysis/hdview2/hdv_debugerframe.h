
#ifndef _HDV_DEBUGERFRAME_H_
#define _HDV_DEBUGERFRAME_H_

// This class is made into a ROOT dictionary ala rootcint.
// Therefore, do NOT include anything Hall-D specific here.
// It is OK to do that in the .cc file, just not here in the 
// header.

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <map>

#include <TGClient.h>
#include <TGButton.h>
#include <TCanvas.h>
#include <TText.h>
#include <TRootEmbeddedCanvas.h>
#include <TTUBE.h>
#include <TNode.h>
#include <TGComboBox.h>
#include <TPolyLine.h>
#include <TEllipse.h>
#include <TMarker.h>
#include <TVector3.h>
#include <TGLabel.h>
#include <TTimer.h>


class hdv_mainframe;
class DKinematicData;

#ifndef __CINT__
#include "hdv_mainframe.h"
#include <PID/DKinematicData.h>
#endif

class hdv_debugerframe:public TGMainFrame {
  
 public:
  hdv_debugerframe(hdv_mainframe *hdvmf, const TGWindow *p, UInt_t w, UInt_t h);
  virtual ~hdv_debugerframe(){};
  
  void DoDone(void);
  void UpdateTrackLabels();
  TGHorizontalFrame *topframe;
  TGHorizontalFrame *botframe;
  TGGroupFrame *hitdrawopts ;
  TGTextButton *done ;

  Int_t GetNTrCand(void) {return NTrCand;}
  void SetNTrCand(Int_t d) { NTrCand = d;}
  void SetTrackCandidates(vector<const DKinematicData*> d) {TrackCandidates=d;}

 private:

  Int_t NTrCand;
  vector<const DKinematicData*> TrackCandidates;

  map<string, TGVerticalFrame *> tf;
  map<string, vector<TGLabel*> > candlabs;
	       
  hdv_mainframe *hdvmf;
  map<string, TGCheckButton*> checkbuttons;
  
  ClassDef(hdv_debugerframe,1)
    };

// The following line is supposed to avoid the warning messages about:
// "dereferencing type-punned pointer will break strict-aliasing rules"
#ifdef __CINT__
#pragma link C++ class hdv_debugerframe+;
#endif



#endif //_HDV_DEBUGERFRAME_H_
