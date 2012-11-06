
#ifndef _HDV_OPTIONSFRAME_H_
#define _HDV_OPTIONSFRAME_H_

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
#ifndef __CINT__
#include "hdv_mainframe.h"
#endif

class hdv_optionsframe:public TGMainFrame {

	public:
		hdv_optionsframe(hdv_mainframe *hdvmf, const TGWindow *p, UInt_t w, UInt_t h);
		virtual ~hdv_optionsframe(){};
		
		void DoDone(void);
		
	private:
	
		hdv_mainframe *hdvmf;
		map<string, TGCheckButton*> checkbuttons;
		
	ClassDef(hdv_optionsframe,1)
};

// The following line is supposed to avoid the warning messages about:
// "dereferencing type-punned pointer will break strict-aliasing rules"
#ifdef __CINT__
#pragma link C++ class hdv_optionsframe+;
#endif



#endif //_HDV_OPTIONSFRAME_H_
