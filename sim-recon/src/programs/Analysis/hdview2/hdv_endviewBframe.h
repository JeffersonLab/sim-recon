
#ifndef _HDV_ENDVIEWBFRAME_H_
#define _HDV_ENDVIEWBFRAME_H_

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

class hdv_endviewBframe:public TGMainFrame {

	public:
		hdv_endviewBframe(hdv_mainframe *hdvmf, const TGWindow *p, UInt_t w, UInt_t h);
		virtual ~hdv_endviewBframe(){};
		
		void DoDismiss(void);

		void SetRange(double xlo, double ylo, double xhi, double yhi);
		void DrawObjects(vector<TObject*> &graphics_endB);
		
	private:
	
		hdv_mainframe *hdvmf;
		TRootEmbeddedCanvas *ecanvas;
		
	ClassDef(hdv_endviewBframe,1)
};

// The following line is supposed to avoid the warning messages about:
// "dereferencing type-punned pointer will break strict-aliasing rules"
#ifdef __CINT__
#pragma link C++ class hdv_endviewBframe+;
#endif



#endif //_HDV_ENDVIEWBFRAME_H_
