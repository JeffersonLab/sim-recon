
#ifndef _HDV_MAINFRAME_H_
#define _HDV_MAINFRAME_H_

#include "hdview.h"

#include <TGClient.h>
#include <TGButton.h>
#include <TCanvas.h>
#include <TRootEmbeddedCanvas.h>
#include <TTUBE.h>
#include <TNode.h>

class hdv_mainframe:public TGMainFrame {

	public:
		hdv_mainframe(const TGWindow *p, UInt_t w, UInt_t h);
		~hdv_mainframe(){};
		
		Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);

	private:
		TGTextButton *quit, *next;
		TGLayoutHints *fLayout;
		TRootEmbeddedCanvas *emcanvas;
};


#endif //_HDV_MAINFRAME_H_
