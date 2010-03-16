
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <fstream>
using namespace std;

#include <pthread.h>

#include <TRACKING/DMCThrown.h>
#include "hdv_optionsframe.h"
#include "hdview2.h"
#include "MyProcessor.h"
#include "FDC/DFDCGeometry.h"
#include "FCAL/DFCALGeometry.h"
#include "DVector2.h"
#include "HDGEOMETRY/DGeometry.h"
#include "PID/DPhoton.h"

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
hdv_optionsframe::hdv_optionsframe(hdv_mainframe *hdvmf, const TGWindow *p, UInt_t w, UInt_t h):TGMainFrame(p,w,h)
{
	this->hdvmf = hdvmf;

	// First, define all of the of the graphics objects. Below that, make all
	// of the connections to the methods so these things will work!

	// The main GUI window is divided into three sections, top, middle, and bottom.
	// Create those frames here.
	TGLayoutHints *hints = new TGLayoutHints(kLHintsNormal|kLHintsExpandX|kLHintsExpandY, 5,5,5,5);
	TGLayoutHints *lhints = new TGLayoutHints(kLHintsNormal, 2,2,2,2);
	TGLayoutHints *rhints = new TGLayoutHints(kLHintsCenterY|kLHintsRight, 2,2,2,2);
	TGLayoutHints *chints = new TGLayoutHints(kLHintsCenterY|kLHintsCenterX, 2,2,2,2);
	TGLayoutHints *bhints = new TGLayoutHints(kLHintsBottom|kLHintsCenterX, 2,2,2,2);
	TGLayoutHints *xhints = new TGLayoutHints(kLHintsNormal|kLHintsExpandX, 2,2,2,2);
	TGLayoutHints *yhints = new TGLayoutHints(kLHintsNormal|kLHintsExpandY, 2,2,2,2);
	TGLayoutHints *dhints = new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 0,0,0,0);
	TGLayoutHints *ehints = new TGLayoutHints(kLHintsNormal, 2,2,0,0);
	TGLayoutHints *thints = new TGLayoutHints(kLHintsTop|kLHintsCenterX, 2,2,0,0);
	TGLayoutHints *lxhints = new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 2,2,0,0);
	TGHorizontalFrame *topframe = new TGHorizontalFrame(this, w, h);
	TGHorizontalFrame *botframe = new TGHorizontalFrame(this, w, h);
	AddFrame(topframe, lhints);
	AddFrame(botframe, chints);

	TGGroupFrame *hitdrawopts = new TGGroupFrame(topframe, "Extra Drawing Options", kVerticalFrame);
	topframe->AddFrame(hitdrawopts, lhints);

	checkbuttons["fdcintersection"]				= new TGCheckButton(hitdrawopts,	"FDC Intersection");
	TGHorizontal3DLine *separator0				= new TGHorizontal3DLine(hitdrawopts, 250);
	checkbuttons["bcaltruth"]						= new TGCheckButton(hitdrawopts,	"BCALTruth");
	checkbuttons["thrown_charged_bcal"]			= new TGCheckButton(hitdrawopts,	"Draw thrown charged track projections on BCAL");
	checkbuttons["recon_charged_bcal"]			= new TGCheckButton(hitdrawopts,	"Draw reconstructed charged track projections on BCAL");
	checkbuttons["thrown_photons_bcal"]			= new TGCheckButton(hitdrawopts,	"Draw thrown photon projections on BCAL");
	checkbuttons["recon_photons_bcal"]			= new TGCheckButton(hitdrawopts,	"Draw reconstructed photon projections on BCAL");
	TGHorizontal3DLine *separator1				= new TGHorizontal3DLine(hitdrawopts, 250);
	checkbuttons["fcaltruth"]						= new TGCheckButton(hitdrawopts,	"FCALTruth");
	checkbuttons["thrown_charged_fcal"]			= new TGCheckButton(hitdrawopts,	"Draw thrown charged track projections on FCAL");
	checkbuttons["recon_charged_fcal"]			= new TGCheckButton(hitdrawopts,	"Draw reconstructed charged track projections on FCAL");
	checkbuttons["thrown_photons_fcal"]			= new TGCheckButton(hitdrawopts,	"Draw thrown photon projections on FCAL");
	checkbuttons["recon_photons_fcal"]			= new TGCheckButton(hitdrawopts,	"Draw reconstructed photon projections on FCAL");
	TGHorizontal3DLine *separator2				= new TGHorizontal3DLine(hitdrawopts, 250);
	checkbuttons["recon_photons_track_match"]	= new TGCheckButton(hitdrawopts,	"Draw reconstructed photons matched to charged tracks");
	TGHorizontal3DLine *separator3				= new TGHorizontal3DLine(hitdrawopts, 250);
	checkbuttons["trajectories_photon"]			= new TGCheckButton(hitdrawopts,	"When drawing DMCTrajectoryPoint, draw photon tracks");
	checkbuttons["trajectories_electron"]		= new TGCheckButton(hitdrawopts,	"When drawing DMCTrajectoryPoint, draw electron tracks");
	checkbuttons["trajectories_positron"]		= new TGCheckButton(hitdrawopts,	"When drawing DMCTrajectoryPoint, draw positron tracks");
	checkbuttons["trajectories_proton"]			= new TGCheckButton(hitdrawopts,	"When drawing DMCTrajectoryPoint, draw proton tracks");
	checkbuttons["trajectories_neutron"]		= new TGCheckButton(hitdrawopts,	"When drawing DMCTrajectoryPoint, draw neutron tracks");
	checkbuttons["trajectories_piplus"]			= new TGCheckButton(hitdrawopts,	"When drawing DMCTrajectoryPoint, draw piplus tracks");
	checkbuttons["trajectories_piminus"]		= new TGCheckButton(hitdrawopts,	"When drawing DMCTrajectoryPoint, draw piminus tracks");
	checkbuttons["trajectories_other"]			= new TGCheckButton(hitdrawopts,	"When drawing DMCTrajectoryPoint, draw all other tracks");

	hitdrawopts->AddFrame(checkbuttons["fdcintersection"], lhints);
	hitdrawopts->AddFrame(separator0, chints);
	hitdrawopts->AddFrame(checkbuttons["bcaltruth"], lhints);
	hitdrawopts->AddFrame(checkbuttons["thrown_charged_bcal"], lhints);
	hitdrawopts->AddFrame(checkbuttons["recon_charged_bcal"], lhints);
	hitdrawopts->AddFrame(checkbuttons["thrown_photons_bcal"], lhints);
	hitdrawopts->AddFrame(checkbuttons["recon_photons_bcal"], lhints);
	hitdrawopts->AddFrame(separator1, chints);
	hitdrawopts->AddFrame(checkbuttons["fcaltruth"], lhints);
	hitdrawopts->AddFrame(checkbuttons["thrown_charged_fcal"], lhints);
	hitdrawopts->AddFrame(checkbuttons["recon_charged_fcal"], lhints);
	hitdrawopts->AddFrame(checkbuttons["thrown_photons_fcal"], lhints);
	hitdrawopts->AddFrame(checkbuttons["recon_photons_fcal"], lhints);
	hitdrawopts->AddFrame(separator2, chints);
	hitdrawopts->AddFrame(checkbuttons["recon_photons_track_match"], lhints);
	hitdrawopts->AddFrame(separator3, chints);
	hitdrawopts->AddFrame(checkbuttons["trajectories_photon"], lhints);
	hitdrawopts->AddFrame(checkbuttons["trajectories_electron"], lhints);
	hitdrawopts->AddFrame(checkbuttons["trajectories_positron"], lhints);
	hitdrawopts->AddFrame(checkbuttons["trajectories_proton"], lhints);
	hitdrawopts->AddFrame(checkbuttons["trajectories_neutron"], lhints);
	hitdrawopts->AddFrame(checkbuttons["trajectories_piplus"], lhints);
	hitdrawopts->AddFrame(checkbuttons["trajectories_piminus"], lhints);
	hitdrawopts->AddFrame(checkbuttons["trajectories_other"], lhints);

	//========== Done Button ===========
	TGTextButton *done = new TGTextButton(botframe,	"Done");
	botframe->AddFrame(done, chints);

	//&&&&&&&&&&&&&&&& Connections
	map<string, TGCheckButton*>::iterator iter = checkbuttons.begin();
	for(; iter!=checkbuttons.end(); iter++){
		iter->second->Connect("Clicked()","hdv_mainframe", hdvmf, "DoMyRedraw()");
	}

	// Add out checkbuttons to the list kept in hdv_mainframe
	hdvmf->AddCheckButtons(checkbuttons);

	// Finish up and map the window
	SetWindowName("Hall-D Event Viewer Options");
	SetIconName("HDView");

	done->Connect("Clicked()","hdv_optionsframe", this, "DoDone()");

	MapSubwindows();
	Resize(GetDefaultSize());
	//MapWindow();
	//LowerWindow();
}

//-------------------
// DoDone
//-------------------
void hdv_optionsframe::DoDone(void)
{
	//LowerWindow();
	UnmapWindow();
}


