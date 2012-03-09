
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <fstream>
using namespace std;

#include <pthread.h>

#include <TRACKING/DMCThrown.h>
#include "hdv_fulllistframe.h"
#include "hdview2.h"
#include "MyProcessor.h"
#include "FDC/DFDCGeometry.h"
#include "FCAL/DFCALGeometry.h"
#include "DVector2.h"
#include "HDGEOMETRY/DGeometry.h"
#include <PID/DNeutralParticle.h>
#include <PID/DParticleSet.h>
#include <PID/DPhysicsEvent.h>
#include <PID/DTwoGammaFit.h>
#include <TRACKING/DTrackTimeBased.h>
#include <TRACKING/DTrackWireBased.h>
#include <TRACKING/DMCThrown.h>

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
#include <TGListView.h>
#include <TGTextView.h>

extern JApplication *japp;
//TGeoVolume *MOTHER = NULL;
//TGeoCombiTrans *MotherRotTrans = NULL;

extern int GO;


//-------------------
// Constructor
//-------------------
hdv_fulllistframe::hdv_fulllistframe(hdv_mainframe *hdvmf, const TGWindow *p, UInt_t w, UInt_t h):TGMainFrame(p,w,h)
{
	this->hdvmf = hdvmf;

	// First, define all of the of the graphics objects. Below that, make all
	// of the connections to the methods so these things will work!

	// The main GUI window is divided into three sections, top, middle, and bottom.
	// Create those frames here.
	TGLayoutHints *lhints = new TGLayoutHints(kLHintsNormal, 2,2,2,2);
	TGLayoutHints *chints = new TGLayoutHints(kLHintsCenterY|kLHintsCenterX, 2,2,2,2);
	TGLayoutHints *xhints = new TGLayoutHints(kLHintsNormal|kLHintsExpandX, 2,2,2,2);
	TGHorizontalFrame *topframe = new TGHorizontalFrame(this, w, h);
	TGHorizontalFrame *botframe = new TGHorizontalFrame(this, w, h);
	AddFrame(topframe, lhints);
	AddFrame(botframe, chints);

	TGGroupFrame *trackinfo = new TGGroupFrame(topframe, "Track Info", kHorizontalFrame);
	topframe->AddFrame(trackinfo, xhints);

		//------ Track Info ------
		throwninfo = new TGGroupFrame(trackinfo, "Thrown", kHorizontalFrame);
		reconinfo = new TGGroupFrame(trackinfo, "Reconstructed", kHorizontalFrame);
		trackinfo->AddFrame(throwninfo, lhints);
		trackinfo->AddFrame(reconinfo, lhints);
			
			// Column names
			vector<string> colnames;
			colnames.push_back("trk");
			colnames.push_back("type");
			colnames.push_back("p");
			colnames.push_back("theta");
			colnames.push_back("phi");
			colnames.push_back("z");
			colnames.push_back("chisq/Ndof");
			colnames.push_back("Ndof");
			colnames.push_back("FOM");
			colnames.push_back("cand");

			TGTextView *tview = new TGTextView(throwninfo, 700, 400);
			throwninfo->AddFrame(tview);
			
			tview->AddLine("Here is a line");
			//tview->AdjustWidth();
#if 0
			
			// Create TGListView objects to hold the thrown the thrown and reconstructed values
			TGListView *thrownview = new TGListView(throwninfo, 700, 400);
			throwninfo->AddFrame(thrownview, hints);

			TGLVContainer *lvcontainer = new TGLVContainer(thrownview);

			thrownview->SetHeaders(6);
			for(unsigned int i=0; i<colnames.size(); i++){
				if(i<6)thrownview->SetHeader(colnames[i].c_str(), 2, 1, i);
			}
			
			TGLVEntry *item = new TGLVEntry(lvcontainer, "One item", "TFolder");
			lvcontainer->AddItem(item);

			// Create a vertical frame for each column and insert the label as the first item
			for(unsigned int i=0; i<colnames.size(); i++){
				// create frames
				tf[colnames[i]] = new TGVerticalFrame(throwninfo);
				rf[colnames[i]] = new TGVerticalFrame(reconinfo);
				throwninfo->AddFrame(tf[colnames[i]], bhints);
				reconinfo->AddFrame(rf[colnames[i]], bhints);

				// create column labels
				string lab = colnames[i]+":";
				TGLabel *tl = new TGLabel(tf[colnames[i]], lab.c_str());
				TGLabel *rl = new TGLabel(rf[colnames[i]], lab.c_str());
				if(i<6)tf[colnames[i]]->AddFrame(tl, chints);
				rf[colnames[i]]->AddFrame(rl, chints);
				vector<TGLabel*> tv;
				vector<TGLabel*> rv;
				tv.push_back(tl);
				rv.push_back(rl);

				// Record the label object pointers for later use
				thrownlabs[colnames[i]] = tv;
				reconlabs[colnames[i]] = rv;
			}

			// Reconstruction factory and full list button
			TGVerticalFrame *vf = new TGVerticalFrame(reconinfo);
			reconinfo->AddFrame(vf, yhints);
			reconfactory = new TGComboBox(vf, "DTrackCandidate:", 0);
			reconfactory->Resize(160,20);
			vf->AddFrame(reconfactory, thints);
#endif

	//========== Close Button ===========
	TGTextButton *close = new TGTextButton(botframe,	"Close");
	botframe->AddFrame(close, chints);

	//&&&&&&&&&&&&&&&& Connections
	close->Connect("Clicked()","hdv_fulllistframe", this, "DoClose()");

	// Finish up and map the window
	SetWindowName("Hall-D Event Viewer Full Particle Listing");
	SetIconName("HDView");

	MapSubwindows();
	Resize(GetDefaultSize());
	
	tview->Update();
}

//-------------------
// DoClose
//-------------------
void hdv_fulllistframe::DoClose(void)
{
	UnmapWindow();
}

//-------------------
// UpdateTrackLabels
//-------------------
void hdv_fulllistframe::UpdateTrackLabels(vector<const DMCThrown*> &throwns, vector<const DKinematicData*> &trks)
{
#if 0
_DBG__;
	//if(!IsMapped())return; // don't bother updating unless we can actually be seen
_DBG__;

	// Make sure there are enough labels for all of the thrown particles
	map<string, vector<TGLabel*> >::iterator iter = thrownlabs.begin();
	for(; iter!=thrownlabs.end(); iter++){
_DBG__;
		vector<TGLabel*> &labels = iter->second;
		unsigned int max_rows = throwns.size();
		if(max_rows>20)max_rows=20;
		for(unsigned int row = labels.size(); row<max_rows+1; row++){
			stringstream ss;
			ss<<row;
			TGLabel *lab = new TGLabel(tf[iter->first], iter==thrownlabs.begin() ? ss.str().c_str():"--------");
			tf[iter->first]->AddFrame(lab, new TGLayoutHints(kLHintsTop|kLHintsCenterX, 2,2,0,0));
			labels.push_back(lab);
_DBG_<<"Adding row:"<<row<<endl;
		}
	}

	// Loop over thrown particles and fill in labels
	int ii=0;
	for(unsigned int i=0; i<throwns.size(); i++){
		const DMCThrown *trk = throwns[i];
		//if(trk->type==0)continue;
		int row = thrownlabs["trk"].size()-(ii++)-1;
		if(row<1)break;
		
		stringstream trkno, type, p, theta, phi, z;
		trkno<<setprecision(4)<<i+1;
		thrownlabs["trk"][row]->SetText(trkno.str().c_str());

		thrownlabs["type"][row]->SetText(ParticleType((Particle_t)trk->type));

		p<<setprecision(4)<<trk->momentum().Mag();
		thrownlabs["p"][row]->SetText(p.str().c_str());

		theta<<setprecision(4)<<trk->momentum().Theta()*TMath::RadToDeg();
		thrownlabs["theta"][row]->SetText(theta.str().c_str());

		double myphi = trk->momentum().Phi();
		if(myphi<0.0)myphi+=2.0*M_PI;
		phi<<setprecision(4)<<myphi;
		thrownlabs["phi"][row]->SetText(phi.str().c_str());

		z<<setprecision(4)<<trk->position().Z();
		thrownlabs["z"][row]->SetText(z.str().c_str());
	}

	// Loop over tracks and fill in labels
	for(unsigned int i=0; i<trks.size(); i++){
		const DKinematicData *trk = trks[i];
		int row = reconlabs["trk"].size()-i-1;
		if(row<1)break;
		
		stringstream trkno, type, p, theta, phi, z, chisq_per_dof, Ndof;
		stringstream fom;
		trkno<<setprecision(4)<<i+1;
		reconlabs["trk"][row]->SetText(trkno.str().c_str());
		
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
		reconlabs["type"][row]->SetText(type.str().c_str());

		p<<setprecision(4)<<trk->momentum().Mag();
		reconlabs["p"][row]->SetText(p.str().c_str());

		theta<<setprecision(4)<<trk->momentum().Theta()*TMath::RadToDeg();
		reconlabs["theta"][row]->SetText(theta.str().c_str());

		double myphi = trk->momentum().Phi();
		if(myphi<0.0)myphi+=2.0*M_PI;
		phi<<setprecision(4)<<myphi;
		reconlabs["phi"][row]->SetText(phi.str().c_str());

		z<<setprecision(4)<<trk->position().Z();
		reconlabs["z"][row]->SetText(z.str().c_str());

		// Get chisq and Ndof for DTrackTimeBased or DTrackWireBased objects
		const DTrackTimeBased *timetrack=dynamic_cast<const DTrackTimeBased*>(trk);
		const DTrackWireBased *track=dynamic_cast<const DTrackWireBased*>(trk);	
		const DTwoGammaFit *twogammafit=dynamic_cast<const DTwoGammaFit*>(trk);	
	
		if(timetrack){
			chisq_per_dof<<setprecision(4)<<timetrack->chisq/timetrack->Ndof;
			Ndof<<timetrack->Ndof;
			fom << timetrack->FOM;
		}else if(track){
			chisq_per_dof<<setprecision(4)<<track->chisq/track->Ndof;
			Ndof<<track->Ndof;
			fom << "N/A";
		}else if(twogammafit){
			chisq_per_dof<<setprecision(4)<<twogammafit->getChi2();
			Ndof<<twogammafit->getNdf();
			fom << twogammafit->getProb();
		}else{
			chisq_per_dof<<"N/A";
			Ndof<<"N/A";
			fom << "N/A";
		}
		reconlabs["chisq/Ndof"][row]->SetText(chisq_per_dof.str().c_str());
		reconlabs["Ndof"][row]->SetText(Ndof.str().c_str());
		reconlabs["FOM"][row]->SetText(fom.str().c_str());
	}
#endif

	MapSubwindows();
	Resize(GetDefaultSize());
}


