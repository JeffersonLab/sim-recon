
#include "hdview.h"

#include <TVector3.h>
#include <TMarker.h>
#include <TTUBE.h>
#include <TEllipse.h>

TEllipse *ellipse[10];
int Nellipse = 0;

TMarker *cdcMarkers[1000];
int NcdcMarkers=0;

//-------------------
// hdv_getevent
//-------------------
derror_t hdv_getevent(void)
{
	// Read in next event. The MyProcessor object's event
	// method will automatically be called which will convert
	// the hddm to more succienct variables used in
	// hdv_drawevent() below.
	derror_t err = eventloop->OneEvent();
	if(err!=NOERROR)return err;
	
	cout<<"Event Number:"<<myproc->eventNo<<endl;
	cout<<"Ncdchits = "<<myproc->Ncdchits<<endl;
	cout<<endl;

	// Do this here so it is doesn't slow down "orbit"
	// Draw an ellipse for each track
	for(int i=0;i<Nellipse;i++)delete ellipse[i];
	Nellipse = 0;
	s_Cdc_track_t *cdc_track = myproc->hddm->cdc_tracks->in;
	for(int i=0;i<myproc->hddm->cdc_tracks->mult;i++, cdc_track++){
		float B=-2.0*0.61; // The 0.61 is empirical
		float hbarc = 197.326;
		float r0 = -cdc_track->p*sin(cdc_track->theta)*hbarc/B/fabs(cdc_track->q);
		float phi = cdc_track->phi - cdc_track->q*M_PI_2;
		float x0 = r0*cos(phi);
		float y0 = r0*sin(phi);

		ellipse[Nellipse] = new TEllipse(x0, y0, r0, r0);
		ellipse[Nellipse++]->Draw();
		maincanvas->Update();
		if(Nellipse>=10)break;
	}

	return NOERROR;
}

//-------------------
// hdv_drawevent
//-------------------
derror_t hdv_drawevent(void)
{
	// --------- CDC ----------

	// Delete old markers
	for(int i=0;i<NcdcMarkers;i++)delete cdcMarkers[i];
	NcdcMarkers = 0;

	// Create markers for all CDC hits
	TVector3 *v = myproc->cdchits;
	for(int i=0;i<myproc->Ncdchits;i++ ,v++){
		//cdcMarkers[NcdcMarkers] = new TMarker(0.5+v->x()/500.0, 0.5+v->y()/300.0, 20);
		cdcMarkers[NcdcMarkers] = new TMarker(v->x(), v->y(), 20);
		cdcMarkers[NcdcMarkers]->SetMarkerColor(myproc->cdchit_tracks[i]);
		cdcMarkers[NcdcMarkers]->SetMarkerSize(exp(v->z()/500.0)/2.0);
		cdcMarkers[NcdcMarkers++]->Draw();
	}
	

	// Update canvas
	maincanvas->Update();

	return NOERROR;
}

//-------------------
// RotatePoints
//-------------------
void RotatePoints(int how)
{
	TVector3 *v = myproc->cdchits;
	float angle = 0.1;
	
	for(int i=0;i<NcdcMarkers;i++, v++){
		switch(how){
			case 3: v->RotateX(0.1);	break;
			case 4: v->RotateY(0.1);	break;
			case 5: v->RotateZ(0.1);	break;
		}
	}
	
	hdv_drawevent();
}

//-------------------
// Orbit
//-------------------
void Orbit(void)
{
	
	float angle = 0.01;
	float tilt = 0.4;
	
	TVector3 *v;
	for(int j=0;j<tilt/angle;j++){
		v = myproc->cdchits;
		for(int i=0;i<NcdcMarkers;i++, v++){
			v->RotateX(angle);
		}
		hdv_drawevent();
	}
	
	for(int j=0;j<2.0*3.14159265/angle;j++){
		v = myproc->cdchits;
		for(int i=0;i<NcdcMarkers;i++, v++){
			v->RotateX(-tilt);
			v->RotateY(angle);
			v->RotateX(tilt);
		}
		hdv_drawevent();
	}
	
	for(int j=0;j<tilt/angle;j++){
		v = myproc->cdchits;
		for(int i=0;i<NcdcMarkers;i++, v++){
			v->RotateX(-angle);
		}
		hdv_drawevent();
	}
}


