
#include "hdview.h"

#include <TVector3.h>
#include <TMarker.h>
#include <TTUBE.h>


TMarker *cdcMarkers[1000];
int NcdcMarkers=0;

//-------------------
// hdv_getevent
//-------------------
derror_t hdv_getevent(void)
{
	// Read in next event. 
	derror_t err = eventloop->OneEvent();
	if(err!=NOERROR)return err;
		

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


