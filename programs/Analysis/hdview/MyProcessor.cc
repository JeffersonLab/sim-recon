// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>

#include "hdview.h"
#include "hdv_mainframe.h"
#include "MyProcessor.h"
#include "hddm.h"

extern TCanvas *maincanvas;
extern DEventLoop *eventloop;


//------------------------------------------------------------------
// evnt 
//------------------------------------------------------------------
derror_t MyProcessor::evnt(int eventnumber)
{

	// Delete objects from previous event
	for(int i=0;i<Nellipse;i++)delete ellipse[i];
	for(int i=0;i<Nhelix;i++)delete helix[i];
	for(int i=0;i<Ncdctracks;i++)delete cdchits3D[i];
	Nellipse = Nhelix = Ncdctracks = 0;
	
	// Draw Detector
	extern TGeoVolume *MOTHER;
	MOTHER->Draw();

	// Need to get Rotation matrix of detector. For some reason, ROOT
	// uses two independant rotation classes.
	extern TGeoCombiTrans *MotherRotTrans;
	double theta1,phi1,theta2,phi2,theta3,phi3;
	MotherRotTrans->GetRotation()->GetAngles(theta1,phi1,theta2,phi2,theta3,phi3);
	TRotation rot;
	TVector3 xx,yy,zz;
	xx.SetPtThetaPhi(1.0, theta1, phi1);
	yy.SetPtThetaPhi(1.0, theta2, phi2);
	zz.SetPtThetaPhi(1.0, theta3, phi3);
	rot.SetXAxis(xx);
	rot.SetYAxis(yy);
	rot.SetZAxis(zz);
	TLorentzRotation lrot(rot);

	// Create 3D polymarkers from CDChits
	CDChit_t *CDChit = hddm->CDChits->CDChit;
	CDCtrack_t *CDCtrack = hddm->CDCtracks->CDCtrack;
	int last_track = -1;
	for(Ncdchits=0; Ncdchits<hddm->CDChits->nrows; Ncdchits++, CDChit++){
		
		TVector3 *v = &CDChit->pos;
		if(CDChit->track != last_track){
		
			// create new polymarker for track
			last_track = CDChit->track;
			cdchits3D[Ncdctracks] = new TPolyMarker3D();
			
			// create new helix for track
			if(Nhelix<10){
				// To find the vertex z position, we must extrapolate from
				// the current point back to the beamline along the helix.
				// First, find dphi/dz. Then find deltaphi between the beamline
				// and the current cdchit.
				CDCtrack_t *trk = &CDCtrack[Ncdctracks];
				float r0 = sqrt(trk->x0*trk->x0 + trk->y0*trk->y0);
				float dphidz = r0*tan(trk->dir.Theta());
				TVector3 center(trk->x0,trk->y0, 0);
				TVector3 pt(v->x()-trk->x0, v->y()-trk->y0, 0);
				float dphi = pt.Phi() - center.Phi();
				float z0 = 70. - (dphi/dphidz);
				z0=70.0;
				cout <<" center.Phi()="<<center.Phi()<<" v->Phi()="<<v->Phi()<<" ";
				cout << "z0="<<z0<<" dphidz="<<dphidz<<" dphi="<<dphi<<endl;
			
				TLorentzVector p = trk->p;
				float B=-2.0*0.593/197.326; // empirical
				helix[Nhelix] = new THelix(0, 0, z0, p.X(), p.Y(), p.Z(), CDCtrack[Ncdctracks].q*B);
				helix[Nhelix]->SetRange(z0,400);
				helix[Nhelix]->SetLineColor(Ncdctracks+2);
				helix[Nhelix++]->Draw();
			}
			Ncdctracks++;
		}
		
		// Add hit to polymarker
		cdchits3D[Ncdctracks-1]->SetNextPoint(v->x(), v->y(), v->z());
		//cout<<__FILE__<<":"<<__LINE__<<" x="<<v->x()<<" y="<<v->y()<<" z="<<v->z()<<endl;
	}
	for(int i=0;i<Ncdctracks;i++){
		cout<<__FILE__<<":"<<__LINE__<<" nhits="<<cdchits3D[i]->GetLastPoint()<<endl;
		cdchits3D[i]->SetMarkerColor(i+2);
		cdchits3D[i]->SetMarkerSize(0.4);
		cdchits3D[i]->SetMarkerStyle(20);
		cdchits3D[i]->Draw();
	}
	
	cout<<"Next Event"<<endl;
	//cout<<"\t Ncdctracks="<<hddm->CDCtracks->nrows<<endl;
	
	return NOERROR;
}

