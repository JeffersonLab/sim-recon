// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include "MyProcessor.h"
#include "hddm_s.h"
#include "cdc_fithelix.h"

int qsort_cdc(const void* arg1, const void* arg2);

//------------------------------------------------------------------
// evnt 
//------------------------------------------------------------------
derror_t MyProcessor::evnt(int eventnumber)
{
	eventNo = eventnumber;
	Ncdchits = 0;

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm->physicsEvents;
	if(!PE)return NOERROR;
	for(int i=0; i<PE->mult; i++){
	
		s_Rings_t *rings=NULL;
		s_HitView_t *HV = PE->in[i].hitView;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->centralDC)
				rings = PE->in[i].hitView->centralDC->rings;
		if(rings){
			for(int j=0;j<rings->mult;j++){
				float radius = rings->in[j].radius;
				s_Straws_t *straws = rings->in[j].straws;
				if(straws){
					for(int k=0;k<straws->mult;k++){
						float phim = straws->in[k].phim;
						s_CdcPoints_t *cdcPoints = straws->in[k].cdcPoints;
						if(cdcPoints){						
							for(int m=0;m<cdcPoints->mult;m++){
								float	r = cdcPoints->in[m].r;
								float phi = cdcPoints->in[m].phi;
								float z = cdcPoints->in[m].z;
								if(Ncdchits>=1000)continue;
								cdchit_tracks[Ncdchits] = cdcPoints->in[m].track;
								cdchits[Ncdchits++].SetXYZ(r*cos(phi), r*sin(phi), z);								
							}
						}
					}
				}
			}
		}
		
		// Tracks are interwoven at this point. Sort them by track number
		// and z value. Use an index array since this is easiest for now.
		int index[Ncdchits];
		for(int j=0;j<Ncdchits;j++)index[j] = j;
		qsort(index, Ncdchits, sizeof(int), qsort_cdc);
		
		// Re-arrange actual arrays according to index
		TVector3 tmpv[Ncdchits];
		int tmpi[Ncdchits];
		for(int j=0;j<Ncdchits;j++){
			int k = index[j];
			tmpv[j] = cdchits[k];
			tmpi[j] = cdchit_tracks[k];
		}
		for(int j=0;j<Ncdchits;j++){
			cdchits[j] = tmpv[j];
			cdchit_tracks[j] = tmpi[j];
		}

		// Fit CDC tracks
		int track_start = 0;
		for(int j=1;j<Ncdchits;j++){
			if(cdchit_tracks[j]!= cdchit_tracks[j-1]){
				cout<<__FILE__<<":"<<__LINE__<<" track="<<cdchit_tracks[j-1]<<" ";
				cdc_fithelix(&cdchits[track_start], j-track_start);
				track_start = j;
			}
		}
		if(Ncdchits-track_start > 1){
			cout<<__FILE__<<":"<<__LINE__<<" track="<<cdchit_tracks[Ncdchits-1]<<" ";
			cdc_fithelix(&cdchits[track_start], Ncdchits-track_start);
		}
		
		s_Reactions_t *reactions = PE->in[i].reactions;
		if(reactions){
			for(int j=0;j<reactions->mult;j++){
				s_Vertices_t *vertices = reactions->in[j].vertices;
				if(!vertices)continue;
				for(int k=0;k<vertices->mult;k++){
					s_Products_t *products = vertices->in[j].products;
					if(!products)continue;
					for(int m=0;m<products->mult;m++){
						if(products->in[m].momentum){
						}
					}
				}
			}
		}
		
	}

	return NOERROR;
}

//------------------------------------------------------------------
// qsort_cdc 
//------------------------------------------------------------------
int qsort_cdc(const void* arg1, const void* arg2)
{
	int *a = (int*)arg1;
	int *b = (int*)arg2;
	
	extern MyProcessor *myproc;
	
	if(myproc->cdchit_tracks[*a] != myproc->cdchit_tracks[*b]){
		return myproc->cdchit_tracks[*a] - myproc->cdchit_tracks[*b];
	}
	
	if(myproc->cdchits[*a].z() == myproc->cdchits[*b].z())return 0;
	return myproc->cdchits[*b].z() > myproc->cdchits[*a].z() ? 1:-1;
}

