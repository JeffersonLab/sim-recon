// Author: David Lawrence  June 24, 2004
//
//
// DEventSourceFile methods
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "DEventSourceFile.h"

//----------------
// Constructor
//----------------
DEventSourceFile::DEventSourceFile(int narg, char *argv[]):DEventSource(narg,argv)
{
	/// Constructor for DEventSourceFile object
	fin = NULL;
	hddm_s = NULL;
	
}

//----------------
// Destructor
//----------------
DEventSourceFile::~DEventSourceFile()
{
	if(hddm_s)flush_s_HDDM(hddm_s, fin);
	if(fin)close_s_HDDM(fin);
}

//----------------
// Open
//----------------
derror_t DEventSourceFile::Open(char *source)
{
	/// Implementation of DEventSource virtual function
	/// Support only HDDM formatted files for now.
	Close();
	fin = open_s_HDDM(source);
	
	return fin==NULL ? ERROR_OPENING_EVENT_SOURCE:NOERROR;
}

//----------------
// Close
//----------------
derror_t DEventSourceFile::Close(void)
{
	/// Implementation of DEventSource virtual function
	/// Support only HDDM formatted files for now.
	if(fin)close_s_HDDM(fin);
	fin = NULL;
	
	return NOERROR;
}

//----------------
// GetEvent
//----------------
derror_t DEventSourceFile::GetEvent(void)
{
	/// Implementation of DEventSource virtual function
	/// Support only HDDM formatted files for now.
	if(hddm_s)flush_s_HDDM(hddm_s, fin);
	hddm_s = read_s_HDDM(fin);
	if(!hddm_s)return NO_MORE_EVENTS_IN_SOURCE;
	
	//----- Convert historic hddm format to hddm_containers_t -----

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(PE){
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
									
									hddm->CDChits->nrows++;
									hddm->CDChits->Grow();
									CDChit_t *CDChit = &hddm->CDChits->CDChit[hddm->CDChits->nrows-1];
									
									CDChit->pos.SetXYZ(r*cos(phi), r*sin(phi), z);
									CDChit->t = 0.0;
									CDChit->track = cdcPoints->in[m].track;
								}
							}
						}
					}
				}
			}
		}
	}
	
	return NOERROR;
}



