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
	
	// Fill in hddm_containers pointers to classic hddm tree
	// (Number of rows in all containers should now be zero (set in DEventSource.cc))
	
	
	
	//----- Convert historic hddm format to hddm_containers_t -----
	// centralDC
	// forwardDC
	// startCntr
	// barrelEMcal
	// Cerenkov
	// forwardTOF
	// forwardEMcal

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(PE){
		hddm->PhysicsEvents->Set(PE->mult, (void**)&PE->in);
		for(int i=0; i<PE->mult; i++){
		
			

			// ------------ CdcPoints, Hits --------------
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
							
							//----- CDC "Hits"
							s_Hits_t *hits = straws->in[k].hits;
							if(hits){

							}							
							//----- CDC "Points"
							s_CdcPoints_t *cdcPoints = straws->in[k].cdcPoints;
							if(cdcPoints){

							}
						}
					}
				}
			}
		}
	}
	
	return NOERROR;
}



