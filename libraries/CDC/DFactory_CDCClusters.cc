// $Id$

#include "DEvent.h"
#include "DFactory_CDCClusters.h"


//-------------------
// evnt
//-------------------
derror_t DFactory_CDCClusters::evnt(int eventnumber)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	// I assume that the cdcPoints array matches up
	// one to one with the hits array. The hits array
	// should be one to one with the entries in CDCHits
	CDCHit_t *cdchit = (CDCHit_t*)event->Get("CDCHits")->first();
	
	int ntracks=0;
	int tracknumber[100];
	
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
						s_Hits_t *hits = straws->in[k].hits;
						s_CdcPoints_t *cdcPoints = straws->in[k].cdcPoints;
						if(hits){
							
							for(int m=0;m<hits->mult;m++){
								
								// Loop over existing clusters to see if this track
								// already exists
								int clusterindex=-1;
								for(int n=0;n<ntracks;n++){
									if(tracknumber[n] == cdcPoints->in[m].track){
										clusterindex = n;
										break;
									}
								}
								
								// If this doesn't have a cluster, then add one
								if(clusterindex<0){
									CDCCluster_t *cdccluster = (CDCCluster_t*)_data->Add();
									cdccluster->nhits = 0;
									tracknumber[ntracks] = cdcPoints->in[m].track;
									clusterindex = ntracks++;
								}
								
								// Get pointer to cluster row
								CDCCluster_t *cdccluster = (CDCCluster_t*)_data->index(clusterindex);
								
								// Add this hit to the cluster
								cdccluster->hits[cdccluster->nhits++] = cdchit++;
							}
						}
					}
				}
			}
		}
	}

	
	return NOERROR;
}


