// $Id$

#include "DEvent.h"
#include "DFactory_CDCClusters.h"
#include "DQuickFit.h"


//-------------------
// evnt
//-------------------
derror_t DFactory_CDCClusters::evnt(int eventnumber)
{
	// Use the MC cheat codes in place of cluster finding for now.
	derror_t err = CopyFromMCCheatCodes();

	// Do a quick fit to the clusters
	CDCCluster_t *cdccluster = (CDCCluster_t*)_data->first();
	for(int i=0; i<_data->nrows; i++, cdccluster++){
		if(cdccluster->nhits<2)continue;

		DQuickFit *qf = new DQuickFit();
		CDCHit_t **hits = cdccluster->hits;
		for(int j=0;j<cdccluster->nhits; j++, hits++){
			qf->AddHit((*hits)->radius, (*hits)->phim);
		}
		
		// Do the fit and copy the results into the cluster
		qf->FitCircle();
		qf->CopyToFitParms(&cdccluster->fit);
				
		delete qf;
	}

	return err;
}

//-------------------
// CopyFromMCCheatCodes
//-------------------
derror_t DFactory_CDCClusters::CopyFromMCCheatCodes(void)
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

//------------
// Print
//------------
derror_t DFactory_CDCClusters::Print(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(!_data)return NOERROR;
	if(_data->nrows<=0)return NOERROR; // don't print anything if we have no data!

	cout<<name<<endl;
	cout<<"------------------------------------------------------------------------------"<<endl;
	cout<<"row: radius(cm): phim(rad):   dE(MeV):   t(ns):     x0:    y0: p_trans(GeV/c):"<<endl;
	cout<<endl;
	
	CDCCluster_t *cdccluster = (CDCCluster_t*)_data->first();
	for(int i=0; i<_data->nrows; i++, cdccluster++){
		char str[80];
		memset(str,' ',80);
		str[79] = 0;

		char num[32];
		sprintf(num, "%d", i);
		strncpy(&str[3-strlen(num)], num, strlen(num));

		CDCHit_t **cdchit = cdccluster->hits;
		for(int j=0;j<cdccluster->nhits; j++, cdchit++){
			if(j!=0){
				memset(str,' ',80);
				str[79] = 0;
			}
			sprintf(num, "%3.1f", (*cdchit)->radius);
			strncpy(&str[15-strlen(num)], num, strlen(num));
			sprintf(num, "%1.3f", (*cdchit)->phim);
			strncpy(&str[26-strlen(num)], num, strlen(num));
			sprintf(num, "%2.3f", (*cdchit)->dE);
			strncpy(&str[37-strlen(num)], num, strlen(num));
			sprintf(num, "%4.0f", (*cdchit)->t);
			strncpy(&str[46-strlen(num)], num, strlen(num));
			if(j==0){
				sprintf(num, "%3.2f", cdccluster->fit.x0);
				strncpy(&str[54-strlen(num)], num, strlen(num));
				sprintf(num, "%3.2f", cdccluster->fit.y0);
				strncpy(&str[61-strlen(num)], num, strlen(num));
				sprintf(num, "%2.3f", cdccluster->fit.p_trans);
				strncpy(&str[77-strlen(num)], num, strlen(num));
			}
			cout<<str<<endl;
		}
	}
	cout<<endl;
}


