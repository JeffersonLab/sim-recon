// $Id$


#include "DEvent.h"
#include "DFactory_FDCClusters.h"


//-------------------
// evnt
//-------------------
derror_t DFactory_FDCClusters::evnt(int eventnumber)
{

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return NOERROR;
	
	// Since the FDCHits data have anode and cathode hits intermixed,
	// and cheat codes are only available for the anode hits, we'll
	// assume each cdcPoint has a one to one correspondence with the
	// anode hits only. In other words, we'll loop over cdcPoints and
	// assume each one matches the next hit in the CDCHits container
	// (we skip over the cathode hits).
	FDCHit_t *fdchit = (FDCHit_t*)event->Get("FDCHits")->first();
	
	int ntracks=0;
	int tracknumber[100];
	
	for(int i=0; i<PE->mult; i++){
		s_Chambers_t *chambers = NULL;
		s_HitView_t *HV = PE->in[i].hitView;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->forwardDC)
				chambers = PE->in[i].hitView->forwardDC->chambers;
		if(!chambers)continue;
		
		for(int j=0;j<chambers->mult;j++){

			s_AnodePlanes_t *anodeplanes = chambers->in[j].anodePlanes;
			if(anodeplanes){
			
				for(int k=0;k<anodeplanes->mult;k++){
					s_Wires_t *wires = anodeplanes->in[k].wires;
					if(!wires)continue;
				
					for(int m=0;m<wires->mult;m++){
						s_Hits_t *hits = wires->in[m].hits;
						s_FdcPoints_t *fdcPoints = wires->in[k].fdcPoints;
						if(!hits)continue;
						for(int n=0;n<hits->mult;n++){
						
							// Loop over existing clusters to see if this track
							// already exists
							int clusterindex=-1;
							for(int p=0;p<ntracks;p++){
								if(tracknumber[p] == fdcPoints->in[n].track){
									clusterindex = p;
									break;
								}
							}
								
							// If this doesn't have a cluster, then add one
							if(clusterindex<0){
								FDCCluster_t *fdccluster = (FDCCluster_t*)_data->Add();
								fdccluster->nhits = 0;
								tracknumber[ntracks] = fdcPoints->in[n].track;
								clusterindex = ntracks++;
							}
								
							// Get pointer to cluster row
							FDCCluster_t *fdccluster = (FDCCluster_t*)_data->index(clusterindex);
								
							// Add this hit to the cluster
							fdccluster->hits[fdccluster->nhits++] = fdchit++;
						}
					}
				}
			}

			s_CathodePlanes_t *cathodeplanes = chambers->in[j].cathodePlanes;
			if(cathodeplanes){
			
				for(int k=0;k<cathodeplanes->mult;k++){
					float tau = cathodeplanes->in[k].tau;
					float z = cathodeplanes->in[k].z;
					s_Strips_t *strips = cathodeplanes->in[k].strips;
					if(!strips)continue;
				
					for(int m=0;m<strips->mult;m++){
						// Just skip cathode hits
						fdchit++;
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
derror_t DFactory_FDCClusters::Print(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(!_data)return NOERROR;
	if(_data->nrows<=0)return NOERROR; // don't print anything if we have no data!

	cout<<name<<endl;
	cout<<"---------------------------------------"<<endl;
	cout<<"row: layer: module: tau(rad):    z(cm):  u(cm):  dE(MeV):   t(ns):   type:"<<endl;
	cout<<endl;
	
	FDCCluster_t *fdccluster = (FDCCluster_t*)_data->first();
	for(int i=0; i<_data->nrows; i++, fdccluster++){
		char str[80];
		memset(str,' ',80);
		str[79] = 0;

		char num[32];
		sprintf(num, "%d", i);
		strncpy(&str[3-strlen(num)], num, strlen(num));

		FDCHit_t **fdchit = fdccluster->hits;
		for(int j=0;j<fdccluster->nhits; j++, fdchit++){
			if(j!=0){
				memset(str,' ',80);
				str[79] = 0;
			}
			sprintf(num, "%d", (*fdchit)->layer);
			strncpy(&str[10-strlen(num)], num, strlen(num));
			sprintf(num, "%d", (*fdchit)->module);
			strncpy(&str[18-strlen(num)], num, strlen(num));
			sprintf(num, "%3.1f", (*fdchit)->tau);
			strncpy(&str[28-strlen(num)], num, strlen(num));
			sprintf(num, "%3.1f", (*fdchit)->z);
			strncpy(&str[38-strlen(num)], num, strlen(num));
			sprintf(num, "%2.3f", (*fdchit)->u);
			strncpy(&str[46-strlen(num)], num, strlen(num));
			if(!(*fdchit)->type){
				sprintf(num, "%1.3f", (*fdchit)->dE*1000.0);
				strncpy(&str[56-strlen(num)], num, strlen(num));
				sprintf(num, "%4.0f", (*fdchit)->t);
				strncpy(&str[65-strlen(num)], num, strlen(num));
			}
			sprintf(num, "%s", (*fdchit)->type ? "cathode":"anode");
			strncpy(&str[73-strlen(num)], num, strlen(num));

			cout<<str<<endl;
		}
	}
	cout<<endl;
}


