// $Id: DFCALShower_factory.cc 2001 2006-11-28 11:09:47Z mikornic $
//
//    File: DFCALShower_factory.cc
// Created: Tue Nov 28 11:57:50 EST 2006
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#include <iostream>
#include <fstream>
#include <math.h>
#include <DVector3.h>
using namespace std;

#include <JANA/JEvent.h>
using namespace jana;

#include "DFCALCluster_factory.h"
#include "DFCALHit.h"

// Used to sort hits by Energy
bool FCALHitsSort_C(const DFCALHit* const &thit1, const DFCALHit* const &thit2) {
	return thit1->E > thit2->E;
}


int ambiguous_events = 0;

//----------------
// Constructor
//----------------
DFCALCluster_factory::DFCALCluster_factory()
{
	// Set defaults
        MIN_CLUSTER_BLOCK_COUNT = 1;
        MIN_CLUSTER_SEED_ENERGY = 0.035; // GeV

	gPARMS->SetDefaultParameter("FCAL:MIN_CLUSTER_BLOCK_COUNT", MIN_CLUSTER_BLOCK_COUNT);
	gPARMS->SetDefaultParameter("FCAL:MIN_CLUSTER_SEED_ENERGY", MIN_CLUSTER_SEED_ENERGY);

}

//------------------
// brun
//------------------
jerror_t DFCALCluster_factory::brun(JEventLoop *eventLoop, int runnumber)
{
	
	vector<const DFCALGeometry*> fcalGeoms;
	eventLoop->Get(fcalGeoms);
	if(fcalGeoms.size()<1){
		_DBG_<<"fcalGeoms.size()<1 !!!"<<endl;
		throw JException("fcalGeoms.size()<1 !!!");
	}
	fcalGeom = fcalGeoms[0];

	return NOERROR;
}

//------------------
// evnt
//    Implementation of UConn LGD clusterizer (M. Kornicer)
//------------------
jerror_t DFCALCluster_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{

	vector<const DFCALHit*> fcalhits;
	eventLoop->Get(fcalhits);
	
	
	// Sort hits by energy
	sort(fcalhits.begin(), fcalhits.end(), FCALHitsSort_C);

	// fill user's hit list
        int nhits = 0;
        DFCALCluster::userhits_t* hits = (DFCALCluster::userhits_t*) malloc(sizeof(DFCALCluster::userhits_t)*FCAL_USER_HITS_MAX);

// Fill the structure that used to be used by clusterizers in Radphi 
	for (vector<const DFCALHit*>::const_iterator hit  = fcalhits.begin(); 
                                                     hit != fcalhits.end(); hit++ ) {
           if ( (**hit).E <  1e-6 ) continue;
           hits->hit[nhits].id = (**hit).id;
           hits->hit[nhits].x = (**hit).x;
           hits->hit[nhits].y = (**hit).y;
           hits->hit[nhits].E = (**hit).E; // adjust a hit energy either in the hit or photon factory
           hits->hit[nhits].t = (**hit).t;
           nhits++;
           
           if (nhits >= (int) FCAL_USER_HITS_MAX)  { 
              cout << "ERROR: FCALCluster_factory: number of hits " << nhits << " larger than " << FCAL_USER_HITS_MAX << endl;
              break;
           }

        }
        hits->nhits = nhits;

//        vector<int>  hitUsed(nhits, 0); // number of clusters that use a hit
                                          // or -1 for seed
        int hitUsed[nhits]; 
        for ( int i = 0; i < nhits; i++ ) { hitUsed[i] = 0; }
 
        const unsigned int max = 999;
	DFCALCluster* clusterList[max];
	unsigned int clusterCount = 0;
	int iter;
	for ( iter=0; iter < 99; iter++ ) {

          // 1. At beginning of iteration, recompute info for all clusters.
          //    If something changed, return all hits to the pool and repeat.

	   bool something_changed = false;
	   for ( unsigned int c = 0; c < clusterCount; c++ ) {
              //cout << " --------- Update iteration " << iter << endl;
      	      something_changed |= clusterList[c]->update( hits );
           }
      	   if (something_changed) {
              for ( unsigned int c = 0; c < clusterCount; c++ ) {
                  clusterList[c]->resetClusterHits();
              }
              // reset hits in factory also:
              for ( int h = 0; h < nhits; h++ ) hitUsed[h] = 0;
           }
           else if (iter > 0) {
              break;
           }
           
	   // 2. Look for blocks with energy large enough to require formation
	   //    of a new cluster, and assign them as cluster seeds.

	   for ( int ih = 0; ih < hits->nhits; ih++ ) {
	      double energy = hits->hit[ih].E;
              //cout << "hit: " << ih <<  " E: " << energy << endl;
	      if (energy < MIN_CLUSTER_SEED_ENERGY)
		 break;
	      double totalAllowed = 0;
	      for ( unsigned int c = 0; c < clusterCount; c++ ) {
		 totalAllowed += clusterList[c]->getEallowed(ih);
                 //cout << " totalAlowed from clust " << c <<  " is " << totalAllowed << endl;
                 
	      }
	      if (energy > totalAllowed) {
		 clusterList[clusterCount] = new DFCALCluster( hits->nhits );
                 hitUsed[ih] = -1;
		 clusterList[clusterCount]->addHit(ih,1.);
		 clusterList[clusterCount]->update( hits );
		 ++clusterCount;
	      }
	      else if (iter > 0) {
		 for ( unsigned int c = 0; c < clusterCount; c++ ) {
                    int nh_clust = clusterList[c]->getHits();
                    //cout << " Nhits " << nh_clust << " from clust " << c << " ? " << endl;
		    if ( nh_clust )
		       continue;
		    totalAllowed -= clusterList[c]->getEallowed(ih);
                    //cout << " else totalAlowed from " << c <<  " E: " << totalAllowed << endl;
		    if (energy > totalAllowed) {
                       if ( clusterList[c]->getHits() == 0  ) {
                          hitUsed[ih] = -1;
                       }
                       else {
                          ++hitUsed[ih];
                       }
		       clusterList[c]->addHit(ih,1.);
		       break;
		    }
	         }
	      }
	   }


	   // 3. Share all non-seed blocks among seeded clusters, where
	   //    any cluster shares a block if it expects at least 1 KeV in it.

	   for ( int ih = 0; ih < hits->nhits; ih++ ) {
              if ( hitUsed[ih]  < 0) // cannot share seed 
		 continue;
	      double totalExpected = 0;
              //cout << " Share hit: " << ih <<  " E: " << hits->hit[ih].E << endl;
	      for ( unsigned int c = 0; c < clusterCount; c++ ) {
		 if (clusterList[c]->getHits() > 0) {
		    totalExpected += clusterList[c]->getEexpected(ih);
		 }
	      }
              //cout << " totExpected " << totalExpected ;
	      for ( unsigned int c = 0; c < clusterCount; c++ ) {
		 if (clusterList[c]->getHits() > 0) {
		    double expected = clusterList[c]->getEexpected(ih);
                    //cout << " expected " << expected ;
		    if (expected > 1e-6) {
		       clusterList[c]->addHit(ih,expected/totalExpected);
                       ++hitUsed[ih];
		    }
		 }
	      }
	   }

        }

	if (iter == 99) {
	   ++ambiguous_events;
	}

        for ( unsigned int c = 0; c < clusterCount; c++) {
           unsigned int blockCount = clusterList[c]->getHits();
           //cout << " Blocks " << blockCount << endl;
	   if (blockCount < MIN_CLUSTER_BLOCK_COUNT) {
              delete clusterList[c];
	      continue;
	   }
	   else {

              clusterList[c]->saveHits( hits );

              _data.push_back( clusterList[c] );

	   }
        }
  

        if (hits) {
           free(hits);
           hits=0;
        }

	return NOERROR;

}

