// $Id: DFCALShower_factory.cc 2001 2006-11-28 11:09:47Z mikornic $
//
//    File: DFCALShower_factory.cc
// Created: Tue Nov 28 11:57:50 EST 2006
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#include <math.h>
#include <DVector3.h>

#include "DFCALCluster_factory.h"
#include "DFCALCluster.h"
#include "DFCALHit.h"
#include "JANA/JEvent.h"

#ifndef SQR
# define SQR(x) (x)*(x)
#endif

userhits_t* hits = NULL;

/* uncomment below to enable log-weighting method for cluster centroid
   (imported from IU clusterizer, it was default method for Radphi) - MK */
// #define LOG2_WEIGHTING 

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
        MIN_CLUSTER_BLOCK_COUNT = 2;
        MIN_CLUSTER_SEED_ENERGY = 0.035; // GeV

	gPARMS->SetDefaultParameter("FCAL:MIN_CLUSTER_BLOCK_COUNT", MIN_CLUSTER_BLOCK_COUNT);
	gPARMS->SetDefaultParameter("FCAL:MIN_CLUSTER_SEED_ENERGY", MIN_CLUSTER_SEED_ENERGY);

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
        if (hits == 0) {
           hits = (userhits_t*) malloc(sizeof(userhits_t)*FCAL_USER_HITS_MAX);
        }

// Fill the structure that used to be used by clusterizers in Radphi 
	for (vector<const DFCALHit*>::const_iterator hit  = fcalhits.begin(); 
                                                     hit != fcalhits.end(); hit++ ) {
           hits->hit[nhits].col = (**hit).column;
           hits->hit[nhits].row = (**hit).row;
           hits->hit[nhits].x = (**hit).x;
           hits->hit[nhits].y = (**hit).y;
           hits->hit[nhits].E = (**hit).E; // adjust a hit energy either in the hit or photon factory
           hits->hit[nhits].t = (**hit).t;
           nhits++;
           
           if (nhits >= (int) FCAL_USER_HITS_MAX)  { 
              cout << "nhits gt nhmax" << endl;
              break;
           }

        }
        hits->nhits = nhits;
        
        const unsigned int max = 999;
	DFCALCluster::setHitlist(hits);
	DFCALCluster* clusterList[max];
	unsigned int clusterCount = 0;
	int iter;
	for ( iter=0; iter < 99; iter++ ) {

          // 1. At beginning of iteration, recompute info for all clusters.
          //    If something changed, return all hits to the pool and repeat.

	   bool something_changed = false;
	   for ( unsigned int c = 0; c < clusterCount; c++ ) {
      	      something_changed |= clusterList[c]->update();
           }
      	   if (something_changed) {
              for ( unsigned int c = 0; c < clusterCount; c++ ) {
                  clusterList[c]->resetHits();
              }
           }
           else if (iter > 0) {
              break;
           }
           
	   // 2. Look for blocks with energy large enough to require formation
	   //    of a new cluster, and assign them as cluster seeds.

	   for ( int ih = 0; ih < hits->nhits; ih++ ) {
	      double energy = hits->hit[ih].E;
	      if (energy < MIN_CLUSTER_SEED_ENERGY)
		 break;
	      double totalAllowed = 0;
	      for ( unsigned int c = 0; c < clusterCount; c++ ) {
		 totalAllowed += clusterList[c]->getEallowed(ih);
	      }
	      if (energy > totalAllowed) {
		 clusterList[clusterCount] = new DFCALCluster();
		 clusterList[clusterCount]->addHit(ih,1.);
		 clusterList[clusterCount]->update();
		 ++clusterCount;
	      }
	      else if (iter > 0) {
		 for ( unsigned int c = 0; c < clusterCount; c++ ) {
		    if (clusterList[c]->getHits())
		       continue;
		    totalAllowed -= clusterList[c]->getEallowed(ih);
		    if (energy > totalAllowed) {
		       clusterList[c]->addHit(ih,1.);
		       break;
		    }
	         }
	      }
	   }


	   // 3. Share all non-seed blocks among seeded clusters, where
	   //    any cluster shares a block if it expects at least 1 KeV in it.

	   for ( int ih = 0; ih < hits->nhits; ih++ ) {
	      if (DFCALCluster::getUsed(ih) < 0)
		 continue;
	      double totalExpected = 0;
	      for ( unsigned int c = 0; c < clusterCount; c++ ) {
		 if (clusterList[c]->getHits() > 0) {
		    totalExpected += clusterList[c]->getEexpected(ih);
		 }
	      }
	      for ( unsigned int c = 0; c < clusterCount; c++ ) {
		 if (clusterList[c]->getHits() > 0) {
		    double expected = clusterList[c]->getEexpected(ih);
		    if (expected > 1e-6) {
		       clusterList[c]->addHit(ih,expected/totalExpected);
		    }
		 }
	      }
	   }

        }

	if (iter == 99) {
	   ++ambiguous_events;
	}

        for ( unsigned int c = 0; c < clusterCount; c++) {
           int hitlist[FCAL_USER_HITS_MAX];
           unsigned int blockCount = clusterList[c]->getHits(hitlist,hits->nhits);
	   if (blockCount < MIN_CLUSTER_BLOCK_COUNT) {
              delete clusterList[c];
	      continue;
	   }
	   else {
              _data.push_back( clusterList[c] );
	   }
        }
  

        if (hits) {
           free(hits);
           hits=0;
        }

	return NOERROR;

}


//------------------
// toString
//------------------
const string DFCALCluster_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("index:   x(cm):   y(cm):   E(GeV):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DFCALCluster *fcalclust = _data[i];

		printnewrow();
		printcol("%d",	i);
		printcol("%3.1f", fcalclust->getCentroid().x());
		printcol("%3.1f", fcalclust->getCentroid().y());
		printcol("%2.3f", fcalclust->getEnergy());

		printrow();
	}

	return _table;
}

const userhits_t* DFCALCluster::fHitlist = 0;
//const DFCALHit* DFCALCluster::fHitlist = 0;
int* DFCALCluster::fHitused = 0;

DFCALCluster::DFCALCluster()
{
   fEnergy = 0;
   fEmax = 0;
   fCentroid.SetX(0);
   fCentroid.SetY(0);
   fCentroid.SetZ(0);
   fRMS = 0;
   fRMS_u = 0;
   fRMS_v = 0;
   fNhits = 0;

   if (fHitlist && fHitlist->nhits > 0) {
      fHit = new int[fHitlist->nhits];
      fHitf = new double[fHitlist->nhits];
      fEallowed = new double[fHitlist->nhits];
      fEexpected = new double[fHitlist->nhits];
   }
   else {
      fHit = 0;
      fHitf = 0;
      fEallowed = 0;
      fEexpected = 0;
   }
}

DFCALCluster::~DFCALCluster()
{
   if (fHit) delete [] fHit;
   if (fHitf) delete [] fHitf;
   if (fEallowed) delete [] fEallowed;
   if (fEexpected) delete [] fEexpected;
}

 void DFCALCluster::setHitlist( const userhits_t* const hits)
//void DFCALCluster::setHitlist( const DFCALHits* const hits, const )
{
   fHitlist = hits;
   if (fHitused) {
      delete [] fHitused;
      fHitused = 0;
   }
   if (hits->nhits > 0) {
      fHitused = new int[hits->nhits];
      for ( int h = 0; h < hits->nhits; h++ ) {
         fHitused[h] = 0;
      }
   }
}

/* not needed to solve memory leak
void DFCALCluster::unsetHitlist()
{
   if (fHitused) {
      delete [] fHitused;
      fHitused = 0;
   }
   if (fHitlist) {
       delete [] fHitlist;
       fHitlist = 0;
   }   
}
*/

int DFCALCluster::addHit(const int ihit, const double frac)
{
   if (ihit >= 0 && ihit < fHitlist->nhits && fNhits < fHitlist->nhits) {
      fHit[fNhits] = ihit;
      fHitf[fNhits] = frac;
      if (fNhits == 0) {
        fHitused[ihit] = -1;    // special used code for cluster seeds
      }
      else {
        ++fHitused[ihit];       // otherwise, just count owning clusters
      }
      ++fNhits;
      return 0;
   }
   else {
      return 1;
   }
}

void DFCALCluster::resetHits()
{
   if (fNhits) {
      fHitused[fHit[0]] = 0;
      for (int h = 1; h < fNhits; h++) {
         --fHitused[fHit[h]];
      }
      fNhits = 0;
   }
}

bool DFCALCluster::update()
{
   double energy = 0;
   for ( int h = 0; h < fNhits; h++ ) {
      int ih = fHit[h];
      double frac = fHitf[h];
      energy += fHitlist->hit[ih].E*frac;
   }
   double Emax=0;
   if (fNhits > 0) Emax = fHitlist->hit[fHit[0]].E;

   DVector3 centroid;
   centroid.SetX(0);
   centroid.SetY(0);
   centroid.SetZ(0);
   double xc=0;
   double yc=0;
#ifdef LOG2_WEIGHTING
  /* This complicated centroid algorithm was copied from
   * Scott Teige's lgdClusterIU.c code -- don't ask [rtj]
   */
   double weight = 0.0;
   double weightSum = 0.0;
   double centerWeight = 0.0;
   double neighborMaxWeight = 0.0;
   double logFraction = exp(-0.23*(energy));
   double currentOffset = log(energy)/7.0 + 3.7;
   for (int h = 0; h < fNhits; h++) {
      int ih = fHit[h];
      double frac = fHitf[h];

   /*
    * Find the weight - offset determines minimum energy of block
    * to be used in weighting, since negative weights are thrown out
    */

      weight = currentOffset + log(fHitlist->hit[ih].E*frac/energy);
      if (h == 0) {
         centerWeight = weight;
      }
      else {
         neighborMaxWeight =
               (neighborMaxWeight < weight)? weight:neighborMaxWeight;
      }
      if (weight > 0) {
         xc  += fHitlist->hit[ih].x*weight;
         yc  += fHitlist->hit[ih].y*weight;
         weightSum += weight;
      }
   }
   /*
    * Now patch up the center block's weight if it's got a neighbor
    * in the cluster that had positive weight
    */
   if (neighborMaxWeight > 0) {
      xc += (logFraction-1)*(centerWeight-neighborMaxWeight)
                             *fHitlist->hit[0].x;
      yc += (logFraction-1)*(centerWeight-neighborMaxWeight)
                             *fHitlist->hit[0].y;
      weightSum += (logFraction-1)*(centerWeight-neighborMaxWeight);
   }
   centroid.SetX(xc/weightSum);
   centroid.SetY(yc/weightSum);
#else
   for (int h = 0; h < fNhits; h++) {
      int ih = fHit[h];
      double frac = fHitf[h];
      xc += fHitlist->hit[ih].x*(fHitlist->hit[ih].E*frac);
      yc += fHitlist->hit[ih].y*(fHitlist->hit[ih].E*frac);

   }
   centroid.SetX(xc/energy);
   centroid.SetY(yc/energy);
#endif

   double MOM1x = 0;
   double MOM1y = 0;
   double MOM2 = 0;
   double MOM1u = 0;
   double MOM2u = 0;
   double MOM1v = 0;
   double MOM2v = 0;
   for ( int h = 0; h < fNhits; h++ ) {
      int ih = fHit[h];
      double frac = fHitf[h];
      double x = fHitlist->hit[ih].x;
      double y = fHitlist->hit[ih].y;
/*      RMS += fHitlist->hit[ih].E*frac
             *(SQR(x-centroid.x) + SQR(y-centroid.y));*/
      MOM1x += fHitlist->hit[ih].E*frac*x;
      MOM1y += fHitlist->hit[ih].E*frac*y;
      MOM2 += fHitlist->hit[ih].E*frac*(SQR(x)+SQR(y));
//      double u0 = sqrt(SQR(centroid.x) + SQR(centroid.y));
//      double v0 = 0;
      double phi = atan2(centroid.y(),centroid.x());
      double u = x*cos(phi) + y*sin(phi);
      double v =-x*sin(phi) + y*cos(phi);
      MOM1u += fHitlist->hit[ih].E*frac*u;
      MOM1v += fHitlist->hit[ih].E*frac*v;
      MOM2u += fHitlist->hit[ih].E*frac*SQR(u);
      MOM2v += fHitlist->hit[ih].E*frac*SQR(v);
   }
//   fRMS = sqrt(RMS)/energy;
   fRMS = sqrt(energy*MOM2-SQR(MOM1x)-SQR(MOM1y))/(energy);
   fRMS_u = sqrt(energy*MOM2u - SQR(MOM1u))/(energy);
   fRMS_v = sqrt(energy*MOM2v - SQR(MOM1v))/(energy);

   bool something_changed = false;
   if (fabs(energy-fEnergy) > 0.001) {
      fEnergy = energy;
      something_changed = true;
   }
   if (fabs(Emax-fEmax) > 0.001) {
      fEmax = Emax;
      something_changed = true;
   }
   if (fabs(centroid.x()-fCentroid.x()) > 0.1 ||
       fabs(centroid.y()-fCentroid.y()) > 0.1) {
      fCentroid = centroid;
      something_changed = true;
   }
   if (something_changed) {
      for (int ih = 0; ih < fHitlist->nhits; ih++) {
         shower_profile(ih,fEallowed[ih],fEexpected[ih]);
      }
   }
   return something_changed;
}

#define MOLIER_RADIUS 3.696
#define MAX_SHOWER_RADIUS 25

void DFCALCluster::shower_profile(const int ihit,
                                double& Eallowed, double& Eexpected) const
{
   Eallowed = Eexpected = 0;
   if (fEnergy == 0) return;
   double x = fHitlist->hit[ihit].x;
   double y = fHitlist->hit[ihit].y;
   double dist = sqrt(SQR(x - fCentroid.x()) + SQR(y - fCentroid.y()));
   if (dist > MAX_SHOWER_RADIUS) return;
   double theta = atan2((double)sqrt(SQR(fCentroid.x()) + SQR(fCentroid.y())),(double)FCAL_Zmid);
   double phi = atan2(fCentroid.y(),fCentroid.x());
   double u0 = sqrt(SQR(fCentroid.x())+SQR(fCentroid.y()));
   double v0 = 0;
   double u = x*cos(phi) + y*sin(phi);
   double v =-x*sin(phi) + y*cos(phi);
   double vVar = SQR(MOLIER_RADIUS);
   double uVar = vVar+SQR(SQR(8*theta));
   double vTail = 4.5+0.9*log(fEnergy+0.05);
   double uTail = vTail+SQR(10*theta);
   double core = exp(-0.5*SQR(SQR(u-u0)/uVar + SQR(v-v0)/vVar));
   double tail = exp(-sqrt(SQR((u-u0)/uTail)+SQR((v-v0)/vTail)));
   Eexpected = fEnergy*core;
   Eallowed = 2*fEmax*core + (0.2+0.5*log(fEmax+1.))*tail;

   if ((dist <= 4.) && (Eallowed < fEmax) ) {
      std::cerr << "Warning: lgdClusterUC Eallowed value out of range!\n";
      Eallowed = fEmax;
   }
}




