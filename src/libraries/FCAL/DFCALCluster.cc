// 
// FCALCluster member functions
//
#include <math.h>
#include "DFCALCluster.h"
#include "DFCALGeometry.h"

#ifndef SQR
# define SQR(x) (x)*(x)
#endif

DFCALCluster::DFCALCluster( const int nhits )
{
   fEnergy = 0;
   fEmax = 0;
   fTime = 0;
   fTimeMaxE = 0;
   fTimeEWeight = 0;
   fCentroid.SetXYZ( 0., 0., 0.);
   fRMS = 0;
   fRMS_t = 0;
   fRMS_x = 0;
   fRMS_y = 0;
   fRMS_u = 0;
   fRMS_v = 0;
   fNhits = 0;
   m_nFcalHits = nhits;

   if ( nhits > 0) {
      fHit = new int[ nhits ];
      fHitf = new double[ nhits ];
      fEallowed = new double[ nhits ];
      fEexpected = new double[ nhits ];
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
   if (fHit)
      delete [] fHit;
   if (fHitf)
      delete [] fHitf;
   if (fEallowed)
      delete [] fEallowed;
   if (fEexpected)
      delete [] fEexpected;
}


void DFCALCluster::saveHits( const userhits_t* const hits )
{
   
   for ( int i=0; i < fNhits; i++) {
      DFCALClusterHit_t h;
      JObject::oid_t id = getHitID( hits, i ) ;
      if ( id != 0 ) {  
         h.id = (JObject::oid_t) id;
	 h.ch = getHitCh( hits, i );
         h.E = getHitE( hits, i ) ;
         h.x = getHitX( hits, i ) ;
         h.y = getHitY( hits, i ) ;
         h.t = getHitT( hits, i ) ;
	 h.intOverPeak = getHitIntOverPeak( hits, i );
         my_hits.push_back(h);
      }
      else {
         static uint32_t Nwarns=0;
         if (++Nwarns < 100)
            std::cout << "Warning: DFCALCluster : corrupted cluster hit "
                      << i << std::endl;
         if (Nwarns == 100)
            std::cout << "Last warning!!! (further warnings supressed)"
                      << std::endl;
      }
   }
}


int DFCALCluster::addHit(const int ihit, const double frac)
{
   if (ihit >= 0 ) {
      fHit[fNhits] = ihit;
      fHitf[fNhits] = frac;
      ++fNhits;

      return 0;
   }
   else {
      return 1;
   }
}

void DFCALCluster::resetClusterHits()
{
   if (fNhits) {
      fNhits = 0;
   }
}

bool DFCALCluster::update( const userhits_t* const hitList,
			   double fcalFaceZ )
{

   double energy = 0;
   double t2EWeight = 0, tEWeight = 0;
   for ( int h = 0; h < fNhits; h++ ) {
      int ih = fHit[h];
      double frac = fHitf[h];
      double hitEnergy = hitList->hit[ih].E*frac;

      energy += hitEnergy;
      
      t2EWeight += hitEnergy * hitList->hit[ih].t * hitList->hit[ih].t;
      tEWeight += hitEnergy * hitList->hit[ih].t;
   }

   tEWeight /= energy;
   t2EWeight /= energy;

   double eMax=0, timeMax=0;
   if (fNhits > 0) { 
       eMax = hitList->hit[fHit[0]].E;
       timeMax = hitList->hit[fHit[0]].t;
   }

   DVector3 centroid;
   centroid.SetXYZ(0., 0., fcalFaceZ );
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

      weight = currentOffset + log(hitList->hit[ih].E*frac/energy);
      if (h == 0) {
         centerWeight = weight;
      }
      else {
         neighborMaxWeight =
               (neighborMaxWeight < weight)? weight:neighborMaxWeight;
      }
      if (weight > 0) {
         xc  += hitList->hit[ih].x*weight;
         yc  += hitList->hit[ih].y*weight;
         weightSum += weight;
      }
   }
   /*
    * Now patch up the center block's weight if it's got a neighbor
    * in the cluster that had positive weight
    */
   if (neighborMaxWeight > 0) {
      xc += (logFraction-1)*(centerWeight-neighborMaxWeight)
                             *hitList->hit[0].x;
      yc += (logFraction-1)*(centerWeight-neighborMaxWeight)
                             *hitList->hit[0].y;
      weightSum += (logFraction-1)*(centerWeight-neighborMaxWeight);
   }
   centroid.SetX(xc/weightSum);
   centroid.SetY(yc/weightSum);
#else
   for (int h = 0; h < fNhits; h++) {
      int ih = fHit[h];
      double frac = fHitf[h];
      xc += hitList->hit[ih].x*(hitList->hit[ih].E*frac);
      yc += hitList->hit[ih].y*(hitList->hit[ih].E*frac);

   }
   centroid.SetX(xc/energy);
   centroid.SetY(yc/energy);
#endif

   double MOM1x = 0;
   double MOM2x = 0;
   double MOM1y = 0;
   double MOM2y = 0;
   double MOM1u = 0;
   double MOM2u = 0;
   double MOM1v = 0;
   double MOM2v = 0;
   for ( int h = 0; h < fNhits; h++ ) {
      int ih = fHit[h];
      double frac = fHitf[h];
      double x = hitList->hit[ih].x;
      double y = hitList->hit[ih].y;

      MOM1x += hitList->hit[ih].E*frac*x;
      MOM1y += hitList->hit[ih].E*frac*y;
      MOM2x += hitList->hit[ih].E*frac*SQR(x);
      MOM2y += hitList->hit[ih].E*frac*SQR(y);

      double phi = atan2( centroid.y() , centroid.x() );
      double u = x*cos(phi) + y*sin(phi);
      double v =-x*sin(phi) + y*cos(phi);
      MOM1u += hitList->hit[ih].E*frac*u;
      MOM1v += hitList->hit[ih].E*frac*v;
      MOM2u += hitList->hit[ih].E*frac*SQR(u);
      MOM2v += hitList->hit[ih].E*frac*SQR(v);
   }

   bool something_changed = false;
   if (fabs(energy-fEnergy) > 0.001) {
      fEnergy = energy;
      something_changed = true;
   }
   if (fabs(eMax-fEmax) > 0.001) {
      fEmax = eMax;
      fTimeMaxE = timeMax;
      something_changed = true;
   }
   if (fabs(centroid.x()-fCentroid.x()) > 0.1 ||
       fabs(centroid.y()-fCentroid.y()) > 0.1) {
      fCentroid = centroid;
      something_changed = true;
   }

   if (something_changed) {

      fTime = timeMax;
      fTimeEWeight = tEWeight;
      fRMS_t = sqrt( t2EWeight - ( tEWeight * tEWeight ) );

      fRMS = sqrt(energy*(MOM2x+MOM2y)-SQR(MOM1x)-SQR(MOM1y))/(energy);
      fRMS_x = sqrt(energy*MOM2x - SQR(MOM1x))/(energy);
      fRMS_y = sqrt(energy*MOM2y - SQR(MOM1y))/(energy);
      fRMS_u = sqrt(energy*MOM2u - SQR(MOM1u))/(energy);
      fRMS_v = sqrt(energy*MOM2v - SQR(MOM1v))/(energy);

      for (int ih = 0; ih < hitList->nhits; ih++) {
	shower_profile( hitList, ih,fEallowed[ih],fEexpected[ih],
			fcalFaceZ+0.5*DFCALGeometry::blockLength());
      }
   }

   return something_changed;
}

void DFCALCluster::shower_profile( const userhits_t* const hitList, 
                                   const int ihit,
                                   double& Eallowed, double& Eexpected,
				   double fcalMidplaneZ) const
{

   //std::cout << " Run profile for hit " << ihit; 
   Eallowed = Eexpected = 0;
   if (fEnergy == 0)
      return;
   double x = hitList->hit[ihit].x;
   double y = hitList->hit[ihit].y;
   double dist = sqrt(SQR(x - fCentroid.x()) + SQR(y - fCentroid.y()));
   if (dist > MAX_SHOWER_RADIUS)
      return;
   double theta = atan2((double)sqrt(SQR(fCentroid.x()) + SQR(fCentroid.y())), fcalMidplaneZ);
   double phi = atan2( fCentroid.y(), fCentroid.x() );
   double u0 = sqrt(SQR(fCentroid.x())+SQR(fCentroid.y()));
   double v0 = 0;
   double u = x*cos(phi) + y*sin(phi);
   double v =-x*sin(phi) + y*cos(phi);
   double vVar = SQR(MOLIERE_RADIUS);
   double uVar = vVar+SQR(SQR(8*theta));
   double vTail = 4.5+0.9*log(fEnergy+0.05);
   double uTail = vTail+SQR(10*theta);
   double core = exp(-0.5*SQR(SQR(u-u0)/uVar + SQR(v-v0)/vVar));
   double tail = exp(-sqrt(SQR((u-u0)/uTail)+SQR((v-v0)/vTail)));
   Eexpected = fEnergy*core;
   Eallowed = 2*fEmax*core + (0.2+0.5*log(fEmax+1.))*tail;

   if ((dist <= 4.) && (Eallowed < fEmax) ) {
      std::cerr << "Warning: FCAL cluster Eallowed value out of range!\n";
      Eallowed = fEmax;
   }
}
