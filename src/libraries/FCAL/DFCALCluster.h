// $Id: DFCALShower.h 1899 2006-07-13 16:29:56Z davidl $
//
//    File: DFCALShower.h
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DFCALCluster_
#define _DFCALCluster_

#include <DVector3.h>
#include "DFCALHit.h"
using namespace std;

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

#define FCAL_USER_HITS_MAX 2800
#define MOLIERE_RADIUS 3.696
#define MAX_SHOWER_RADIUS 25

class DFCALCluster : public JObject {
   public:
      JOBJECT_PUBLIC(DFCALCluster);
                
      DFCALCluster( const int nhits );
      ~DFCALCluster();

   typedef struct {
      oid_t id;
      int ch;
      float x;
      float y;
      float E;
      float t;
     // this may be pruned at a later stage to reduce size
     // but it is useful for development now
      float intOverPeak;
   } userhit_t;

   typedef struct {
      int nhits;
      userhit_t hit[1];
   } userhits_t;

   typedef struct {
      oid_t id;
      int ch;
      float x;
      float y;
      float E;
      float t;
      float intOverPeak;
   } DFCALClusterHit_t;

   void saveHits( const userhits_t* const hit );

   double getEexpected(const int ihit) const;
   double getEallowed(const int ihit) const;
   double getEnergy() const;
   double getEmax() const;
   double getTime() const;
   double getTimeMaxE() const;
   double getTimeEWeight() const;
   // the centroid returned by this function is in 
   // a frame where z=0 is at the center of the target
   DVector3 getCentroid() const;

   double getRMS() const;
   double getRMS_t() const;
   double getRMS_x() const;
   double getRMS_y() const;
   double getRMS_u() const;
   double getRMS_v() const;

   int getHits() const; // get number of hits owned by a cluster
   int addHit(const int ihit, const double frac);
   void resetClusterHits();
   bool update( const userhits_t* const hitList, double fcalFaceZ );

// get hits that form a cluster after clustering is finished
   inline const vector<DFCALClusterHit_t> GetHits() const { return my_hits; }
   inline uint32_t GetNHits(void) const { return my_hits.size(); }

   void toStrings(vector<pair<string,string> > &items) const {
      AddString(items, "x(cm)", "%3.1f", getCentroid().x());
      AddString(items, "y(cm)", "%3.1f", getCentroid().y());
      AddString(items, "z(cm)", "%3.1f", getCentroid().z());
      AddString(items, "E(GeV)", "%2.3f", getEnergy());
      AddString(items, "t(ns)", "%2.3f", getTime());
   }

   private:

   void shower_profile( const userhits_t* const hitList, 
                        const int ihit,
                        double& Eallowed, double& Eexpected, 
			double fcalMidplaneZ ) const ;

   // internal parsers of properties for a hit belonging to a cluster 
   oid_t  getHitID( const userhits_t* const hitList, const int ihit) const;
   int    getHitCh( const userhits_t* const hitList, const int ihit) const;
   double getHitX( const userhits_t* const hitList, const int ihit) const;
   double getHitY( const userhits_t* const hitList, const int ihit) const;
   double getHitT( const userhits_t* const hitList, const int ihit) const;
   double getHitIntOverPeak( const userhits_t* const hitList, const int ihit) const;
   double getHitE( const userhits_t* const hitList, const int ihit) const;  // hit energy owned by cluster
   double getHitEhit( const userhits_t* const hitList, const int ihit) const; // energy in a FCAL block

   double fEnergy;        // total cluster energy (GeV) or 0 if stale
   double fTime;          // cluster time(ns) set equivalent to fTimeMaxE below
   double fTimeMaxE;      // cluster time(ns) set by first (max E) block, for now
   double fTimeEWeight;   // cluster time(ns) energy weighted average
   double fEmax;          // energy in the first block of the cluster
   DVector3 fCentroid;    // cluster centroid position (cm) [z=0 at target center]
   double fRMS;           // cluster r.m.s. size (cm)
   double fRMS_t;         // r.m.s. of energy weighted hits (ns)
   double fRMS_x;         // cluster r.m.s. size along X-axis (cm)
   double fRMS_y;         // cluster r.m.s. size along Y-axis (cm)
   double fRMS_u;         // cluster r.m.s. size in radial direction (cm)
   double fRMS_v;         // cluster r.m.s. size in azimuth direction (cm)
   int m_nFcalHits;       // total number of hits to work with
                          //   need to rename other member data  
   int fNhits;            // number of hits owned by this cluster
   int *fHit;             // index list of hits owned by this cluster
   double *fHitf;         // list of hit fractions owned by this cluster
   double *fEexpected;    // expected energy of hit by cluster (GeV)
   double *fEallowed;     // allowed energy of hit by cluster (GeV)

   // container for hits that form a cluster to be used after clustering is done
   vector<DFCALClusterHit_t> my_hits; 

};

inline double DFCALCluster::getEexpected(const int ihit) const
{
   if ( ihit >= 0 && ihit < m_nFcalHits )
      return fEexpected[ ihit ];
   else
      return 0;
}

inline double DFCALCluster::getEallowed(const int ihit) const
{

   if ( ihit >= 0 && ihit < m_nFcalHits ) 
      return fEallowed[ ihit ];
   else
      return 0;
}

inline double DFCALCluster::getEnergy() const
{
   return fEnergy;
}

inline double DFCALCluster::getEmax() const
{
   return fEmax;
}

inline double DFCALCluster::getTime() const
{
   return fTime;
}

inline double DFCALCluster::getTimeMaxE() const
{
   return fTimeMaxE;
}

inline double DFCALCluster::getTimeEWeight() const
{
   return fTimeEWeight;
}

inline DVector3 DFCALCluster::getCentroid() const
{
   return fCentroid;
}

inline double DFCALCluster::getRMS() const
{
   return fRMS;
}

inline double DFCALCluster::getRMS_t() const
{
   return fRMS_t;
}

inline double DFCALCluster::getRMS_x() const
{
   return fRMS_x;
}

inline double DFCALCluster::getRMS_y() const
{
   return fRMS_y;
}

inline double DFCALCluster::getRMS_u() const
{
   return fRMS_u;
}

inline double DFCALCluster::getRMS_v() const
{
   return fRMS_v;
}

inline int DFCALCluster::getHits() const
{
   return fNhits;
}

inline JObject::oid_t DFCALCluster::getHitID(const userhits_t* const hitList, const int ihit ) const
{
   if ( ihit >= 0  && ihit < fNhits && hitList && ihit < hitList->nhits ) {
     return hitList->hit[ fHit[ ihit ] ].id;
   }
   else {
     return 0;
   }
}

inline int DFCALCluster::getHitCh(const userhits_t* const hitList, const int ihit ) const
{
   if ( ihit >= 0  && ihit < fNhits && hitList && ihit < hitList->nhits ) {
     return hitList->hit[ fHit[ ihit ] ].ch;
   }
   else {
     return 0;
   }
}

inline double DFCALCluster::getHitX(const userhits_t* const hitList, const int ihit ) const
{
   if ( ihit >= 0  && ihit < fNhits && hitList && ihit < hitList->nhits ) {
     return  hitList->hit[ fHit[ ihit ] ].x;
   }
   else {
     return 0.;
   }
}

inline double DFCALCluster::getHitY(const userhits_t* const hitList, const int ihit ) const
{
   if ( ihit >= 0  && ihit < fNhits && hitList && ihit < hitList->nhits ) {
     return  hitList->hit[ fHit[ ihit ] ].y;
   }
   else {
     return 0.;
   }
}

inline double DFCALCluster::getHitT(const userhits_t* const hitList, const int ihit ) const
{
   if ( ihit >= 0  && ihit < fNhits && hitList && ihit < hitList->nhits ) { 
     return  hitList->hit[ fHit[ ihit ] ].t;
   }
   else {
     return 0.;
   }
}

inline double DFCALCluster::getHitIntOverPeak(const userhits_t* const hitList, const int ihit ) const
{
   if ( ihit >= 0  && ihit < fNhits && hitList && ihit < hitList->nhits ) { 
     return  hitList->hit[ fHit[ ihit ] ].intOverPeak;
   }
   else {
     return 0.;
   }
}

inline double DFCALCluster::getHitE(const userhits_t* const hitList, const int ihit ) const
{
   if ( ihit >= 0  && ihit < fNhits && hitList && ihit < hitList->nhits ) {
     return fHitf[ ihit ] * hitList->hit[ fHit[ ihit ] ].E ;
   }
   else {
     return -1.;
   }
}

inline double DFCALCluster::getHitEhit(const userhits_t* const hitList, const int ihit ) const
{
   if ( ihit >= 0  && ihit < fNhits && hitList && ihit < hitList->nhits ) {
     return hitList->hit[ fHit[ ihit ] ].E ;
   }
   else {
     return -1.;
   }
}

#endif // _DFCALCluster_

