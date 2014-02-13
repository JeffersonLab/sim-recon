// $Id: DFCALShower.h 1899 2006-07-13 16:29:56Z davidl $
//
//    File: DFCALShower.h
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DFCALCluster_
#define _DFCALCluster_

#include <DVector3.h>
//#include "DFCALHit.h"
using namespace std;

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

#define FCAL_USER_HITS_MAX 2800


/* Define FCAL mid plane for the cluster angle determination 
   (it should be the front face of the FCAL for real clusters).
   Make sure the same vertex from the Geant control.in file is used here
   and elsewhere in analysis - MK.
*/
//static float FCAL_Zmid = FCAL_Zmin - Shower_Vertex_Z +FCAL_Zlen/2.0;

class DFCALCluster:public JObject{
	public:
		JOBJECT_PUBLIC(DFCALCluster);
                
                DFCALCluster( const int nhits );
		~DFCALCluster();

// hits structures taken from radphi

typedef struct {
           oid_t id;
	   float x;
	   float y;
	   float E;
	   float t;
	} userhit_t;

typedef struct {
	   int nhits;
	   userhit_t hit[1];
	} userhits_t;

// to replace c-like array pointers with C++ vectors 
	typedef struct {
           oid_t id;
	   float x;
	   float y;
	   float E;
	   float t;
	} DFCALClusterHit_t;

		void saveHits( const userhits_t* const hit );

		double getEexpected(const int ihit) const;
		double getEallowed(const int ihit) const;
		double getEnergy() const;
		double getEmax() const;
		double getTime() const;
		DVector3 getCentroid() const;
		double getRMS() const;
		double getRMS_x() const;
		double getRMS_y() const;
		double getRMS_u() const;
		double getRMS_v() const;
		int getHits() const; // get number of hits owned by a cluster
		int addHit(const int ihit, const double frac);
		void resetClusterHits();
		bool update( const userhits_t* const hitList );
// get hits that form a cluster after clustering is finished
                inline const vector<DFCALClusterHit_t> GetHits() const { return my_hits; }

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "x(cm)", "%3.1f", getCentroid().x());
			AddString(items, "y(cm)", "%3.1f", getCentroid().y());
			AddString(items, "E(GeV)", "%2.3f", getEnergy());
			AddString(items, "t(ns)", "%2.3f", getTime());
		}

	private:

                //const double VERTEX_Z = 65.;
		void shower_profile( const userhits_t* const hitList, 
                                     const int ihit,
		                     double& Eallowed, double& Eexpected) const ;
// internal parsers of properties for a hit belonging to a cluster 
                int getHitID( const userhits_t* const hitList, const int ihit) const;
                double getHitX( const userhits_t* const hitList, const int ihit) const;
                double getHitY( const userhits_t* const hitList, const int ihit) const;
                double getHitT( const userhits_t* const hitList, const int ihit) const;
                double getHitE( const userhits_t* const hitList, const int ihit) const;  // hit energy owned by cluster
                double getHitEhit( const userhits_t* const hitList, const int ihit) const; // energy in a FCAL block

		double fEnergy;              // total cluster energy (GeV) or 0 if stale
		double fTime;                // cluster time(ns) set by first (max E) block, for now
		double fEmax;                // energy in the first block of the cluster
		DVector3 fCentroid;         // cluster centroid position (cm)
		double fRMS;                 // cluster r.m.s. size (cm)
		double fRMS_x;               // cluster r.m.s. size along X-axis (cm)
		double fRMS_y;               // cluster r.m.s. size along Y-axis (cm)
		double fRMS_u;               // cluster r.m.s. size in radial direction (cm)
		double fRMS_v;               // cluster r.m.s. size in azimuth direction (cm)
                int m_nFcalHits;             // total number of hits to work with
					     //	need to rename other member data  
		int fNhits;                  // number of hits owned by this cluster
		int *fHit;                   // index list of hits owned by this cluster
		double *fHitf;               // list of hit fractions owned by this cluster
		double *fEexpected;          // expected energy of hit by cluster (GeV)
		double *fEallowed;           // allowed energy of hit by cluster (GeV)
                vector<DFCALClusterHit_t> my_hits; // container for hits that form a cluster
						   // to be used after clustering is done

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

inline DVector3 DFCALCluster::getCentroid() const
{
   return fCentroid;
}

inline double DFCALCluster::getRMS() const
{
   return fRMS;
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

inline int DFCALCluster::getHitID(const userhits_t* const hitList, const int ihit ) const
{
   if ( ihit >= 0  && ihit < fNhits && hitList && ihit < hitList->nhits ) {
     return (int) hitList->hit[ fHit[ ihit ] ].id;
   }
   else {
     return -1;
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

