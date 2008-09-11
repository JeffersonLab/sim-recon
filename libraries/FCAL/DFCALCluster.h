// $Id: DFCALShower.h 1899 2006-07-13 16:29:56Z davidl $
//
//    File: DFCALShower.h
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DFCALCluster_
#define _DFCALCluster_

//#include <ntypes.h>
//#include "userhits.h"
#include <DVector3.h>
//#include "DFCALHit.h"

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

#define FCAL_USER_HITS_MAX 2800

// hits structures taken from radphi
	typedef struct {
	   float x;
	   float y;
	   float E;
	   float t;
	} userhit_t;

	typedef struct {
	   int nhits;
	   userhit_t hit[1];
	} userhits_t;


/* Define FCAL mid plane for the cluster angle determination 
   (it should be the front face of the FCAL for real clusters).
   Make sure the same vertex from the Geant control.in file is used here
   and elsewhere in analysis - MK.
*/
static float FCAL_Zmin = 625.3;
static float FCAL_Zlen = 45.0;
static float Shower_Vertex_Z = 65.0; 
static float FCAL_Zmid = FCAL_Zmin - Shower_Vertex_Z +FCAL_Zlen/2.0;

class DFCALCluster:public JObject{
	public:
		JOBJECT_PUBLIC(DFCALCluster);
                
                DFCALCluster();
		~DFCALCluster();

		static void setHitlist(const userhits_t* const hits);

		double getEexpected(const int ihit) const;
		double getEallowed(const int ihit) const;
		double getEnergy() const;
		double getEmax() const;
		DVector3 getCentroid() const;
		double getRMS() const;
		double getRMS_x() const;
		double getRMS_y() const;
		double getRMS_u() const;
		double getRMS_v() const;
		int getHits() const;
		int getHits(int hitlist[], int nhits) const;
		static int getUsed(const int ihit);
		int addHit(const int ihit, const double frac);
		void resetHits();
		bool update();

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "x(cm)", "%3.1f", getCentroid().x());
			AddString(items, "y(cm)", "%3.1f", getCentroid().y());
			AddString(items, "E(GeV)", "%2.3f", getEnergy());
		}

	private:

		void shower_profile(const int ihit,
		                    double& Eallowed, double& Eexpected) const ;
		double fEnergy;              // total cluster energy (GeV) or 0 if stale
		double fEmax;                // energy in the first block of the cluster
		DVector3 fCentroid;         // cluster centroid position (cm)
		double fRMS;                 // cluster r.m.s. size (cm)
		double fRMS_x;               // cluster r.m.s. size along X-axis (cm)
		double fRMS_y;               // cluster r.m.s. size along Y-axis (cm)
		double fRMS_u;               // cluster r.m.s. size in radial direction (cm)
		double fRMS_v;               // cluster r.m.s. size in azimuth direction (cm)
		int fNhits;                  // number of hits owned by this cluster
		int *fHit;                   // index list of hits owned by this cluster
		double *fHitf;               // list of hit fractions owned by this cluster
		static const userhits_t* fHitlist;   // pointer to user's hit list
		static int *fHitused;        // number of clusters that use hit,
					     // or -1 if it is a cluster seed
		double *fEexpected;          // expected energy of hit by cluster (GeV)
		double *fEallowed;           // allowed energy of hit by cluster (GeV)

};

inline double DFCALCluster::getEexpected(const int ihit) const
{

   if ( (ihit >= 0) && (fHitlist) && (ihit < fHitlist->nhits) ) 
      return fEexpected[ihit];
   else
      return 0;
}

inline double DFCALCluster::getEallowed(const int ihit) const
{

   if (ihit >= 0 && fHitlist && ihit < fHitlist->nhits) 
      return fEallowed[ihit];
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

inline int DFCALCluster::getHits(int hitlist[], const int nhits) const
{
   int ih;
   for (ih = 0; (ih < fNhits) && (ih < nhits); ih++) {
      hitlist[ih] = fHit[ih];
   }
   return ih;
}

inline int DFCALCluster::getUsed(const int ihit)
{

   if (ihit >= 0 && fHitlist && ihit < fHitlist->nhits) 
      return fHitused[ihit];
   else
      return 0;
}

#endif // _DFCALCluster_

