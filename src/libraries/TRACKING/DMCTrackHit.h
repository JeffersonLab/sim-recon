// $Id$
//
//    File: DMCTrackHit.h
// Created: Mon Apr  4 08:18:07 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DMCTrackHit_
#define _DMCTrackHit_

#include <JANA/JObject.h>
using namespace jana;

#include "GlueX.h"

class DMCTrackHit:public JObject{
 public:
  JOBJECT_PUBLIC(DMCTrackHit);
  
  float r,phi,z;	///< coordinates of hit in cm and rad
  int track;		///< Track number
  int itrack;		///< MC track index
  int primary;	        ///< primary track=1    not primary track=0
  int ptype;            /// particle type  
  DetectorSystem_t system;///< 1=CDC 2=FDC 4=BCAL 8=TOF 16=Cherenkov 32=FCAL 64=UPV

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "r(cm)", "%3.1f", r);
			AddString(items, "phi(rad)", "%1.3f", phi);
			AddString(items, "z(cm)", "%3.1f", z);
			AddString(items, "track", "%d", track);
			AddString(items, "itrack", "%d", itrack);
			AddString(items, "primary", "%d", primary);
			AddString(items, "ptype", "%d", ptype);
			AddString(items, "system", "%s", SystemName(system));
		}
};

#endif // _DMCTrackHit_

