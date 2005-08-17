// $Id$
//
//    File: DMCCheatHit.h
// Created: Mon Apr  4 08:18:07 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DMCCheatHit_
#define _DMCCheatHit_

#include "DObject.h"
#include "DFactory.h"

class DMCCheatHit:public DObject{
	public:
		HDCLASSDEF(DMCCheatHit);
		
		float r,phi,z;	///< coordinates of hit in cm and rad
		int track;		///< Track number
		int primary;	///< primary track=1    not primary track=0
		int system;		///< 1=CDC 2=FDC 3=BCAL 4=TOF 5=Cherenkov 6=FCAL 7=UPV

		inline bool operator<(const DMCCheatHit &mccheathit) const{
			if(track < mccheathit.track)return true;
			if(track > mccheathit.track)return false;
			if(z < mccheathit.z)return true;
			return false;
		}
};

#endif // _DMCCheatHit_

