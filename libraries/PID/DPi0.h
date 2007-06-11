//
//    File: DPhoton.h
// Created: Tue Apr 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#ifndef _DPi0_
#define _DPi0_

#include <TLorentzVector.h>
#include "DPhoton.h"

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DPi0:public DKinematicData {
	public:
		HDCLASSDEF(DPi0);
                
                DPi0();
		~DPi0();

                oid_t getChildrenID(int child)  const; // JANA id's of children
                unsigned int getChildrenTag(int child) const; // the origin of the 1st photon (FCAL, BCAL)

                void setChildrenID(oid_t id1, oid_t id2 ); 
                void setChildrenTag(unsigned int tag1, unsigned int tag2 ); 

	private:

               unsigned int fTags[2]; 
               oid_t fIDs[2];  

};


// return origin of pi0  (Fcal=0, Bcal=1)
inline unsigned int DPi0::getChildrenTag(int child) const
{
      return fTags[child];
}

// return child ID's
inline oid_t DPi0::getChildrenID(int child) const
{
      return fIDs[child];
}

#endif // _DPi0_

