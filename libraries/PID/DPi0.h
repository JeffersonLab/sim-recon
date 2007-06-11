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

//		DLorentzVector getMom4() const; 
                unsigned int getTag1() const; // the origin of the 1st photon (FCAL, BCAL)
                unsigned int getTag2() const; 
                unsigned int getID1() const; // JANA id's of children
                unsigned int getID2() const; 

// form pi0 candidate from two photons 
//                void setMom4(const DLorentzVector gamma1, const DLorentzVector gamma2); 
// set Pi0 bits regarding detected photons  
                void setChildrenTag(unsigned int tag1, unsigned int tag2 ); 
                void setChildrenID(oid_t id1, oid_t id2 ); 

	private:

               unsigned int fTag1; 
               unsigned int fTag2; 
               oid_t fID1;  
               oid_t fID2; 
//               DLorentzVector fMom4;  // pi0 4-momentum

};


// return origin of pi0  (Fcal=0, Bcal=1)
inline unsigned int DPi0::getTag1() const
{
      return fTag1;
}

inline unsigned int DPi0::getTag2() const
{
      return fTag2;
}

// return child ID's
inline unsigned int DPi0::getID1() const
{
      return fID1;
}
inline unsigned int DPi0::getID2() const
{
      return fID2;
}

/*inline DLorentzVector DPi0::getMom4() const
{
      return fMom4;
}
*/

#endif // _DPi0_

