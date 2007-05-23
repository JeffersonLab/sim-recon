//
//    File: DPIDPhoton.h
// Created: Tue Apr 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#ifndef _DPIDPi0_
#define _DPIDPi0_

#include <TLorentzVector.h>
#include "DPIDPhoton.h"

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DPIDPi0:public JObject{
	public:
		HDCLASSDEF(DPIDPi0);
                
                DPIDPi0();
		~DPIDPi0();

		const TLorentzVector getMom4() const; 
                const unsigned int getTag1() const; // the origin of the 1st photon (FCAL, BCAL)
                const unsigned int getTag2() const; 

		// form pi0 candidate from two photons 
                void setMom4(const TLorentzVector gamma1, const TLorentzVector gamma2); 
                void setOrig(const unsigned int tag1, const unsigned int tag2 ); // set Pi0 bits regarding detected photons  

	private:

               unsigned int fTag1; 
               unsigned int fTag2; 
               TLorentzVector fMom4;  // pi0 4-momentum

};


// return origin of pi0 kids (Fcal=0, Bcal=1)
inline const unsigned int DPIDPi0::getTag1() const
{
      return fTag1;
}

inline const unsigned int DPIDPi0::getTag2() const
{
      return fTag2;
}


inline const TLorentzVector DPIDPi0::getMom4() const
{
      return fMom4;
}


#endif // _DPIDPi0_

