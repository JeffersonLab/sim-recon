//
//    File: DPhoton.h
// Created: Tue Apr 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#ifndef _DPhoton_
#define _DPhoton_

//#include <TLorentzVector.h>
#include <DKinematicData.h>

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DPhoton:public JObject{
	public:
		HDCLASSDEF(DPhoton);
                
                DPhoton();
		~DPhoton();

		TLorentzVector getMom4() const; 
                unsigned int getTag() const;
		const double getDtRT() const; 
                void setMom4(const TLorentzVector gamma);  
                void setTag(const unsigned int tag);  
                void setDtRT(const double dtrt);  
      
	private:

                TLorentzVector fMom4;  // Photon 4-momentum
                unsigned int fTag; //Photon origin (FCAL/BCAL 0/1))
                double fDtRT; //Distance to closest track's RefenceTrajectory

};


inline TLorentzVector DPhoton::getMom4() const
{
      return fMom4;
}

inline unsigned int DPhoton::getTag() const
{
      return fTag;
}

inline const double DPhoton::getDtRT() const
{
      return fDtRT;
}
#endif // _DPhoton_

