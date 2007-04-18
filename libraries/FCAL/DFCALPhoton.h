// $Id: DFCALPhoton.h 1899 2006-07-13 16:29:56Z davidl $
//
//    File: DFCALPhoton.h
// Created: Tue Jan 22 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DFCALPhoton_
#define _DFCALPhoton_

#include <DVector3.h>
#include <TLorentzVector.h>
#include "DFCALCluster.h"

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DFCALPhoton:public JObject{
	public:
		HDCLASSDEF(DFCALPhoton);
                
                DFCALPhoton();
		~DFCALPhoton();

		const double getEnergy() const;  // get photon energy and momentum
		const DVector3 getMom3() const; 
                const TLorentzVector getMom4() const;

		// fix photon energy and momentu from cluster energy and position
                void setEnergy(const double energy);  
                void setMom3(const double energy, const DVector3 mom);    
                void setMom4();  

	private:

                double fEnergy; 
                DVector3 fMom3;  // Photon 3-momentum
                TLorentzVector fMom4;  // Photon 4-momentum

};


inline const DVector3 DFCALPhoton::getMom3() const
{
      return fMom3;
}

inline const double DFCALPhoton::getEnergy() const
{
      return fEnergy;
}


inline const TLorentzVector DFCALPhoton::getMom4() const
{
      return fMom4;
}


#endif // _DFCALPhoton_

