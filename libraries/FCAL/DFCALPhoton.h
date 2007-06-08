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

		// getter functions
			DVector3 getPosition() const; 
			double getEnergy() const;  
			DVector3 getMom3() const; 
			TLorentzVector getMom4() const;

		// fix photon energy and momentu from cluster energy and position
			void fixDepth(const double energy, const DVector3 centroid);  
			void fixEnergy(const double energy);  

                // set photon momentum
			void setMom3(const double energy, const DVector3 pos);    
			void setMom4();  

	private:

			double fEnergy; 
			DVector3 fPosition;  // Photon position in the FCAL
			DVector3 fMom3;  // Photon 3-momentum
			TLorentzVector fMom4;  // Photon 4-momentum
};


inline DVector3 DFCALPhoton::getPosition() const
{
      return fPosition;
}

inline DVector3 DFCALPhoton::getMom3() const
{
      return fMom3;
}

inline double DFCALPhoton::getEnergy() const
{
      return fEnergy;
}


inline TLorentzVector DFCALPhoton::getMom4() const
{
      return fMom4;
}


#endif // _DFCALPhoton_

