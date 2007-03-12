// $Id: DFCALPhoton.h 1899 2006-07-13 16:29:56Z davidl $
//
//    File: DFCALPhoton.h
// Created: Tue Jan 22 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DFCALPhoton_
#define _DFCALPhoton_

#include <TVector3.h>
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
		const TVector3 get3Mom() const; 

		// fix photon energy and momentu from cluster energy and position
                void setEnergy(const double energy);  
                void set3Mom(const double energy, const TVector3 pos);    
// just a placeholder right now
                static void makePhoton(); 

	private:

                double fEnergy; 
                TVector3 fMom;  // Photon 3-momentum
//                TLorentzVector fPos;  // Photon 4-position

};


inline const TVector3 DFCALPhoton::get3Mom() const
{
      return fMom;
}

inline const double DFCALPhoton::getEnergy() const
{
      return fEnergy;
}


/*
inline TLorentzVector DFCALPhoton::get4Position()
{
      return fPos;
}*/


#endif // _DFCALPhoton_

