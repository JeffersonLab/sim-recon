// $Id: DFCALPhoton.h 1899 2006-07-13 16:29:56Z davidl $
//
//    File: DFCALPhoton.h
// Created: Tue Jan 22 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DFCALPhoton_
#define _DFCALPhoton_

#include <math.h>
#include <DVector3.h>
#include <DLorentzVector.h>
#include "DFCALCluster.h"
using namespace std;

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DFCALPhoton:public JObject{
	public:
		JOBJECT_PUBLIC(DFCALPhoton);
                
			DFCALPhoton();
			~DFCALPhoton();

		// getter functions
			DVector3 getPosition() const; 
                        DVector3 getPositionError() const;
			double getEnergy() const;  
			double getTime() const;  
			DVector3 getMom3() const; 
			DLorentzVector getMom4() const;

		// set photon energy and position 
			void setPosition( const DVector3 aPosition );  
			void setEnergy(const double energy);  
			void setTime(const double time);  

                        void setPosError(const double aXerr, const double aYerr, const double aZerr);

                // set photon momentum
			void setMom3(const double energy, const DVector3 pos);    
			void setMom4();  

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "E(GeV)", "%6.2f", getEnergy());
			AddString(items, "Px(GeV)", "%6.2f", getMom3().X());
			AddString(items, "Py(GeV)", "%6.2f", getMom3().Y());
			AddString(items, "Pz(GeV)", "%6.2f", getMom3().Z());
			AddString(items, "X(cm)", "%7.2f", getPosition().X());
			AddString(items, "Y(cm)", "%7.2f", getPosition().Y());
			AddString(items, "Z(cm)", "%7.2f", getPosition().Z());
			AddString(items, "t(ns)", "%7.2f", getTime());
		}

	private:

			double fEnergy; 
			double fTime; 
			DVector3 fPosition;  // Photon position in the FCAL
                        DVector3 fPositionError;  // Errors in X and Y are estimated from
			DVector3 fMom3;  // Photon 3-momentum
			DLorentzVector fMom4;  // Photon 4-momentum
};


inline DVector3 DFCALPhoton::getPosition() const
{
      return fPosition;
}

inline DVector3 DFCALPhoton::getPositionError() const
{
      return fPositionError;
}

inline DVector3 DFCALPhoton::getMom3() const
{
      return fMom3;
}

inline double DFCALPhoton::getEnergy() const
{
      return fEnergy;
}

inline double DFCALPhoton::getTime() const
{
      return fTime;
}


inline DLorentzVector DFCALPhoton::getMom4() const
{
      return fMom4;
}


#endif // _DFCALPhoton_

