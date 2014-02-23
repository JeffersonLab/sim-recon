// $Id: DFCALShower.h 1899 2006-07-13 16:29:56Z davidl $
//
//    File: DFCALShower.h
// Created: Tue Jan 22 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DFCALShower_
#define _DFCALShower_

#include <math.h>
#include <DVector3.h>
#include <DLorentzVector.h>
#include "DFCALCluster.h"
using namespace std;

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DFCALShower:public JObject{
	public:
		JOBJECT_PUBLIC(DFCALShower);
                
			DFCALShower();
			~DFCALShower();

		// getter functions
//  shower position in calorimeter, after depth correction
//  default vertex is in the target center
			DVector3 getPosition() const; 
                        DVector3 getPositionError() const;
			double getEnergy() const;  
			double getTime() const;  

		// set shower energy and position 
			void setPosition( const DVector3 aPosition );  
			void setEnergy(const double energy);  
			void setTime(const double time);  

                        void setPosError(const double aXerr, const double aYerr, const double aZerr);

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "E(GeV)", "%6.2f", getEnergy());
			AddString(items, "X(cm)", "%7.2f", getPosition().X());
			AddString(items, "Y(cm)", "%7.2f", getPosition().Y());
			AddString(items, "Z(cm)", "%7.2f", getPosition().Z());
			AddString(items, "t(ns)", "%7.2f", getTime());
		}

	private:

			double fEnergy; 
			double fTime; 
			DVector3 fPosition;  // Shower position in the FCAL
                        DVector3 fPositionError;  // Errors in X and Y are estimated from
};


inline DVector3 DFCALShower::getPosition() const
{
      return fPosition;
}

inline DVector3 DFCALShower::getPositionError() const
{
      return fPositionError;
}

inline double DFCALShower::getEnergy() const
{
      return fEnergy;
}

inline double DFCALShower::getTime() const
{
      return fTime;
}

#endif // _DFCALShower_

