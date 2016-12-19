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
#include <DMatrix.h>
#include <TMatrixFSym.h>
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

			TMatrixFSym ExyztCovariance;

			float const EErr() const { return sqrt(ExyztCovariance(0,0)); }
			float const xErr() const { return sqrt(ExyztCovariance(1,1)); }
			float const yErr() const { return sqrt(ExyztCovariance(2,2)); }
			float const zErr() const { return sqrt(ExyztCovariance(3,3)); }
			float const tErr() const { return sqrt(ExyztCovariance(4,4)); }
			float const XYcorr() const {
				if (xErr()>0 && yErr()>0) return ExyztCovariance(1,2)/xErr()/yErr();
				else return 0;
			}
			float const XZcorr() const {
				if (xErr()>0 && zErr()>0) return ExyztCovariance(1,3)/xErr()/zErr();
				else return 0;
			}
			float const YZcorr() const {
				if (yErr()>0 && zErr()>0) return ExyztCovariance(2,3)/yErr()/zErr();
				else return 0;
			}
			float const EXcorr() const {
				if (EErr()>0 && xErr()>0) return ExyztCovariance(0,1)/EErr()/xErr();
				else return 0;
			}
			float const EYcorr() const {
				if (EErr()>0 && yErr()>0) return ExyztCovariance(0,2)/EErr()/yErr();
				else return 0;
			}
			float const EZcorr() const {
				if (EErr()>0 && zErr()>0) return ExyztCovariance(0,3)/EErr()/zErr();
				else return 0;
			}
			float const XTcorr() const {
				if (xErr()>0 && tErr()>0) return ExyztCovariance(1,4)/xErr()/tErr();
				else return 0;
			}
			float const YTcorr() const {
				if (yErr()>0 && tErr()>0) return ExyztCovariance(2,4)/yErr()/tErr();
				else return 0;
			}
			float const ZTcorr() const {
				if (zErr()>0 && tErr()>0) return ExyztCovariance(3,4)/zErr()/tErr();
				else return 0;
			}
			float const ETcorr() const {
				if (EErr()>0 && tErr()>0) return ExyztCovariance(0,4)/EErr()/tErr();
				else return 0;
			}



		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "E(GeV)", "%6.2f", getEnergy());
			AddString(items, "X(cm)", "%7.2f", getPosition().X());
			AddString(items, "Y(cm)", "%7.2f", getPosition().Y());
			AddString(items, "Z(cm)", "%7.2f", getPosition().Z());
			AddString(items, "t(ns)", "%7.2f", getTime());
			AddString(items, "dE", "%5.3f", EErr());
			AddString(items, "dx", "%5.3f", xErr());
			AddString(items, "dy", "%5.3f", yErr());
			AddString(items, "dz", "%5.3f", zErr());
			AddString(items, "dt", "%5.3f", tErr());
			AddString(items, "EXcorr", "%5.3f", EXcorr());
			AddString(items, "EYcorr", "%5.3f", EYcorr());
			AddString(items, "EZcorr", "%5.3f", EZcorr());
			AddString(items, "ETcorr", "%5.3f", ETcorr());
			AddString(items, "XYcorr", "%5.3f", XYcorr());
			AddString(items, "XZcorr", "%5.3f", XZcorr());
			AddString(items, "XTcorr", "%5.3f", XTcorr());
			AddString(items, "YZcorr", "%5.3f", YZcorr());
			AddString(items, "YTcorr", "%5.3f", YTcorr());
			AddString(items, "ZTcorr", "%5.3f", ZTcorr());
		}

	private:

			double fEnergy; 
			double fTime; 
			DVector3 fPosition;  // Shower position in the FCAL
};


inline DVector3 DFCALShower::getPosition() const
{
      return fPosition;
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

