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

  //  shower position in calorimeter in global coordinates, after depth correction
  //  for the depth correction, the default vertex is in the target center
  DVector3 getPosition() const; 
  DVector3 getPositionError() const;
  double getEnergy() const;  
  double getTime() const;
  double getDocaTrack() const;
  double getTimeTrack() const;
  double getSumU() const;
  double getSumV() const;
  double getE9E25() const;
  double getE1E9() const;
  int getNumBlocks() const;

  // set shower information
  void setPosition( const DVector3& aPosition );  
  void setEnergy( const double energy );  
  void setTime( const double time );
  void setDocaTrack( const double docaTrack );
  void setTimeTrack( const double tTrack );
  void setSumU( const double sumU );
  void setSumV( const double sumV );
  void setE9E25( const double e9e25 );
  void setE1E9( const double e1e9 );
  void setNumBlocks( const int numBlocks );

  TMatrixFSym ExyztCovariance;

  float EErr() const { return sqrt(ExyztCovariance(0,0)); }
  float xErr() const { return sqrt(ExyztCovariance(1,1)); }
  float yErr() const { return sqrt(ExyztCovariance(2,2)); }
  float zErr() const { return sqrt(ExyztCovariance(3,3)); }
  float tErr() const { return sqrt(ExyztCovariance(4,4)); }
  float XYcorr() const {
    if (xErr()>0 && yErr()>0) return ExyztCovariance(1,2)/xErr()/yErr();
    else return 0;
  }
  float XZcorr() const {
    if (xErr()>0 && zErr()>0) return ExyztCovariance(1,3)/xErr()/zErr();
    else return 0;
  }
  float YZcorr() const {
    if (yErr()>0 && zErr()>0) return ExyztCovariance(2,3)/yErr()/zErr();
    else return 0;
  }
  float EXcorr() const {
    if (EErr()>0 && xErr()>0) return ExyztCovariance(0,1)/EErr()/xErr();
    else return 0;
  }
  float EYcorr() const {
    if (EErr()>0 && yErr()>0) return ExyztCovariance(0,2)/EErr()/yErr();
    else return 0;
  }
  float EZcorr() const {
    if (EErr()>0 && zErr()>0) return ExyztCovariance(0,3)/EErr()/zErr();
    else return 0;
  }
  float XTcorr() const {
    if (xErr()>0 && tErr()>0) return ExyztCovariance(1,4)/xErr()/tErr();
    else return 0;
  }
  float YTcorr() const {
    if (yErr()>0 && tErr()>0) return ExyztCovariance(2,4)/yErr()/tErr();
    else return 0;
  }
  float ZTcorr() const {
    if (zErr()>0 && tErr()>0) return ExyztCovariance(3,4)/zErr()/tErr();
    else return 0;
  }
  float ETcorr() const {
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
    AddString(items, "docaTr","%7.2f", getDocaTrack());
    AddString(items, "timeTr","%7.2f", getTimeTrack());
    AddString(items, "sumU", "%7.2f", getSumU());
    AddString(items, "sumV", "%7.2f", getSumV());
    AddString(items, "E9E25", "%7.2f", getE9E25());
    AddString(items, "E1E9", "%7.2f", getE1E9());
    AddString(items, "numBlocks", "%3.0i",getNumBlocks());
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
  DVector3 fPosition;        // Shower position in the FCAL
  double fTimeTr;
  double fDocaTr;
  double fSumU;
  double fSumV;
  double fE9E25;
  double fE1E9;
  int iNumBlocks;
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

inline double DFCALShower::getDocaTrack() const { return fDocaTr; }
inline double DFCALShower::getTimeTrack() const { return fTimeTr; }
inline double DFCALShower::getSumU() const { return fSumU; }
inline double DFCALShower::getSumV() const { return fSumV; }
inline double DFCALShower::getE9E25() const { return fE9E25; }
inline double DFCALShower::getE1E9() const { return fE1E9; }
inline int DFCALShower::getNumBlocks() const { return iNumBlocks; }

#endif // _DFCALShower_

