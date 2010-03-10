//
//    File: DPhoton.h
// Created: Tue Apr 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#ifndef _DPhoton_
#define _DPhoton_

//#include <TLorentzVector.h>
#include "DKinematicData.h"

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

//class DPhoton: public JObject {
class DPhoton: public DKinematicData {
  public:
    JOBJECT_PUBLIC(DPhoton);


    enum { // kDefaultTag = 0 ,
           kDefaultDistance = 1000
    } ;

    enum PhotonTag {  kDefaultTag, kFcal, kBcal, kCharge }; 
    

    DPhoton();
    DPhoton( const oid_t id );
    ~DPhoton();

    // getters             
    // Introducing detection point in calorimeter to avoid confusion between 
    // measured cluster position and vertex.
    // The last can be obtained by position() method of DKinemtaicData.
    DVector3 getPositionCal() const; 

    //unsigned int getTag() const; 
    PhotonTag getTag() const; 
    double getDtRT() const; 
    double getdThetaCharge() const; 
    double getTime() const; 

    // setters
    void setTag( PhotonTag aTag );  
    void setTime( double aTime );  
    void setDtRT( double aDtRT );  
    void setdThetaCharge( double adTheta );  
    void setPositionCal( const DVector3& aPosition );
    void makeErrorMatrix( const DMatrixDSym& aSigmas );  

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "E(GeV)", "%5.2f", energy());
			AddString(items, "Px(GeV/c)", "%5.2f", momentum().X());
			AddString(items, "Py(GeV/c)", "%5.2f", momentum().Y());
			AddString(items, "Pz(GeV/c)", "%5.2f", momentum().Z());
			AddString(items, "X(cm)", "%7.2f", position().X());
			AddString(items, "Y(cm)", "%7.2f", position().Y());
			AddString(items, "Z(cm)", "%7.2f", position().Z());
			AddString(items, "t(ns)", "%7.2f", getTime());
			AddString(items, "Tag", "%5i", getTag());
		}

  private:

    //unsigned int fTag; //Photon origin (FCAL/BCAL/charged 1/2/3))
    PhotonTag fTag; 
    double fTime; 
    double fDtRT; //Distance to closest track's ReferenceTrajectory
    double fdThetaCharge; //Distance to closest generated charge in polar angle
    DVector3 fPositionCal; // position in calorimeter

};


inline DVector3 DPhoton::getPositionCal() const
{
  return fPositionCal;
}


inline DPhoton::PhotonTag DPhoton::getTag() const
{
  return fTag;
}

inline double DPhoton::getTime() const
{
  return fTime;
}

inline double DPhoton::getDtRT() const
{
  return fDtRT;
}

inline double DPhoton::getdThetaCharge() const
{
  return fdThetaCharge;
}
#endif // _DPhoton_

