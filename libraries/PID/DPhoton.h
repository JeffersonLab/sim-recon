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
    HDCLASSDEF(DPhoton);


    enum { kDefaultTag = 0 ,
           kDefaultDistance = 1000
    } ;
    

    DPhoton();
    DPhoton( const oid_t id );
    ~DPhoton();

    // getters             
    // Introducing detection point in calorimeter to avoid confusion between 
    // measured cluster position and vertex, 
    // which is position() in terms of DKinemtaicData.
    DVector3 getPositionCal() const; 

    unsigned int getTag() const; 
    double getDtRT() const; 
    double getdThetaCharge() const; 

    // setters
    void setTag(unsigned int tag);  
    void setDtRT(double aDtRT);  
    void setdThetaCharge(double adTheta);  
    void setPositionCal( const DVector3& aPosition );
    void makeErrorMatrix( const DMatrixDSym& aSigmas );  

  private:

    unsigned int fTag; //Photon origin (FCAL/BCAL/charged 1/2/3))
    double fDtRT; //Distance to closest track's ReferenceTrajectory
    double fdThetaCharge; //Distance to closest generated charge in polar angle
    DVector3 fPositionCal; // position in calorimeter

};


inline DVector3 DPhoton::getPositionCal() const
{
  return fPositionCal;
}


inline unsigned int DPhoton::getTag() const
{
  return fTag;
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

