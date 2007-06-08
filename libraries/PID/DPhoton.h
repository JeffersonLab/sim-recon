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

//class DPhoton: public JObject {
class DPhoton: public DKinematicData {
	public:
		HDCLASSDEF(DPhoton);
                
                DPhoton();
		~DPhoton();
                               
		DVector3 getPosition() const; 
		DVector3 getMomentum() const; 
		double getEnergy() const; 
		DLorentzVector getMom4() const; 
		DMatrixDSym* getErrorMatrix() const; 
                unsigned int getTag() const;
		double getDtRT() const; 
                void setPosition(const DVector3 aPosition);  
                void setMomentum(const DVector3 aMom);  
                void setEnergy(const double aEnergy);  
                void setTag(const unsigned int tag);  
                void setDtRT(const double aDtRT);  
                void setErrorMatrix(const DVector3 aPosition, const DVector3 aVertex, const double aEnergy, DMatrixDSym* aSigmas );  
      
	private:

                double fEnergy;  // Photon energy
                DVector3 fPosition;  // Photon position
                DVector3 fVertex;  // Photon vertex (set to zero for the moment)
                DVector3 fMomentum;  // Photon 3-momentum
                //DLorentzVector fMom4;  // Photon 4-momentum
		DMatrixDSym* fErrorMatrix; 
                unsigned int fTag; //Photon origin (FCAL/BCAL 0/1))
                double fDtRT; //Distance to closest track's RefenceTrajectory

};

inline double DPhoton::getEnergy() const
{
      return fEnergy;
}

inline DVector3 DPhoton::getPosition() const
{
      return fPosition;
}

inline DVector3 DPhoton::getMomentum() const
{
      return fMomentum;
}

inline DLorentzVector DPhoton::getMom4() const
{
      return DLorentzVector( fMomentum, fEnergy );
}

inline unsigned int DPhoton::getTag() const
{
      return fTag;
}

inline double DPhoton::getDtRT() const
{
      return fDtRT;
}

inline DMatrixDSym* DPhoton::getErrorMatrix() const
{
      return fErrorMatrix;
}
#endif // _DPhoton_

