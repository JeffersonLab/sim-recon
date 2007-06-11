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
                
		// getters             
//		DVector3 getPosition() const; 
//		DVector3 getMomentum() const; 
//		double getEnergy() const; 
//		DLorentzVector getMom4() const; 
//		DMatrixDSym* getErrorMatrix() const; 
                oid_t getID() const;  // returns JANA object ID
                unsigned int getTag() const; 
		double getDtRT() const; 
//		DVector3 getVertex() const; // this is position now

		// setters
//                void setEnergy(const double aEnergy);  
//                void setPosition(const DVector3 aPosition);  
//                void setMomentum(const DVector3 aMom);  
//                void setVertex(const DVector3& aVertex);  
                void setTag(unsigned int tag);  
                void setDtRT(double aDtRT);  
                void makeErrorMatrix( const DMatrixDSym& aSigmas );  
      
	private:

//                double fEnergy;  // Photon energy
//                DVector3 fPosition;  // Photon position
//                DVector3 fMomentum;  // Photon 3-momentum
//		DMatrixDSym fErrorMatrix; 
//                DVector3 fVertex;  // Photon vertex (set to 0,0,0 for the moment)
                unsigned int fTag; //Photon origin (FCAL/BCAL 0/1))
                double fDtRT; //Distance to closest track's RefenceTrajectory

};

/*inline double DPhoton::getEnergy() const
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
*/

inline oid_t DPhoton::getID() const
{
      return id;
}


/*inline DVector3 DPhoton::getVertex() const
{
      return fVertex;
}*/

inline unsigned int DPhoton::getTag() const
{
      return fTag;
}

inline double DPhoton::getDtRT() const
{
      return fDtRT;
}

#endif // _DPhoton_

