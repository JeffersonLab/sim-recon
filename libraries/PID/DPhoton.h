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
// Intrioducing detection point in calorimeter to avoid confusion between 
// measured cluster position and vertex, 
// which is position() in terms of DKinemtaicData.
		DVector3 positionCal() const; 

                oid_t getID() const;  // returns JANA object ID
                unsigned int getTag() const; 
		double getDtRT() const; 

		// setters
                void setTag(unsigned int tag);  
                void setDtRT(double aDtRT);  
                void setPositionCal( const DVector3& aPosition );
                void makeErrorMatrix( const DMatrixDSym& aSigmas );  
      
	private:

//                DVector3 fVertex;  // Photon vertex (set to 0,0,0 for the moment)
                unsigned int fTag; //Photon origin (FCAL/BCAL 0/1))
                double fDtRT; //Distance to closest track's RefenceTrajectory
                DVector3 fPositionCal; // position in calorimeter

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

