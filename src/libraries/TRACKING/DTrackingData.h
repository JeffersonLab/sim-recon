#ifndef _DTrackingData_
#define _DTrackingData_

#include <memory>
#include "TMatrixFSym.h"

#include <PID/DKinematicData.h>

using namespace std;
using namespace jana;

class DTrackingData : public DKinematicData
{
	public:

		// constructors and destructor
		DTrackingData(void);
		DTrackingData(const DTrackingData& locSourceData, bool locShareTrackingFlag = false, bool locShareKinematicsFlag = false);
		DTrackingData(const DKinematicData& locSourceData, bool locShareKinematicsFlag = false);
		virtual ~DTrackingData(void) {};
		void Reset(void);

		//Assignment operator
		DTrackingData& operator=(const DTrackingData& locSourceData);

		//GETTERS
		const TMatrixFSym* TrackingErrorMatrix(void) const{return dTrackingInfo->m_TrackingErrorMatrix;}
		bool forwardParmFlag(void) const{return dTrackingInfo->m_use_forward_parameters;}
		void TrackingStateVector(double aVec[5]) const;

		//SETTERS
		void setForwardParmFlag(bool aFlag){dTrackingInfo->m_use_forward_parameters = aFlag;}
		void setTrackingErrorMatrix(const TMatrixFSym* aMatrix){dTrackingInfo->m_TrackingErrorMatrix = aMatrix;}
		void setTrackingStateVector(double a1, double a2, double a3, double a4, double a5);

	private:

		struct DTrackingInfo
		{
			// NONE OF THIS DEPENDS ON THE KINEMATIC FIT
			// so, this can be stored as a pointer & shared between multiple objects instead of stored separately for each
			const TMatrixFSym *m_TrackingErrorMatrix = nullptr;  // order is q/pt,phi,tanl,D,z
			bool m_use_forward_parameters = false; // Flag indicating the use of the forward parameterization (x,y,tx,ty,q/p)
			double m_TrackingStateVector[5] = {0.0, 0.0, 0.0, 0.0, 0.0}; // order is q/pt,phi,tanl,D,z
		};

		//memory of object in shared_ptr is managed automatically: deleted automatically when no references are left
		shared_ptr<DTrackingInfo> dTrackingInfo;
};

/************************************************************** CONSTRUCTORS & OPERATORS ***************************************************************/

inline DTrackingData::DTrackingData(void) : dTrackingInfo(make_shared<DTrackingInfo>()){}

inline DTrackingData::DTrackingData(const DTrackingData& locSourceData, bool locShareTrackingFlag, bool locShareKinematicsFlag) :
DKinematicData(locSourceData, locShareKinematicsFlag)
{
	//Default is NOT to share: create a new, independent copy of the input data (tracked separately from input so it can be modified)
	dTrackingInfo = locShareTrackingFlag ? locSourceData.dTrackingInfo : std::make_shared<DTrackingInfo>(*(locSourceData.dTrackingInfo));
}

inline DTrackingData::DTrackingData(const DKinematicData& locSourceData, bool locShareKinematicsFlag) :
DKinematicData(locSourceData, locShareKinematicsFlag), dTrackingInfo(make_shared<DTrackingInfo>()) {}

inline DTrackingData& DTrackingData::operator=(const DTrackingData& locSourceData)
{
	//Replace current data with a new, independent copy of the input data: tracked separately from input so it can be modified
	dTrackingInfo = make_shared<DTrackingInfo>(*(locSourceData.dTrackingInfo));
	return *this;
}

/********************************************************************** GETTERS ************************************************************************/

inline void DTrackingData::TrackingStateVector(double aVec[5]) const
{
	for (unsigned int i = 0; i < 5; ++i)
		aVec[i] = dTrackingInfo->m_TrackingStateVector[i];
}

/********************************************************************** SETTERS ************************************************************************/

inline void DTrackingData::setTrackingStateVector(double a1, double a2, double a3, double a4, double a5)
{
	dTrackingInfo->m_TrackingStateVector[0]=a1;
	dTrackingInfo->m_TrackingStateVector[1]=a2;
	dTrackingInfo->m_TrackingStateVector[2]=a3;
	dTrackingInfo->m_TrackingStateVector[3]=a4;
	dTrackingInfo->m_TrackingStateVector[4]=a5;
}

inline void DTrackingData::Reset(void)
{
	DKinematicData::Reset();
	dTrackingInfo = make_shared<DTrackingInfo>(); //not safe to reset individually, since you don't know what it's shared with
}

#endif /* _DTrackingData_ */
