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

		//GETTERS
		const TMatrixFSym* TrackingErrorMatrix(void) const{return m_TrackingErrorMatrix;}
		bool forwardParmFlag(void) const{return m_use_forward_parameters;}
		void TrackingStateVector(double aVec[5]) const;
		double t0(void) const{return dt0;}
		double t0_err(void) const{return dt0_err;}
		DetectorSystem_t t0_detector(void) const{return dt0_detector;}

		//SETTERS
		void setForwardParmFlag(bool aFlag){m_use_forward_parameters = aFlag;}
		void setTrackingErrorMatrix(const TMatrixFSym* aMatrix){m_TrackingErrorMatrix = aMatrix;}
		void setTrackingStateVector(double a1, double a2, double a3, double a4, double a5);
		void setT0(double at0, double at0_err, DetectorSystem_t at0_detector);

	private:

		const TMatrixFSym *m_TrackingErrorMatrix = nullptr;  // order is q/pt,phi,tanl,D,z
		bool m_use_forward_parameters = false; // Flag indicating the use of the forward parameterization (x,y,tx,ty,q/p)
		double m_TrackingStateVector[5] = {0.0, 0.0, 0.0, 0.0, 0.0}; // order is q/pt,phi,tanl,D,z

		double dt0 = 0.0;
		double dt0_err = 0.0;
		DetectorSystem_t dt0_detector = SYS_NULL;
};

/********************************************************************** GETTERS ************************************************************************/

inline void DTrackingData::TrackingStateVector(double aVec[5]) const
{
	for (unsigned int i = 0; i < 5; ++i)
		aVec[i] = m_TrackingStateVector[i];
}

/********************************************************************** SETTERS ************************************************************************/

inline void DTrackingData::setTrackingStateVector(double a1, double a2, double a3, double a4, double a5)
{
	m_TrackingStateVector[0]=a1;
	m_TrackingStateVector[1]=a2;
	m_TrackingStateVector[2]=a3;
	m_TrackingStateVector[3]=a4;
	m_TrackingStateVector[4]=a5;
}

inline void DTrackingData::setT0(double at0, double at0_err, DetectorSystem_t at0_detector)
{
	dt0 = at0;
	dt0_err = at0_err;
	dt0_detector = at0_detector;
}

#endif /* _DTrackingData_ */
