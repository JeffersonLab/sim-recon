#include "PID/DKinematicData.h"

void DKinematicData::Reset(void)
{
	setPID(Unknown);
	setMomentum(DVector3());
	setPosition(DVector3());
	setTime(0.0);

	setForwardParmFlag(false);
	setTrackingStateVector(0.0, 0.0, 0.0, 0.0, 0.0);

	//Other classes are responsible for managing the error matrix memory!
	dKinematics->m_errorMatrix = nullptr;
	dNonKinematics->m_TrackingErrorMatrix = nullptr;

	ClearAssociatedObjects();
}
