#include "PID/DKinematicData.h"

void DKinematicData::Reset(void)
{
	setPID(Unknown);
	setMass(0.0);
	setCharge(0.0);

	setMomentum(DVector3());
	setPosition(DVector3());

	setTime(0.0);
	setT0(0.0, 0.0, SYS_NULL);
	setT1(0.0, 0.0, SYS_NULL);

	setPathLength(0.0, 0.0);
	setdEdx(0.0);
	setForwardParmFlag(false);
	setTrackingStateVector(0.0, 0.0, 0.0, 0.0, 0.0);

	//Other classes are responsible for managing the error matrix memory!
	m_errorMatrix = nullptr;
	m_TrackingErrorMatrix = nullptr;

	ClearAssociatedObjects();
}
