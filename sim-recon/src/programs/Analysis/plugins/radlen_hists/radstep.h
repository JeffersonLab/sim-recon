// $Id$
//
//    File: radstep.h
// Created: Fri Apr 17 13:27:38 EDT 2009
// Creator: davidl
//

#ifndef _radstep_
#define _radstep_

#include <TObject.h>
#include <TVector3.h>


class radstep:public TObject{

	public:
		TVector3 pthrown;
		TVector3 pos;
		TVector3 B;
		double s;
		double stot;
		double Xo;
		double ix_over_Xo;
		double iB_cross_p_dl; // integral of magnitude of B x p_hat dl up to and including this step in Tesla-cm
		double iB_dl; // integral of magnitude of B times the step dl in Tesla-cm

	private:
		ClassDef(radstep,1);

};

// Notes:
//
// The iB_cross_p_dl value is the running sum over steps of the magnitude of the cross
// product of the B-field vector and the unit vector in the momentum direction
// times the step size. This is the relevant quantity for assessing the analyzing
// power of the field for the current trajectory.
//
// The iB_dl value is just the running sum over steps of the magnitude of the B-field
// times the step size. This gives and idea of how much field is "seen" by a particle
// even the component parallel to the particle direction. This is more typical for
// specifying a magnet's strength. Therefore, this may be most relevant for particles
// going straight down the beamline.
// 

#endif // _radstep_

