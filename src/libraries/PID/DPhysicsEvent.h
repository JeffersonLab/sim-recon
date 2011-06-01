// $Id$
//
//    File: DPhysicsEvent.h
// Created: Wed Aug  4 10:37:55 EDT 2010
// Creator: davidl (on Darwin eleanor.jlab.org 10.4.0 i386)
//

#ifndef _DPhysicsEvent_
#define _DPhysicsEvent_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

#include <PID/DParticleSet.h>

/// Objects of this class are intended to hold collections of
/// reconstructed particles believed to have come from a single
/// physics event. Since it is possible to have more than one
/// physics event in a single DAQ event, this allows an analysis
/// to loop over all physics in a single DAQ event.
///
/// For the charged particles, the DChargedTrack object from
/// which it came is added as an associated object. This, in
/// principle, will allow one to get at the other fit hypotheses
/// for a given track in case they want to try a different
/// configuration.

class DPhysicsEvent:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DPhysicsEvent);

		vector<const DParticleSet*>particle_sets;

		void toStrings(vector<pair<string,string> > &items)const{			
			AddString(items, "Num particle sets", "%d", particle_sets.size());
		}
};

#endif // _DPhysicsEvent_

