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

#include <PID/DPhoton.h>
#include <PID/DVertex.h>
#include <TRACKING/DTrackTimeBased.h>

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

		const DVertex* vertex;					// vertex
		vector<const DPhoton*> photon;			// photons
		vector<const DTrackTimeBased*> pip;		// pi+
		vector<const DTrackTimeBased*> pim;		// pi-
		vector<const DTrackTimeBased*> proton;	// proton
		vector<const DTrackTimeBased*> Kp;		// K+
		vector<const DTrackTimeBased*> Km;		// K-
		vector<const DTrackTimeBased*> otherp;	// other positively charged tracks (positrons?)
		vector<const DTrackTimeBased*> otherm;	// other positively charged tracks (anti-protons?)

		// There is too much info to fit on a single line here so
		// we limit toStrings to just saying how manyod each type
		// of particle is here. For anything else, the id values
		// for each particle should be kept at least.
		void toStrings(vector<pair<string,string> > &items)const{
			
			//AddString(items, "x", "%3.2f", vertex->x.X());
			//AddString(items, "y", "%3.2f", vertex->x.Y());
			//AddString(items, "z", "%3.2f", vertex->x.Z());
			AddString(items, "Nphoton",      "%d", photon.size());
			AddString(items, "Npi_plus",     "%d", pip.size());
			AddString(items, "Npi_minus",    "%d", pim.size());
			AddString(items, "Nproton",      "%d", proton.size());
			AddString(items, "NK_plus",      "%d", Kp.size());
			AddString(items, "NK_minus",     "%d", Km.size());
			AddString(items, "Nother_plus",  "%d", otherp.size());
			AddString(items, "Nother_minus", "%d", otherm.size());
		}
};

#endif // _DPhysicsEvent_

