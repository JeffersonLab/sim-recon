// $Id$
//
//    File: DChargedTrackHypothesis_factory.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DChargedTrackHypothesis_
#define _DChargedTrackHypothesis_

#include <vector>
#include <JANA/JObject.h>
#include <TRACKING/DTrackTimeBased.h>
#include <particleType.h>

using namespace std;

class DChargedTrackHypothesis : public jana::JObject {
	public:
		JOBJECT_PUBLIC(DChargedTrackHypothesis);

		const DTrackTimeBased* dTrackTimeBased;
		Particle_t dPID;
		float dProjectedTime; //Time at the track position in the DTrackTimeBased object, calculated from matching to either the FCAL, BCAL, or TOF
		float dPathLength; //Path length from the track position in the DTrackTimeBased object to the matched hit in either the FCAL, BCAL, or TOF
		float dFlightTime; //The amount of time that the track took to traverse the dPathLength
		DetectorSystem_t dMatchedTimeDetector;

		unsigned int dNDF_Timing;
		float dChiSq_Timing;

		unsigned int dNDF_DCdEdx;
		float dChiSq_DCdEdx;

		unsigned int dNDF;
		float dChiSq;
		float dFOM;

		float mass() const{return ParticleMass(dPID);} //this may be different than the value in DTrackTimeBased!!
		float charge() const{return ParticleCharge(dPID);} //this may be different than the value in DTrackTimeBased!!
		DVector3 momentum() const{return dTrackTimeBased->momentum();}
		DVector3 position() const{return dTrackTimeBased->position();}
		float energy() const{ //this may be different than the value in DTrackTimeBased!!
			float locMomentum = momentum().Mag();
			float locMass = mass();
			return sqrt(locMomentum*locMomentum + locMass*locMass);
		}
		float beta() const{return momentum().Mag()/energy();}

		void toStrings(vector<pair<string,string> > &items) const{
			AddString(items, "PID", "%d", int(dPID));

			AddString(items, "q", "%+1.0f", charge());
			AddString(items, "x(cm)", "%3.1f", position().x());
			AddString(items, "y(cm)", "%3.1f", position().y());
			AddString(items, "z(cm)", "%3.1f", position().z());
			AddString(items, "E(GeV)", "%2.3f", energy());
			AddString(items, "t(ns)", "%2.3f", dTrackTimeBased->t0());
			AddString(items, "p(GeV/c)", "%2.3f", momentum().Mag());
			AddString(items, "theta(deg)", "%2.3f", momentum().Theta()*180.0/M_PI);
			AddString(items, "phi(deg)", "%2.3f", momentum().Phi()*180.0/M_PI);

			AddString(items, "T_Proj", "%3.5f", dProjectedTime);
			AddString(items, "Path", "%3.2f", dPathLength);
			AddString(items, "TOF", "%3.5f", dFlightTime);
			AddString(items, "PID_ChiSq", "%f", dChiSq);
			AddString(items, "PID_FOM", "%f", dFOM);
		}

};

#endif // _DChargedTrackHypothesis_

