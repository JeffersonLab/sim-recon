#ifndef _DMCThrownMatching_
#define _DMCThrownMatching_

#include <map>
#include <deque>

#include "JANA/JObject.h"
#include "TRACKING/DMCThrown.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DChargedTrack.h"
#include "PID/DNeutralParticleHypothesis.h"
#include "PID/DNeutralParticle.h"

using namespace std;
using namespace jana;

class DMCThrownMatching : public JObject
{
	//uses measured tracks for matching, not kinfit ones
	public:

		JOBJECT_PUBLIC(DMCThrownMatching);

		//SETTERS
		inline void Set_ChargedHypoToThrownMap(const map<const DChargedTrackHypothesis*, const DMCThrown*>& locChargedHypoToThrownMap){dChargedHypoToThrownMap = locChargedHypoToThrownMap;}
		inline void Set_ThrownToChargedHypoMap(const map<const DMCThrown*, deque<const DChargedTrackHypothesis*> >& locThrownToChargedHypoMap){dThrownToChargedHypoMap = locThrownToChargedHypoMap;}

		inline void Set_ChargedToThrownMap(const map<const DChargedTrack*, const DMCThrown*>& locChargedToThrownMap){dChargedToThrownMap = locChargedToThrownMap;}
		inline void Set_ThrownToChargedMap(const map<const DMCThrown*, const DChargedTrack*>& locThrownToChargedMap){dThrownToChargedMap = locThrownToChargedMap;}

		inline void Set_NeutralHypoToThrownMap(const map<const DNeutralParticleHypothesis*, const DMCThrown*>& locNeutralHypoToThrownMap){dNeutralHypoToThrownMap = locNeutralHypoToThrownMap;}
		inline void Set_ThrownToNeutralHypoMap(const map<const DMCThrown*, deque<const DNeutralParticleHypothesis*> >& locThrownToNeutralHypoMap){dThrownToNeutralHypoMap = locThrownToNeutralHypoMap;}

		inline void Set_NeutralToThrownMap(const map<const DNeutralParticle*, const DMCThrown*>& locNeutralToThrownMap){dNeutralToThrownMap = locNeutralToThrownMap;}
		inline void Set_ThrownToNeutralMap(const map<const DMCThrown*, const DNeutralParticle*>& locThrownToNeutralMap){dThrownToNeutralMap = locThrownToNeutralMap;}

		//GETTERS: INDIVIDUAL PARTICLES
		void Get_MatchingChargedHypotheses(const DMCThrown* locInputMCThrown, deque<const DChargedTrackHypothesis*>& locMatchingChargedHypotheses) const;
		const DChargedTrack* Get_MatchingChargedTrack(const DMCThrown* locInputMCThrown) const;

		void Get_MatchingNeutralHypotheses(const DMCThrown* locInputMCThrown, deque<const DNeutralParticleHypothesis*>& locMatchingNeutralHypotheses) const;
		const DNeutralParticle* Get_MatchingNeutralParticle(const DMCThrown* locInputMCThrown) const;

		const DMCThrown* Get_MatchingMCThrown(const DChargedTrackHypothesis* locChargedTrackHypothesis) const;
		const DMCThrown* Get_MatchingMCThrown(const DChargedTrack* locChargedTrack) const;

		const DMCThrown* Get_MatchingMCThrown(const DNeutralParticleHypothesis* locNeutralParticleHypothesis) const;
		const DMCThrown* Get_MatchingMCThrown(const DNeutralParticle* locNeutralParticle) const;

		//GETTERS: WHOLE MAPS
		inline void Get_ChargedHypoToThrownMap(map<const DChargedTrackHypothesis*, const DMCThrown*>& locChargedHypoToThrownMap) const{locChargedHypoToThrownMap = dChargedHypoToThrownMap;}
		inline void Get_ThrownToChargedHypoMap(map<const DMCThrown*, deque<const DChargedTrackHypothesis*> >& locThrownToChargedHypoMap) const{locThrownToChargedHypoMap = dThrownToChargedHypoMap;}

		inline void Get_ChargedToThrownMap(map<const DChargedTrack*, const DMCThrown*>& locChargedToThrownMap) const{locChargedToThrownMap = dChargedToThrownMap;}
		inline void Get_ThrownToChargedMap(map<const DMCThrown*, const DChargedTrack*>& locThrownToChargedMap) const{locThrownToChargedMap = dThrownToChargedMap;}

		inline void Get_NeutralHypoToThrownMap(map<const DNeutralParticleHypothesis*, const DMCThrown*>& locNeutralHypoToThrownMap) const{locNeutralHypoToThrownMap = dNeutralHypoToThrownMap;}
		inline void Get_ThrownToNeutralHypoMap(map<const DMCThrown*, deque<const DNeutralParticleHypothesis*> >& locThrownToNeutralHypoMap) const{locThrownToNeutralHypoMap = dThrownToNeutralHypoMap;}

		inline void Get_NeutralToThrownMap(map<const DNeutralParticle*, const DMCThrown*>& locNeutralToThrownMap) const{locNeutralToThrownMap = dNeutralToThrownMap;}
		inline void Get_ThrownToNeutralMap(map<const DMCThrown*, const DNeutralParticle*>& locThrownToNeutralMap) const{locThrownToNeutralMap = dThrownToNeutralMap;}

	private:

		map<const DChargedTrackHypothesis*, const DMCThrown*> dChargedHypoToThrownMap;
		map<const DMCThrown*, deque<const DChargedTrackHypothesis*> > dThrownToChargedHypoMap;

		map<const DChargedTrack*, const DMCThrown*> dChargedToThrownMap;
		map<const DMCThrown*, const DChargedTrack*> dThrownToChargedMap;

		map<const DNeutralParticleHypothesis*, const DMCThrown*> dNeutralHypoToThrownMap;
		map<const DMCThrown*, deque<const DNeutralParticleHypothesis*> > dThrownToNeutralHypoMap;

		map<const DNeutralParticle*, const DMCThrown*> dNeutralToThrownMap;
		map<const DMCThrown*, const DNeutralParticle*> dThrownToNeutralMap;
};


#endif // _DMCThrownMatching_

