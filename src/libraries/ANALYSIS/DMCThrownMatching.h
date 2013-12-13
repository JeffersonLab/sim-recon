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
#include "PID/DNeutralShower.h"

#include "TOF/DTOFPoint.h"
#include "TOF/DTOFTruth.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALTruthShower.h"
#include "FCAL/DFCALShower.h"
#include "FCAL/DFCALTruthShower.h"

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

		inline void Set_TOFPointToTruthMap(const map<const DTOFPoint*, const DTOFTruth*>& locTOFPointToTruthMap){dTOFPointToTruthMap = locTOFPointToTruthMap;}
		inline void Set_TOFTruthToPointMap(const map<const DTOFTruth*, const DTOFPoint*>& locTOFTruthToPointMap){dTOFTruthToPointMap = locTOFTruthToPointMap;}

		inline void Set_BCALShowerToTruthMap(map<const DBCALShower*, const DBCALTruthShower*>& locBCALShowerToTruthMap){dBCALShowerToTruthMap = locBCALShowerToTruthMap;}
		inline void Set_BCALTruthToShowerMap(map<const DBCALTruthShower*, const DBCALShower*>& locBCALTruthToShowerMap){dBCALTruthToShowerMap = locBCALTruthToShowerMap;}

		inline void Set_FCALShowerToTruthMap(map<const DFCALShower*, const DFCALTruthShower*>& locFCALShowerToTruthMap){dFCALShowerToTruthMap = locFCALShowerToTruthMap;}
		inline void Set_FCALTruthToShowerMap(map<const DFCALTruthShower*, const DFCALShower*>& locFCALTruthToShowerMap){dFCALTruthToShowerMap = locFCALTruthToShowerMap;}

		//GETTERS: INDIVIDUAL PARTICLES

		//the below two functions return the hypothesis with PID = MC PID. if not available, returns one with best PID FOM
		const DChargedTrackHypothesis* Get_MatchingChargedHypothesis(const DMCThrown* locInputMCThrown) const;
		const DNeutralParticleHypothesis* Get_MatchingNeutralHypothesis(const DMCThrown* locInputMCThrown) const;

		void Get_MatchingChargedHypotheses(const DMCThrown* locInputMCThrown, deque<const DChargedTrackHypothesis*>& locMatchingChargedHypotheses) const;
		const DChargedTrack* Get_MatchingChargedTrack(const DMCThrown* locInputMCThrown) const;

		void Get_MatchingNeutralHypotheses(const DMCThrown* locInputMCThrown, deque<const DNeutralParticleHypothesis*>& locMatchingNeutralHypotheses) const;
		const DNeutralParticle* Get_MatchingNeutralParticle(const DMCThrown* locInputMCThrown) const;

		const DMCThrown* Get_MatchingMCThrown(const DChargedTrackHypothesis* locChargedTrackHypothesis) const;
		const DMCThrown* Get_MatchingMCThrown(const DChargedTrack* locChargedTrack) const;

		const DMCThrown* Get_MatchingMCThrown(const DNeutralParticleHypothesis* locNeutralParticleHypothesis) const;
		const DMCThrown* Get_MatchingMCThrown(const DNeutralParticle* locNeutralParticle) const;

		//GETTERS: INDIVIDUAL HITS
		const DTOFPoint* Get_MatchingTOFPoint(const DTOFTruth* locTOFTruth) const;
		const DTOFTruth* Get_MatchingTOFTruth(const DTOFPoint* locTOFPoint) const;

		const DBCALShower* Get_MatchingBCALShower(const DBCALTruthShower* locBCALTruthShower) const;
		const DBCALTruthShower* Get_MatchingBCALTruthShower(const DBCALShower* locBCALShower) const;

		const DFCALShower* Get_MatchingFCALShower(const DFCALTruthShower* locFCALTruthShower) const;
		const DFCALTruthShower* Get_MatchingFCALTruthShower(const DFCALShower* locFCALShower) const;

		//GETTERS: WHOLE MAPS
		inline void Get_ChargedHypoToThrownMap(map<const DChargedTrackHypothesis*, const DMCThrown*>& locChargedHypoToThrownMap) const{locChargedHypoToThrownMap = dChargedHypoToThrownMap;}
		inline void Get_ThrownToChargedHypoMap(map<const DMCThrown*, deque<const DChargedTrackHypothesis*> >& locThrownToChargedHypoMap) const{locThrownToChargedHypoMap = dThrownToChargedHypoMap;}

		inline void Get_ChargedToThrownMap(map<const DChargedTrack*, const DMCThrown*>& locChargedToThrownMap) const{locChargedToThrownMap = dChargedToThrownMap;}
		inline void Get_ThrownToChargedMap(map<const DMCThrown*, const DChargedTrack*>& locThrownToChargedMap) const{locThrownToChargedMap = dThrownToChargedMap;}

		inline void Get_NeutralHypoToThrownMap(map<const DNeutralParticleHypothesis*, const DMCThrown*>& locNeutralHypoToThrownMap) const{locNeutralHypoToThrownMap = dNeutralHypoToThrownMap;}
		inline void Get_ThrownToNeutralHypoMap(map<const DMCThrown*, deque<const DNeutralParticleHypothesis*> >& locThrownToNeutralHypoMap) const{locThrownToNeutralHypoMap = dThrownToNeutralHypoMap;}

		inline void Get_NeutralToThrownMap(map<const DNeutralParticle*, const DMCThrown*>& locNeutralToThrownMap) const{locNeutralToThrownMap = dNeutralToThrownMap;}
		inline void Get_ThrownToNeutralMap(map<const DMCThrown*, const DNeutralParticle*>& locThrownToNeutralMap) const{locThrownToNeutralMap = dThrownToNeutralMap;}

		inline void Get_TOFPointToTruthMap(map<const DTOFPoint*, const DTOFTruth*>& locTOFPointToTruthMap) const{locTOFPointToTruthMap = dTOFPointToTruthMap;}
		inline void Get_TOFTruthToPointMap(map<const DTOFTruth*, const DTOFPoint*>& locTOFTruthToPointMap) const{locTOFTruthToPointMap = dTOFTruthToPointMap;}

		inline void Get_BCALShowerToTruthMap(map<const DBCALShower*, const DBCALTruthShower*>& locBCALShowerToTruthMap) const{locBCALShowerToTruthMap = dBCALShowerToTruthMap;}
		inline void Get_BCALTruthToShowerMap(map<const DBCALTruthShower*, const DBCALShower*>& locBCALTruthToShowerMap) const{locBCALTruthToShowerMap = dBCALTruthToShowerMap;}

		inline void Get_FCALShowerToTruthMap(map<const DFCALShower*, const DFCALTruthShower*>& locFCALShowerToTruthMap) const{locFCALShowerToTruthMap = dFCALShowerToTruthMap;}
		inline void Get_FCALTruthToShowerMap(map<const DFCALTruthShower*, const DFCALShower*>& locFCALTruthToShowerMap) const{locFCALTruthToShowerMap = dFCALTruthToShowerMap;}

	private:

		map<const DChargedTrackHypothesis*, const DMCThrown*> dChargedHypoToThrownMap;
		map<const DMCThrown*, deque<const DChargedTrackHypothesis*> > dThrownToChargedHypoMap;

		map<const DChargedTrack*, const DMCThrown*> dChargedToThrownMap;
		map<const DMCThrown*, const DChargedTrack*> dThrownToChargedMap;

		map<const DNeutralParticleHypothesis*, const DMCThrown*> dNeutralHypoToThrownMap;
		map<const DMCThrown*, deque<const DNeutralParticleHypothesis*> > dThrownToNeutralHypoMap;

		map<const DNeutralParticle*, const DMCThrown*> dNeutralToThrownMap;
		map<const DMCThrown*, const DNeutralParticle*> dThrownToNeutralMap;

		map<const DTOFPoint*, const DTOFTruth*> dTOFPointToTruthMap;
		map<const DTOFTruth*, const DTOFPoint*> dTOFTruthToPointMap;

		map<const DBCALShower*, const DBCALTruthShower*> dBCALShowerToTruthMap;
		map<const DBCALTruthShower*, const DBCALShower*> dBCALTruthToShowerMap;

		map<const DFCALShower*, const DFCALTruthShower*> dFCALShowerToTruthMap;
		map<const DFCALTruthShower*, const DFCALShower*> dFCALTruthToShowerMap;
};

#endif // _DMCThrownMatching_

