#include <deque>
#include <set>
#include <unordered_map>
#include <utility>
#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>

#include <JANA/JEventLoop.h>

#include <particleType.h>
#include <PID/DChargedTrackHypothesis.h>
#include <PID/DChargedTrack.h>
#include <PID/DNeutralShower.h>
#include <ANALYSIS/DReaction.h>
#include <ANALYSIS/DKinFitUtils_GlueX.h>

using namespace std;


class DVertexInfo
{
	public:
		DVertexInfo(vector<Particle_t> locDetectedPIDs) : dDetectedPIDs(locDetectedPIDs)
		{
			std::sort(dDetectedPIDs.begin(), dDetectedPIDs.end()); //so that < operator works!
		}

		void Add_DecayingPIDVertex(Particle_t locDecayingPID, shared_ptr<DVertexInfo>& locConstrainingVertex)
		{
			dDecayingPIDVertices.emplace_back(locDecayingPID, locConstrainingVertex);
			//expensive to sort for every new entry; however, this will probably only be called once or maybe twice (if at all)
			std::sort(dDecayingPIDVertices.begin(), dDecayingPIDVertices.end());  //so that < operator works!
		}

		vector<Particle_t> Get_DetectedPIDs(void){return dDetectedPIDs;}
		vector<pair<Particle_t, shared_ptr<DVertexInfo> > > Get_DecayingPIDVertices(void){return dDecayingPIDVertices;}

		bool operator<(const DVertexInfo& locOtherInfo) const;

	private:

		vector<Particle_t> dDetectedPIDs; //at vertex, charged, & detected //sorted, so no longer in order from DReactionStep's

		//e.g. if Lambda is needed at this vertex, but its position is constrained at the other vertex, then this points to the other vertex
		vector<pair<Particle_t, shared_ptr<DVertexInfo> > > dDecayingPIDVertices;
};

struct DVertexInfoComparer
{
	bool operator() (const DVertexInfo*& lhs, const DVertexInfo*& rhs) const{return *lhs < *rhs;}
};

class DComboVertex
{
	public:
		DComboVertex(shared_ptr<DVertexInfo>& locVertexInfo, vector<const DChargedTrackHypothesis*>& locComboHypos) :
			dVertexInfo(locVertexInfo), dComboHypos(locComboHypos){}

		void Add_DecayingPIDVertex(Particle_t locDecayingPID, shared_ptr<DComboVertex*>& locComboVertex)
		{
			dDecayingPIDVertices.emplace_back(locDecayingPID, locComboVertex);
			//is this still necessary
			std::sort(dDecayingPIDVertices.begin(), dDecayingPIDVertices.end());  //so that < operator works!
		}

		void Set_Vertex(TVector3 locVertex){dVertex = locVertex;}
		TVector3 Get_Vertex(void) const{return dVertex;}

		vector<const DChargedTrackHypothesis*> Get_ComboHypos(void) const{return dComboHypos;}

	private:
		shared_ptr<DVertexInfo> dVertexInfo;
		vector<const DChargedTrackHypothesis*> dComboHypos; //detected //sorted, so no longer in order from DVertexInfo::dDetectedPIDs
		TVector3 dVertex;

		//e.g. if Lambda is needed at this vertex, but its position is constrained at the other vertex, then this points to the other vertex
		vector<pair<Particle_t, shared_ptr<DComboVertex*> > > dDecayingPIDVertices;
};






class DVertexCreator
{
	DVertexCreator(JEventLoop* locEventLoop);

	private:
		DKinFitUtils_GlueX* dKinFitUtils;
		vector<Particle_t> dTrackingPIDs;

		//CONTROL INFORMATION
		bool dUseSigmaForRFSelectionFlag; //if false, all votes created equal
		string dShowerSelectionTag;

		//EXPERIMENT INFORMATION
		double dTargetCenterZ;
		double dBeamBunchPeriod;

		//STANDARDIZED CUTS
		map<pair<Particle_t, DetectorSystem_t>, double> dPIDTimeCutMap;

		//Store DVertexInfos
		auto DVertexInfoComparer = [](const DVertexInfo*& lhs, const DVertexInfo*& rhs) -> bool {return *lhs < *rhs;};
		auto dAllVertexInfos_Set = set<shared_ptr<DVertexInfo> >(DVertexInfoComparer); //for finding duplicate vertices
		vector<shared_ptr<DVertexInfo> > dAllVertexInfos_Vector; //for creation (stored in dependency order)
		vector<pair<shared_ptr<DVertexInfo>, vector<const DReaction*> > > dProductionVertexInfos;
		unordered_map<shared_ptr<DVertexInfo>, unordered_map<const DReaction*, vector<const DReactionStep*> > > dReactionVertexMap;

		//created combo vertices
		unordered_map<shared_ptr<DVertexInfo>, vector<shared_ptr<DComboVertex> > > dComboVertices;

		//EVENT INFORMATION //so that they don't have to be passed around everywhere
		const DESSkimData* dESSkimData;
		vector<const DChargedTrack*> dChargedTracks;
		vector<const DNeutralShower*> dNeutralShowers;
		const DEventRFBunch* dEventRFBunch; //reaction-independent version

		//maps to facilitate hypo searching for combo building
		unordered_map<Particle_t, vector<const DChargedTrackHypothesis*> > dTrackMap_ByPID;
		unordered_map<const DChargedTrackHypothesis*, unordered_map<Particle_t, vector<const DChargedTrackHypothesis*> > > dTrackMap_DeltaTMatch; //sorted by PID

		//when looking for a given PID (first key) for a combo, and you've already selected a few hypos for that combo (vector key)
			//this stores the vector of valid hypos for that PID (not from the same track as before, in time with other tracks)
			//this value vector is the intersection of the combo-hypo vectors from dTrackMap_DeltaTMatch
			//once it's needed, the intersection is stored here, so that if it's needed again later, it doesn't have to be recomputed for other DReactions
		unordered_map<Particle_t, unordered_map<vector<const DChargedTrackHypothesis*>, vector<const DChargedTrackHypothesis*> > > dTrackSearchVectors;



		//UTILITY FUNCTIONS
		int Calc_RFBunchShift(double locTimeToStep, double locTimeToStepTo) const; //returns integer shift
		unordered_map<int, double> Calc_NeutralRFDeltaTs(const DNeutralShower* locNeutralShower, const TVector3& locVertex, double locRFTime) const;

		DVertexCreator(void);
};

inline int DVertexCreator::Calc_RFBunchShift(double locTimeToStep, double locTimeToStepTo) const
{
	double locDeltaT = locTimeToStepTo - locTimeToStep;
	return (locDeltaT > 0.0) ? int(locDeltaT/dBeamBunchPeriod + 0.5) : int(locDeltaT/dBeamBunchPeriod - 0.5);
}

class DParticleComboInfo
{
	private:
		const DReaction* dReaction; //can get whether beam at prod vertex or not
		vector<DVertexInfo*> dVertexInfos; //in dependence-order
		size_t dNumPhotons;
};

inline bool DVertexInfo::operator<(const DVertexInfo& locOtherInfo) const
{
	//This is used to determine whether the charged-hypo combo-ing, and thus resulting vertex position, will be identical or not
	//So, don't compare all information, just that used to compute the vertex
	if(dBeamAtVertexFlag < locOtherInfo.dBeamAtVertexFlag)
		return true;
	if(dBeamAtVertexFlag > locOtherInfo.dBeamAtVertexFlag)
		return false;

	if(dDetectedPIDs < locOtherInfo.dDetectedPIDs)
		return true;
	if(dDetectedPIDs > locOtherInfo.dDetectedPIDs)
		return false;

	if(dDecayingPIDVertices < locOtherInfo.dDecayingPIDVertices)
		return true;

	return false; //either greater than or equivalent: not less-than
}
