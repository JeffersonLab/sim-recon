#ifndef DPhotonComboer_h
#define DPhotonComboer_h

#include <map>
#include <set>

#include "particleType.h"
#include "PID/DNeutralShower.h"

using namespace std;

namespace DAnalysis
{

class DPhotonComboer
{

private:

	//VERTEX-DEPENDENT PHOTON INFORMATION
	//For every 10cm in vertex-z, calculate the photon p4 & time for placing mass & delta-t cuts
	//The z-range extends from the upstream end of the target - 5cm to the downstream end + 15cm
	//so for a 30-cm-long target, it's a range of 50cm: 5bins, evaluated at the center of each bin
	float dPhotonVertexZBinWidth;
	float dPhotonVertexZRangeLow;
	int dNumPhotonVertexZBins;

	//num fully-photon decays needed at production vertex
	//first PID is decay PID, second is product PID, int is # of that product PID
	using DParticleDecay = pair<Particle_t, map<Particle_t, int> >; //if first PID is unknown then is direct #photons (not a decay)

	//int below is # of that decaying particle
	auto DParticleDecayComparer = [](const pair<DParticleDecay, int>& lhs, const pair<DParticleDecay, int>& rhs) -> bool
	{
		//check if decay is the same
		if(lhs.first == rhs.first)
			return lhs.second < rhs.second; //yes, sort by #of type

		//check if decaying PID is the same
		if(lhs.first.first == rhs.first.first)
			return lhs.first.second < rhs.first.second; //yes, sort by decay products (order doesn't really matter)

		//decaying PID different: sort by mass //will automatically get the decaying dependence correct!!
		return (ParticleMass(lhs.first.first) < ParticleMass(rhs.first.first));
	};
	auto dProductionVertexNeutralsSet = set<pair<DParticleDecay, int> >(DParticleDecayComparer);

};

inline int DVertexCreator::Get_PhotonVertexZBin(double locVertexZ) const
{
	//given some vertex-z, what bin am I in?
	int locPhotonVertexZBin = int((locVertexZ - dPhotonVertexZRangeLow)/dPhotonVertexZBinWidth);
	if(locPhotonVertexZBin < 0)
		return 0;
	else if(locPhotonVertexZBin >= dNumPhotonVertexZBins)
		return dNumPhotonVertexZBins - 1;
	return locPhotonVertexZBin;
}

inline DVector3 DVertexCreator::Calc_PhotonP3Time(const DNeutralShower* locNeutralShower, const TVector3& locVertex, double& locVertexTime) const
{
	DVector3 locPath = locNeutralShower->dSpacetimeVertex.Vect() - locVertex;
	double locPathLength = locPath.Mag();
	locVertexTime = locNeutralShower->dSpacetimeVertex.T() - locPathLength/29.9792458;
	return locNeutralShower->dEnergy*locPath.Unit();
}

inline double DVertexCreator::Calc_DeltaTError(const DNeutralShower* locNeutralShower, const shared_ptr<const DKinematicData*>& locKinematicData) const
{
	float& locZError = dPhotonVertexZBinWidth/2.0; //evaluated at center of bin
	double locTheta = locKinematicData->momentum().Theta();
	double locR = locNeutralShower->dSpacetimeVertex.Vect().Perp();
	double locPathError = locR*(1.0/sin(locTheta) - sqrt(1.0 + pow(1.0/tan(locTheta) - locZError/locR, 2))) - locZError;
	return locPathError/SPEED_OF_LIGHT;
}

} //end DAnalysis namespace

#endif // DPhotonComboer_h
