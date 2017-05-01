#include "ANALYSIS/DVertexCreator.h"

namespace DAnalysis
{

DVertexCreator::DVertexCreator(JEventLoop* locEventLoop)
{

	//GET TRACKING HYPOTHESES
	vector<int> locHypotheses = {Positron, PiPlus, KPlus, Proton, Electron, PiMinus, KMinus, Antiproton};
	ostringstream locMassStream;
	for(size_t loc_i = 0; loc_i < locHypotheses.size(); ++loc_i)
	{
		locMassStream << locHypotheses[loc_i];
		if(loc_i != (locHypotheses.size() - 1))
			locMassStream << ",";
	}
	string HYPOTHESES = locMassStream.str();
	gPARMS->SetDefaultParameter("TRKFIT:HYPOTHESES", HYPOTHESES);

	// Parse MASS_HYPOTHESES strings to make list of Particle_t's
	locHypotheses.clear();
	SplitString(HYPOTHESES, locHypotheses, ",");
	for(size_t loc_i = 0; loc_i < locHypotheses.size(); ++loc_i)
		dTrackingPIDs.push_back(Particle_t(locHypotheses[loc_i]));
	std::sort(dTrackingPIDs.begin(), dTrackingPIDs.end()); //so that can search later

}

void DVertexCreator::Do_All(JEventLoop* locEventLoop, const vector<const DReaction*>& locReactions)
{

	//Make initial delta-t cuts to select RF bunches for every neutral shower
	//Still need to place tighter cuts, calculate chisq on vert-by-vert basis


}

} //end DAnalysis namespace

