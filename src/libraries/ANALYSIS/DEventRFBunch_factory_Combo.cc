// $Id$
//
//    File: DEventRFBunch_factory_Combo.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt
//

#include "DEventRFBunch_factory_Combo.h"

//------------------
// init
//------------------
jerror_t DEventRFBunch_factory_Combo::init(void)
{
	dRFBunchFrequency = 2.004;
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventRFBunch_factory_Combo::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	// Get Target parameters from XML
	DApplication *locApplication = dynamic_cast<DApplication*> (locEventLoop->GetJApplication());
	DGeometry *locGeometry = locApplication ? locApplication->GetDGeometry(runnumber):NULL;
	dTargetCenterZ = 65.0;
	dTargetLength = 30.0;
	dTargetRadius = 1.5; //FIX: grab from database!!!
	if(locGeometry)
	{
		locGeometry->GetTargetZ(dTargetCenterZ);
		locGeometry->GetTargetLength(dTargetLength);
	}

	locEventLoop->GetSingle(dParticleID);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventRFBunch_factory_Combo::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
 	vector<const DParticleComboBlueprint*> locParticleComboBlueprints;
	locEventLoop->Get(locParticleComboBlueprints);

	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches);

/*
//the below disables this routine until everything is working
if (locEventRFBunches.size() > 0) {
   DEventRFBunch* locEventRFBunch = new DEventRFBunch(*locEventRFBunches[0]);
   for(size_t loc_i = 0; loc_i < locParticleComboBlueprints.size(); ++loc_i)
	locEventRFBunch->AddAssociatedObject(locParticleComboBlueprints[loc_i]);
   _data.push_back(locEventRFBunch);
}
return NOERROR;
*/

	//use all particles in the combo to find the vertex
		//ignore detached vertices

 	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector, "Combo");

	map<int, DEventRFBunch*> locComboRFBunchMap;

	for(size_t loc_i = 0; loc_i < locParticleComboBlueprints.size(); ++loc_i)
	{
		const DParticleComboBlueprint* locParticleComboBlueprint = locParticleComboBlueprints[loc_i];
		double locWeightedAverageStartTime_Numerator = 0.0, locWeightedAverageStartTime_Denominator = 0.0;
		double locWeightedAverageVertexZ_Numerator = 0.0, locWeightedAverageVertexZ_Denominator = 0.0;
		for(size_t loc_j = 0; loc_j < locParticleComboBlueprint->Get_NumParticleComboBlueprintSteps(); ++loc_j)
		{
			const DParticleComboBlueprintStep* locParticleComboBlueprintStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_j);
			for(size_t loc_k = 0; loc_k < locParticleComboBlueprintStep->Get_NumFinalParticleSourceObjects(); ++loc_k)
			{
				const JObject* locSourceObject = locParticleComboBlueprintStep->Get_FinalParticle_SourceObject(loc_k);
				if(locSourceObject == NULL)
					continue; //missing or decaying

				//charged tracks
				const DChargedTrack* locChargedTrack = dynamic_cast<const DChargedTrack*>(locSourceObject);
				if(locChargedTrack == NULL)
					continue; //will do neutrals later

				Particle_t locPID = locParticleComboBlueprintStep->Get_FinalParticleID(loc_k);
				const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(locPID);
				if(locChargedTrackHypothesis != NULL)
				{
					double locStartTime = locChargedTrackHypothesis->time();
					double locStartTimeVariance = (locChargedTrackHypothesis->errorMatrix())(6, 6);
					if(locStartTimeVariance > 0.0)
					{
						locWeightedAverageStartTime_Numerator += locStartTime/locStartTimeVariance;
						locWeightedAverageStartTime_Denominator += 1.0/locStartTimeVariance;
					}
					double locVertexZ = locChargedTrackHypothesis->position().Z();
					double locVertexZVariance = (locChargedTrackHypothesis->errorMatrix())(5, 5);
					locWeightedAverageVertexZ_Numerator += locVertexZ/locVertexZVariance;
					locWeightedAverageVertexZ_Denominator += 1.0/locVertexZVariance;
					continue;
				}

				//get from locTrackTimeBasedVector
				for(size_t loc_l = 0; loc_l < locTrackTimeBasedVector.size(); ++loc_l)
				{
					const DTrackTimeBased* locTrackTimeBased = locTrackTimeBasedVector[loc_l];
					if(locTrackTimeBased->PID() != locPID)
						continue;
					const DChargedTrack* locTimeBasedSourceChargedTrack = NULL;
					locTrackTimeBased->GetSingleT(locTimeBasedSourceChargedTrack);
					if(locTimeBasedSourceChargedTrack != locChargedTrack)
						continue;

					double locStartTime, locStartTimeVariance;
					Get_StartTime(locEventLoop, locTrackTimeBased, locStartTime, locStartTimeVariance);
					if(locStartTimeVariance > 0.0)
					{
						locWeightedAverageStartTime_Numerator += locStartTime/locStartTimeVariance;
						locWeightedAverageStartTime_Denominator += 1.0/locStartTimeVariance;
					}
					double locVertexZ = locTrackTimeBased->position().Z();
					double locVertexZVariance = (locTrackTimeBased->errorMatrix())(5, 5);
					locWeightedAverageVertexZ_Numerator += locVertexZ/locVertexZVariance;
					locWeightedAverageVertexZ_Denominator += 1.0/locVertexZVariance;
					break;
				}
				continue;
			}
		}

		double locVertexZ = (locWeightedAverageVertexZ_Denominator > 0.0) ? locWeightedAverageVertexZ_Numerator/locWeightedAverageVertexZ_Denominator : dTargetCenterZ;
		DVector3 locVertex(0.0, 0.0, locVertexZ);

		//neutrals
		for(size_t loc_j = 0; loc_j < locParticleComboBlueprint->Get_NumParticleComboBlueprintSteps(); ++loc_j)
		{
			const DParticleComboBlueprintStep* locParticleComboBlueprintStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_j);
			for(size_t loc_k = 0; loc_k < locParticleComboBlueprintStep->Get_NumFinalParticleSourceObjects(); ++loc_k)
			{
				const JObject* locSourceObject = locParticleComboBlueprintStep->Get_FinalParticle_SourceObject(loc_k);
				if(locSourceObject == NULL)
					continue; //missing or decaying

				const DNeutralShower* locNeutralShower = dynamic_cast<const DNeutralShower*>(locSourceObject);
				if(locNeutralShower == NULL)
					continue; //charged, already done
				Particle_t locPID = locParticleComboBlueprintStep->Get_FinalParticleID(loc_k);
				if(locPID != Gamma)
					continue; //other neutrals (e.g. neutron) can't be used to pick the time: their momentum is defined by the time

				double locStartTime, locStartTimeVariance;
				Calc_StartTime(locNeutralShower, locVertex, locStartTime, locStartTimeVariance);
				if(locStartTimeVariance > 0.0)
				{
					locWeightedAverageStartTime_Numerator += locStartTime/locStartTimeVariance;
					locWeightedAverageStartTime_Denominator += 1.0/locStartTimeVariance;
				}
			}
		}

		if(!(locWeightedAverageStartTime_Denominator > 0.0))
		{
			//no timing information somehow: use the pre-existing value
			DEventRFBunch* locEventRFBunch = new DEventRFBunch(*locEventRFBunches[0]);
			locEventRFBunch->dMatchedToTracksFlag = false;
			locEventRFBunch->AddAssociatedObject(locParticleComboBlueprint);
			_data.push_back(locEventRFBunch);
			continue;
		}

		// Find # RF Bunch Shifts
		double locStartTime = locWeightedAverageStartTime_Numerator/locWeightedAverageStartTime_Denominator;
		double locPropagatedRFTime = locEventRFBunches[0]->dTime + (locVertexZ - dTargetCenterZ)/SPEED_OF_LIGHT;

		int locNumBunchShifts = 0;
		while((locPropagatedRFTime - locStartTime) > (0.5*dRFBunchFrequency))
		{
			--locNumBunchShifts;
			locPropagatedRFTime -= dRFBunchFrequency;
		}
		while((locPropagatedRFTime - locStartTime) < (-0.5*dRFBunchFrequency))
		{
			++locNumBunchShifts;
			locPropagatedRFTime += dRFBunchFrequency;
		}

		// Create new RF Bunch if doesn't already exist
		if(locComboRFBunchMap.find(locNumBunchShifts) != locComboRFBunchMap.end()) //already created, don't recreate identical object!
			locComboRFBunchMap[locNumBunchShifts]->AddAssociatedObject(locParticleComboBlueprint);
		else
		{
			DEventRFBunch* locEventRFBunch = new DEventRFBunch();
			locEventRFBunch->dMatchedToTracksFlag = true;
			locEventRFBunch->dTime = locEventRFBunches[0]->dTime + (double)(locNumBunchShifts)*dRFBunchFrequency;
			locEventRFBunch->dTimeVariance = locEventRFBunches[0]->dTimeVariance;
			locEventRFBunch->AddAssociatedObject(locParticleComboBlueprint);
			_data.push_back(locEventRFBunch);
			locComboRFBunchMap[locNumBunchShifts] = locEventRFBunch;
		}
	}

	return NOERROR;
}

void DEventRFBunch_factory_Combo::Get_StartTime(JEventLoop* locEventLoop, const DTrackTimeBased* locTrackTimeBased, double& locStartTime, double& locStartTimeVariance)
{
	vector<const DTOFPoint*> locTOFPoints;
	locEventLoop->Get(locTOFPoints);

	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locFCALShowers);

	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locSCHits);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	// Use time-based tracking time as initial guess
	locStartTime = 0.0;
	locStartTimeVariance = 0.0;

	//BCAL
	DShowerMatchParams locBCALShowerMatchParams;
	if(dParticleID->Get_BestBCALMatchParams(locTrackTimeBased, locDetectorMatches, locBCALShowerMatchParams))
	{
		const DBCALShower* locBCALShower = dynamic_cast<const DBCALShower*>(locBCALShowerMatchParams.dShowerObject);
		locStartTime = locBCALShower->t - locBCALShowerMatchParams.dFlightTime;
//		locStartTimeVariance = sqrt(locShowerMatchParams.dFlightTimeVariance) - locBCALShower->dCovarianceMatrix(4, 4); //uncomment when ready!!
//		locStartTimeVariance *= locStartTimeVariance; //uncomment when ready!!
		locStartTimeVariance = 0.5*0.5;
		return;
	}

	//TOF
	DTOFHitMatchParams locTOFHitMatchParams;
	if(dParticleID->Get_BestTOFMatchParams(locTrackTimeBased, locDetectorMatches, locTOFHitMatchParams))
	{
		const DTOFPoint* locTOFPoint = locTOFHitMatchParams.dTOFPoint;
		locStartTime = locTOFPoint->t - locTOFHitMatchParams.dFlightTime;
//		locStartTimeVariance = sqrt(locTOFHitMatchParams.dFlightTimeVariance) - locTOFPoints[loc_i]->tErr; //uncomment when ready!
//		locStartTimeVariance *= locStartTimeVariance; //uncomment when ready!
		locStartTimeVariance = 0.1*0.1;
		return;
	}

	//FCAL
	DShowerMatchParams locFCALShowerMatchParams;
	if(dParticleID->Get_BestFCALMatchParams(locTrackTimeBased, locDetectorMatches, locFCALShowerMatchParams))
	{
		const DFCALShower* locFCALShower = dynamic_cast<const DFCALShower*>(locFCALShowerMatchParams.dShowerObject);
		locStartTime = locFCALShower->getTime() - locFCALShowerMatchParams.dFlightTime;
//		locStartTimeVariance = sqrt(locShowerMatchParams.dFlightTimeVariance) - sqrt(locFCALShowers[loc_i]->dCovarianceMatrix(4, 4)); //uncomment when ready!
//		locStartTimeVariance *= locStartTimeVariance; //uncomment when ready!
		locStartTimeVariance = 0.5*0.5;
		return;
	}
}

void DEventRFBunch_factory_Combo::Calc_StartTime(const DNeutralShower* locNeutralShower, DVector3 locVertex, double& locStartTime, double& locStartTimeVariance)
{
	locStartTime = 0.0;
	locStartTimeVariance = 0.0;

	//doesn't work for neutrons!!
	double locHitTime = locNeutralShower->dSpacetimeVertex.T();
	DVector3 locHitPoint = locNeutralShower->dSpacetimeVertex.Vect();

	// Calculate DNeutralParticleHypothesis Quantities (projected time at vertex for given id, etc.)
	DVector3 locPath = locHitPoint - locVertex;
	double locPathLength = locPath.Mag();
	if(!(locPathLength > 0.0))
		return;

	double locFlightTime = locPathLength/29.9792458;
	locStartTime = locHitTime - locFlightTime;

	locStartTimeVariance = Calc_StartTimeVariance(locNeutralShower, locPath);
}

double DEventRFBunch_factory_Combo::Calc_StartTimeVariance(const DNeutralShower* locNeutralShower, const DVector3& locPathVector)
{
	//build 8x8 matrix: 5x5 shower, 3x3 vertex position
	DMatrixDSym locShowerPlusVertCovariance(8);
	for(unsigned int loc_l = 0; loc_l < 5; ++loc_l) //shower: e, x, y, z, t
	{
		for(unsigned int loc_m = 0; loc_m < 5; ++loc_m)
			locShowerPlusVertCovariance(loc_l, loc_m) = locNeutralShower->dCovarianceMatrix(loc_l, loc_m);
	}
	locShowerPlusVertCovariance(5, 5) = 0.25*dTargetRadius*dTargetRadius/12.0; //vertex position x
	locShowerPlusVertCovariance(6, 6) = 0.25*dTargetRadius*dTargetRadius/12.0; //vertex position y
	locShowerPlusVertCovariance(7, 7) = dTargetLength*dTargetLength/12.0; //vertex position z

	DVector3 locDeltaX = -1.0*locPathVector; //defined oppositely in document! //delta_x defined here as common_vertex - hit_point
	DVector3 locUnitDeltaXOverC = (1.0/(SPEED_OF_LIGHT*locDeltaX.Mag()))*locDeltaX;

	//build transform matrix
	DMatrix locTransformMatrix(1, 8);

	locTransformMatrix(0, 0) = 0.0; //partial deriv of t wrst shower-e //beta defined = 1
	locTransformMatrix(0, 1) = locUnitDeltaXOverC.X(); //partial deriv of t wrst shower-x
	locTransformMatrix(0, 2) = locUnitDeltaXOverC.Y(); //partial deriv of t wrst shower-y
	locTransformMatrix(0, 3) = locUnitDeltaXOverC.Z(); //partial deriv of t wrst shower-z
	locTransformMatrix(0, 4) = 1.0; //partial deriv of t wrst shower-t
	locTransformMatrix(0, 5) = -1.0*locTransformMatrix(0, 1); //partial deriv of t wrst vert-x
	locTransformMatrix(0, 6) = -1.0*locTransformMatrix(0, 2); //partial deriv of t wrst vert-y
	locTransformMatrix(0, 7) = -1.0*locTransformMatrix(0, 3); //partial deriv of t wrst vert-z

	//convert
	DMatrixDSym locParticleCovariance(1);
	locParticleCovariance = locShowerPlusVertCovariance.Similarity(locTransformMatrix);
	return locParticleCovariance(0, 0);
}

//------------------
// erun
//------------------
jerror_t DEventRFBunch_factory_Combo::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventRFBunch_factory_Combo::fini(void)
{
	return NOERROR;
}


