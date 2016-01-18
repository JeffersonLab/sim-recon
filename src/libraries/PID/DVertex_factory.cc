// $Id$
//
//    File: DVertex_factory.cc
// Created: Tue Apr  6 17:01:54 EDT 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#include "DVertex_factory.h"

//------------------
// init
//------------------
jerror_t DVertex_factory::init(void)
{
	dKinFitDebugLevel = 0;
	dMinTrackingFOM = 5.73303E-7;
	dNoKinematicFitFlag = false;
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DVertex_factory::brun(jana::JEventLoop* locEventLoop, int32_t runnumber)
{
	dKinFitUtils = new DKinFitUtils_GlueX(locEventLoop);
	dKinFitter = new DKinFitter(dKinFitUtils);

	// Get Target parameters from XML
	dTargetZCenter = 65.0;
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());
	locGeometry->GetTargetZ(dTargetZCenter);

	gPARMS->SetDefaultParameter("VERTEX:NO_KINFIT_FLAG", dNoKinematicFitFlag);
	gPARMS->SetDefaultParameter("KINFIT:DEBUGLEVEL", dKinFitDebugLevel);

	dKinFitter->Set_DebugLevel(dKinFitDebugLevel);

	locEventLoop->GetSingle(dAnalysisUtilities);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DVertex_factory::evnt(JEventLoop* locEventLoop, uint64_t eventnumber)
{
	//preferentially (kinematic fit):
		//use tracks with a matched hit & good tracking FOM
		//if no good tracks (or none with matched hits), use all tracks
		//if no tracks, use target center

	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	//select the best DTrackTimeBased for each track: use best tracking FOM
	map<JObject::oid_t, const DTrackTimeBased*> locBestTrackTimeBasedMap; //lowest tracking chisq/ndf for each candidate id
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		JObject::oid_t locCandidateID = locTrackTimeBasedVector[loc_i]->candidateid;
		if(locBestTrackTimeBasedMap.find(locCandidateID) == locBestTrackTimeBasedMap.end())
			locBestTrackTimeBasedMap[locCandidateID] = locTrackTimeBasedVector[loc_i];
		else if(locTrackTimeBasedVector[loc_i]->FOM > locBestTrackTimeBasedMap[locCandidateID]->FOM)
			locBestTrackTimeBasedMap[locCandidateID] = locTrackTimeBasedVector[loc_i];
	}

	//separate the tracks based on high/low tracking FOM & has hit-match
	map<JObject::oid_t, const DTrackTimeBased*>::iterator locIterator;
	deque<const DTrackTimeBased*> locTrackTimeBasedVector_OnePerTrack, locTrackTimeBasedVector_OnePerTrack_Good;
	for(locIterator = locBestTrackTimeBasedMap.begin(); locIterator != locBestTrackTimeBasedMap.end(); ++locIterator)
	{
		const DTrackTimeBased* locTrackTimeBased = locIterator->second;
		locTrackTimeBasedVector_OnePerTrack.push_back(locTrackTimeBased);
		if((locTrackTimeBased->FOM >= dMinTrackingFOM) && locDetectorMatches->Get_IsMatchedToHit(locTrackTimeBased))
			locTrackTimeBasedVector_OnePerTrack_Good.push_back(locTrackTimeBased);
	}

	deque<const DTrackTimeBased*> locTrackTimeBasedVectorToUse = (locTrackTimeBasedVector_OnePerTrack_Good.size() >= 2) ? locTrackTimeBasedVector_OnePerTrack_Good : locTrackTimeBasedVector_OnePerTrack;

	//handle cases of no/one track
	if(locTrackTimeBasedVectorToUse.empty())
		return Create_Vertex(DVector3(0.0, 0.0, dTargetZCenter), locEventRFBunch->dTime);
	if(locTrackTimeBasedVectorToUse.size() == 1)
		return Create_Vertex(locTrackTimeBasedVectorToUse[0]->position(), locEventRFBunch->dTime);

	// first calculate a rough vertex
	DVector3 locRoughPosition = dAnalysisUtilities->Calc_CrudeVertex(locTrackTimeBasedVectorToUse);

	// if only want rough guess, save it and exit
	if(dNoKinematicFitFlag)
		return Create_Vertex(locRoughPosition, locEventRFBunch->dTime);

	//prepare for kinematic fit
	dKinFitter->Reset_NewEvent();
	TVector3 locTRoughPosition(locRoughPosition.X(), locRoughPosition.Y(), locRoughPosition.Z());

	// create particles for kinematic fit
	set<DKinFitParticle*> locKinFitParticles;
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVectorToUse.size(); ++loc_i)
		locKinFitParticles.insert(dKinFitUtils->Make_DetectedParticle(locTrackTimeBasedVectorToUse[loc_i]));

	// create vertex constraint
	set<DKinFitParticle*> locNoConstrainParticles;
	DKinFitConstraint_Vertex* locVertexConstraint = dKinFitUtils->Make_VertexConstraint(locKinFitParticles, locNoConstrainParticles, locTRoughPosition);

	dKinFitter->Add_Constraint(locVertexConstraint);

	if(!dKinFitter->Fit_Reaction()) //if fit fails to converge: use rough results
		return Create_Vertex(locRoughPosition, locEventRFBunch->dTime);

	DKinFitConstraint_Vertex* locResultVertexConstraint = dynamic_cast<DKinFitConstraint_Vertex*>(*dKinFitter->Get_KinFitConstraints().begin());

	//save kinfit results
	TVector3 locFitVertex = locResultVertexConstraint->Get_CommonVertex();
	DVector3 locTFitVertex(locFitVertex.X(), locFitVertex.Y(), locFitVertex.Z());
	unsigned int locKinFitNDF = dKinFitter->Get_NDF();
	double locKinFitChiSq = dKinFitter->Get_ChiSq();
	Create_Vertex(locTFitVertex, locEventRFBunch->dTime, locKinFitNDF, locKinFitChiSq);
	dKinFitter->Get_Pulls(_data.back()->dKinFitPulls);

	return NOERROR;
}

jerror_t DVertex_factory::Create_Vertex(DVector3 locPosition, double locRFTime, unsigned int locKinFitNDF, double locKinFitChiSq)
{
	double locTime = locRFTime + (locPosition.Z() - dTargetZCenter)/29.9792458;

	DVertex* locVertex = new DVertex();
	locVertex->dSpacetimeVertex = DLorentzVector(locPosition, locTime);
	locVertex->dKinFitNDF = locKinFitNDF;
	locVertex->dKinFitChiSq = locKinFitChiSq;

	_data.push_back(locVertex);
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DVertex_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DVertex_factory::fini(void)
{
	return NOERROR;
}

