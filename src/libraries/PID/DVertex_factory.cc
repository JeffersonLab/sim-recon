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
	dTargetLength = 30.0;
	dTargetRadius = 1.5; //FIX: grab from database!!!
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());
	locGeometry->GetTargetZ(dTargetZCenter);
	locGeometry->GetTargetLength(dTargetLength);

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
		return Create_Vertex_NoTracks(locEventRFBunch);
	if(locTrackTimeBasedVectorToUse.size() == 1)
		return Create_Vertex_OneTrack(locTrackTimeBasedVectorToUse[0], locEventRFBunch);

	// first calculate a rough vertex
	DVector3 locRoughPosition = dAnalysisUtilities->Calc_CrudeVertex(locTrackTimeBasedVectorToUse);

	// if only want rough guess, save it and exit
	if(dNoKinematicFitFlag)
		return Create_Vertex_Rough(locRoughPosition, locEventRFBunch);

	//prepare for kinematic fit
	dKinFitUtils->Reset_NewEvent(locEventLoop->GetJEvent().GetEventNumber());
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

	// fit, save, and return
	if(!dKinFitter->Fit_Reaction()) //if fit fails to converge: use rough results
		return Create_Vertex_Rough(locRoughPosition, locEventRFBunch);

	//save kinfit results
	return Create_Vertex_KinFit(locEventRFBunch);
}

jerror_t DVertex_factory::Create_Vertex_NoTracks(const DEventRFBunch* locEventRFBunch)
{
	DVertex* locVertex = new DVertex();
	locVertex->dSpacetimeVertex = DLorentzVector(DVector3(0.0, 0.0, dTargetZCenter), locEventRFBunch->dTime);
	locVertex->dKinFitNDF = 0;
	locVertex->dKinFitChiSq = 0.0;

	//error matrix
	locVertex->dCovarianceMatrix.ResizeTo(4, 4);
	locVertex->dCovarianceMatrix.Zero();
	locVertex->dCovarianceMatrix(0, 0) = dTargetRadius*dTargetRadius/12.0; //x variance //should instead use beam spot size
	locVertex->dCovarianceMatrix(1, 1) = dTargetRadius*dTargetRadius/12.0; //y variance //should instead use beam spot size
	locVertex->dCovarianceMatrix(2, 2) = dTargetLength*dTargetLength/12.0; //z variance
	locVertex->dCovarianceMatrix(3, 3) = locEventRFBunch->dTimeVariance; //t variance

	_data.push_back(locVertex);
	return NOERROR;
}

jerror_t DVertex_factory::Create_Vertex_OneTrack(const DTrackTimeBased* locTrackTimeBased, const DEventRFBunch* locEventRFBunch)
{
	DVector3 locPosition = locTrackTimeBased->position();
	double locTime = locEventRFBunch->dTime + (locPosition.Z() - dTargetZCenter)/29.9792458;

	DVertex* locVertex = new DVertex();
	locVertex->dSpacetimeVertex = DLorentzVector(locPosition, locTime);
	locVertex->dKinFitNDF = 0;
	locVertex->dKinFitChiSq = 0.0;

	//error matrix
    const TMatrixFSym& locTrackErrorMatrix = *(locTrackTimeBased->errorMatrix());
	locVertex->dCovarianceMatrix.ResizeTo(4, 4);
	locVertex->dCovarianceMatrix.Zero();
	for(size_t loc_i = 0; loc_i < 3; ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < 3; ++loc_j)
			locVertex->dCovarianceMatrix(loc_i, loc_j) = locTrackErrorMatrix(loc_i + 3, loc_j + 3);
	}
	locVertex->dCovarianceMatrix(3, 3) = locEventRFBunch->dTimeVariance; //t variance

	_data.push_back(locVertex);
	return NOERROR;
}

jerror_t DVertex_factory::Create_Vertex_Rough(DVector3 locPosition, const DEventRFBunch* locEventRFBunch)
{
	double locTime = locEventRFBunch->dTime + (locPosition.Z() - dTargetZCenter)/29.9792458;

	DVertex* locVertex = new DVertex();
	locVertex->dSpacetimeVertex = DLorentzVector(locPosition, locTime);
	locVertex->dKinFitNDF = 0;
	locVertex->dKinFitChiSq = 0.0;

	//error matrix //too lazy to compute properly right now ... need to hack DAnalysisUtilities::Calc_DOCA()
	locVertex->dCovarianceMatrix.ResizeTo(4, 4);
	locVertex->dCovarianceMatrix.Zero();
	locVertex->dCovarianceMatrix(0, 0) = 0.65; //x variance //from monitoring plots of vertex
	locVertex->dCovarianceMatrix(1, 1) = 0.65; //y variance //from monitoring plots of vertex
	locVertex->dCovarianceMatrix(2, 2) = 1.5; //z variance //a guess, semi-guarding against the worst case scenario //ugh
	locVertex->dCovarianceMatrix(3, 3) = locEventRFBunch->dTimeVariance; //t variance

	_data.push_back(locVertex);
	return NOERROR;
}

jerror_t DVertex_factory::Create_Vertex_KinFit(const DEventRFBunch* locEventRFBunch)
{
	DKinFitConstraint_Vertex* locResultVertexConstraint = dynamic_cast<DKinFitConstraint_Vertex*>(*dKinFitter->Get_KinFitConstraints().begin());

	TVector3 locFitVertex = locResultVertexConstraint->Get_CommonVertex();
	DVector3 locDFitVertex(locFitVertex.X(), locFitVertex.Y(), locFitVertex.Z());
	double locTime = locEventRFBunch->dTime + (locDFitVertex.Z() - dTargetZCenter)/29.9792458;

	DVertex* locVertex = new DVertex();
	locVertex->dSpacetimeVertex = DLorentzVector(locDFitVertex, locTime);
	locVertex->dKinFitNDF = dKinFitter->Get_NDF();
	locVertex->dKinFitChiSq = dKinFitter->Get_ChiSq();

	//error matrix
	const TMatrixDSym& locMatrixDSym = dKinFitter->Get_VXi();
	locVertex->dCovarianceMatrix.ResizeTo(4, 4);
	locVertex->dCovarianceMatrix.Zero();
	for(size_t loc_i = 0; loc_i < 3; ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < 3; ++loc_j)
			locVertex->dCovarianceMatrix(loc_i, loc_j) = locMatrixDSym(loc_i, loc_j);
	}
	locVertex->dCovarianceMatrix(3, 3) = locEventRFBunch->dTimeVariance; //t variance

	//Particle Maps & Pulls
	//Build pulls from this:
	map<DKinFitParticle*, map<DKinFitPullType, double> > locPulls_KinFitParticle;
	dKinFitter->Get_Pulls(locPulls_KinFitParticle);

	//By looping over the pulls:
	map<DKinFitParticle*, map<DKinFitPullType, double> >::iterator locMapIterator = locPulls_KinFitParticle.begin();
	for(; locMapIterator != locPulls_KinFitParticle.end(); ++locMapIterator)
	{
		DKinFitParticle* locOutputKinFitParticle = locMapIterator->first;
		DKinFitParticle* locInputKinFitParticle = dKinFitUtils->Get_InputKinFitParticle(locOutputKinFitParticle);

		const JObject* locSourceJObject = dKinFitUtils->Get_SourceJObject(locInputKinFitParticle);
		locVertex->dKinFitPulls[locSourceJObject] = locMapIterator->second;
	}

	_data.push_back(locVertex);

	return NOERROR;
}
