// $Id$
//
//    File: DChargedTrackHypothesis_factory_KinFit.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DChargedTrackHypothesis_factory_KinFit.h"

//------------------
// init
//------------------
jerror_t DChargedTrackHypothesis_factory_KinFit::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DChargedTrackHypothesis_factory_KinFit::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DChargedTrackHypothesis_factory_KinFit::evnt(jana::JEventLoop* locEventLoop, int eventnumber)
{
 	vector<const DKinFitResults*> locKinFitResultsVector;
	locEventLoop->Get(locKinFitResultsVector);

	const DParticleCombo* locParticleCombo;
	const DParticleComboStep* locParticleComboStep;
	const DChargedTrack* locChargedTrack;
	const DChargedTrackHypothesis* locChargedTrackHypothesis;
	DChargedTrackHypothesis* locNewChargedTrackHypothesis;
	const DKinFitParticle* locKinFitParticle;

	map<const DKinFitParticle*, DChargedTrackHypothesis*> locKinFitParticleMap;
	map<DChargedTrackHypothesis*, deque<const DParticleCombo*> > locChargedParticleComboMap;

	deque<deque<const DKinFitParticle*> > locFinalKinFitParticles;
	for(size_t loc_i = 0; loc_i < locKinFitResultsVector.size(); ++loc_i)
	{
		locParticleCombo = locKinFitResultsVector[loc_i]->Get_ParticleCombo();
		locKinFitResultsVector[loc_i]->Get_FinalKinFitParticles(locFinalKinFitParticles);

		//loop over previous kinfitresults objects, see if the kinfit results are identical: if so, don't need to correct new tracks for the results
		bool locMatchFlag = false;
		for(size_t loc_j = 0; loc_j < loc_i; ++loc_j)
		{
			if(!locParticleCombo->Will_KinFitBeIdentical(locKinFitResultsVector[loc_j]->Get_ParticleCombo()))
				continue;
			//kinfit results are identical: setup the maps so that the charged tracks are copied instead of created anew
			for(size_t loc_k = 0; loc_k < locFinalKinFitParticles.size(); ++loc_k)
			{
				for(size_t loc_l = 0; loc_l < locFinalKinFitParticles[loc_k].size(); ++loc_l)
				{
					locKinFitParticle = locFinalKinFitParticles[loc_k][loc_l];
					if(locKinFitParticle == NULL)
						continue; //e.g. a decaying resonance particle not involved in the kinfit
					if((locKinFitParticle->Get_KinFitParticleType() != d_DetectedParticle) || (locKinFitParticle->Get_Charge() == 0))
						continue;
					locNewChargedTrackHypothesis = locKinFitParticleMap[locKinFitParticle];
					locChargedParticleComboMap[locNewChargedTrackHypothesis].push_back(locParticleCombo);
				}
			}
			locMatchFlag = true;
			break;
		}
		if(locMatchFlag)
			continue;

		for(size_t loc_j = 0; loc_j < locParticleCombo->Get_NumParticleComboSteps(); ++loc_j)
		{
			locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_j);
			for(size_t loc_k = 0; loc_k < locParticleComboStep->Get_NumFinalParticles(); ++loc_k)
			{
				if(!locParticleComboStep->Is_FinalParticleDetected(loc_k))
					continue;
				if(!locParticleComboStep->Is_FinalParticleCharged(loc_k))
					continue;
				locKinFitParticle = locFinalKinFitParticles[loc_j][loc_k];
				locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_k));
				locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticleComboStep->Get_FinalParticle(loc_k));
				locNewChargedTrackHypothesis = Build_ChargedTrackHypothesis(locChargedTrackHypothesis, locKinFitParticle, locChargedTrack, locParticleCombo);
				locKinFitParticleMap[locKinFitParticle] = locNewChargedTrackHypothesis;
				locChargedParticleComboMap[locNewChargedTrackHypothesis] = deque<const DParticleCombo*>(1, locParticleCombo);
			}
		}
	}

	//now set the particle combos as associated objects of the charged tracks, and save the tracks //this marks which combos they originated from
	map<DChargedTrackHypothesis*, deque<const DParticleCombo*> >::iterator locIterator;
	for(locIterator = locChargedParticleComboMap.begin(); locIterator != locChargedParticleComboMap.end(); ++locIterator)
	{
		locNewChargedTrackHypothesis = locIterator->first;
		deque<const DParticleCombo*>& locParticleCombos = locIterator->second;
		for(size_t loc_i = 0; loc_i < locParticleCombos.size(); ++loc_i)
			locNewChargedTrackHypothesis->AddAssociatedObject(locParticleCombos[loc_i]);
		_data.push_back(locNewChargedTrackHypothesis);
	}

	return NOERROR;
}

DChargedTrackHypothesis* DChargedTrackHypothesis_factory_KinFit::Build_ChargedTrackHypothesis(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DKinFitParticle* locKinFitParticle, const DChargedTrack* locChargedTrack, const DParticleCombo* locParticleCombo)
{
	DChargedTrackHypothesis* locNewChargedTrackHypothesis = new DChargedTrackHypothesis(*locChargedTrackHypothesis);
	locNewChargedTrackHypothesis->AddAssociatedObject(locChargedTrackHypothesis);
	locNewChargedTrackHypothesis->AddAssociatedObject(locChargedTrack);

 	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locChargedTrackHypothesis->GetT(locTrackTimeBasedVector);
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
		locNewChargedTrackHypothesis->AddAssociatedObject(locTrackTimeBasedVector[loc_i]);

	vector<const DTOFPoint*> locTOFPoints;
	locChargedTrackHypothesis->GetT(locTOFPoints);
	for(size_t loc_i = 0; loc_i < locTOFPoints.size(); ++loc_i)
		locNewChargedTrackHypothesis->AddAssociatedObject(locTOFPoints[loc_i]);

	vector<const DBCALShower*> locBCALShowers;
	locChargedTrackHypothesis->GetT(locBCALShowers);
	for(size_t loc_i = 0; loc_i < locBCALShowers.size(); ++loc_i)
		locNewChargedTrackHypothesis->AddAssociatedObject(locBCALShowers[loc_i]);

	vector<const DFCALShower*> locFCALShowers;
	locChargedTrackHypothesis->GetT(locFCALShowers);
	for(size_t loc_i = 0; loc_i < locFCALShowers.size(); ++loc_i)
		locNewChargedTrackHypothesis->AddAssociatedObject(locFCALShowers[loc_i]);

	vector<const DSCHit*> locSCHits;
	locChargedTrackHypothesis->GetT(locSCHits);
	for(size_t loc_i = 0; loc_i < locSCHits.size(); ++loc_i)
	  locNewChargedTrackHypothesis->AddAssociatedObject(locSCHits[loc_i]);

	locNewChargedTrackHypothesis->setMomentum(DVector3(locKinFitParticle->Get_Momentum().X(),locKinFitParticle->Get_Momentum().Y(),locKinFitParticle->Get_Momentum().Z()));
	locNewChargedTrackHypothesis->setPosition(DVector3(locKinFitParticle->Get_Position().X(),locKinFitParticle->Get_Position().Y(),locKinFitParticle->Get_Position().Z()));
	locNewChargedTrackHypothesis->setTime(locKinFitParticle->Get_Time());
	locNewChargedTrackHypothesis->setErrorMatrix(*locKinFitParticle->Get_CovarianceMatrix());

	double locPathLength = locNewChargedTrackHypothesis->pathLength() - locKinFitParticle->Get_PathLength();
	double locPathLengthUncertainty_Orig = locNewChargedTrackHypothesis->pathLength_err();
	double locPathLengthUncertainty_KinFit = locKinFitParticle->Get_PathLengthUncertainty();
	double locPathLengthUncertainty = sqrt(locPathLengthUncertainty_Orig*locPathLengthUncertainty_Orig + locPathLengthUncertainty_KinFit*locPathLengthUncertainty_KinFit);
	locNewChargedTrackHypothesis->setPathLength(locPathLength, locPathLengthUncertainty);

	//don't recompute EITHER dedx chisq OR timing chisq: after kinfit timing auto lined up! (even to RF) only orig info and kinfit FOM matters
	return locNewChargedTrackHypothesis;
}

//------------------
// erun
//------------------
jerror_t DChargedTrackHypothesis_factory_KinFit::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DChargedTrackHypothesis_factory_KinFit::fini(void)
{
	return NOERROR;
}


