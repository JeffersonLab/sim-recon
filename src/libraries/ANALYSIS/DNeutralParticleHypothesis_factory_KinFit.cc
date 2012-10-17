// $Id$
//
//    File: DNeutralParticleHypothesis_factory_KinFit.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DNeutralParticleHypothesis_factory_KinFit.h"

//------------------
// init
//------------------
jerror_t DNeutralParticleHypothesis_factory_KinFit::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DNeutralParticleHypothesis_factory_KinFit::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DNeutralParticleHypothesis_factory_KinFit::evnt(jana::JEventLoop* locEventLoop, int eventnumber)
{
 	vector<const DKinFitResults*> locKinFitResultsVector;
	locEventLoop->Get(locKinFitResultsVector);

	const DParticleCombo* locParticleCombo;
	const DParticleComboStep* locParticleComboStep;
	const DNeutralShower* locNeutralShower;
	const DNeutralParticleHypothesis* locNeutralParticleHypothesis;
	DNeutralParticleHypothesis* locNewNeutralParticleHypothesis;
	const DKinFitParticle* locKinFitParticle;

	map<const DKinFitParticle*, DNeutralParticleHypothesis*> locKinFitParticleMap;
	map<DNeutralParticleHypothesis*, deque<const DParticleCombo*> > locNeutralParticleComboMap;

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
					if((locKinFitParticle->Get_KinFitParticleType() != d_DetectedParticle) || (locKinFitParticle->Get_Charge() != 0))
						continue;
					locNewNeutralParticleHypothesis = locKinFitParticleMap[locKinFitParticle];
					locNeutralParticleComboMap[locNewNeutralParticleHypothesis].push_back(locParticleCombo);
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
				if(!locParticleComboStep->Is_FinalParticleNeutral(loc_k))
					continue;
				locKinFitParticle = locFinalKinFitParticles[loc_j][loc_k];
				locNeutralShower = static_cast<const DNeutralShower*>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_k));
				locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticleComboStep->Get_FinalParticle(loc_k));
				locNewNeutralParticleHypothesis = Build_NeutralParticleHypothesis(locNeutralParticleHypothesis, locKinFitParticle, locNeutralShower, locParticleCombo);
				locKinFitParticleMap[locKinFitParticle] = locNewNeutralParticleHypothesis;
				locNeutralParticleComboMap[locNewNeutralParticleHypothesis] = deque<const DParticleCombo*>(1, locParticleCombo);
			}
		}
	}

	//now set the particle combos as associated objects of the neutral tracks, and save the tracks //this marks which combos they originated from
	map<DNeutralParticleHypothesis*, deque<const DParticleCombo*> >::iterator locIterator;
	for(locIterator = locNeutralParticleComboMap.begin(); locIterator != locNeutralParticleComboMap.end(); ++locIterator)
	{
		locNewNeutralParticleHypothesis = locIterator->first;
		deque<const DParticleCombo*>& locParticleCombos = locIterator->second;
		for(size_t loc_i = 0; loc_i < locParticleCombos.size(); ++loc_i)
			locNewNeutralParticleHypothesis->AddAssociatedObject(locParticleCombos[loc_i]);
		_data.push_back(locNewNeutralParticleHypothesis);
	}

	return NOERROR;
}

DNeutralParticleHypothesis* DNeutralParticleHypothesis_factory_KinFit::Build_NeutralParticleHypothesis(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, const DKinFitParticle* locKinFitParticle, const DNeutralShower* locNeutralShower, const DParticleCombo* locParticleCombo)
{
	DNeutralParticleHypothesis* locNewNeutralParticleHypothesis = new DNeutralParticleHypothesis(*locNeutralParticleHypothesis);
	locNewNeutralParticleHypothesis->AddAssociatedObject(locNeutralParticleHypothesis);
	locNewNeutralParticleHypothesis->AddAssociatedObject(locNeutralShower);

	vector<const DBCALShower*> locBCALShowers;
	locNeutralParticleHypothesis->GetT(locBCALShowers);
	for(size_t loc_i = 0; loc_i < locBCALShowers.size(); ++loc_i)
		locNewNeutralParticleHypothesis->AddAssociatedObject(locBCALShowers[loc_i]);

	vector<const DFCALShower*> locFCALShowers;
	locNeutralParticleHypothesis->GetT(locFCALShowers);
	for(size_t loc_i = 0; loc_i < locFCALShowers.size(); ++loc_i)
		locNewNeutralParticleHypothesis->AddAssociatedObject(locFCALShowers[loc_i]);

	locNewNeutralParticleHypothesis->setMomentum(DVector3(locKinFitParticle->Get_Momentum().X(),locKinFitParticle->Get_Momentum().Y(),locKinFitParticle->Get_Momentum().Z()));
	locNewNeutralParticleHypothesis->setPosition(DVector3(locKinFitParticle->Get_CommonVertex().X(),locKinFitParticle->Get_CommonVertex().Y(),locKinFitParticle->Get_CommonVertex().Z()));
	locNewNeutralParticleHypothesis->setTime(locKinFitParticle->Get_CommonTime());
	locNewNeutralParticleHypothesis->setErrorMatrix(*locKinFitParticle->Get_CovarianceMatrix());

	double locPathLength = locNewNeutralParticleHypothesis->pathLength() - locKinFitParticle->Get_PathLength();
	double locPathLengthUncertainty_Orig = locNewNeutralParticleHypothesis->pathLength_err();
	double locPathLengthUncertainty_KinFit = locKinFitParticle->Get_PathLengthUncertainty();
	double locPathLengthUncertainty = sqrt(locPathLengthUncertainty_Orig*locPathLengthUncertainty_Orig + locPathLengthUncertainty_KinFit*locPathLengthUncertainty_KinFit);
	locNewNeutralParticleHypothesis->setPathLength(locPathLength, locPathLengthUncertainty);

	//don't recompute EITHER dedx chisq OR timing chisq: after kinfit timing auto lined up! (even to RF) only orig info and kinfit FOM matters
	return locNewNeutralParticleHypothesis;
}

//------------------
// erun
//------------------
jerror_t DNeutralParticleHypothesis_factory_KinFit::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DNeutralParticleHypothesis_factory_KinFit::fini(void)
{
	return NOERROR;
}


