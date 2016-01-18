// $Id$
//
//    File: DChargedTrackHypothesis_factory_KinFit.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifdef VTRACE
#include "vt_user.h"
#endif

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
jerror_t DChargedTrackHypothesis_factory_KinFit::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	locEventLoop->GetSingle(dPIDAlgorithm);
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DChargedTrackHypothesis_factory_KinFit::evnt(jana::JEventLoop* locEventLoop, uint64_t eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DChargedTrackHypothesis_factory_KinFit::evnt()");
#endif

 	vector<const DKinFitResults*> locKinFitResultsVector;
	locEventLoop->Get(locKinFitResultsVector);

	map<DKinFitParticle*, DChargedTrackHypothesis*> locNewObjectMap;

	for(size_t loc_i = 0; loc_i < locKinFitResultsVector.size(); ++loc_i)
	{
		map<const DParticleCombo*, const DKinFitChain*> locParticleComboMap;
		locKinFitResultsVector[loc_i]->Get_ParticleComboMap(locParticleComboMap);

		map<const DParticleCombo*, const DKinFitChain*>::iterator locComboIterator = locParticleComboMap.begin();
		for(; locComboIterator != locParticleComboMap.end(); ++locComboIterator)
		{
			const DParticleCombo* locParticleCombo = locComboIterator->first;
			for(size_t loc_j = 0; loc_j < locParticleCombo->Get_NumParticleComboSteps(); ++loc_j)
			{
				const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_j);
				for(size_t loc_k = 0; loc_k < locParticleComboStep->Get_NumFinalParticles(); ++loc_k)
				{
					if(!locParticleComboStep->Is_FinalParticleDetected(loc_k))
						continue;
					if(!locParticleComboStep->Is_FinalParticleCharged(loc_k))
						continue;

					const DKinematicData* locMeasuredParticle = locParticleComboStep->Get_FinalParticle_Measured(loc_k);
					DKinFitParticle* locKinFitParticle = locKinFitResultsVector[loc_i]->Get_OutputKinFitParticle(locMeasuredParticle);
					if(locKinFitParticle == NULL)
						continue; //not used in fit

					map<DKinFitParticle*, DChargedTrackHypothesis*>::iterator locNewHypoIterator = locNewObjectMap.find(locKinFitParticle);
					if(locNewHypoIterator != locNewObjectMap.end())
					{
						locNewHypoIterator->second->AddAssociatedObject(locParticleCombo);
						continue; //new particle already created for this kinfit particle
					}

					const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_k));
					const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locParticleComboStep->Get_FinalParticle(loc_k));

					DChargedTrackHypothesis* locNewChargedTrackHypothesis = Build_ChargedTrackHypothesis(locChargedTrackHypothesis, locKinFitParticle, locChargedTrack, locParticleCombo);
					locNewObjectMap[locKinFitParticle] = locNewChargedTrackHypothesis;
					locNewChargedTrackHypothesis->AddAssociatedObject(locParticleCombo);

					_data.push_back(locNewChargedTrackHypothesis);
				}
			}
		}
	}

	return NOERROR;
}

DChargedTrackHypothesis* DChargedTrackHypothesis_factory_KinFit::Build_ChargedTrackHypothesis(const DChargedTrackHypothesis* locChargedTrackHypothesis, DKinFitParticle* locKinFitParticle, const DChargedTrack* locChargedTrack, const DParticleCombo* locParticleCombo)
{
	DChargedTrackHypothesis* locNewChargedTrackHypothesis = new DChargedTrackHypothesis(*locChargedTrackHypothesis);
	locNewChargedTrackHypothesis->AddAssociatedObject(locChargedTrackHypothesis);
	locNewChargedTrackHypothesis->AddAssociatedObject(locChargedTrack);

 	vector<const JObject*> locObjects;
	locChargedTrackHypothesis->GetT(locObjects);
	for(size_t loc_i = 0; loc_i < locObjects.size(); ++loc_i)
		locNewChargedTrackHypothesis->AddAssociatedObject(locObjects[loc_i]);

	locNewChargedTrackHypothesis->setMomentum(DVector3(locKinFitParticle->Get_Momentum().X(),locKinFitParticle->Get_Momentum().Y(),locKinFitParticle->Get_Momentum().Z()));
	locNewChargedTrackHypothesis->setPosition(DVector3(locKinFitParticle->Get_Position().X(),locKinFitParticle->Get_Position().Y(),locKinFitParticle->Get_Position().Z()));
	locNewChargedTrackHypothesis->setTime(locKinFitParticle->Get_Time());
	locNewChargedTrackHypothesis->setErrorMatrix(*locKinFitParticle->Get_CovarianceMatrix());

	double locPathLength = locNewChargedTrackHypothesis->pathLength() - locKinFitParticle->Get_PathLength();
	double locPathLengthUncertainty_Orig = locNewChargedTrackHypothesis->pathLength_err();
	double locPathLengthUncertainty_KinFit = locKinFitParticle->Get_PathLengthUncertainty();
	double locPathLengthUncertainty = sqrt(locPathLengthUncertainty_Orig*locPathLengthUncertainty_Orig + locPathLengthUncertainty_KinFit*locPathLengthUncertainty_KinFit);
	locNewChargedTrackHypothesis->setPathLength(locPathLength, locPathLengthUncertainty);

	dPIDAlgorithm->Calc_ChargedPIDFOM(locNewChargedTrackHypothesis, locParticleCombo->Get_EventRFBunch());

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


