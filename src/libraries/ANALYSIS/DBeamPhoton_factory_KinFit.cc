// $Id$
//
//    File: DBeamPhoton_factory_KinFit.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DBeamPhoton_factory_KinFit.h"

//------------------
// init
//------------------
jerror_t DBeamPhoton_factory_KinFit::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DBeamPhoton_factory_KinFit::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBeamPhoton_factory_KinFit::evnt(jana::JEventLoop* locEventLoop, int eventnumber)
{
 	vector<const DKinFitResults*> locKinFitResultsVector;
	locEventLoop->Get(locKinFitResultsVector);

	const DParticleCombo* locParticleCombo;
	deque<deque<const DKinFitParticle*> > locInitialKinFitParticles;
	const DKinFitParticle* locKinFitParticle;
	const DBeamPhoton* locBeamPhoton;
	DBeamPhoton* locNewBeamPhoton;

	map<const DKinFitParticle*, DBeamPhoton*> locKinFitParticleMap;
	map<DBeamPhoton*, deque<const DParticleCombo*> > locBeamParticleComboMap;

	for(size_t loc_i = 0; loc_i < locKinFitResultsVector.size(); ++loc_i)
	{
		locParticleCombo = locKinFitResultsVector[loc_i]->Get_ParticleCombo();
		if(!locParticleCombo->Get_ParticleComboStep(0)->Is_InitialParticleDetected())
			continue;
		locKinFitResultsVector[loc_i]->Get_InitialKinFitParticles(locInitialKinFitParticles);
		locKinFitParticle = locInitialKinFitParticles[0][0];

		//loop over previous kinfitresults objects, see if the kinfit results are identical: if so, don't need to correct new tracks for the results
		bool locMatchFlag = false;
		for(size_t loc_j = 0; loc_j < loc_i; ++loc_j)
		{
			if(!locParticleCombo->Will_KinFitBeIdentical(locKinFitResultsVector[loc_j]->Get_ParticleCombo()))
				continue;
			//kinfit results are identical: setup the maps so that the charged tracks are copied instead of created anew
			locNewBeamPhoton = locKinFitParticleMap[locKinFitParticle];
			locBeamParticleComboMap[locNewBeamPhoton].push_back(locParticleCombo);
			locMatchFlag = true;
			break;
		}
		if(locMatchFlag)
			continue;

		locBeamPhoton = static_cast<const DBeamPhoton*>(locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticle_Measured());
		locNewBeamPhoton = Build_BeamPhoton(locBeamPhoton, locKinFitParticle, locParticleCombo);
		locKinFitParticleMap[locKinFitParticle] = locNewBeamPhoton;
		locBeamParticleComboMap[locNewBeamPhoton] = deque<const DParticleCombo*>(1, locParticleCombo);
	}

	//now set the particle combos as associated objects of the beam photons, and save the tracks //this marks which combos they originated from
	map<DBeamPhoton*, deque<const DParticleCombo*> >::iterator locIterator;
	for(locIterator = locBeamParticleComboMap.begin(); locIterator != locBeamParticleComboMap.end(); ++locIterator)
	{
		locNewBeamPhoton = locIterator->first;
		deque<const DParticleCombo*>& locParticleCombos = locIterator->second;
		for(size_t loc_i = 0; loc_i < locParticleCombos.size(); ++loc_i)
			locNewBeamPhoton->AddAssociatedObject(locParticleCombos[loc_i]);
		_data.push_back(locNewBeamPhoton);
	}

	return NOERROR;
}

DBeamPhoton* DBeamPhoton_factory_KinFit::Build_BeamPhoton(const DBeamPhoton* locBeamPhoton, const DKinFitParticle* locKinFitParticle, const DParticleCombo* locParticleCombo)
{
	DBeamPhoton* locNewBeamPhoton = new DBeamPhoton(*locBeamPhoton);
	locNewBeamPhoton->AddAssociatedObject(locBeamPhoton);

	locNewBeamPhoton->setMomentum(DVector3(locKinFitParticle->Get_Momentum().X(),locKinFitParticle->Get_Momentum().Y(),locKinFitParticle->Get_Momentum().Z()));
	locNewBeamPhoton->setPosition(DVector3(locKinFitParticle->Get_Position().X(),locKinFitParticle->Get_Position().Y(),locKinFitParticle->Get_Position().Z()));
	locNewBeamPhoton->setTime(locKinFitParticle->Get_Time());
	locNewBeamPhoton->setErrorMatrix(*locKinFitParticle->Get_CovarianceMatrix());
	double locPathLength = locBeamPhoton->pathLength() - locKinFitParticle->Get_PathLength(); //locKinFitParticle->Get_PathLength() = (X_common - X_track).Dot(UnitP)
	double locPathLengthUncertainty_Orig = locBeamPhoton->pathLength_err();
	double locPathLengthUncertainty_KinFit = locKinFitParticle->Get_PathLengthUncertainty();
	double locPathLengthUncertainty = sqrt(locPathLengthUncertainty_Orig*locPathLengthUncertainty_Orig + locPathLengthUncertainty_KinFit*locPathLengthUncertainty_KinFit);
	locNewBeamPhoton->setPathLength(locPathLength, locPathLengthUncertainty);
	return locNewBeamPhoton;
}

//------------------
// erun
//------------------
jerror_t DBeamPhoton_factory_KinFit::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DBeamPhoton_factory_KinFit::fini(void)
{
	return NOERROR;
}


