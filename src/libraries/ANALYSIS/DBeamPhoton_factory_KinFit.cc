// $Id$
//
//    File: DBeamPhoton_factory_KinFit.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifdef VTRACE
#include "vt_user.h"
#endif

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
#ifdef VTRACE
	VT_TRACER("DBeamPhoton_factory_KinFit::evnt()");
#endif

 	vector<const DKinFitResults*> locKinFitResultsVector;
	locEventLoop->Get(locKinFitResultsVector);

	map<const DKinFitParticle*, DBeamPhoton*> locKinFitParticleMap;

	for(size_t loc_i = 0; loc_i < locKinFitResultsVector.size(); ++loc_i)
	{
		set<const DParticleCombo*> locParticleCombos;
		locKinFitResultsVector[loc_i]->Get_ParticleCombos(locParticleCombos);

		set<const DParticleCombo*>::iterator locComboIterator = locParticleCombos.begin();
		const DParticleCombo* locParticleCombo = *locComboIterator;

		if(!locParticleCombo->Get_ParticleComboStep(0)->Is_InitialParticleDetected())
			continue;

		map<const DKinematicData*, const DKinFitParticle*> locReverseParticleMapping;
		locKinFitResultsVector[loc_i]->Get_ReverseParticleMapping(locReverseParticleMapping);

		const DKinematicData* locInitialParticle = locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticle_Measured();
		const DKinFitParticle* locKinFitParticle = locReverseParticleMapping[locInitialParticle];

		map<const DKinFitParticle*, DBeamPhoton*>::iterator locNewPhotonIterator = locKinFitParticleMap.find(locKinFitParticle);
		if(locNewPhotonIterator != locKinFitParticleMap.end())
		{
			for(locComboIterator = locParticleCombos.begin(); locComboIterator != locParticleCombos.end(); ++locComboIterator)
				locNewPhotonIterator->second->AddAssociatedObject(*locComboIterator);
			continue; //new particle already created for this kinfit particle
		}

		const DBeamPhoton* locBeamPhoton = static_cast<const DBeamPhoton*>(locInitialParticle);
		DBeamPhoton* locNewBeamPhoton = Build_BeamPhoton(locBeamPhoton, locKinFitParticle, locParticleCombo);
		locKinFitParticleMap[locKinFitParticle] = locNewBeamPhoton;

		for(locComboIterator = locParticleCombos.begin(); locComboIterator != locParticleCombos.end(); ++locComboIterator)
			locNewBeamPhoton->AddAssociatedObject(*locComboIterator);

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


