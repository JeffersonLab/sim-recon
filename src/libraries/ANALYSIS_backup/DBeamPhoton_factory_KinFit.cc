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
jerror_t DBeamPhoton_factory_KinFit::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBeamPhoton_factory_KinFit::evnt(jana::JEventLoop* locEventLoop, uint64_t eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DBeamPhoton_factory_KinFit::evnt()");
#endif

 	vector<const DKinFitResults*> locKinFitResultsVector;
	locEventLoop->Get(locKinFitResultsVector);

	for(size_t loc_i = 0; loc_i < locKinFitResultsVector.size(); ++loc_i)
	{
		map<const DParticleCombo*, const DKinFitChain*> locParticleComboMap;
		locKinFitResultsVector[loc_i]->Get_ParticleComboMap(locParticleComboMap);
		set<DKinFitParticle*> locOutputKinFitParticles = locKinFitResultsVector[loc_i]->Get_OutputKinFitParticles();

		map<DKinFitParticle*, DBeamPhoton*> locNewObjectMap;
		map<const DParticleCombo*, const DKinFitChain*>::iterator locComboIterator = locParticleComboMap.begin();
		for(; locComboIterator != locParticleComboMap.end(); ++locComboIterator)
		{
			const DParticleCombo* locParticleCombo = locComboIterator->first;

			if(!locParticleCombo->Get_ParticleComboStep(0)->Is_InitialParticleDetected())
				continue;

			const DKinematicData* locInitialParticle = locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticle_Measured();
			DKinFitParticle* locKinFitParticle = locKinFitResultsVector[loc_i]->Get_OutputKinFitParticle(locInitialParticle);
			if(locKinFitParticle == NULL)
				continue; //should be impossible
			if(locOutputKinFitParticles.find(locKinFitParticle) == locOutputKinFitParticles.end())
				continue; //not used in fit

			map<DKinFitParticle*, DBeamPhoton*>::iterator locNewPhotonIterator = locNewObjectMap.find(locKinFitParticle);
			if(locNewPhotonIterator != locNewObjectMap.end())
			{
				locNewPhotonIterator->second->AddAssociatedObject(locParticleCombo);
				continue; //new particle already created for this kinfit particle
			}

			const DBeamPhoton* locBeamPhoton = static_cast<const DBeamPhoton*>(locInitialParticle);
			DBeamPhoton* locNewBeamPhoton = Build_BeamPhoton(locBeamPhoton, locKinFitParticle, locParticleCombo);
			locNewObjectMap[locKinFitParticle] = locNewBeamPhoton;
			locNewBeamPhoton->AddAssociatedObject(locParticleCombo);

			_data.push_back(locNewBeamPhoton);
		}
	}

	return NOERROR;
}

DBeamPhoton* DBeamPhoton_factory_KinFit::Build_BeamPhoton(const DBeamPhoton* locBeamPhoton, DKinFitParticle* locKinFitParticle, const DParticleCombo* locParticleCombo)
{
	DBeamPhoton* locNewBeamPhoton = new DBeamPhoton(*locBeamPhoton);
	locNewBeamPhoton->AddAssociatedObject(locBeamPhoton);

	locNewBeamPhoton->setMomentum(DVector3(locKinFitParticle->Get_Momentum().X(),locKinFitParticle->Get_Momentum().Y(),locKinFitParticle->Get_Momentum().Z()));
	locNewBeamPhoton->setPosition(DVector3(locKinFitParticle->Get_Position().X(),locKinFitParticle->Get_Position().Y(),locKinFitParticle->Get_Position().Z()));
	locNewBeamPhoton->setTime(locKinFitParticle->Get_Time());
	locNewBeamPhoton->setErrorMatrix(locKinFitParticle->Get_CovarianceMatrix());
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


