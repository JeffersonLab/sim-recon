#ifndef _DKinFitResults_
#define _DKinFitResults_

#include <string>
#include <deque>

#include "TMatrixDSym.h"

#include "PID/DKinematicData.h"
#include "KINFITTER/DKinFitChain.h"
#include "KINFITTER/DKinFitParticle.h"
#include "KINFITTER/DKinFitConstraint.h"

using namespace std;
using namespace jana;

class DParticleCombo;

class DKinFitResults : public JObject
{
	public:
		DKinFitResults(void){}
		JOBJECT_PUBLIC(DKinFitResults);

		/************************************************************ SET FIT INFORMATION ***********************************************************/

		void Set_KinFitType(DKinFitType locKinFitType){dKinFitType = locKinFitType;}

		void Set_NumConstraints(unsigned int locNumConstraints){dNumConstraints = locNumConstraints;}
		void Set_NumUnknowns(unsigned int locNumUnknowns){dNumUnknowns = locNumUnknowns;}

		void Set_NDF(unsigned int locNDF){dNDF = locNDF;}
		void Set_ChiSq(double locChiSq){dChiSq = locChiSq;}
		void Set_ConfidenceLevel(double locConfidenceLevel){dConfidenceLevel = locConfidenceLevel;}

		void Set_VXi(const TMatrixDSym& locVXi){dVXi = locVXi;}
		//void Set_VEta(const TMatrixDSym* locVEta){dVEta = locVEta;}
		//void Set_VXi(const TMatrixDSym* locVXi){dVXi = locVXi;}
		//void Set_V(const TMatrixDSym* locV){dV = locV;}

		void Set_Pulls(const map<const JObject*, map<DKinFitPullType, double> >& locPulls){dPulls = locPulls;}

		/************************************************************ GET FIT INFORMATION ***********************************************************/

		DKinFitType Get_KinFitType(void) const{return dKinFitType;}

		unsigned int Get_NumConstraints(void) const{return dNumConstraints;}
		unsigned int Get_NumUnknowns(void) const{return dNumUnknowns;}

		unsigned int Get_NDF(void) const{return dNDF;}
		double Get_ChiSq(void) const{return dChiSq;}
		double Get_ConfidenceLevel(void) const{return dConfidenceLevel;}

		const TMatrixDSym& Get_VXi(void) const{return dVXi;}
		//const TMatrixDSym* Get_VEta(void) const{return dVEta;}
		//const TMatrixDSym* Get_VXi(void) const{return dVXi;}
		//const TMatrixDSym* Get_V(void) const{return dV;}

		void Get_Pulls(map<const JObject*, map<DKinFitPullType, double> >& locPulls) const{locPulls = dPulls;}

		/************************************************** SET PARTICLES, COMBOS, AND CONSTRAINTS **************************************************/

		void Add_ParticleCombo(const DParticleCombo* locParticleCombo, const DKinFitChain* locKinFitChain);
		void Add_OutputKinFitParticles(const set<DKinFitParticle*>& locOutputKinFitParticles);
		void Add_KinFitConstraints(const set<DKinFitConstraint*>& locKinFitConstraints);

		void Add_ParticleMapping_SourceToOutput(const JObject* locSourceJObject, DKinFitParticle* locOutputKinFitParticle);

		/************************************************** GET PARTICLES, COMBOS, AND CONSTRAINTS **************************************************/

		void Get_ParticleComboMap(map<const DParticleCombo*, const DKinFitChain*>& locParticleComboMap) const{locParticleComboMap = dParticleComboMap;}
		set<DKinFitParticle*> Get_OutputKinFitParticles(void) const;
		set<DKinFitParticle*> Get_OutputKinFitParticles(DKinFitParticleType locKinFitParticleType) const;
		set<const DKinFitConstraint*> Get_OutputKinFitConstraints(void) const{return dKinFitConstraints;}

		//Source: JObject from DParticleCombo
		//Output: DKinFitParticle's containing the fit results (if not included in fit, is still the INPUT object)
		DKinFitParticle* Get_OutputKinFitParticle(const JObject* locSourceObject) const;

	private:

		DKinFitType dKinFitType;

		unsigned int dNumConstraints;
		unsigned int dNumUnknowns;

		double dConfidenceLevel;
		double dChiSq;
		unsigned int dNDF;

		map<const JObject*, map<DKinFitPullType, double> > dPulls; //JObject is the MEASURED particle (or shower!)

		TMatrixDSym dVXi;
		//const TMatrixDSym* dVXi; //covariance matrix of dXi, the unmeasured parameters in the fit
		//const TMatrixDSym* dVEta; //covariance matrix of dEta, the measured parameters in the fit
		//const TMatrixDSym* dV; //full covariance matrix: dVEta at top-left and dVXi at bottom-right (+ the eta, xi covariance)

		//OUTPUT PARTICLES AND CONSTRAINTS
		map<DKinFitParticleType, set<DKinFitParticle*> > dOutputKinFitParticles; //does not include particles not used in the constraints!
		set<const DKinFitConstraint*> dKinFitConstraints;

		//PARTICLE MAPS
		//Source: JObject from DParticleCombo
		//Output: DKinFitParticle's containing the fit results
		map<const JObject*, DKinFitParticle*> dParticleMap_SourceToOutput;

		//multiple combos may have the same kinfit result, and different DKinFitChain's
		//chain contains output kinfit particles (if a particle not in a constraint, is the input particle)
		map<const DParticleCombo*, const DKinFitChain*> dParticleComboMap;
};

/****************************************************** SET PARTICLES, COMBOS, AND CONSTRAINTS ******************************************************/

inline void DKinFitResults::Add_ParticleCombo(const DParticleCombo* locParticleCombo, const DKinFitChain* locKinFitChain)
{
	dParticleComboMap[locParticleCombo] = locKinFitChain;
}

inline void DKinFitResults::Add_OutputKinFitParticles(const set<DKinFitParticle*>& locOutputKinFitParticles)
{
	set<DKinFitParticle*>::const_iterator locIterator = locOutputKinFitParticles.begin();
	for(; locIterator != locOutputKinFitParticles.end(); ++locIterator)
		dOutputKinFitParticles[(*locIterator)->Get_KinFitParticleType()].insert(*locIterator);
}

inline void DKinFitResults::Add_KinFitConstraints(const set<DKinFitConstraint*>& locKinFitConstraints)
{
	dKinFitConstraints.insert(locKinFitConstraints.begin(), locKinFitConstraints.end());
}

inline void DKinFitResults::Add_ParticleMapping_SourceToOutput(const JObject* locSourceJObject, DKinFitParticle* locOutputKinFitParticle)
{
	dParticleMap_SourceToOutput[locSourceJObject] = locOutputKinFitParticle;
}

/****************************************************** GET PARTICLES, COMBOS, AND CONSTRAINTS ******************************************************/

inline DKinFitParticle* DKinFitResults::Get_OutputKinFitParticle(const JObject* locSourceObject) const
{
	map<const JObject*, DKinFitParticle*>::const_iterator locIterator = dParticleMap_SourceToOutput.find(locSourceObject);
	return ((locIterator != dParticleMap_SourceToOutput.end()) ? locIterator->second : NULL);
}

inline set<DKinFitParticle*> DKinFitResults::Get_OutputKinFitParticles(void) const
{
	set<DKinFitParticle*> locOutputKinFitParticles;
	map<DKinFitParticleType, set<DKinFitParticle*> >::const_iterator locIterator = dOutputKinFitParticles.begin();
	for(; locIterator != dOutputKinFitParticles.end(); ++locIterator)
		locOutputKinFitParticles.insert(locIterator->second.begin(), locIterator->second.end());
	return locOutputKinFitParticles;
}

inline set<DKinFitParticle*> DKinFitResults::Get_OutputKinFitParticles(DKinFitParticleType locKinFitParticleType) const
{
	map<DKinFitParticleType, set<DKinFitParticle*> >::const_iterator locIterator = dOutputKinFitParticles.find(locKinFitParticleType);
	return ((locIterator != dOutputKinFitParticles.end()) ? locIterator->second : set<DKinFitParticle*>());
}

#endif // _DKinFitResults_

