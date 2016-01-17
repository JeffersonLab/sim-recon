#ifndef _DKinFitResults_
#define _DKinFitResults_

#include <string>
#include <deque>

#include "PID/DKinematicData.h"
#include "KINFITTER/DKinFitChain.h"
#include "KINFITTER/DKinFitParticle.h"
#include "KINFITTER/DKinFitConstraint.h"
#include "ANALYSIS/DParticleCombo.h"

using namespace std;
using namespace jana;

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

		void Set_VEta(const TMatrixDSym* locVEta){dVEta = locVEta;}
		void Set_VXi(const TMatrixDSym* locVXi){dVXi = locVXi;}
		void Set_V(const TMatrixDSym* locV){dV = locV;}

		void Set_Pulls(const map<const JObject*, map<DKinFitPullType, double> >& locPulls){dPulls = locPulls;}

		/************************************************************ GET FIT INFORMATION ***********************************************************/

		DKinFitType Get_KinFitType(void) const{return dKinFitType;}

		unsigned int Get_NumConstraints(void) const{return dNumConstraints;}
		unsigned int Get_NumUnknowns(void) const{return dNumUnknowns;}

		unsigned int Get_NDF(void) const{return dNDF;}
		double Get_ChiSq(void) const{return dChiSq;}
		double Get_ConfidenceLevel(void) const{return dConfidenceLevel;}

		const TMatrixDSym* Get_VEta(void) const{return dVEta;}
		const TMatrixDSym* Get_VXi(void) const{return dVXi;}
		const TMatrixDSym* Get_V(void) const{return dV;}

		void Get_Pulls(map<const JObject*, map<DKinFitPullType, double> >& locPulls) const{locPulls = dPulls;}

		/************************************************** SET PARTICLES, COMBOS, AND CONSTRAINTS **************************************************/

		void Add_ParticleCombo(const DParticleCombo* locParticleCombo, const DKinFitChain* locKinFitChain);
		void Add_OutputKinFitParticles(const set<DKinFitParticle*>& locOutputKinFitParticles);
		void Add_KinFitConstraints(const set<DKinFitParticle*>& locKinFitConstraints);

		void Add_ParticleMapping_SourceToInput(const JObject* locSourceJObject, const DKinFitParticle* locInputKinFitParticle);
		void Add_ParticleMapping_InputToOutput(const DKinFitParticle* locInputKinFitParticle, const DKinFitParticle* locOutputKinFitParticle);

		/************************************************** GET PARTICLES, COMBOS, AND CONSTRAINTS **************************************************/

		void Get_ParticleComboMap(map<const DParticleCombo*, const DKinFitChain*>& locParticleComboMap) const{locParticleComboMap = dParticleComboMap;}
		set<const DKinFitParticle*> Get_OutputKinFitParticles(DKinFitParticleType locKinFitParticleType) const;
		set<const DKinFitConstraint*> Get_OutputKinFitConstraints(void) const{return dKinFitConstraints;}

		//Source: JObject from DParticleCombo
		//Input: DKinFitParticle used to create the kinematic fit constraints
		//Output: DKinFitParticle's containing the fit results
		const DKinFitParticle* Get_InputKinFitParticle(const JObject* locSourceObject) const;
		const DKinFitParticle* Get_OutputKinFitParticle(const DKinFitParticle* locInputKinFitParticle) const;
		const DKinFitParticle* Get_OutputKinFitParticle(const JObject* locSourceObject) const;

	private:

		DKinFitType dKinFitType;

		unsigned int dNumConstraints;
		unsigned int dNumUnknowns;

		double dConfidenceLevel;
		double dChiSq;
		unsigned int dNDF;

		map<const JObject*, map<DKinFitPullType, double> > dPulls; //JObject is the MEASURED particle (or shower!)

		const TMatrixDSym* dVXi; //covariance matrix of dXi, the unmeasured parameters in the fit
		const TMatrixDSym* dVEta; //covariance matrix of dEta, the measured parameters in the fit
		const TMatrixDSym* dV; //full covariance matrix: dVEta at top-left and dVXi at bottom-right (+ the eta, xi covariance)

		//OUTPUT PARTICLES AND CONSTRAINTS
		map<DKinFitParticleType, set<const DKinFitParticle*> > dOutputKinFitParticles;
		set<const DKinFitConstraint*> dKinFitConstraints;

		//PARTICLE MAPS
		//Source: JObject from DParticleCombo
		//Input: DKinFitParticle used to create the kinematic fit constraints
		//Output: DKinFitParticle's containing the fit results
		map<const JObject*, const DKinFitParticle*> dParticleMap_SourceToInput;
		map<const DKinFitParticle*, const DKinFitParticle*> dParticleMap_InputToOutput;

		//multiple combos may have the same kinfit result, and different DKinFitChain's
		map<const DParticleCombo*, const DKinFitChain*> dParticleComboMap; //chain contains INPUT kinfit particles!
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

inline void DKinFitResults::Add_KinFitConstraints(const set<DKinFitParticle*>& locKinFitConstraints)
{
	set<DKinFitParticle*>::const_iterator locIterator = locKinFitConstraints.begin();
	for(; locIterator != locKinFitConstraints.end(); ++locIterator)
		dKinFitConstraints.insert(*locIterator);
}

inline void DKinFitResults::Add_ParticleMapping_SourceToInput(const JObject* locSourceJObject, const DKinFitParticle* locInputKinFitParticle)
{
	dParticleMap_SourceToInput[locSourceJObject] = locInputKinFitParticle;
}

inline void DKinFitResults::Add_ParticleMapping_InputToOutput(const DKinFitParticle* locInputKinFitParticle, const DKinFitParticle* locOutputKinFitParticle)
{
	dParticleMap_InputToOutput[locInputKinFitParticle] = locOutputKinFitParticle;
}

/****************************************************** GET PARTICLES, COMBOS, AND CONSTRAINTS ******************************************************/

inline const DKinFitParticle* DKinFitResults::Get_InputKinFitParticle(const JObject* locSourceObject) const
{
	map<const JObject*, const DKinFitParticle*>::const_iterator locIterator = dParticleMap_SourceToInput.find(locSourceObject);
	return (locIterator != dParticleMap_SourceToInput.end()) ? locIterator->second : NULL;
}

inline const DKinFitParticle* DKinFitResults::Get_OutputKinFitParticle(const DKinFitParticle* locInputKinFitParticle) const
{
	map<const DKinFitParticle*, const DKinFitParticle*>::const_iterator locIterator = dParticleMap_InputToOutput.find(locInputKinFitParticle);
	return (locIterator != dParticleMap_InputToOutput.end()) ? locIterator->second : NULL;
}

inline const DKinFitParticle* DKinFitResults::Get_OutputKinFitParticle(const JObject* locSourceObject) const
{
	DKinFitParticle* locInputKinFitParticle = Get_InputKinFitParticle(locSourceObject);
	return Get_OutputKinFitParticle(locInputKinFitParticle);
}

inline set<const DKinFitParticle*> DKinFitResults::Get_OutputKinFitParticles(DKinFitParticleType locKinFitParticleType) const
{
	map<DKinFitParticleType, set<const DKinFitParticle*> >::const_iterator locIterator = dOutputKinFitParticles.find(locKinFitParticleType);
	return ((locIterator != dOutputKinFitParticles.end()) ? locIterator->second : set<const DKinFitParticle*>();
}

#endif // _DKinFitResults_

