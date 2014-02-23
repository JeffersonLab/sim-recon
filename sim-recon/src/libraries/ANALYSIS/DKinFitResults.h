#ifndef _DKinFitResults_
#define _DKinFitResults_

#include <string>
#include <deque>

#include "PID/DKinematicData.h"
#include "ANALYSIS/DKinFitParticle.h"
#include "ANALYSIS/DKinFitConstraints.h"

using namespace std;
using namespace jana;

class DParticleCombo;

enum DKinFitType
{
	d_NoFit = 0, 
	d_P4Fit, //also includes invariant mass constraints
	d_VertexFit,
	d_SpacetimeFit,
	d_P4AndVertexFit, //also includes invariant mass constraints
	d_P4AndSpacetimeFit //also includes invariant mass constraints
};

class DKinFitResults : public JObject
{
	public:
		DKinFitResults(void) : dMissingParticle(NULL) {}
		JOBJECT_PUBLIC(DKinFitResults);

		inline void Set_KinFitType(DKinFitType locKinFitType){dKinFitType = locKinFitType;}
		inline DKinFitType Get_KinFitType(void) const{return dKinFitType;}

		inline void Set_ConfidenceLevel(double locConfidenceLevel){dConfidenceLevel = locConfidenceLevel;}
		inline double Get_ConfidenceLevel(void) const{return dConfidenceLevel;}

		inline void Set_ChiSq(double locChiSq){dChiSq = locChiSq;}
		inline double Get_ChiSq(void) const{return dChiSq;}

		inline void Set_NDF(unsigned int locNDF){dNDF = locNDF;}
		inline unsigned int Get_NDF(void) const{return dNDF;}

		inline void Set_ParticleCombo(const DParticleCombo* locParticleCombo){dParticleCombo = locParticleCombo;}
		inline const DParticleCombo* Get_ParticleCombo(void) const{return dParticleCombo;}

		inline void Set_Pulls(const map<const DKinematicData*, map<DKinFitPullType, double> >& locPulls){dPulls = locPulls;}
		inline void Get_Pulls(map<const DKinematicData*, map<DKinFitPullType, double> >& locPulls) const{locPulls = dPulls;}

		inline void Set_KinFitConstraints(const deque<const DKinFitConstraint*>& locKinFitConstraints){dKinFitConstraints = locKinFitConstraints;}
		inline void Get_KinFitConstraints(deque<const DKinFitConstraint*>& locKinFitConstraints) const{locKinFitConstraints = dKinFitConstraints;}

		inline void Set_ParticleMapping(const map<const DKinFitParticle*, const DKinematicData*>& locParticleMapping){dParticleMapping = locParticleMapping;}
		inline void Get_ParticleMapping(map<const DKinFitParticle*, const DKinematicData*>& locParticleMapping) const{locParticleMapping = dParticleMapping;}

		inline void Set_ReverseParticleMapping(const map<const DKinematicData*, const DKinFitParticle*>& locReverseParticleMapping){dReverseParticleMapping = locReverseParticleMapping;}
		inline void Get_ReverseParticleMapping(map<const DKinematicData*, const DKinFitParticle*>& locReverseParticleMapping) const{locReverseParticleMapping = dReverseParticleMapping;}

		inline void Set_MissingParticle(const DKinFitParticle* locMissingParticle){dMissingParticle = locMissingParticle;}
		inline const DKinFitParticle* Get_MissingParticle(void) const{return dMissingParticle;}

		inline void Add_DecayingParticle(pair<Particle_t, deque<const DKinematicData*> > locDecayPair, const DKinFitParticle* locKinFitParticle){dDecayingParticles[locDecayPair] = locKinFitParticle;}
		inline void Set_DecayingParticles(const map<pair<Particle_t, deque<const DKinematicData*> >, const DKinFitParticle*>& locDecayingParticles){dDecayingParticles = locDecayingParticles;}
		inline void Get_DecayingParticles(map<pair<Particle_t, deque<const DKinematicData*> >, const DKinFitParticle*>& locDecayingParticles) const{locDecayingParticles = dDecayingParticles;}

		inline void Set_VEta(const TMatrixDSym* locVEta){dVEta = locVEta;}
		inline const TMatrixDSym* Get_VEta(void) const{return dVEta;}

		inline void Set_VXi(const TMatrixDSym* locVXi){dVXi = locVXi;}
		inline const TMatrixDSym* Get_VXi(void) const{return dVXi;}

	private:
		double dConfidenceLevel;
		double dChiSq;
		unsigned int dNDF;
		map<const DKinematicData*, map<DKinFitPullType, double> > dPulls; //DKinematicData is the MEASURED particle

		const TMatrixDSym* dVXi; //covariance matrix of dXi, the unmeasured parameters in the fit
		const TMatrixDSym* dVEta; //covariance matrix of dEta, the measured parameters in the fit

		//PARTICLES: DETECTED (including beam)
		map<const DKinFitParticle*, const DKinematicData*> dParticleMapping; //from output -> source (MEASURED) //is NULL for missing/decaying/target particles
		map<const DKinematicData*, const DKinFitParticle*> dReverseParticleMapping; //from source (MEASURED) -> output //decaying/missing/target particles aren't present

		//PARTICLES: NON-DETECTED
		map<pair<Particle_t, deque<const DKinematicData*> >, const DKinFitParticle*> dDecayingParticles; //key is PID of decaying particle + all final-state decay products
		const DKinFitParticle* dMissingParticle; //NULL if none

		const DParticleCombo* dParticleCombo;
		DKinFitType dKinFitType;

		deque<const DKinFitConstraint*> dKinFitConstraints;
};

#endif // _DKinFitResults_

