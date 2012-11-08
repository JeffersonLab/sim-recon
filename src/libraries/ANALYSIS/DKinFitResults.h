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
	d_P4Fit, 
	d_VertexFit,
	d_SpacetimeFit,
	d_P4AndVertexFit, 
	d_P4AndSpacetimeFit
};

class DKinFitResults : public JObject
{
	public:
		JOBJECT_PUBLIC(DKinFitResults);

		inline void Set_KinFitName(string locKinFitName){dKinFitName = locKinFitName;}
		inline string Get_KinFitName(void) const{return dKinFitName;}

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

		inline void Get_InitialKinFitParticles(deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles) const{locInitialKinFitParticles = dInitialKinFitParticles;}
		void Set_InitialKinFitParticles(const deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles){dInitialKinFitParticles = locInitialKinFitParticles;}

		inline void Get_FinalKinFitParticles(deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles) const{locFinalKinFitParticles = dFinalKinFitParticles;}
		void Set_FinalKinFitParticles(const deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles){dFinalKinFitParticles = locFinalKinFitParticles;}

		inline void Set_VEta(const TMatrixDSym* locVEta){dVEta = locVEta;}
		inline const TMatrixDSym* Get_VEta(void) const{return dVEta;}

		inline void Set_VXi(const TMatrixDSym* locVXi){dVXi = locVXi;}
		inline const TMatrixDSym* Get_VXi(void) const{return dVXi;}

	private:
		string dKinFitName; //is DReaction name
		double dConfidenceLevel;
		double dChiSq;
		unsigned int dNDF;
		const TMatrixDSym* dVXi; //covariance matrix of dXi, the unmeasured parameters in the fit
		const TMatrixDSym* dVEta; //covariance matrix of dEta, the measured parameters in the fit

		map<const DKinematicData*, map<DKinFitPullType, double> > dPulls; //DKinematicData is the MEASURED particle
		map<const DKinFitParticle*, const DKinematicData*> dParticleMapping; //from output -> source (MEASURED) //is NULL for missing/decaying/target particles

		const DParticleCombo* dParticleCombo;
		DKinFitType dKinFitType;

		//these particles correspond to the dParticleCombo particles (1st dim is step index, 2nd is particle index
			//for initial: if 2nd dimension is size 2, then 2nd object is the target
			//note that these particles can be NULL: if not involved in kinfit (e.g. p4 kinfit, excluded or resonance particle)
		deque<deque<const DKinFitParticle*> > dInitialKinFitParticles;
		deque<deque<const DKinFitParticle*> > dFinalKinFitParticles;

		deque<const DKinFitConstraint*> dKinFitConstraints;
};

#endif // _DKinFitResults_

