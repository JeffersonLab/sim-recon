#ifndef _DKinFitter_GlueX_
#define _DKinFitter_GlueX_

#include "TVector3.h"
#include "TLorentzVector.h"

#include "particleType.h"
#include "HDGEOMETRY/DMagneticFieldMap.h"

#include "PID/DBeamPhoton.h"
#include "PID/DNeutralShower.h"
#include "PID/DKinematicData.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DNeutralParticleHypothesis.h"

#include "ANALYSIS/DKinFitter.h"

using namespace std;

class DKinFitter_GlueX : public DKinFitter
{
	public:
		DKinFitter_GlueX(void);

		bool Get_IsBFieldNearBeamline(void) const;
		void Set_BField(const DMagneticFieldMap* locMagneticFieldMap);

		void Reset_NewFit(void);
		void Reset_NewEvent(void);

		const DKinFitParticle* Make_BeamParticle(const DBeamPhoton* locBeamPhoton);
		const DKinFitParticle* Make_DetectedParticle(const DChargedTrackHypothesis* locChargedTrackHypothesis);
		const DKinFitParticle* Make_DetectedParticle(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, bool locWillOnlyBeUsedInP4FitFlag = false); //if true uses particle, false uses associated shower
		const DKinFitParticle* Make_DecayingParticle(Particle_t locPID);
		const DKinFitParticle* Make_MissingParticle(Particle_t locPID);
		const DKinFitParticle* Make_TargetParticle(Particle_t locPID);

		bool Fit_Reaction(void);

		bool Propagate_TrackInfoToCommonVertex(DKinematicData* locKinematicData, const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi);

		void Get_ParticleMapping_InputToSource(map<const DKinFitParticle*, const DKinematicData*>& locParticleMapping) const{locParticleMapping = dParticleMapping_InputToSource;}
		void Get_ParticleMapping_OutputToSource(map<const DKinFitParticle*, const DKinematicData*>& locParticleMapping) const{locParticleMapping = dParticleMapping_OutputToSource;}
		inline void Get_Pulls(map<const DKinematicData*, map<DKinFitPullType, double> >& locPulls) const{locPulls = dPulls;} //key is source data (NULL for rf-bunch), 2nd key is param type

	private:
		TVector3 Get_BField(const TVector3& locPosition) const;

		bool dIsBFieldNearBeamline;
		const DMagneticFieldMap* dMagneticFieldMap;
		map<const DKinFitParticle*, const DKinematicData*> dParticleMapping_InputToSource; //source data is NULL for decaying/target/missing objects!!
		map<const DKinFitParticle*, const DKinematicData*> dParticleMapping_OutputToSource; //source data is NULL for decaying/target/missing objects!!
		//cannot map source to output: source is NULL often!! //User should could keep track. 
		map<const DKinematicData*, map<DKinFitPullType, double> > dPulls; //input source data to pulls

};

#endif // _DKinFitter_GlueX_

