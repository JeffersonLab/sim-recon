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
		inline void Set_TargetCenterZ(double locTargetCenterZ){dTargetCenterZ = locTargetCenterZ;}; //for checking field

		void Reset_NewFit(void);
		void Reset_NewEvent(void);

		const DKinFitParticle* Make_BeamParticle(const DBeamPhoton* locBeamPhoton);
		const DKinFitParticle* Make_DetectedParticle(const DChargedTrackHypothesis* locChargedTrackHypothesis);
		const DKinFitParticle* Make_DetectedShower(const DNeutralParticleHypothesis* locNeutralParticleHypothesis); //DO NOT call this unless the neutral is also in a vertex fit!
		const DKinFitParticle* Make_DetectedParticle(const DNeutralParticleHypothesis* locNeutralParticleHypothesis);
		const DKinFitParticle* Make_DecayingParticle(Particle_t locPID);
		const DKinFitParticle* Make_MissingParticle(Particle_t locPID);
		const DKinFitParticle* Make_TargetParticle(Particle_t locPID);
		using DKinFitter::Make_DetectedParticle; //this is necessary because the above declaration hides the base class function, which is needed by DKinFitResults_factory

		bool Fit_Reaction(void);

		bool Propagate_TrackInfoToCommonVertex(DKinematicData* locKinematicData, const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi);

		inline void Get_ParticleMapping_InputToSource(map<const DKinFitParticle*, const DKinematicData*>& locParticleMapping) const{locParticleMapping = dParticleMapping_InputToSource;}
		inline void Get_ParticleMapping_OutputToSource(map<const DKinFitParticle*, const DKinematicData*>& locParticleMapping) const{locParticleMapping = dParticleMapping_OutputToSource;}
		const DKinematicData* Get_Source_FromInput(const DKinFitParticle* locKinFitParticle) const;
		const DKinematicData* Get_Source_FromOutput(const DKinFitParticle* locKinFitParticle) const;
		inline void Get_Pulls(map<const DKinematicData*, map<DKinFitPullType, double> >& locPulls) const{locPulls = dPulls;} //key is source data (NULL for rf-bunch), 2nd key is param type

	private:
		TVector3 Get_BField(const TVector3& locPosition) const;

		bool dIsBFieldNearBeamline;
		double dTargetCenterZ;
		const DMagneticFieldMap* dMagneticFieldMap;
		map<const DKinFitParticle*, const DKinematicData*> dParticleMapping_InputToSource; //source data is NULL for decaying/target/missing objects!!
		map<const DKinFitParticle*, const DKinematicData*> dParticleMapping_OutputToSource; //source data is NULL for decaying/target/missing objects!!
		//cannot map source to output: source is NULL often!! //User should could keep track. 
		map<const DKinematicData*, map<DKinFitPullType, double> > dPulls; //input source data to pulls

};

#endif // _DKinFitter_GlueX_

