// $Id$
//
//    File: DParticle_factory_THROWN.h
// Created: Sat Oct  4 22:04:56 EDT 2008
// Creator: davidl (on Darwin Amelia.local 8.11.1 i386)
//

#ifndef _DParticle_factory_THROWN_
#define _DParticle_factory_THROWN_

#include <JANA/JFactory.h>
#include <TRACKING/DReferenceTrajectory.h>
#include "DParticle.h"

class DParticle_factory_THROWN:public jana::JFactory<DParticle>{
	public:
		DParticle_factory_THROWN();
		~DParticle_factory_THROWN(){};
		const char* Tag(void){return "THROWN";}

		typedef DReferenceTrajectory::swim_step_t swim_step_t;

	private:

		enum fit_type_t{
			kWireBased,
			kTimeBased
		};
		
		class hitInfo{
			public:
				vector<const DCoordinateSystem*> wires;	// Wire definitions
				vector<DVector3> shifts;						// Effective wire shifts due to drift time
				vector<double> errs;								// Errors on drift time (or wire position) measurement
				vector<double> u_dists;							// Distances along the wire (for FDC cathodes)
				vector<double> u_errs;							// Errors on distance along the wire (for FDC cathodes)
				vector<double> all_errs;						// Merging of errs and u_errs so elements correspond to those in chisqv

				void PrintDebug(void){_DBG_<<"sizes: wires="<<wires.size()<<" shifts="<<shifts.size()<<" errs="<<errs.size()<<" u_dists="<<u_dists.size()<<" u_errs="<<u_errs.size()<<" all_errs="<<all_errs.size()<<endl;}
		};

		DCoordinateSystem *target;
		
		int DEBUG_LEVEL;
		double SIGMA_CDC;
		double SIGMA_FDC_ANODE;
		double SIGMA_FDC_CATHODE;
		
		//jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		//jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//jerror_t fini(void);						///< Called after last event of last event source has been processed.

		vector<DReferenceTrajectory*> rt_pool;
		vector<DMatrixDSym*> cov;
		vector<const DCDCTrackHit*> cdchits;
		vector<const DFDCPseudo*> fdchits;
		
		void AddCDCTrackHits(DReferenceTrajectory *rt, vector<const DCDCTrackHit*> &cdctrackhits);
		void AddFDCPseudoHits(DReferenceTrajectory *rt, vector<const DFDCPseudo*> &fdcpseudos);

		double ChiSq(DReferenceTrajectory *rt, hitInfo &hinfo, vector<double> &chisqv, double *chisq_ptr, int *dof_ptr);
		void GetWiresShiftsErrs(fit_type_t fit_type, DReferenceTrajectory *rt, hitInfo &hinfo);

};

#endif // _DParticle_factory_THROWN_

