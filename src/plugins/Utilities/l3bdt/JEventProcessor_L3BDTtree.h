// $Id$
//
//    File: JEventProcessor_L3BDTtree.h
// Created: Wed May 11 22:26:46 EDT 2016
// Creator: davidl (on Linux gluon49.jlab.org 2.6.32-431.20.3.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_L3BDTtree_
#define _JEventProcessor_L3BDTtree_

#include <TTree.h>

#include <JANA/JEventProcessor.h>

#include <DAQ/Df250TriggerTime.h>
#include <DAQ/Df250PulseData.h>
#include <DAQ/Df125TriggerTime.h>
#include <DAQ/Df125PulseIntegral.h>
#include <DAQ/Df125PulseTime.h>
#include <DAQ/Df125PulsePedestal.h>
#include <DAQ/Df125CDCPulse.h>
#include <DAQ/Df125FDCPulse.h>
#include <DAQ/DF1TDCHit.h>
#include <DAQ/DF1TDCTriggerTime.h>
#include <DAQ/DCAEN1290TDCHit.h>
#include <DAQ/DL1Info.h>
#include <DAQ/DEventTag.h>
#include <PID/DVertex.h>

#include <BCAL/DBCALDigiHit.h>
#include <BCAL/DBCALTDCDigiHit.h>
#include <CDC/DCDCDigiHit.h>
#include <FCAL/DFCALDigiHit.h>
#include <FDC/DFDCWireDigiHit.h>
#include <FDC/DFDCCathodeDigiHit.h>
#include <PAIR_SPECTROMETER/DPSCTDCDigiHit.h>
#include <PAIR_SPECTROMETER/DPSCDigiHit.h>
#include <PAIR_SPECTROMETER/DPSDigiHit.h>
#include <START_COUNTER/DSCDigiHit.h>
#include <START_COUNTER/DSCTDCDigiHit.h>
#include <START_COUNTER/DSCHit.h>
#include <TAGGER/DTAGHDigiHit.h>
#include <TAGGER/DTAGHTDCDigiHit.h>
#include <TAGGER/DTAGMDigiHit.h>
#include <TAGGER/DTAGMTDCDigiHit.h>
#include <TOF/DTOFTDCDigiHit.h>
#include <TOF/DTOFDigiHit.h>
#include <TOF/DTOFHit.h>
#include <TOF/DTOFPoint.h>
#include <TPOL/DTPOLRingDigiHit.h>
#include <TPOL/DTPOLSectorDigiHit.h>

#include <BCAL/DBCALPoint.h>
#include <BCAL/DBCALCluster.h>
#include <FCAL/DFCALCluster.h>
#include <PID/DNeutralParticle.h>
#include <PID/DNeutralShower.h>
#include <PID/DChargedTrack.h>
#include <PID/DBeamPhoton.h>
#include <PID/DEventRFBunch.h>
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DTrackWireBased.h>
#include <TRACKING/DTrackTimeBased.h>

#define MyTypes(X) \
	X(Df250TriggerTime) \
	X(Df250PulseData) \
	X(Df125TriggerTime) \
	X(Df125PulseIntegral) \
	X(Df125PulseTime) \
	X(Df125PulsePedestal) \
	X(Df125CDCPulse) \
	X(Df125FDCPulse) \
	X(DF1TDCHit) \
	X(DF1TDCTriggerTime) \
	X(DCAEN1290TDCHit) \
	X(DL1Info) \
	X(DEventTag) \
	X(DVertex) \
	\
	X(DBCALDigiHit) \
	X(DBCALTDCDigiHit) \
	X(DCDCDigiHit) \
	X(DFCALDigiHit) \
	X(DFDCWireDigiHit) \
	X(DFDCCathodeDigiHit) \
	X(DPSCTDCDigiHit) \
	X(DPSCDigiHit) \
	X(DPSDigiHit) \
	X(DSCDigiHit) \
	X(DSCTDCDigiHit) \
	X(DTAGHDigiHit) \
	X(DTAGHTDCDigiHit) \
	X(DTAGMDigiHit) \
	X(DTAGMTDCDigiHit) \
	X(DTOFTDCDigiHit) \
	X(DTOFDigiHit) \
	X(DTPOLRingDigiHit) \
	X(DTPOLSectorDigiHit) \
	\
	X(DBCALPoint) \
	X(DBCALCluster) \
	X(DFCALCluster) \
	X(DNeutralParticle) \
	X(DNeutralShower) \
	X(DChargedTrack) \
	X(DBeamPhoton) \
	X(DEventRFBunch) \
	X(DSCHit) \
	X(DTOFHit) \
	X(DTOFPoint) \
	X(DTrackCandidate) \
	X(DTrackWireBased) \
	X(DTrackTimeBased)

#define MyIntBranchTypes(X) \
	X(trig_mask) \
	X(fp_trig_mask) \
	\
	X(NCDC_superlayer1) \
	X(NCDC_superlayer2) \
	X(NCDC_superlayer3) \
	X(NCDC_superlayer4) \
	X(NCDC_superlayer5) \
	X(NFDCwires_package1) \
	X(NFDCwires_package2) \
	X(NFDCwires_package3) \
	X(NFDCwires_package4) \
	X(NFDCCathodes_package1) \
	X(NFDCCathodes_package2) \
	X(NFDCCathodes_package3) \
	X(NFDCCathodes_package4) \
	X(NTOF_half_length) \
	X(NTOF_half_width) \
	\
	X(Nbeam_photons_coherent) \
	X(Nbeam_photons_3_4) \
	X(Nbeam_photons_4_5) \
	X(Nbeam_photons_5_6) \
	X(Nbeam_photons_6_7) \
	X(Nbeam_photons_7_8) \
	X(Nbeam_photons_8_9) \
	X(Nbeam_photons_9_10) \
	X(Nbeam_photons_10_11) \
	X(Nbeam_photons_11_12)

#define MyFloatBranchTypes(X) \
	X(Esc_tot) \
	X(Etof_tot) \
	X(Ebcal_points) \
	X(Ebcal_clusters) \
	X(Efcal_clusters) \
	X(Rfcal_max) \
	X(Rfcal_min) \
	\
	X(Ptot_candidates) \
	\
	X(Evisible_neutrals) \
	X(Evisible_tracks) \
	X(Evisible_charged_Kaons) \
	X(Evisible_charged_pions) \
	X(Evisible_protons) \
	X(Evisible)

class JEventProcessor_L3BDTtree:public jana::JEventProcessor{
	public:
		JEventProcessor_L3BDTtree();
		~JEventProcessor_L3BDTtree();
		const char* className(void){return "JEventProcessor_L3BDTtree";}
		
		class bdt_params_t{
			public:
				// Include number of objects for all JANA types
				#define Nobjs(A) Int_t N##A;
				MyTypes(Nobjs)

				// Include sorted int types
				#define Intobjs(A) Int_t A;
				MyIntBranchTypes(Intobjs)

				// Include sorted float types
				#define Floatobjs(A) Int_t A;
				MyFloatBranchTypes(Floatobjs)
		};
		
		TTree *l3tree;
		bdt_params_t bdt;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_L3BDTtree_

