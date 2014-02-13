// $Id$
//
//    File: DEventProcessor_bcalfcaltof_res_tree.h
// Created: Mon Apr  3 11:38:03 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#ifndef _DEventProcessor_bcalfcaltof_res_tree_
#define _DEventProcessor_bcalfcaltof_res_tree_

#include <JANA/JEventProcessor.h>
using namespace jana;

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TTree.h>
#include <DVector3.h>
#include <particleType.h>

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <BCAL/DBCALShower.h>
#include <FCAL/DFCALShower.h>
#include <TOF/DTOFPoint.h>
#include <TOF/DTOFTruth.h>
#include <TOF/DTOFHit.h>
#include <HDGEOMETRY/DRootGeom.h>
#include <BCALMCComparison.h>
#include <FCALMCComparison.h>
#include <TOFMCComparison.h>
#include <TRACKING/DTrackFitter.h>

class DEventProcessor_bcalfcaltof_res_tree:public JEventProcessor{
	public:
		DEventProcessor_bcalfcaltof_res_tree(){};
		~DEventProcessor_bcalfcaltof_res_tree(){};
		const char* className(void){return "DEventProcessor_bcalfcaltof_res_tree";}

		void Convert_Coordinates_BCALToLab(float locBCALR, float locBCALPhi, float locBCALZ, DVector3& locLabVertex);
		void Convert_Coordinates_LabToBCAL(const DVector3& locLabVertex, float& locBCALR, float& locBCALPhi, float& locBCALZ);
		double Calc_MostProbableTOFdE(const DVector3 &locMomentum, double mass);
		void Calc_MostProbableTOFdEdx(double p, double mass, double dx, double &dE, double &dEdx);

		float Calc_BCALPathLengthCorrection(float locEnergy);
		float Calc_BCALPathLengthCorrectionZPostE(float locZ);
		float Calc_FCALPathLengthCorrection(float locEnergy);

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		const DRootGeom *dRootGeom;                                 
		bool dBCALStudyFlag;
		bool dFCALStudyFlag;
		bool dTOFStudyFlag;
		bool locSCStudyFlag;
		DTrackFitter *dTrackFitter;

		double dRhoZoverA;
		double dKRhoZoverA;
		double dLnI;

		BCALMCComparison *dBCALMCComparison;
		FCALMCComparison *dFCALMCComparison;
		TOFMCComparison *dTOFMCComparison;

		TTree* dPluginTree_BCALMCComparison;
		TTree* dPluginTree_FCALMCComparison;
		TTree* dPluginTree_TOFMCComparison;
};

#endif // _DEventProcessor_bcalfcaltof_res_tree_

