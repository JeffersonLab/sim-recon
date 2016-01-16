// $Id$
//
//    File: DTrackCandidate_factory_CDCCOSMIC.h
// Created: Sat Jun 28 16:50:07 EDT 2014
// Creator: davidl (on Darwin harriet.local 13.2.0 i386)
//

#ifndef _DTrackCandidate_factory_CDCCOSMIC_
#define _DTrackCandidate_factory_CDCCOSMIC_

#include <TH2.h>
#include <TMinuit.h>
#include <JANA/JFactory.h>
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <HDGEOMETRY/DMagneticFieldMapNoField.h>

class DTrackCandidate_factory_CDCCOSMIC:public jana::JFactory<DTrackCandidate>{
	public:
		DTrackCandidate_factory_CDCCOSMIC():rt(NULL),bfield(NULL),residual_vs_ring(NULL){};
		~DTrackCandidate_factory_CDCCOSMIC(){};
		const char* Tag(void){return "CDCCOSMIC";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

        double CDCDriftDistance(double delta, double t);
        double CDCDriftVariance(double t);
        unsigned int Locate(vector<double>&xx,double x);
        double CDCTrackError(const DCDCWire *, const double *, double *);
        void GetDOCAPhiandZ(const DCDCWire *, DVector3, DVector3, double &, double &);
        DReferenceTrajectory *rt;
        DMagneticFieldMapNoField *bfield;
        DTrackFinder *finder;
        //TMinuit *ptMinuit;
        vector<const DCDCTrackHit *> hits;
        vector<double> cdc_drift_table;
        double cdc_drift_table_min, cdc_drift_table_max;
        // Resolution parameters
        double CDC_RES_PAR1,CDC_RES_PAR2;
        int EXCLUDERING;
        vector<vector<double> >max_sag;
        vector<vector<double> >sag_phi_offset;
        double long_drift_func[3][3];
        double short_drift_func[3][3];

        TH2D *residual_vs_ring;
        TH1D *h_chisq;
        TH1D *h_Ndof;
        TH1D *h_chisq_per_Ndof;

        void CalcChisq(DTrackCandidate *can, vector<const DCDCTrackHit*> &axial_hits, vector<const DCDCTrackHit*> &stereo_hits);
};

#endif // _DTrackCandidate_factory_CDCCOSMIC_

