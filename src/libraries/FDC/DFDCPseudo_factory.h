//******************************************************************
// DFDCPseudo_factory.h: class definition for a factory creating
// pseudopoints from anode hits and cathode clusters.
// Author: Craig Bookwalter
// Date: Apr 2006
// Several revisions made by Simon Taylor, Fall 2006
//******************************************************************
#ifndef DFACTORY_DFDCPSEUDO_H
#define DFACTORY_DFDCPSEUDO_H

#include <JANA/JFactory.h>
#include <JANA/JObject.h>
#include <JANA/JException.h>
#include <JANA/JStreamLog.h>
using namespace jana;

#include "DFDCPseudo.h"
#include "DFDCCathodeCluster.h"
#include "DFDCHit.h"
#include "DFDCGeometry.h"
#include "HDGEOMETRY/DGeometry.h"
#include <TRACKING/DMCTrackHit.h>

#include <DMatrix.h>
#include <DMatrixSIMD.h>
#include <TDecompLU.h>
#include <TH2.h>
#include <TH1.h>

#include <algorithm>
#include <map>
#include <cmath>

typedef struct {
  double pos;
  double q;
  double q_from_pulse_height;
  int numstrips;
  double t; // mean time of strips in peak
  double t_rms; // rms of strips in peak
  unsigned int cluster; // index for cluster from which this centroid was generated
}centroid_t;

///
/// class DFDCPseudo_factory: definition for a JFactory that
/// produces pseudopoints from anode hits and DFDCCathodeClusters.
/// For now, it is purely geometry-based.
/// 
class DFDCPseudo_factory : public JFactory<DFDCPseudo> {
	public:
		
		///
		/// DFDCPseudo_factory::DFDCPseudo_factory():
		/// default constructor -- initializes log file
		///
		DFDCPseudo_factory();
		
		///
		/// DFDCPseudo_factory::~DFDCPseudo_factory():
		/// default destructor -- closes log file
		///
		~DFDCPseudo_factory();	
							

	protected:
		///
		/// DFDCPseudo_factory::evnt():
		/// this is the place that anode hits and DFDCCathodeClusters 
		/// are organized into pseudopoints.
		/// For now, this is done purely by geometry, with no drift
		/// information. See also
		/// DFDCPseudo_factory::makePseudo().
		///
		jerror_t init(void);
		jerror_t evnt(JEventLoop *eventLoop, uint64_t eventNo);
		jerror_t brun(JEventLoop *loop, int32_t runnumber);
		jerror_t erun(void);

		/// 
		/// DFDCPseudo_factory::makePseudo():
		/// performs UV+X matching to create pseudopoints
		///
		void makePseudo( vector<const DFDCHit*>& x,
				 vector<const DFDCCathodeCluster*>& u,
				 vector<const DFDCCathodeCluster*>& v,
				 int layer,
				 vector<const DMCTrackHit*> &mctrackhits);
		
		///
		/// DFDCPseudo_factory::CalcMeanTime()
		/// Calculates mean and RMS time for a cluster of cathode hits
		///
		void CalcMeanTime(const vector<const DFDCHit*>& H, double &t, double &t_rms);
		void CalcMeanTime(vector<const DFDCHit *>::const_iterator peak, double &t, double &t_rms);
		
		///
		/// DFDCPseudo_factory::FindCentroid()
		/// Calculates the centroids of groups of three adjacent strips
		/// containing a peak.
		///
		jerror_t FindCentroid(const vector<const DFDCHit*>& H, 
				      vector<const DFDCHit *>::const_iterator peak,
				      vector<centroid_t> &centroids);
		// Backtracking routine needed by FindCentroid 
		jerror_t FindNewParmVec(const DMatrix3x1 &N,const DMatrix3x1 &X,
					const DMatrix3x1 &F,const DMatrix3x3 &J,
					const DMatrix3x1 &par,
					DMatrix3x1 &newpar);
		
		///
		/// DFDCPseudo_factory::TwoStripCluster()
		/// Calculates the center-of-gravity of two adjacent strips
		///
		jerror_t TwoStripCluster(const vector<const DFDCHit*>& H,
					 vector<const DFDCHit *>::const_iterator peak,
					 vector<centroid_t> &centroids);
		
		///
		/// DFDCPseudo_factory::ThreeStripCluster()
		/// Calculates the center-of-gravity of Three adjacent strips
		///
		jerror_t ThreeStripCluster(const vector<const DFDCHit*>& H,
					 vector<const DFDCHit *>::const_iterator peak,
					 vector<centroid_t> &centroids);
		
 		
	private:		
		vector<vector<DFDCWire*> >fdcwires;
		vector<vector<DFDCCathode*> >fdccathodes;
		vector<double>xshifts;
		vector<double>yshifts;

		double ROUT_FIDUCIAL,RIN_FIDUCIAL;
		double r2_out,r2_in;
		double STRIP_ANODE_TIME_CUT;
		unsigned int MAX_ALLOWED_FDC_HITS;
//		bool DEBUG_HISTS,USE_FDC,MATCH_TRUTH_HITS;
		bool DEBUG_HISTS,USE_FDC;
		double MIDDLE_STRIP_THRESHOLD;
		double FDC_RES_PAR1,FDC_RES_PAR2;

		TH2F *qv_vs_qu, *dtv_vs_dtu;
		TH2F *uv_dt_vs_u,*uv_dt_vs_v,*v_wire_dt_vs_wire;
		TH2F *tv_vs_tu,*u_wire_dt_vs_wire;
		TH2F *Hxy[24],*ut_vs_u,*vt_vs_v;
		TH2F *v_vs_u,*dx_vs_dE;
		TH1F *u_cl_size, *v_cl_size, *u_cl_n, *v_cl_n, *x_dist_2, *x_dist_3, *x_dist_23, *x_dist_33;
		TH1F *d_uv;

//		JStreamLog* _log;
};

#endif // DFACTORY_DFDCPSEUDO_H

