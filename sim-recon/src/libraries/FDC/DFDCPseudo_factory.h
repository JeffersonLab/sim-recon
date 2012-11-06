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
#include <TDecompLU.h>
#include <TH2.h>
#include <TH1.h>

#include <algorithm>
#include <map>
#include <cmath>

typedef struct {
  double pos;
  double q;
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
		jerror_t evnt(JEventLoop *eventLoop, int eventNo);
		jerror_t brun(JEventLoop *loop, int runnumber);

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
		jerror_t FindNewParmVec(DMatrix N,
						       DMatrix X,
						       DMatrix F,
						       DMatrix J,DMatrix par,
						       DMatrix &newpar);
 		
	private:		
		vector<vector<DFDCWire*> >fdcwires;

		double ROUT_FIDUCIAL,RIN_FIDUCIAL;
		double STRIP_ANODE_TIME_CUT;
		unsigned int MAX_ALLOWED_FDC_HITS;
		bool DEBUG_HISTS,USE_FDC;

		TH2F *qa_qc_diff;
		TH2F *qa_vs_qc, *dtv_vs_dtu;
		TH2F *uv_dt_vs_u,*uv_dt_vs_v,*v_wire_dt_vs_wire,*u_wire_dt_vs_wire;

		JStreamLog* _log;
};

#endif // DFACTORY_DFDCPSEUDO_H

