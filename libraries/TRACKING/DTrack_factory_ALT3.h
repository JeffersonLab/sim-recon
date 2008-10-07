#ifndef _DTrack_factory_ALT3_
#define _DTrack_factory_ALT3_

#include <vector>

#include <DMatrix.h>
#include <TH2.h>
#include <TH3.h>
#include <DVector3.h>

#include <JANA/JFactory.h>
#include <JANA/JGeometry.h>
using namespace jana;

#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "HDGEOMETRY/DGeometry.h"
#include "DTrack.h"
#include "DReferenceTrajectory.h"
#include "DCoordinateSystem.h"
#include "FDC/DFDCSegment.h"

using namespace std;

class DTrackCandidate;
class DTrack;
class DCDCTrackHit;
class DFDCPseudo;

class DTrack_factory_ALT3:public JFactory<DTrack>{
 public:
  DTrack_factory_ALT3();
  ~DTrack_factory_ALT3();
  const char* Tag(void){return "ALT3";}

  jerror_t GetPositionAndMomentum(const DFDCSegment *segment,
				  DVector3 &pos, DVector3 &mom);

 private:

  double TOF_MASS;
  double MIN_FDC_HIT_PROB;
  unsigned int MIN_HITS;
  double endplate_z,endplate_dz,endplate_rmin,endplate_rmax;
  vector<double>cdc_half_length;

  jerror_t init(void);
  jerror_t brun(JEventLoop *loop, int runnumber);
  jerror_t evnt(JEventLoop *eventLoop, int eventnumber); 
  jerror_t fini(void);
  jerror_t erun(void);

  vector<const DTrackCandidate*>     trackcandidates;
  vector<const DCDCTrackHit* >       cdctrackhits;
  vector<const DFDCPseudo* >         fdctrackhits;
  vector<const DCDCTrackHit* >       cdchits_on_track;
  vector<const DFDCPseudo* >         fdchits_on_track;
  vector<vector<double> > cdcprobs;
  vector<vector<double> > fdcprobs;

  int eventnumber;
  const DMagneticFieldMap *bfield;
  const DGeometry *dgeom;

  bool DEBUG_HISTS;
  TH2F *cdc_residuals,*fdc_xresiduals,*fdc_yresiduals;
  TH2F *cdc_pulls_histo,*fdc_xpulls,*fdc_ypulls;

};

#endif // _DTrack_factory_ALT3_
