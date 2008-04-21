#ifndef _DTrack_factory_ALT3_
#define _DTrack_factory_ALT3_

#include <vector>

#include <DVector3.h>
#include <DMatrix.h>
#include <JANA/JFactory.h>
#include <JANA/JGeometry.h>
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "DTrack.h"
#include "DReferenceTrajectory.h"
#include "DCoordinateSystem.h"

using namespace std;

class DTrackCandidate;
class DTrack;
class DCDCTrackHit;
class DFDCPseudo;

class DTrack_factory_ALT3:public JFactory<DTrack>{
 public:
  DTrack_factory_ALT3();
  ~DTrack_factory_ALT3();
  const string toString(void);
  const char* Tag(void){return "ALT3";}


 private:

  double TOF_MASS;
  double MIN_FDC_HIT_PROB;
  unsigned int MIN_HITS;

  jerror_t brun(JEventLoop *loop, int runnumber);
  jerror_t evnt(JEventLoop *eventLoop, int eventnumber); 


  vector<const DTrackCandidate*>     trackcandidates;
  vector<const DCDCTrackHit* >       cdctrackhits;
  vector<const DFDCPseudo* >         fdctrackhits;
  vector<const DCDCTrackHit* >       cdchits_on_track;
  vector<const DFDCPseudo* >         fdchits_on_track;
  vector<vector<double> > cdcprobs;
  vector<vector<double> > fdcprobs;

  int eventnumber;
  const DMagneticFieldMap *bfield;



};

#endif // _DTrack_factory_ALT3_
