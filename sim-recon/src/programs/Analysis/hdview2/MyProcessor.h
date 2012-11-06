// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.h
//
/// Example program for a Hall-D analyzer which uses DANA
///

#ifndef _MYPROCESSOR_H_
#define _MYPROCESSOR_H_

#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
#include <JANA/JEvent.h>

#include <GlueX.h>

#include <HDGEOMETRY/DMagneticFieldMap.h>
#include <HDGEOMETRY/DRootGeom.h>
#include <HDGEOMETRY/DGeometry.h>
#include <PID/DNeutralParticle.h>
#include <PID/DKinematicData.h>
#include <DCoordinateSystem.h>
#include <TRACKING/DReferenceTrajectory.h>

class DQuickFit;
class DTrackCandidate_factory;
class DCDCTrackHit;

#include "hdv_mainframe.h"

class hdv_fulllistframe;
class hdv_debugerframe;

#include <TPolyLine.h>
#include <TEllipse.h>
#include <TVector3.h>
#include <TMarker.h>
#include <TFile.h>
#include <TH1.h>

#define MAX_HIT_MARKERS 2000
#define MAX_LINES 100
#define MAX_CIRCLES 100

class MyProcessor;
extern MyProcessor *gMYPROC;

class MyProcessor:public JEventProcessor
{
 public:
  MyProcessor();
  ~MyProcessor();
  
  jerror_t init(void);	///< Called once at program start.
  jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(JEventLoop *eventLoop, int eventnumber);///< Called every event.
  
  //void DrawXY(void);
  //void DrawRPhi(void);
  
  //jerror_t DrawHelicalTrack(DQuickFit *qf, int color);
  //jerror_t DrawStraightTrack(TVector3 p, TVector3 vertex, int color, int style);
  //jerror_t DrawTrack(double q, TVector3 pos, TVector3 mom, int color);
  //jerror_t DrawDetectors(void);
  
  const DMagneticFieldMap *Bfield;
  int eventNo;
  
  enum poly_type{
    kMarker =0,
    kLine   =1
  };
  
  class DGraphicSet{
  public:
  DGraphicSet(Color_t c, poly_type t, double s):color(c),type(t),size(s),marker_style(8){}
    vector<TVector3> points;
    Color_t color;
    poly_type type; // 0=markers, 1=lines
    double size;
    int marker_style;
  };
  vector<DGraphicSet> graphics;
  void FillGraphics(void);
  void UpdateTrackLabels(void);
  
  // Additional graphics that may be appropriate for only certain views
  vector<TObject*> graphics_xyA;
  vector<TObject*> graphics_xyB;
  vector<TObject*> graphics_xz;
  vector<TObject*> graphics_yz;
  
  void GetFactoryNames(vector<string> &facnames);
  void GetFactories(vector<JFactory_base*> &factories);
  unsigned int GetNrows(const string &factory, string tag);
  void GetDReferenceTrajectory(string dataname, string tag, 
			       unsigned int index, DReferenceTrajectory* &rt, vector<const DCDCTrackHit*> &cdchits);
  void GetAllWireHits(vector<pair<const DCoordinateSystem*,double> > &allhits);
  
 private:	
  
  hdv_mainframe *hdvmf;
  hdv_fulllistframe *fulllistmf;
  hdv_debugerframe *debugermf;
  JEventLoop *loop;
  JEvent last_jevent;
  DRootGeom *RootGeom;
  DGeometry *geom;
  string MATERIAL_MAP_MODEL;
  
  map<string, double> photon_track_matching;
  double DELTA_R_FCAL;
  
  void AddKinematicDataTrack(const DKinematicData* kd, int color, double size);
  void GetIntersectionWithCalorimeter(const DKinematicData* kd, DVector3 &pos, DetectorSystem_t &who);
  
  //DTrackCandidate_factory* factory;
  //void DrawTrackXY(const DKinematicData *, int color, float size);
};

extern MyProcessor* gMYPROC;

#endif // _MYPROCESSOR_H_
