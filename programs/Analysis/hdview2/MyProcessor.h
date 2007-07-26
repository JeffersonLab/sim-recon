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
#include "HDGEOMETRY/DMagneticFieldMap.h"
class DQuickFit;
class DTrackCandidate_factory;

#include "hdv_mainframe.h"

#include <TPolyLine.h>
#include <TEllipse.h>
#include <TVector3.h>
#include <TMarker.h>
#include <TFile.h>
#include <TH1.h>

#define MAX_HIT_MARKERS 2000
#define MAX_LINES 100
#define MAX_CIRCLES 100

class MyProcessor:public JEventProcessor
{
	public:
		MyProcessor();
		~MyProcessor();
	
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber);			///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);		///< Called every event.

		jerror_t ConvertToTop(float x, float y, float z, float &X, float &Y);
		jerror_t ConvertToSide(float x, float y, float z, float &X, float &Y);
		jerror_t ConvertToFront(float x, float y, float z, float &X, float &Y);

		jerror_t DrawHelicalTrack(DQuickFit *qf, int color);
		jerror_t DrawStraightTrack(TVector3 p, TVector3 vertex, int color, int style);
		jerror_t DrawTrack(double q, TVector3 pos, TVector3 mom, int color);
		jerror_t DrawDetectors(void);

		const DMagneticFieldMap *Bfield;
		int eventNo;
		vector<TMarker*> markers;
		vector<TEllipse*> circles;
		vector<TPolyLine*> lines;
		vector<TObject*> graphics;
		
	private:
		void ClearEvent(void);
	
		hdv_mainframe *hdvmf;
		int drew_detectors;
		DTrackCandidate_factory* factory;
		string TRACKHIT_SOURCE;
};

extern MyProcessor* gMYPROC;

#endif // _MYPROCESSOR_H_
