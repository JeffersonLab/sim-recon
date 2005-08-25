// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.h
//
/// Example program for a Hall-D analyzer which uses DANA
///

#ifndef _MYPROCESSOR_H_
#define _MYPROCESSOR_H_

#include "DEventProcessor.h"
#include "DMagneticFieldMap.h"
class DQuickFit;
class DFactory_DTrackCandidate;

#include <TPolyLine.h>
#include <TEllipse.h>
#include <TMarker.h>
#include <TFile.h>
#include <TH1.h>

#define MAX_HIT_MARKERS 2000
#define MAX_LINES 100
#define MAX_CIRCLES 100

class MyProcessor:public DEventProcessor
{
	public:
		MyProcessor();
		~MyProcessor();
	
		derror_t init(void);						///< Called once at program start.
		derror_t brun(DEventLoop *eventLoop, int runnumber);			///< Called everytime a new run number is detected.
		derror_t evnt(DEventLoop *eventLoop, int eventnumber);		///< Called every event.

		derror_t ConvertToTop(float x, float y, float z, float &X, float &Y);
		derror_t ConvertToSide(float x, float y, float z, float &X, float &Y);
		derror_t ConvertToFront(float x, float y, float z, float &X, float &Y);

		derror_t DrawHelicalTrack(DQuickFit *qf, int color);
		derror_t DrawTrack(DQuickFit *qf, int color);
		derror_t DrawDetectors(void);

		DMagneticFieldMap *Bfield;
		int eventNo;
		vector<TMarker*> markers;
		vector<TEllipse*> circles;
		vector<TPolyLine*> lines;
		vector<TObject*> graphics;
		
	private:
		void ClearEvent(void);
	
		int drew_detectors;
		DFactory_DTrackCandidate* factory;
		string TRACKHIT_SOURCE;
};

#endif // _MYPROCESSOR_H_
