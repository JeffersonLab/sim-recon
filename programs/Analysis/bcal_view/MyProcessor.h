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
#include "TRACKING/DMagneticFieldMap.h"
class DQuickFit;
class DTrackCandidate_factory;

#include <TPolyLine.h>
#include <TEllipse.h>
#include <TVector3.h>
#include <TMarker.h>
#include <TWbox.h>
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

		jerror_t DrawDetectors(void);

		int eventNo;
		vector<TObject*> graphics;
		
	private:
		void ClearEvent(void);
	
		int drew_detectors;
		
		vector<TWbox*> left_boxes;
		vector<TWbox*> center_boxes;
		vector<TWbox*> right_boxes;
};

#endif // _MYPROCESSOR_H_
