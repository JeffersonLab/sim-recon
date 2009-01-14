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
#include "FDC/DFDCSegment_factory.h"
#include <TH1F.h>
#include <TFile.h>
#include <TNtuple.h>
#include "DTrackLSFitter.h"

class MyProcessor:public JEventProcessor
{
 public:
  MyProcessor();
  ~MyProcessor();
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(JEventLoop *eventLoop, int runnumber);			///< Called everytime a new run number is detected.
  jerror_t evnt(JEventLoop *eventLoop, int eventnumber);		///< Called every event.
  int eventNo;
  jerror_t erun(void);
  jerror_t fini(void);
		
 private:
  
  TFile *rootfile;
  TNtuple *nFitter;

  DTrackLSFitter lsfitter;
};

#endif // _MYPROCESSOR_H_
