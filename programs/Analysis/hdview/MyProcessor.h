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

#include <TMarker.h>
#include <TFile.h>
#include <TH1.h>


class MyProcessor:public DEventProcessor
{
	public:
		MyProcessor();
		~MyProcessor(){};
	
		derror_t init(void){};					///< Called once at program start.
		derror_t brun(int runnumber){};		///< Called everytime a new run number is detected.
		derror_t evnt(int eventnumber);		///< Called every event.
		derror_t erun(void){};					///< Called everytime run number changes, provided brun has been called.
		derror_t fini(void){};					///< Called after last event of last event source has been processed.

		derror_t ConvertToTop(float x, float y, float z, float &X, float &Y);
		derror_t ConvertToSide(float x, float y, float z, float &X, float &Y);
		derror_t ConvertToFront(float x, float y, float z, float &X, float &Y);

		derror_t DrawDetectors(void);

		int eventNo;
		TMarker *hitMarkers[1000];
		int NhitMarkers;
	
	private:
		int drew_detectors;
};

#endif // _MYPROCESSOR_H_
