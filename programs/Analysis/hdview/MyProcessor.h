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

#include <TVector3.h>
#include <TFile.h>
#include <TH1.h>
#include <TEllipse.h>


class MyProcessor:public DEventProcessor
{
	public:
		derror_t init(void){};					///< Called once at program start.
		derror_t brun(int runnumber){};		///< Called everytime a new run number is detected.
		derror_t evnt(int eventnumber);		///< Called every event.
		derror_t erun(void){};					///< Called everytime run number changes, provided brun has been called.
		derror_t fini(void){};					///< Called after last event of last event source has been processed.

		int eventNo;
		TVector3 cdchits[1000];
		int cdchit_tracks[1000];
		int Ncdchits;
		
		TEllipse *ellipse[10];
		int Nellipse;
		
		MyProcessor(void){
			Ncdchits=0;
			Nellipse=0;
		}
};

#endif // _MYPROCESSOR_H_
