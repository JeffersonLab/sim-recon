// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include "MyProcessor.h"

//------------------------------------------------------------------
// init   -Open output file here (e.g. a ROOT file)
//------------------------------------------------------------------
derror_t MyProcessor::init(void)
{
	cout<<__FILE__<<":"<<__LINE__<<endl;
	return NOERROR;
}

//------------------------------------------------------------------
// brun   -Read in calibration constants here
//------------------------------------------------------------------
derror_t MyProcessor::brun(int runnumber)
{
	cout<<__FILE__<<":"<<__LINE__<<endl;
	return NOERROR;
}

//------------------------------------------------------------------
// evnt   -Fill histograms here
//------------------------------------------------------------------
derror_t MyProcessor::evnt(int eventnumber)
{
	//cout<<__FILE__<<":"<<__LINE__<<endl;
	event_loop->PrintRate();

	return NOERROR;
}

//------------------------------------------------------------------
// erun   -Update run-indexed info you want to keep here
//------------------------------------------------------------------
derror_t MyProcessor::erun(void)
{
	cout<<__FILE__<<":"<<__LINE__<<endl;
	return NOERROR;
}

//------------------------------------------------------------------
// fini   -Close output file here
//------------------------------------------------------------------
derror_t MyProcessor::fini(void)
{
	cout<<__FILE__<<":"<<__LINE__<<endl;
	return NOERROR;
}

