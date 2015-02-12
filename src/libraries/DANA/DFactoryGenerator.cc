// $Id$
//
//    File: DFactoryGenerator.cc
// Created: Mon Jul  3 21:46:40 EDT 2006
// Creator: davidl (on Darwin Harriet.local 8.6.0 powerpc)
// Modified:	Oct 3, 2012, Yi Qiang: add CERE
//				Oct 11 2012, Yi Qiang: add RICH

#include "DFactoryGenerator.h"

extern jerror_t BCAL_init(JEventLoop *loop);
extern jerror_t CDC_init(JEventLoop *loop);
extern jerror_t FDC_init(JEventLoop *loop);
extern jerror_t FCAL_init(JEventLoop *loop);
extern jerror_t CCAL_init(JEventLoop *loop);
extern jerror_t START_COUNTER_init(JEventLoop *loop);
extern jerror_t TAGGER_init(JEventLoop *loop);
extern jerror_t TOF_init(JEventLoop *loop);
extern jerror_t TRACKING_init(JEventLoop *loop);
extern jerror_t HDDM_init(JEventLoop *loop);
extern jerror_t PID_init(JEventLoop *loop);
extern jerror_t ANALYSIS_init(JEventLoop *loop);
extern jerror_t DAQ_init(JEventLoop *loop);
extern jerror_t TTAB_init(JEventLoop *loop);
extern jerror_t CERE_init(JEventLoop *loop);
extern jerror_t RICH_init(JEventLoop *loop);
extern jerror_t TRIGGER_init(JEventLoop *loop);
extern jerror_t PAIR_SPECTROMETER_init(JEventLoop *loop);

//---------------------------------
// DFactoryGenerator    (Constructor)
//---------------------------------
DFactoryGenerator::DFactoryGenerator()
{

}

//---------------------------------
// ~DFactoryGenerator    (Destructor)
//---------------------------------
DFactoryGenerator::~DFactoryGenerator()
{

}

//---------------------------------
// GenerateFactories
//---------------------------------
jerror_t DFactoryGenerator::GenerateFactories(JEventLoop *loop)
{
	BCAL_init(loop);
	CDC_init(loop);
	FDC_init(loop);
	FCAL_init(loop);
	CCAL_init(loop);
	START_COUNTER_init(loop);
	TAGGER_init(loop);
	TOF_init(loop);
	TRACKING_init(loop);
	HDDM_init(loop);
	PID_init(loop);
	ANALYSIS_init(loop);
	DAQ_init(loop);
	TTAB_init(loop);
	CERE_init(loop);
	RICH_init(loop);
	TRIGGER_init(loop);
	PAIR_SPECTROMETER_init(loop);
	
	return NOERROR;
}
