// $Id$
//
//    File: DFactoryGenerator.cc
// Created: Mon Jul  3 21:46:40 EDT 2006
// Creator: davidl (on Darwin Harriet.local 8.6.0 powerpc)
//

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
extern jerror_t PID_init(JEventLoop *loop);

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
	PID_init(loop);
	
	return NOERROR;
}
