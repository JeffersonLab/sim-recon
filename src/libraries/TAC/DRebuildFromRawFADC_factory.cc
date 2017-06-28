/*
 * DRebuildFromRawFADC_factory.cc
 *
 *  Created on: May 26, 2017
 *      Author: Hovanes Egiyan
 */

#include <JANA/JEventLoop.h>

#include <TAC/DRebuildFromRawFADC_factory.h>
#include <TAC/DTACHit_factory.h>
#include <TAC/HitRebuilderByFit.h>
#include <TAC/WaveformErfcFunctor.h>

using namespace std;

template<>
jerror_t DRebuildFromRawFADC_factory<DTACHit_factory,
		HitRebuilderByFit<WaveformErfcFunctor>>::init(void) {
	DTACHit_factory::init();
	this->setTimeScaleADC(this->getTimeScaleADC() * 64);          // Time scale for a single FADC hit

	return NOERROR;
}

//template<>
//jerror_t DRebuildFromRawFADC_factory<DTACHit_factory,
//		HitRebuilderByFit<WaveformErfcFunctor>>::brun(
//		jana::JEventLoop *eventLoop, int32_t runnumber) {
//	DTACHit_factory::brun(eventLoop, runnumber);
//
////	// load scale factors
////	map<string, double> shapeParameters;
////
////	double riseTime = 0.8;
////	double decayTime = 2.6;
////
////	if (eventLoop->GetCalib("/TAC/pulse_shape", shapeParameters))
////		jout << "Error loading /TAC/pulse_shape !" << endl;
////
////	if (shapeParameters.find("riseTime") != shapeParameters.end())
////		riseTime = shapeParameters["riseTime"];
////	else
////		jerr << "Unable to get riseTime from /TAC/pulse_shape !" << endl;
////
////	if (shapeParameters.find("decayTime") != shapeParameters.end())
////		decayTime = shapeParameters["decayTime"];
////	else
////		jerr << "Unable to get decayTime from /TAC/pulse_shape !" << endl;
//
////	hitBuilder = new HitRebuilderByFit<WaveformErfcFunctor>(
////			this->getTimeScaleADC(), this->getADCTimeOffset(),
////			this->getTimeBaseADC(), riseTime, decayTime);
//	hitBuilder = new HitRebuilderByFit<WaveformErfcFunctor>(eventLoop);
//
//	return NOERROR;
//}

