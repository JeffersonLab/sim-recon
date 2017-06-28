/*
 * DRebuildFromRawFADC_factory.h
 *
 *  Created on: May 26, 2017
 *      Author: Hovanes Egiyan
 */

#ifndef LIBRARIES_DREBUILDFROMRAWFADC_FACTORY_H_
#define LIBRARIES_DREBUILDFROMRAWFADC_FACTORY_H_

#include <vector>
#include <set>
using namespace std;

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
#include <DANA/ReadWriteLock.h>
#include <TTAB/DTTabUtilities.h>
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df250PulseData.h>

#include "DTACDigiHit.h"
#include "DTACTDCDigiHit.h"
#include "DTACHit.h"

#include <TAC/HitRebuilderTAC.h>

#include <TAC/DTACHit_factory.h>

// We will use this factory to create TAC hits with "REBUILD" tag using the existing
// hits and the raw data. This method will skip the making a new DigiHit
// class and simply replace some of the hit info by reanalyzing the FADC250
// waveform. this class is not supposed to be considered as "another" way
// of getting a TACHit since it depends on the default TACHit building
// factory for the main part of the information contained on the TACHit objects
// produced.

template<typename BASE_FACTORY = DTACHit_factory,
		typename BUILDER = HitRebuilderTAC>
class DRebuildFromRawFADC_factory: public BASE_FACTORY {
protected:
	BUILDER* hitBuilder = nullptr;
//	std::set<const Df250WindowRawData*> rawDataPtrSet;

//	virtual jerror_t rebuildHit(jana::JEventLoop *eventLoop,
//			uint64_t eventNumber, const DTACHit* baseHit);
//	virtual const Df250WindowRawData* getRawData(const DTACHit* baseHit);
//	virtual double getTimeFromRawData(const vector<uint16_t>& samples);
public:
	DRebuildFromRawFADC_factory() :
			BASE_FACTORY() {
	}
	virtual const char* Tag() {
		return (std::string(BASE_FACTORY::Tag()) + BUILDER::getTagString()).c_str();
	}
	virtual ~DRebuildFromRawFADC_factory() {
	}
	virtual jerror_t init(void) override;
	virtual jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber)
			override;          ///< Called everytime a new run number is detected
	virtual jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventNumber)
			override;          ///< Called every event
	virtual jerror_t erun(void) override;          ///< Called everytime run number changes, if brun has been called
	virtual jerror_t fini(void) override;          ///< Called after last event of last event source has been processed
};

using namespace std;

template<typename BASE_FACTORY, typename BUILDER>
jerror_t DRebuildFromRawFADC_factory<BASE_FACTORY, BUILDER>::init(void) {
	BASE_FACTORY::init();
//	this->setTimeScaleADC(this->getTimeScaleADC() * 64);          // Time scale for a single FADC hit
	return NOERROR;
}

template<typename BASE_FACTORY, typename BUILDER>
jerror_t DRebuildFromRawFADC_factory<BASE_FACTORY, BUILDER>::brun(
		jana::JEventLoop *eventLoop, int32_t runnumber) {
	BASE_FACTORY::brun(eventLoop, runnumber);
//	hitBuilder = new BUILDER(this->getTimeScaleADC(), this->getADCTimeOffset(),
//			this->getTimeBaseADC());
	hitBuilder = new BUILDER( eventLoop );

	return NOERROR;
}

// override the event processing method
template<typename BASE_FACTORY, typename BUILDER>
jerror_t DRebuildFromRawFADC_factory<BASE_FACTORY, BUILDER>::evnt(
		jana::JEventLoop* eventLoop, uint64_t eventNumber) {

//	cout << "In DTACHit_Rebuild_factory::evnt" << endl;
//
//	cout << "Data structure has " << _data.size() << " elements" << endl;
	vector<const DTACHit*> baseHitVector;
	eventLoop->Get(baseHitVector, "");

//	// Declare comparison functor for this TAC hits basedon pulse peak
//	static auto compareHits =
//			[](const DTACHit* lhs, const DTACHit* rhs ) ->
//			bool {return( (lhs!=nullptr)&(rhs!=nullptr ) ? fabs(lhs->getPulsePeak()) > fabs(rhs->getPulsePeak() ) : false );};
//
//	std::sort(baseHitVector.begin(), baseHitVector.end(), compareHits);
//
////	rawDataPtrSet = new set<const Df250WindowRawData*>();
//	for (auto baseHit : baseHitVector) {
//		rebuildHit(eventLoop, eventNumber, baseHit);
//		vector<DTACHit*> newHitVec = BUILDER( baseHit );
//	}
//	cout << "Here are all the " << _data.size() << " new hits in the rebuild factory" << endl;
//	for( auto hit : _data ) {
//		vector<pair<string,string>> outputStrings;
//
//		hit->toStrings( outputStrings );
//		cout << endl << "Rebuild hit at " << hit  << endl;
//		for( auto& outputPair :  outputStrings ) {
//			cout << outputPair.first << "  = " << outputPair.second << endl;
//		}
//	}
//	if( rawDataPtrSet != 0 ) {
//		delete rawDataPtrSet;
//	}

//	rawDataPtrSet.clear();

	vector<DTACHit*> newHitVector = (*hitBuilder)(baseHitVector);
	for (auto newHit : newHitVector) {
		this->AppendData(newHit);
//		_data.push_back(newHit);
	}

//	cout << "Done rebuilding TACHit " << endl;

	return NOERROR;
}

// erun
//------------------
template<typename BASE_FACTORY, typename BUILDER>
jerror_t DRebuildFromRawFADC_factory<BASE_FACTORY, BUILDER>::erun(void) {
	if (hitBuilder != nullptr) {
		delete hitBuilder;
		hitBuilder = nullptr;
	}
	BASE_FACTORY::erun();
	return NOERROR;
}

//------------------
// fini
//------------------
template<typename BASE_FACTORY, typename BUILDER>
jerror_t DRebuildFromRawFADC_factory<BASE_FACTORY, BUILDER>::fini(void) {
	BASE_FACTORY::fini();
	return NOERROR;
}


//// Declare specializations
//
//class WaveformErfcFunctor;
//template <typename T> class HitRebuilderByFit;
//
//
//template<>
//jerror_t DRebuildFromRawFADC_factory<DTACHit_factory,
//		HitRebuilderByFit<WaveformErfcFunctor>>::init(void);
//template<>
//jerror_t DRebuildFromRawFADC_factory<DTACHit_factory,
//		HitRebuilderByFit<WaveformErfcFunctor>>::brun(
//		jana::JEventLoop *eventLoop, int32_t runnumber);

#endif /* LIBRARIES_DREBUILDFROMRAWFADC_FACTORY_H_ */
