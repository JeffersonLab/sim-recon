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
	this->setTimeScaleADC( this->getTimeScaleADC() * 64 );          // Time scale for a single FADC hit
	return NOERROR;
}

template<typename BASE_FACTORY, typename BUILDER>
jerror_t DRebuildFromRawFADC_factory<BASE_FACTORY, BUILDER>::brun(
		jana::JEventLoop *eventLoop, int32_t runnumber) {
	BASE_FACTORY::brun(eventLoop, runnumber);
	hitBuilder = new BUILDER(this->getTimeScaleADC(), this->getADCTimeOffset(), this->getTimeBaseADC());
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
	if (hitBuilder != nullptr)
		delete hitBuilder;
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

//template<typename BUILDER>
//jerror_t DTACHit_Rebuild_factory<BUILDER>::rebuildHit(
//		jana::JEventLoop* eventLoop, uint64_t eventNumber,
//		const DTACHit* baseHit) {
//
////	cout << "In DTACHit_Rebuild_factory::rebuildHit for " <<  baseHit << endl;
//
////	set<const JObject*> alreadyChecked;
////	set<const DTACDigiHit*> associatedDigiData;
////	int maxDepth = 1;
////	baseHit->GetAssociatedAncestors(alreadyChecked, maxDepth,
////			associatedDigiData);
////	cout << "Found " << associatedDigiData.size() << " digi ancestors" << endl;
////	for (auto& digiHit : associatedDigiData) {
////		cout << "There was a digi hit at " << hex << showbase << digiHit << dec
////				<< endl;
////	}
//
////	set<const Df250PulseData*> associatedPulseData;
////	alreadyChecked.clear();
////	maxDepth = 2;
////	baseHit->GetAssociatedAncestors(alreadyChecked, maxDepth,
////			associatedPulseData);
////	cout << "Found " << associatedPulseData.size() << " pulse ancestors"
////			<< endl;
////	for (auto& pulseHit : associatedPulseData) {
////		cout << "There was a pulse hit at " << hex << showbase << pulseHit
////				<< dec << endl;
////	}
//
//	auto rawData = getRawData(baseHit);
////	cout << "Raw data for " << baseHit << " is at " << rawData << endl;
//	// Check that no hit has already been made out of this raw data. This means that only
//	// one rebuild hit will be produced out of one window data. The fact that we sort the
//	// hits in the begining should guarantee us that the useful hit is copied. This can be
//	// reprogrammed differently at a later time.
//	if (rawDataPtrSet.find(rawData) == rawDataPtrSet.end()) {
//		auto& rawDataVector = rawData->samples;
//		//	auto rawDataVector = getSampleVector(baseHit);
//
//		double newTime = getTimeFromRawData(rawDataVector);
////		cout << "Received new time " << newTime << endl;
//
////		DTACHit *newHit = new DTACHit(*baseHit);
//		DTACHit* newHit = dynamic_cast<DTACHit*>(baseHit->Clone());
//
////		vector<pair<string,string>> outputStrings;
////
////		baseHit->toStrings( outputStrings );
////		cout << endl << "Old hit at " << baseHit  << endl;
////		for( auto& outputPair :  outputStrings ) {
////			cout << outputPair.first << "  = " << outputPair.second << endl;
////		}
////		cout << endl << "New hit at " << newHit << endl;
////		outputStrings.clear();
////		newHit->toStrings( outputStrings );
////		for( auto& outputPair :  outputStrings ) {
////			cout << outputPair.first << "  = " << outputPair.second << endl;
////		}
//
//		newHit->setTimeFADC(newTime);
//		newHit->setT(newTime);
//
////		cout << "Times are : old = " <<  baseHit->getT() << " new = " << newHit->getT() << endl;
//
////		newHit->AddAssociatedObject(rawData);
//		newHit->AddAssociatedObject(baseHit);
//
//		rawDataPtrSet.insert(rawData);
//
//		_data.push_back(newHit);
//	}
//
//	return NOERROR;
//}

//template<typename BUILDER>
//const Df250WindowRawData* DTACHit_Rebuild_factory<BUILDER>::getRawData(
//		const DTACHit* baseHit) {
////const vector<uint16_t>& DTACHit_Rebuild_factory::getSampleVector(
////		const DTACHit* baseHit) {
//	set<const Df250WindowRawData*> rawDataSet;
//	set<const JObject*> alreadyChecked;
//	int maxDepth = 3;
//	baseHit->GetAssociatedAncestors(alreadyChecked, maxDepth, rawDataSet);
////	cout << "Found " << associatedRawData.size() << " raw ancestors" << endl;
////	for (auto& rawHit : associatedRawData) {
////		cout << "There was a raw hit at " << hex << showbase << rawHit << dec
////				<< endl;
////	}
//	// Make sure there is only one raw data object, otherwise return zero pointer
//	if (rawDataSet.size() != 1) {
//		stringstream errMsg;
//		errMsg << "DTACHit at " << showbase << baseHit
//				<< " did not have Df250WindowRawData ancestors";
//		throw JException(errMsg.str());
//	}
//	return *rawDataSet.begin();
//}

//template<typename BUILDER>
//double DTACHit_Rebuild_factory<BUILDER>::getTimeFromRawData(
//		const vector<uint16_t>& samples) {
//	if (samples.size() < 1)
//		return 0;
//	// Find the maximum by going through the raw data and comparing samples
//	pair<double, double> maxInfo(0, 0);
//	auto maxElementIter = std::max_element(samples.begin(), samples.end());
//	double peakLocation = std::distance(samples.begin(), maxElementIter);
////	cout << "Maximum is at " << peakLocation << " bin" << endl;
//
//	// CAlculated the weighted sum of the peak bin and the neighboring bins. This will
//	// be called the "hit time" and hopefully will consistent for the "fast pulses".
//	double weightedSum = 0;
//	double normalizationFactor = 0;
//
//	auto minIter = maxElementIter - 1;
//	auto maxIter = maxElementIter + 2;
//	if (maxElementIter == samples.begin())
//		minIter = samples.begin();
//	if (maxElementIter == (samples.end() - 1))
//		maxIter = samples.end();
//
////	cout << "minIter is at " << distance( samples.begin(), minIter )
////			<< " , maxIter is at " << distance( samples.begin(), maxIter )
////			<< endl;
//	for (auto& elemIter = minIter; elemIter != maxIter; elemIter++) {
//		double elemNumber = std::distance(samples.begin(), elemIter);
//		weightedSum += elemNumber * (*elemIter);
//		normalizationFactor += (*elemIter);
//	}
////	cout << "weighted sum is " << weightedSum << " norm factor is " << normalizationFactor << endl;
//	// If the normalization factor is non-zero, return the weighted average location
//	// of the tree bins around the highest bin.
//	if (normalizationFactor > 0)
//		peakLocation = weightedSum / normalizationFactor;
//
////	cout << "New peak location is " << peakLocation << endl;
//
//	double pulseTime = timeScaleADC * peakLocation - adcTimeOffset
//			+ timeBaseADC;
//
////	cout << "Rebuild pulse time will be " << pulseTime << endl;
//
//	return pulseTime;
//}

#endif /* LIBRARIES_DREBUILDFROMRAWFADC_FACTORY_H_ */
