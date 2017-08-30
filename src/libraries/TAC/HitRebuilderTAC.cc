/*
 * HitRebuilderTAC.cc
 *
 *  Created on: Jun 7, 2017
 *      Author: Hovanes Egiyan
 */

#include <cmath>
using namespace std;

#include <TAC/HitRebuilderTAC.h>

using namespace std;

jerror_t HitRebuilderTAC::readCCDB( JEventLoop* eventLoop) {
	cout << "In HitRebuilderTAC::readCCDB() , reading calibration constants" << endl;
	// load scale factors
	map<string, double> scaleFactors;
	if (eventLoop->GetCalib("/TAC/digi_scales", scaleFactors))
		jout << "Error loading /TAC/digi_scales !" << endl;
	// t_scale (TAC_ADC_SCALE)
	if (scaleFactors.find("TAC_ADC_TSCALE") != scaleFactors.end())
		timeScaleADC = adcTimeRescaleFactor * scaleFactors["TAC_ADC_TSCALE"];
	else
		jerr << "Unable to get TAC_ADC_TSCALE from /TAC/digi_scales !" << endl;

	// load base time offset
	map<string, double> baseTimeOffsets;
	// t_base (TAC_BASE_TIME_OFFSET)
	if (eventLoop->GetCalib("/TAC/base_time_offset", baseTimeOffsets))
		jout << "Error loading /TAC/base_time_offset !" << endl;
	if (baseTimeOffsets.find("TAC_BASE_TIME_OFFSET") != baseTimeOffsets.end())
		timeBaseADC = baseTimeOffsets["TAC_BASE_TIME_OFFSET"];
	else
		jerr
				<< "Unable to get TAC_BASE_TIME_OFFSET from /TAC/base_time_offset !"
				<< endl;

	// load constant tables
	// adc_time_offsets (adc_timing_offsets)
	if (eventLoop->GetCalib("/TAC/adc_timing_offsets", adcTimeOffset))
		jout << "Error loading /TAC/adc_timing_offsets !" << endl;

	cout << "timeScaleADC is " << timeScaleADC << " , timeBaseADC is " << timeBaseADC <<
			" ,  adcTimeOffset is " << adcTimeOffset << endl;

	return NOERROR;
}


const Df250WindowRawData* HitRebuilderTAC::getRawData(const DTACHit* baseHit) {
//const vector<uint16_t>& DTACHit_Rebuild_factory::getSampleVector(
//		const DTACHit* baseHit) {
	set<const Df250WindowRawData*> rawDataSet;
	set<const JObject*> alreadyChecked;
	int maxDepth = 3;
	baseHit->GetAssociatedAncestors(alreadyChecked, maxDepth, rawDataSet);
//	cout << "Found " << associatedRawData.size() << " raw ancestors" << endl;
//	for (auto& rawHit : associatedRawData) {
//		cout << "There was a raw hit at " << hex << showbase << rawHit << dec
//				<< endl;
//	}
	// Make sure there is only one raw data object, otherwise return zero pointer
	if (rawDataSet.size() != 1) {
		stringstream errMsg;
		errMsg << "DTACHit at " << showbase << baseHit
				<< " did not have Df250WindowRawData ancestors";
		throw JException(errMsg.str());
	}
	return *rawDataSet.begin();
}

double HitRebuilderTAC::getTimeFromRawData(const vector<uint16_t>& samples) {
	if (samples.size() < 1)
		return 0;
	// Find the maximum by going through the raw data and comparing samples
	pair<double, double> maxInfo(0, 0);
	auto maxElementIter = std::max_element(samples.begin(), samples.end());
	double peakLocation = std::distance(samples.begin(), maxElementIter);
	//	cout << "Maximum is at " << peakLocation << " bin" << endl;

	// CAlculated the weighted sum of the peak bin and the neighboring bins. This will
	// be called the "hit time" and hopefully will consistent for the "fast pulses".
	double weightedSum = 0;
	double normalizationFactor = 0;

	auto minIter = maxElementIter - 1;
	auto maxIter = maxElementIter + 3;
	if (maxElementIter == samples.begin())
		minIter = samples.begin();
	if (maxElementIter >= (samples.end() - 2))
		maxIter = samples.end();

	//	cout << "minIter is at " << distance( samples.begin(), minIter )
	//			<< " , maxIter is at " << distance( samples.begin(), maxIter )
	//			<< endl;
	for (auto& elemIter = minIter; elemIter != maxIter; elemIter++) {
		double elemNumber = std::distance(samples.begin(), elemIter);
		weightedSum += (elemNumber + 0.5) * (*elemIter);
		normalizationFactor += (*elemIter);
	}
	//	cout << "weighted sum is " << weightedSum << " norm factor is " << normalizationFactor << endl;
	// If the normalization factor is non-zero, return the weighted average location
	// of the tree bins around the highest bin.
	if (normalizationFactor > 0)
		peakLocation = weightedSum / normalizationFactor;

	//	cout << "New peak location is " << peakLocation << endl;

	double pulseTime = timeScaleADC * peakLocation - adcTimeOffset
			+ timeBaseADC;

	//	cout << "Rebuild pulse time will be " << pulseTime << endl;

	return pulseTime;
}

vector<DTACHit*> HitRebuilderTAC::operator()(
		vector<const DTACHit*>& baseHitVector) {
	// Declare comparison functor for this TAC hits basedon pulse peak
	static auto compareHits =
			[](const DTACHit* lhs, const DTACHit* rhs ) ->
			bool {return( (lhs!=nullptr)&(rhs!=nullptr ) ? fabs(lhs->getPulsePeak()) > fabs(rhs->getPulsePeak() ) : false );};

	// Sort the vector of the hits base on the peak height, with the heighest peak being in the first hit.
	std::sort(baseHitVector.begin(), baseHitVector.end(), compareHits);

	rawDataPtrSet.clear();
	vector<DTACHit*> newHitVector;
	for (auto baseHit : baseHitVector) {
		auto rawData = getRawData(baseHit);
		//	cout << "Raw data for " << baseHit << " is at " << rawData << endl;
		// Check that no hit has already been made out of this raw data. This means that only
		// one rebuild hit will be produced out of one window data. The fact that we sort the
		// hits in the begining should guarantee us that the useful hit is copied. This can be
		// reprogrammed differently at a later time.
		double newTime = baseHit->getT();
		if (rawDataPtrSet.find(rawData) == rawDataPtrSet.end()) {
			auto& rawDataVector = rawData->samples;
			//	auto rawDataVector = getSampleVector(baseHit);

			try {
				newTime = getTimeFromRawData(rawDataVector);
			} catch (const runtime_error& e) {
				newTime = baseHit->getT();
				newTime = -10;
			}
//		cout << "Received new time " << newTime << endl;

//		DTACHit *newHit = new DTACHit(*baseHit);
			DTACHit* newHit = dynamic_cast<DTACHit*>(baseHit->Clone());

//		vector<pair<string,string>> outputStrings;
//
//		baseHit->toStrings( outputStrings );
//		cout << endl << "Old hit at " << baseHit  << endl;
//		for( auto& outputPair :  outputStrings ) {
//			cout << outputPair.first << "  = " << outputPair.second << endl;
//		}
//		cout << endl << "New hit at " << newHit << endl;
//		outputStrings.clear();
//		newHit->toStrings( outputStrings );
//		for( auto& outputPair :  outputStrings ) {
//			cout << outputPair.first << "  = " << outputPair.second << endl;
//		}

			newHit->setTimeFADC(newTime);
			newHit->setT(newTime);

//		cout << "Times are : old = " <<  baseHit->getT() << " new = " << newHit->getT() << endl;

//		newHit->AddAssociatedObject(rawData);
			newHit->AddAssociatedObject(baseHit);

			rawDataPtrSet.insert(rawData);

			newHitVector.push_back(newHit);
		}
	}
	return newHitVector;
}

