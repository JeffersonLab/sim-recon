/*
 * HitRebuilderTAC.h
 *
 *  Created on: Jun 7, 2017
 *      Author: Hovanes Egiyan
 */

#ifndef LIBRARIES_TAC_HITREBUILDERTAC_H_
#define LIBRARIES_TAC_HITREBUILDERTAC_H_

#include <vector>
#include <set>

#include <JANA/JFactory.h>
#include <DAQ/Df250WindowRawData.h>

#include "DTACDigiHit.h"
#include "DTACTDCDigiHit.h"
#include "DTACHit.h"

class HitRebuilderTAC {
protected:
	double timeScaleADC = 4.0;
	double adcTimeOffset = 0;
	double timeBaseADC = -130;

	std::set<const Df250WindowRawData*> rawDataPtrSet;

	static std::string tagString;

	virtual double getTimeFromRawData(const std::vector<uint16_t>& samples);
	virtual const Df250WindowRawData* getRawData(const DTACHit* baseHit);

	static std::string& setTagString( std::string tag ) { return (tagString = tag ); }
public:
	HitRebuilderTAC(double timeScale, double timeOffset, double timeBase) :
			timeScaleADC(timeScale), adcTimeOffset(timeOffset), timeBaseADC(
					timeBase) {
	}

	HitRebuilderTAC(const HitRebuilderTAC& f) :
			timeScaleADC(f.timeScaleADC), adcTimeOffset(f.adcTimeOffset), timeBaseADC(
					f.timeBaseADC) {
	}
	HitRebuilderTAC& operator=(const HitRebuilderTAC& f) {
		if (this != &f) {
			timeScaleADC = f.timeScaleADC;
			adcTimeOffset = f.adcTimeOffset;
			timeBaseADC = f.timeBaseADC;
		}
		return *this;
	}
	virtual ~HitRebuilderTAC() {
		// TODO Auto-generated destructor stub
	}
	virtual std::vector<DTACHit*> operator()(vector<const DTACHit*>& baseHitVector);
	static std::string getTagString() {return tagString;}

	double getADCTimeOffset() const {
		return adcTimeOffset;
	}

	void setADCTimeOffset(double adcTimeOffset = 0) {
		this->adcTimeOffset = adcTimeOffset;
	}

	const std::set<const Df250WindowRawData*>& getRawDataPtrSet() const {
		return rawDataPtrSet;
	}

	void setRawDataPtrSet(
			const std::set<const Df250WindowRawData*>& rawDataPtrSet) {
		this->rawDataPtrSet = rawDataPtrSet;
	}

	double getTimeBaseADC() const {
		return timeBaseADC;
	}

	void setTimeBaseADC(double timeBaseADC = -130) {
		this->timeBaseADC = timeBaseADC;
	}

	double getTimeScaleADC() const {
		return timeScaleADC;
	}

	void setTimeScaleADC(double timeScaleADC = 4.0) {
		this->timeScaleADC = timeScaleADC;
	}
};

#endif /* LIBRARIES_TAC_HITREBUILDERTAC_H_ */
