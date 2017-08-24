/*
 * WaveformSpikeFunctor.h
 *
 *  Created on: Jun 8, 2017
 *      Author: Hovanes Egiyan
 */

#ifndef LIBRARIES_TAC_WAVEFORMSPIKEFUNCTOR_H_
#define LIBRARIES_TAC_WAVEFORMSPIKEFUNCTOR_H_

#include <math.h>

#include <iostream>
#include <string>
#include <functional>

class WaveformSpikeFunctor {
public:
	WaveformSpikeFunctor() {
	}
	virtual ~WaveformSpikeFunctor() {
	}
	static std::string getTagString() {
		return "SPIKE";
	}
	double operator()(const double* xArray, const double* pArray) {
		static const double zeroTolerance = 1.0e-25;

		double x = xArray[0];
		double pedestal = pArray[0];
		double peakPoint = pArray[1];
		double riseTime =  pArray[2];
//		double lineSlope = pArray[2];

		double expAmpl = pArray[3];
		double expLambda = pArray[4];

		double lineSlope = 0;
		if( riseTime > zeroTolerance ) {
			lineSlope = expAmpl  / riseTime;
		}

		double lineIntersect = expAmpl - lineSlope * peakPoint;

		double value = 0;

		if (x > peakPoint) {
			if (fabs(expLambda) > 1.0e-25) {
				value = expAmpl * exp(-(x - peakPoint) / expLambda) + pedestal;
			} else {
				value = pedestal;
			}
		} else {
			value = lineSlope * x + lineIntersect + pedestal;
			if (value < pedestal)
				value = pedestal;
		}
		return value;
	}
};

#endif /* LIBRARIES_TAC_WAVEFORMSPIKEFUNCTOR_H_ */
