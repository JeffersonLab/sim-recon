/*
 * WaveformErfcFunctor.h
 *
 *  Created on: Jun 9, 2017
 *      Author: Hovanes Egiyan
 */

#ifndef LIBRARIES_TAC_WAVEFORMERFCFUNCTOR_H_
#define LIBRARIES_TAC_WAVEFORMERFCFUNCTOR_H_

#include <math.h>

#include <iostream>
#include <string>
#include <functional>

#include <TMath.h>

class WaveformErfcFunctor {
public:
	WaveformErfcFunctor() {
	}
	virtual ~WaveformErfcFunctor() {
	}
	static std::string getTagString() {
		return "ERFC";
	}
	double operator()(const double* xArray, const double* pArray) {
		static const double zeroTolerance = 1.0e-25;
		double x = xArray[0];
		double pedestal = pArray[0];
		double peakPoint = pArray[1];
		double riseTime = pArray[2];

		double expAmpl = pArray[3];
		double expLambda = pArray[4];

		double sigma = riseTime / 2.0 ;
//		if (lineSlope > zeroTolerance)
//			sigma = expAmpl / lineSlope;
		double mu = 0;
		if (expLambda > zeroTolerance)
			mu = 1. / expLambda;

		double value = 0;

		if (sigma > zeroTolerance) {
			double variablePart = expAmpl * exp(-mu * (x - peakPoint))
					* TMath::Erfc(
							(mu * sigma * sigma - (x - peakPoint))
									/ (sqrt(2) * sigma));
			if (variablePart < 0)
				variablePart = 0;
			value = pedestal + variablePart;
		}
		return value;
	}
};

#endif /* LIBRARIES_TAC_WAVEFORMERFCFUNCTOR_H_ */
