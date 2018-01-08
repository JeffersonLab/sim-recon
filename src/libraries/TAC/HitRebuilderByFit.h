/*
 * HitRebuilderByFit.h
 * Class that rebuild the TAC time using 5-parameter fit.
 * the meaning of the parameters are fixed since this class
 * needs to provide the initial values for the parameters.
 * The template argument is the functor that will calculate the
 * value of the pulse distribution for a given set of five parameters.
 *  Created on: Jun 8, 2017
 *      Author: Hovanes Egiyan
 */

#ifndef LIBRARIES_TAC_HITREBUILDERBYFIT_H_
#define LIBRARIES_TAC_HITREBUILDERBYFIT_H_

#include <vector>
#include <set>
#include <string>
#include <functional>
#include <iostream>
#include <stdexcept>

#include <TH1D.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TCanvas.h>

#include <JANA/JFactory.h>
#include <DAQ/Df250WindowRawData.h>

#include "DTACDigiHit.h"
#include "DTACTDCDigiHit.h"
#include "DTACHit.h"

#include <TAC/HitRebuilderTAC.h>

template<typename F>
class HitRebuilderByFit: public HitRebuilderTAC {
protected:
	std::function<double(const double *data, const Double_t *param)> waveFunction =
			F();
	double adcPedestal = 100.0;
	double riseTime = 0.8;
	double decayTime = 2.6;
	TF1* fitFun = nullptr;
public:
	HitRebuilderByFit( jana::JEventLoop* eventLoop) :
			HitRebuilderTAC(eventLoop) {

		HitRebuilderByFit::readCCDB( eventLoop );
		// Create ROOT function for the fit using the functor
		std::string funName = this->getTagString() + ":fitFun";
		fitFun = new TF1(funName.c_str(), waveFunction, 5, 70, 5);
		fitFun->Print();
		std::cout << "ROOT function pointer is " << fitFun << std::endl;
	}

	HitRebuilderByFit(const HitRebuilderByFit& b) :
			HitRebuilderTAC(b), riseTime(b.riseTime), decayTime(b.decayTime) {
		fitFun = dynamic_cast<TF1*>(b.fitFun->Clone());
	}
	HitRebuilderByFit& operator=(const HitRebuilderByFit& b) {
		if (&b != this) {
			*(dynamic_cast<HitRebuilderTAC*>(this)) =
					*(dynamic_cast<HitRebuilderTAC*>(&b));
		}
		return *this;
	}
	virtual ~HitRebuilderByFit() {
		std::cout << "In ~HitRebuilderByFit" << std::endl;
		if (fitFun != nullptr)
			delete fitFun;
	}

	jerror_t readCCDB(jana::JEventLoop *eventLoop);

	virtual double getTimeFromRawData(const vector<uint16_t>& samples) override;

	static std::string getTagString() {
		return "REBUILD_" + F::getTagString();
	}

	double getDecayTime() const {
		return decayTime;
	}

	void setDecayTime(double decayTime = 2.6) {
		this->decayTime = decayTime;
	}

	double getRiseTime() const {
		return riseTime;
	}

	void setRiseTime(double riseTime = 0.8) {
		this->riseTime = riseTime;
	}

	const TF1* getFitFun() {
		return fitFun;
	}

	const std::function<double(const double* data, const Double_t* param)>& getWaveFunction() const {
		return waveFunction;
	}
};

template<typename F>
inline jerror_t HitRebuilderByFit<F>::readCCDB(jana::JEventLoop* eventLoop) {
	std::cout << "In HitRebuilderByFit::readCCDB() , reading calibration constants" << std::endl;

	// First re-read the constants for the base class. This is not a virtual method since it is
	// usually being called from the constructor.
	HitRebuilderTAC::readCCDB( eventLoop );

	// a_pedestals (pedestals)
	if (eventLoop->GetCalib("/TAC/pedestals", adcPedestal))
		jout << "Error loading /TAC/pedestals !" << std::endl;


	map<string, double> pulseShapeParameter;
	if (eventLoop->GetCalib("/TAC/pulse_shape", pulseShapeParameter))
		jout << "Error loading /TAC/pulse_shape !" << std::endl;

	// riseTime (riseTime)
	if (pulseShapeParameter.find("riseTime") != pulseShapeParameter.end())
		riseTime = pulseShapeParameter["riseTime"];
	else
		jerr << "Unable to get riseTime from /TAC/pulse_shape !" << std::endl;
	// decayTime (decayTime)
	if (pulseShapeParameter.find("decayTime") != pulseShapeParameter.end())
		decayTime = pulseShapeParameter["decayTime"];
	else
		jerr << "Unable to get decayTime from /TAC/pulse_shape !" << std::endl;


	std::cout << "riseTime is " << riseTime << " , decayTime is " << decayTime <<
			" , adcPedestal is " << adcPedestal << std::endl;

	return NOERROR;
}

template<typename F>
double HitRebuilderByFit<F>::getTimeFromRawData(
		const vector<uint16_t>& samples) {

// Create the histogram for fitting
	string histoName = this->getTagString() + ":sample";
	string histoTitle = string(
			" FADC samples for TAC for builder " + getTagString());
	TH1D sampleHisto(histoName.c_str(), histoTitle.c_str(), samples.size(), 0.0,
			double(samples.size()));
	int iBin = 1;
	for (auto& sample : samples) {
		sampleHisto.SetBinContent(iBin, double(sample));
		sampleHisto.SetBinError(iBin, 1.0);
		iBin++;
	}
	// Find the maximum by going through the raw data and comparing samples
	pair<double, double> maxInfo(0, 0);
	auto maxElementIter = std::max_element(samples.begin(), samples.end());
	double peakLocation = std::distance(samples.begin(), maxElementIter);

//	double pedestal = 101.0;
	double amplitude = *maxElementIter - adcPedestal;
	fitFun->SetParNames("Pedestal", "PeakLocation", "riseTime", "Amplitude",
			"Lambda");

//	std::cout << "Pedestal is " << adcPedestal << " peak is at " << peakLocation
//			<< " amplitude is at " << amplitude << " lambda is " << lambda
//			<< std::endl;
	fitFun->SetParameters(adcPedestal, peakLocation, amplitude, amplitude,
			decayTime);

	//	std::cout << "Set parameters" << std::endl;
	double xMin = fitFun->GetXmin();
	double xMax = peakLocation + 2 * decayTime;
//	std::cout << "Function peak is at " << fitFun->Eval(peakLocation)
//			<< std::cout;
//
//	std::cout << "Will try to fit from " << xMin << " to " << xMax
//			<< " using TF1 at " << fitFun << std::endl;

	fitFun->FixParameter(0, adcPedestal);
	fitFun->FixParameter(2, riseTime);
	fitFun->FixParameter(4, decayTime);
	TFitResultPtr fitResultPtr = sampleHisto.Fit(fitFun, "WWQ0S", "", xMin,
			xMax);
	double newPeak = fitFun->GetParameter(1);

	double pulseTime = timeScaleADC * newPeak - adcTimeOffset + timeBaseADC;

	if (!fitResultPtr->IsValid()) {
		std::string errMsg = "Fit failed in " + this->getTagString();
		throw std::runtime_error(errMsg);
		//		std::cout << "Failed fit " << std::endl;
//		std::cout << "	Convergence status is " << fitResultPtr->Status() << std::endl;
////		TCanvas c();
////		sampleHisto.Draw();
////		usleep(2000000);
	}

//	int fitConvergenceStatus = fitResultPtr.Get()->Status();
//
//	if( fitConvergenceStatus < 2 )
//		std::cout << "Convergence status is " << fitConvergenceStatus << std::endl;
//
//	std::cout << "Convergence status is " << fitConvergenceStatus << std::endl;
	return pulseTime;
}

#endif /* LIBRARIES_TAC_HITREBUILDERBYFIT_H_ */
