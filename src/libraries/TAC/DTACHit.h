/*
 * DTACHit.h
 *
 *  Created on: Mar 24, 2017
 *      Author: hovanes
 */

#ifndef LIBRARIES_TAC_DTACHIT_H_
#define LIBRARIES_TAC_DTACHIT_H_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DTACHit: public jana::JObject {
protected:
	double E = 0;
	double T = 0;

	double integral = 0;
	double pulsePeak = 0;
	double timeTDC = 0;
	double timeFADC = 0;
	double npeFADC = 0;
	bool fadcPresent = false;
	bool tdcPresent = false;


public:
	JOBJECT_PUBLIC(DTACHit)

	DTACHit() {
	}
	virtual ~DTACHit() {
	}

	virtual void toStrings(vector<pair<string, string> > &items) const override {
		AddString(items, "E(MeV)", "%2.3f", E * 1000.0);
		AddString(items, "t(ns)", "%2.3f", T);
		AddString(items, "time_tdc(ns)", "%f", timeTDC);
		AddString(items, "time_fadc(ns)", "%f", timeFADC);
		AddString(items, "integral", "%f", integral);
		AddString(items, "pulse_peak", "%f", pulsePeak);
		AddString(items, "npe_fadc", "%f", npeFADC);
		AddString(items, "has_fADC", "%d", (int) fadcPresent);
		AddString(items, "has_TDC", "%d", (int) tdcPresent);
	}

	double getE() const {
		return E;
	}

	void setE(double e = 0) {
		E = e;
	}

	double getIntegral() const {
		return integral;
	}

	void setIntegral(double integral = 0) {
		this->integral = integral;
	}

	double getNpeFADC() const {
		return npeFADC;
	}

	void setNpeFADC(double npeFadc = 0) {
		npeFADC = npeFadc;
	}

	double getPulsePeak() const {
		return pulsePeak;
	}

	void setPulsePeak(double pulsePeak = 0) {
		this->pulsePeak = pulsePeak;
	}

	double getT() const {
		return T;
	}

	void setT(double t = 0) {
		T = t;
	}

	double getTimeFADC() const {
		return timeFADC;
	}

	void setTimeFADC(double timeFadc = 0) {
		timeFADC = timeFadc;
	}

	double getTimeTDC() const {
		return timeTDC;
	}

	void setTimeTDC(double timeTdc = 0) {
		timeTDC = timeTdc;
	}

	bool isFADCPresent() const {
		return fadcPresent;
	}

	void setFADCPresent(bool fadcPresent = false) {
		this->fadcPresent = fadcPresent;
	}

	bool isTDCPresent() const {
		return tdcPresent;
	}

	void setTDCPresent(bool tdcPresent = false) {
		this->tdcPresent = tdcPresent;
	}
};

#endif /* LIBRARIES_TAC_DTACHIT_H_ */
