/*
 * DTACDigiHit.h
 *
 *  Created on: Mar 24, 2017
 *      Author: Hovanes Egiyan
 */

#ifndef LIBRARIES_TAC_DTACDIGIHIT_H_
#define LIBRARIES_TAC_DTACDIGIHIT_H_

#include <JANA/JObject.h>

class DTACDigiHit: public jana::JObject {
public:
	JOBJECT_PUBLIC(DTACDigiHit)

	uint32_t pulse_integral;          ///< identified pulse integral as returned by FPGA algorithm
	uint32_t pulse_time;          ///< identified pulse time as returned by FPGA algorithm
	uint32_t pedestal;          ///< pedestal info used by FPGA (if any)
	uint32_t QF;             ///< Quality Factor from FPGA algorithms
	uint32_t nsamples_integral;          ///< number of samples used in integral
	uint32_t nsamples_pedestal;          ///< number of samples used in pedestal
	uint32_t pulse_peak;           ///<  maximum sample in pulse

	uint32_t datasource;          ///<  0=window raw data, 1=old(pre-Fall16) firmware, 2=Df250PulseData

	// This method is used primarily for pretty printing
	// the second argument to AddString is printf style format
	void toStrings(vector<pair<string, string> > &items) const {
		AddString(items, "pulse_integral", "%d", pulse_integral);
		AddString(items, "pulse_peak", "%d", pulse_peak);
		AddString(items, "pulse_time", "%d", pulse_time);
		AddString(items, "pedestal", "%d", pedestal);
		AddString(items, "QF", "%d", QF);
		AddString(items, "nsamples_integral", "%d", nsamples_integral);
		AddString(items, "nsamples_pedestal", "%d", nsamples_pedestal);
	}


	uint32_t getDataSource() const {
		return datasource;
	}

	void setDataSource(uint32_t datasource) {
		this->datasource = datasource;
	}

	uint32_t getNsamplesIntegral() const {
		return nsamples_integral;
	}

	void setNsamplesIntegral(uint32_t nsamplesIntegral) {
		nsamples_integral = nsamplesIntegral;
	}

	uint32_t getNsamplesPedestal() const {
		return nsamples_pedestal;
	}

	void setNsamplesPedestal(uint32_t nsamplesPedestal) {
		nsamples_pedestal = nsamplesPedestal;
	}

	uint32_t getPedestal() const {
		return pedestal;
	}

	void setPedestal(uint32_t pedestal) {
		this->pedestal = pedestal;
	}

	uint32_t getPulseIntegral() const {
		return pulse_integral;
	}

	void setPulseIntegral(uint32_t pulseIntegral) {
		pulse_integral = pulseIntegral;
	}

	uint32_t getPulseTime() const {
		return pulse_time;
	}

	void setPulseTime(uint32_t pulseTime) {
		pulse_time = pulseTime;
	}

	uint32_t getQF() const {
		return QF;
	}

	void setQF(uint32_t qf) {
		QF = qf;
	}

	uint32_t getPulsePeak() const {
		return pulse_peak;
	}

	void setPulsePeak(uint32_t pulsePeak) {
		pulse_peak = pulsePeak;
	}
};

#endif /* LIBRARIES_TAC_DTACDIGIHIT_H_ */
