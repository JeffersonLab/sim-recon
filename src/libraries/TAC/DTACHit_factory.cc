/*
 * DTACHIT_factory.cc
 *
 *  Created on: Mar 24, 2017
 *      Author: hovanes
 */
/*
 * DTACHIT_factory.cc
 *
 *  Created on: Mar 24, 2017
 *      Author: Hovanes Egiyan
 */

#include "DTACHit_factory.h"

// The parameter below can be overwritten from the command line
// Define the static data members
bool DTACHit_factory::checkErrorsOnFADC = true;
// Flag to use time-walk corrections
bool DTACHit_factory::useTimeWalkCorrections = false;
// TDC time window in ns
double DTACHit_factory::timeWindowTDC = 50;
// Maximum difference between ADc and TDC time to form a single hit
double DTACHit_factory::timeDifferencInADCandTDC = 100;

// End JANA command line parameter values

// List of announced runs, starts empty
set<int> DTACHit_factory::announcedRuns;

//------------------
// init
//------------------
jerror_t DTACHit_factory::init(void) {
	// Set default configuration parameters

//	gPARMS->SetDefaultParameter("TAC:CHECK_FADC_ERRORS", checkErrorsOnFADC,
//			"Set to 1 to reject hits with fADC250 errors, ser to 0 to keep these hits");
//
//	gPARMS->SetDefaultParameter("TAC:HIT_TIME_WINDOW", timeWindowTDC,
//			"Time window of trigger corrected TDC time in which a hit in"
//					" in the TDC will match to a hit in the fADC to form an TAC hit");
//
//	gPARMS->SetDefaultParameter("TAC:DELTA_T_ADC_TDC_MAX",
//			timeDifferencInADCandTDC,
//			"Maximum difference in ns between a (calibrated) fADC time and"
//					" CAEN TDC time for them to be matched in a single hit");
//
//	gPARMS->SetDefaultParameter("TAC:USE_TIMEWALK_CORRECTION",
//			useTimeWalkCorrections,
//			"Flag to decide if time-walk corrections should be applied.");

	cout << "In DTACHit_factory::init" << endl;
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTACHit_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber) {

	// Only print messages for one thread whenever run number change
	bool printMessages = false;
	if (addRun(runnumber)) {
		printMessages = true;
	}
	/// Read in calibration constants
	if (printMessages)
		jout << "In DTACHit_factory, loading constants..." << endl;

//	// load scale factors
//	map<string, double> scaleFactors;
//	// a_scale (TAC_ADC_SCALE)
//	if (eventLoop->GetCalib("/TAC/digi_scales", scaleFactors))
//		jout << "Error loading /TAC/digi_scales !" << endl;
//	if (scaleFactors.find("TAC_ADC_ASCALE") != scaleFactors.end())
//		energyScale = scaleFactors["TAC_ADC_ASCALE"];
//	else
//		jerr << "Unable to get TAC_ADC_ASCALE from /TAC/digi_scales !" << endl;
//	// t_scale (TAC_ADC_SCALE)
//	if (scaleFactors.find("TAC_ADC_TSCALE") != scaleFactors.end())
//		timeScaleADC = scaleFactors["TAC_ADC_TSCALE"];
//	else
//		jerr << "Unable to get TAC_ADC_TSCALE from /TAC/digi_scales !" << endl;
//
//	// load base time offset
//	map<string, double> baseTimeOffsets;
//	// t_base (TAC_BASE_TIME_OFFSET)
//	if (eventLoop->GetCalib("/TAC/base_time_offset", baseTimeOffsets))
//		jout << "Error loading /TAC/base_time_offset !" << endl;
//	if (baseTimeOffsets.find("TAC_BASE_TIME_OFFSET") != baseTimeOffsets.end())
//		timeBaseADC = baseTimeOffsets["TAC_BASE_TIME_OFFSET"];
//	else
//		jerr
//				<< "Unable to get TAC_BASE_TIME_OFFSET from /TAC/base_time_offset !"
//				<< endl;
//	// t_tdc_base (TAC_TDC_BASE_TIME_OFFSET)
//	if (baseTimeOffsets.find("TAC_TDC_BASE_TIME_OFFSET")
//			!= baseTimeOffsets.end())
//		timeBaseTDC = baseTimeOffsets["TAC_TDC_BASE_TIME_OFFSET"];
//	else
//		jerr
//				<< "Unable to get TAC_BASE_TIME_OFFSET from /TAC/base_time_offset !"
//				<< endl;
//
//	// load constant tables
//	// a_gains (gains)
//	if (eventLoop->GetCalib("/TAC/gains", energyGain))
//		jout << "Error loading /TAC/gains !" << endl;
//	// a_pedestals (pedestals)
//	if (eventLoop->GetCalib("/TAC/pedestals", adcPedestal))
//		jout << "Error loading /TAC/pedestals !" << endl;
//	// adc_time_offsets (adc_timing_offsets)
//	if (eventLoop->GetCalib("/TAC/adc_timing_offsets", adcTimeOffset))
//		jout << "Error loading /TAC/adc_timing_offsets !" << endl;
//	// tdc_time_offsets (tdc_timing_offsets)
//	if (eventLoop->GetCalib("/TAC/tdc_timing_offsets", tdcTimeOffsets))
//		jout << "Error loading /TAC/tdc_timing_offsets !" << endl;
//	// timewalk_parameters (timewalk_parms)
//	if (eventLoop->GetCalib("TAC/timewalk_parms_v2", timeWalkParameters))
//		jout << "Error loading /TAC/timewalk_parms_v2 !" << endl;

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTACHit_factory::evnt(JEventLoop *eventLoop, uint64_t eventNumber) {
//	cout << "Building basic DTACHit objects" << endl;
	makeFADCHits(eventLoop, eventNumber);
	makeTDCHits(eventLoop, eventNumber);

	return NOERROR;
}

void DTACHit_factory::makeFADCHits(JEventLoop *eventLoop,
		uint64_t eventNumber) {
	const DTTabUtilities* locTTabUtilities = nullptr;
	eventLoop->GetSingle(locTTabUtilities);

	// Get TAC hits
	vector<const DTACDigiHit*> tacDigiHits;
	eventLoop->Get(tacDigiHits);

	for (auto& tacDigiHit : tacDigiHits) {
		// Throw away hits with firmware errors (post-summer 2016 firmware)
		if (errorCheckIsNeededForFADC()
				&& !locTTabUtilities->CheckFADC250_NoErrors(
						tacDigiHit->getQF()))
			continue;

		// Initialize pedestal to one found in CCDB, but override it
		// with one found in event if is available (?)
		// For now, only keep events with a correct pedestal
		double pedestal = adcPedestal;
		double nsamples_integral = tacDigiHit->getNsamplesIntegral();
		double nsamples_pedestal = tacDigiHit->getNsamplesPedestal();

		// nsamples_pedestal should always be positive for valid data - err on the side of caution for now
		if (nsamples_pedestal == 0) {
			jerr << "DSCDigiHit with nsamples_pedestal == 0 !   Event = "
					<< eventNumber << endl;
			continue;
		}

		// tacDigiHit->pedestal is the sum of "nsamples_pedestal" samples
		// Calculate the average pedestal per sample
		if ((tacDigiHit->getPedestal() > 0)
				&& locTTabUtilities->CheckFADC250_PedestalOK(tacDigiHit->QF)
				&& nsamples_pedestal > 0) {
			pedestal = (double) tacDigiHit->getPedestal() / nsamples_pedestal;
		}

		// Subtract pedestal from pulse peak
		if (tacDigiHit->getPulseTime() == 0 || tacDigiHit->getPedestal() == 0
				|| tacDigiHit->getPulsePeak() == 0)
			continue;
		double pulse_peak = tacDigiHit->getPulsePeak() - pedestal;

		// Subtract pedestal from pulse integral
		double pulseIntegral = (double) tacDigiHit->getPulseIntegral()
				- pedestal * nsamples_integral;

		double pulseTime = (double) tacDigiHit->getPulseTime();

		DTACHit *tacHit = new DTACHit;

		tacHit->setE(energyScale * energyGain * pulseIntegral);

		tacHit->setTimeFADC(
				timeScaleADC * pulseTime - adcTimeOffset + timeBaseADC);
		tacHit->setT( tacHit->getTimeFADC() );
		tacHit->setTimeTDC(numeric_limits<double>::quiet_NaN());

		tacHit->setTDCPresent(false);
		tacHit->setFADCPresent(true);

		tacHit->setTimeTDC(tacHit->getTimeFADC());          // set time from fADC in case no TDC tacHit
		tacHit->setPulsePeak(pulse_peak);

//		vector<pair<string,string>> outputStrings;
//		tacHit->toStrings( outputStrings );
//		cout << endl << "Newly build TAC hit " << endl;
//		for( auto& outputPair :  outputStrings ) {
//			cout << outputPair.first << "  = " << outputPair.second << endl;
//		}

		tacHit->AddAssociatedObject(tacDigiHit);

		_data.push_back(tacHit);
	}

	return;
}

void DTACHit_factory::makeTDCHits(JEventLoop *eventLoop, uint64_t eventnumber) {
	// Next, loop over TDC hits, matching them to the
	// existing fADC hits where possible and updating
	// their time information. If no match is found, then
	// create a new hit with just the TDC info.
	const DTTabUtilities* locTTabUtilities = nullptr;
	eventLoop->GetSingle(locTTabUtilities);

	vector<const DTACTDCDigiHit*> tacTDCDigiHits;
	eventLoop->Get(tacTDCDigiHits);

	for (auto& tacTDCDigiHit : tacTDCDigiHits) {

		// Get the the time and corrected it for the TDC base
		double tdcTime = locTTabUtilities->Convert_DigiTimeToNs_CAEN1290TDC(
				tacTDCDigiHit) - timeBaseTDC;

		// Look for existing hits to see if there is a match
		//   or create new one if there is no match
		// Require that the trigger corrected TDC time fall within
		//   a reasonable time window so that when a hit is associated with
		//   a hit in the TDC and not the ADC it is a "decent" TDC hit
		if (fabs(tdcTime) < timeWindowTDC) {
			DTACHit *tacHitCombined = findMatch(tdcTime);
			if (tacHitCombined == nullptr) {
				// If there was no matching hit create a new one without FADC info
				// push it into the list of hits for this factory.
				tacHitCombined = new DTACHit();
				tacHitCombined->setE(0.0);
				tacHitCombined->setTimeFADC(
						numeric_limits<double>::quiet_NaN());
				tacHitCombined->setFADCPresent(false);
				_data.push_back(tacHitCombined);
			}

			tacHitCombined->setTDCPresent(true);
			tacHitCombined->setT(tdcTime);

			// Check if the time-walk corrections are needed and adjust time time.
			if (useTimeWalkCorrections && (tacHitCombined->getE() > 0.0)) {
				// Correct for time walk
				// The correction is the form t=t_tdc- C1 (A^C2 - A0^C2)
				// double A  = tacHitCombined->dE;
				// double C1 = timewalk_parameters[id][1];
				// double C2 = timewalk_parameters[id][2];
				// double A0 = timewalk_parameters[id][3];
				// T -= C1*(pow(A,C2) - pow(A0,C2));

				// Correct for timewalk using pulse peak instead of pulse integral
				double A = tacHitCombined->getPulsePeak();
//				double C0 = timeWalkParameters[0];
				double C1 = timeWalkParameters[1];
				double C2 = timeWalkParameters[2];
				double A_THRESH = timeWalkParameters[3];
				double A0 = timeWalkParameters[4];
				// do correction
				tdcTime -= C1
						* (pow(A / A_THRESH, C2) - pow(A0 / A_THRESH, C2));
			}
			tacHitCombined->setT(tdcTime);

			tacHitCombined->AddAssociatedObject(tacTDCDigiHit);
		}
	}
	return;
}

//------------------
// Find a matching ADC hit and return pointer if found. It is assumed that no FADC hits have
// been inserted into the vector of the TAC hits, so looping through
// the vector will only take us through the ADC hits if anything at all.
// The ADC hit s are supposed to be taken care earlier but have NaN
// for TDC time information.
//------------------
DTACHit* DTACHit_factory::findMatch(double tdcTime) {
	DTACHit *bestMatchedHit = nullptr;

	// Loop over existing hits (from fADC) and look for a match
	// in both the sector and the time.
	for (auto& adcHit : _data) {
		// only match to fADC hits, not bachelor TDC hits
		if (!isfinite(adcHit->getTimeFADC()))
			continue;

		double deltaT = fabs(adcHit->getTimeFADC() - tdcTime);
		if (deltaT > timeDifferencInADCandTDC)
			continue;

		if (bestMatchedHit != nullptr) {
			if (deltaT < fabs(bestMatchedHit->getTimeFADC() - tdcTime))
				bestMatchedHit = adcHit;
		} else
			bestMatchedHit = adcHit;
	}
	return bestMatchedHit;
}

//------------------
// Reset_Data()
//------------------
void DTACHit_factory::Reset_Data(void) {
	// Clear _data vector
	_data.clear();
}

//------------------
// erun
//------------------
jerror_t DTACHit_factory::erun(void) {
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTACHit_factory::fini(void) {
	return NOERROR;
}

//std::string DTACHit_factory::SetTag(std::string tag) {
//
//	char* newString = new char[tag.size()+1];
//	std::copy( tag.begin(), tag.end(), newString );
//	newString[tag.size()] = '\0';
//	return string(tag_str);
//}
