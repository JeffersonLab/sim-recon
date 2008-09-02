// $Id$
//
//    File: DTrackFitter.h
// Created: Mon Sep  1 10:30:04 EDT 2008
// Creator: davidl (on Darwin Amelia.local 8.11.1 i386)
//

#ifndef _DTrackFitter_
#define _DTrackFitter_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

#include <PID/DKinematicData.h>
#include <HDGEOMETRY/DMagnetidFieldMap.h>


class DTrackFitter:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DTrackFitter);
		
		enum fit_type_t{
			kWireBased,
			kTimeBased
		};
		
		enum fit_status_t{
			kFitNotDone,
			kFitSuccess,
			kFitFailed
		};
		
		// Constructor and destructor
		DTrackFitter(JEventLoop *loop);	// require JEventLoop in constructor
		virtual ~DTrackFitter();
		
		// Hit accessor methods
		AddHit(const DCDCTrackHit* cdchit);
		AddHits(vector<const DCDCTrackHit*> cdchits);
		AddHit(const DFDCPseudo* fdchit);
		AddHits(vector<const DFDCPseudo*> fdchits);
		const vector<const DCDCTrackHit*>& GetCDCInputHits(void){return cdchits;}
		const vector<const DFDCPseudo*>&   GetFDCInputHits(void){return fdchits;}
		const vector<const DCDCTrackHit*>& GetCDCFitHits(void){return cdchits_used_in_fit;}
		const vector<const DFDCPseudo*>&   GetFDCFitHits(void){return fdchits_used_in_fit;}
		
		// Fit parameter accessor methods
		DKinematicData& GetInputParameters(void){return input_params;}
		DKinematicData& GetFitParameters(void){return fit_params;}
		fit_type_t GetFitType(void){return fit_type;}
		void SetFitType(fit_type_t type){fit_type=type;}
		
		// Wrappers
		fit_status_t FitTrack(const DVector3 &pos, const DVector3 &mom, double q);
		fit_status_t FitTrack(const DKinematicData &starting_params);
		
		//---- The following needs to be supplied by the subclass ----
		virtual fit_status_t FitTrack(void)=0;

	protected:

		// The following should be used as inputs by FitTrack(void)
		vector<const DCDCTrackHit*> cdchits;	//< Hits in the CDC
		vector<const DFDCPseudo*> fdchits;		//< Hits in the FDC
		DKinematicData input_params;				//< Starting parameters for the fit
		fit_type_t fit_type;							//< kWireBased or kTimeBased
		const DMagneticFieldMap *bfield;			//< Magnetic field map for current event (acquired through loop)
		JEventLoop *loop;								//< Pointer to JEventLoop object handling the current event

		// The following should be set as outputs by FitTrack(void)
		DKinematicData fit_params;									//< Results of last fit
		double chisq;													//< Chi-sq of final track fit (not the chisq/dof!)
		int dof;															//< Number of degrees of freedom for final fit parameters
		vector<const DCDCTrackHit*> cdchits_used_in_fit;	//< The CDC hits actually used in the fit
		vector<const DFDCPseudo*> fdchits_used_in_fit;		//< The FDC hits actually used in the fit

	private:
		
		// Prohibit default constructor
		DTrackFitter();
		
};

#endif // _DTrackFitter_

