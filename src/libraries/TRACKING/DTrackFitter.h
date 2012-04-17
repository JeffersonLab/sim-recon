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

#include <prof_time.h>
#include <DANA/DApplication.h>
#include <PID/DKinematicData.h>
#include <HDGEOMETRY/DMagneticFieldMap.h>
#include "HDGEOMETRY/DLorentzMapCalibDB.h"
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>
#include <TRACKING/DReferenceTrajectory.h>
using namespace std;

#define NaN std::numeric_limits<double>::quiet_NaN()


class DReferenceTrajectory;
class DGeometry;

//////////////////////////////////////////////////////////////////////////////////
/// The DTrackFitter class is a base class for different charged track
/// fitting algorithms. It does not actually fit the track itself, but
/// provides the interface and some common support features most algorthims
/// will need to implement.
///
/// The reason this is needed (instead of just using the mechanism already
/// built into JANA) is due to the nature of charged track fitting.
/// Specifically, tracks are usually fit first to the wire positions
/// and then to the drift times. The algorithm for both is (at least
/// usually) the same. However, we want to separate the wire-based and
/// time-based fitting into 2 distinct stages allowing easy access to the 
/// wire-based fits.
///
/// There were a few options on how to handle this within the JANA framework
/// but it was decided passing DTrackFitter objects through the framework
/// was the best way to address it. Sub-classes of DTrackFitter will
/// implement the actual algorithms, but JANA will only see these
/// objects as pointers to the DTrackFitter base class. Only one
/// DTrackFitterXXX object will exist for each thread (i.e. each JEventLoop).
/// As such, the state of that object will likely be overwritten
/// many times in a single event and it's internal data never used
/// by anything outside of the TRACKING package. Also, the factories that
/// produce the DTrackFitterXXX objects will make them as persistent
/// and will turn off the the WRITE_TO_OUTPUT bit by default.
//////////////////////////////////////////////////////////////////////////////////

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
			kFitFailed,
			kFitNoImprovement
		};
		
		class dedx_t{
		public:
		  dedx_t(double dE,double dx, double p):dE(dE),dx(dx),p(p){}
		    double dE; // energy loss in layer
		    double dx; // path length in layer
		    double p;  // momentum at this dE/dx measurement

		};

		class pull_t{
		public:
		  pull_t(double resi, double err, double s=0.0):resi(resi),err(err),s(s){}
		    double resi;	// residual of measurement
		    double err;		// estimated error of measurement
			 double s;
		};
		
		// Constructor and destructor
		DTrackFitter(JEventLoop *loop);	// require JEventLoop in constructor
		virtual ~DTrackFitter();
		
		void Reset(void);
		
		// Hit accessor methods
		void AddHit(const DCDCTrackHit* cdchit);
		void AddHits(vector<const DCDCTrackHit*> cdchits);
		void AddHit(const DFDCPseudo* fdchit);
		void AddHits(vector<const DFDCPseudo*> fdchits);
		const vector<const DCDCTrackHit*>& GetCDCInputHits(void) const {return cdchits;}
		const vector<const DFDCPseudo*>&   GetFDCInputHits(void) const {return fdchits;}
		const vector<const DCDCTrackHit*>& GetCDCFitHits(void) const {return cdchits_used_in_fit;}
		const vector<const DFDCPseudo*>&   GetFDCFitHits(void) const {return fdchits_used_in_fit;}
		
		// Fit parameter accessor methods
		const DKinematicData& GetInputParameters(void) const {return input_params;}
		const DKinematicData& GetFitParameters(void) const {return fit_params;}
		double GetChisq(void) const {return chisq;}
		int GetNdof(void) const {return Ndof;}
		vector<pull_t>& GetPulls(void){return pulls;}
		fit_type_t GetFitType(void) const {return fit_type;}
		const DMagneticFieldMap* GetDMagneticFieldMap(void) const {return bfield;}

		void SetFitType(fit_type_t type){fit_type=type;}
		void SetInputParameters(const DKinematicData &starting_params){input_params=starting_params;}
		
		// Wrappers
		fit_status_t FitTrack(const DVector3 &pos, const DVector3 &mom, double q, double mass,double t0=NaN,DetectorSystem_t t0_det=SYS_NULL);
		fit_status_t FitTrack(const DKinematicData &starting_params);
		
		// Methods that actually do something
		fit_status_t 
		  FindHitsAndFitTrack(const DKinematicData &starting_params, 
				      DReferenceTrajectory *rt, 
				      JEventLoop *loop, double mass=-1.0,
				      double t0=NaN,
				      DetectorSystem_t t0_det=SYS_NULL); ///< mass<0 means get it from starting_params
		jerror_t CorrectForELoss(const DKinematicData &starting_params, DReferenceTrajectory *rt, DVector3 &pos, DVector3 &mom, double mass);
		double CalcDensityEffect(double p,double mass,double density,
					 double Z_over_A,double I);  
		double CalcDensityEffect(double p,double mass,
					 double rho_Z_over_A,double LnI);
		double CalcDensityEffect(double betagamma,
					 double rho_Z_over_A,double LnI);
		
		void GetProfilingTimes(std::map<std::string, prof_time::time_diffs> &my_prof_times) const;
		
		//---- The following need to be supplied by the subclass ----
		virtual string Name(void) const =0;
		virtual fit_status_t FitTrack(void)=0;
		virtual double ChiSq(fit_type_t fit_type, DReferenceTrajectory *rt, double *chisq_ptr=NULL, int *dof_ptr=NULL, vector<pull_t> *pulls_ptr=NULL)=0;

	protected:
	

		// The following should be used as inputs by FitTrack(void)
		vector<const DCDCTrackHit*> cdchits;	//< Hits in the CDC
		vector<const DFDCPseudo*> fdchits;		//< Hits in the FDC
		DKinematicData input_params;				//< Starting parameters for the fit
		fit_type_t fit_type;							//< kWireBased or kTimeBased
		const DMagneticFieldMap *bfield;			//< Magnetic field map for current event (acquired through loop)
		const DLorentzDeflections *lorentz_def;//< Correction to FDC cathodes due to Lorentz force
		const DGeometry *geom;						//< DGeometry pointer used to access materials through calibDB maps for eloss
		const DRootGeom *RootGeom;					//< ROOT geometry used for accessing material for MULS, energy loss
		JEventLoop *loop;								//< Pointer to JEventLoop object handling the current event

		// The following should be set as outputs by FitTrack(void)
		DKinematicData fit_params;									//< Results of last fit
		double chisq;													//< Chi-sq of final track fit (not the chisq/dof!)
		int Ndof;														//< Number of degrees of freedom for final track
		vector<pull_t> pulls;										//< pull_t objects for each contribution to chisq (assuming no correlations)
		fit_status_t fit_status;									//< Status of values in fit_params (kFitSuccess, kFitFailed, ...)
		vector<const DCDCTrackHit*> cdchits_used_in_fit;	//< The CDC hits actually used in the fit
		vector<const DFDCPseudo*> fdchits_used_in_fit;		//< The FDC hits actually used in the fit


		bool CORRECT_FOR_ELOSS;

	private:
		int DEBUG_LEVEL;			

		// Prohibit default constructor
		DTrackFitter();
		
		// gas material properties
		double mKRhoZoverAGas,mRhoZoverAGas,mLnIGas;
		
};

#endif // _DTrackFitter_

