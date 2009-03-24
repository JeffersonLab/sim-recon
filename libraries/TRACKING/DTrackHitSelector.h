// $Id$
//
//    File: DTrackHitSelector.h
// Created: Thu Feb  5 13:34:58 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#ifndef _DTrackHitSelector_
#define _DTrackHitSelector_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

#include <TRACKING/DTrackFitter.h>



class DReferenceTrajectory;
class DCDCTrackHit;
class DFDCPseudo;

/// The DTrackHitSelector class is a base class for algorithms that 
/// will select hits from the drift chamber systems that are likely
/// to belong to a specified trajectory. This class doesn't actually
/// do the hit selection itself, it just provides a standard API so
/// multiple algorithms can be written. It is done this way since at
/// this point in time, we expect at least a couple of algorithms may
/// be tried.    Feb. 6, 2009  DL

class DTrackHitSelector:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DTrackHitSelector);

		DTrackHitSelector(JEventLoop *loop);
		DTrackHitSelector(){};

		enum fit_type_t{
			kWireBased = DTrackFitter::kWireBased, // ensure compatibility with DTrackFitter
			kTimeBased = DTrackFitter::kTimeBased, // ensure compatibility with DTrackFitter
			kHelical
		};
		
		virtual void GetCDCHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DCDCTrackHit*> &cdchits_in, vector<const DCDCTrackHit*> &cdchits_out) const =0;
		virtual void GetFDCHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DFDCPseudo*> &fdchits_in, vector<const DFDCPseudo*> &fdchits_out) const =0;

		void GetCDCHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DCDCTrackHit*> &cdchits_in, DTrackFitter *fitter) const;
		void GetFDCHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DFDCPseudo*> &fdchits_in, DTrackFitter *fitter) const;
		void GetAllHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DCDCTrackHit*> &cdchits_in, const vector<const DFDCPseudo*> &fdchits_in, DTrackFitter *fitter) const;

	protected:
	
		JEventLoop *loop;
};

#endif // _DTrackHitSelector_

