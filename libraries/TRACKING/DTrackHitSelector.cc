// $Id$
//
//    File: DTrackHitSelector.cc
// Created: Thu Feb  5 13:34:58 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#include "DTrackHitSelector.h"


// The methods here are just convenience methods that delegate
// to the virtual methods that get implemented in the subclass.

// Note that we sort the hits in the same way as is done in 
// DTrack_factory just so it is easier to compare the lists
// of hits from that and other factories using this hit
// selector. These routines are defined in DTrack_factory.cc
extern bool CDCSortByRincreasing(const DCDCTrackHit* const &hit1, const DCDCTrackHit* const &hit2);
extern bool FDCSortByZincreasing(const DFDCPseudo* const &hit1, const DFDCPseudo* const &hit2);

//---------------------
// DTrackHitSelector  (Constructor)
//---------------------
DTrackHitSelector::DTrackHitSelector(JEventLoop *loop)
{
	this->loop = loop;
}

//---------------------
// GetCDCHits
//---------------------
void DTrackHitSelector::GetCDCHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DCDCTrackHit*> &cdchits_in, DTrackFitter *fitter) const
{
	/// Get all hits from the CDC and add them to the specified DTrackFitter object

	vector<const DCDCTrackHit*> cdchits_out;
	GetCDCHits(fit_type, rt, cdchits_in, cdchits_out);
	sort(cdchits_out.begin(), cdchits_out.end(), CDCSortByRincreasing);
	for(unsigned int i=0; i<cdchits_out.size(); i++)fitter->AddHit(cdchits_out[i]);
}

//---------------------
// GetFDCHits
//---------------------
void DTrackHitSelector::GetFDCHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DFDCPseudo*> &fdchits_in, DTrackFitter *fitter) const
{
	/// Get all hits from the FDC and add them to the specified DTrackFitter object

	vector<const DFDCPseudo*> fdchits_out;
	GetFDCHits(fit_type, rt, fdchits_in, fdchits_out);
	sort(fdchits_out.begin(), fdchits_out.end(), FDCSortByZincreasing);
	for(unsigned int i=0; i<fdchits_out.size(); i++)fitter->AddHit(fdchits_out[i]);
}

//---------------------
// GetAllHits
//---------------------
void DTrackHitSelector::GetAllHits(fit_type_t fit_type, DReferenceTrajectory *rt, const vector<const DCDCTrackHit*> &cdchits_in, const vector<const DFDCPseudo*> &fdchits_in, DTrackFitter *fitter) const
{
	/// Get all hits from both CDC and FDC and add them to the specified DTrackFitter object
	GetCDCHits(fit_type, rt, cdchits_in, fitter);
	GetFDCHits(fit_type, rt, fdchits_in, fitter);
}


