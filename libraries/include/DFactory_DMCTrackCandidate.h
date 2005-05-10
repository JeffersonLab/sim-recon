// $Id$
//
//    File: DFactory_DMCTrackCandidate.h
// Created: Sun Apr  3 12:38:16 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DMCTrackCandidate_
#define _DFactory_DMCTrackCandidate_

#include <TH1.h>
#include <TH2.h>
#include <TEllipse.h>
#include <TMarker.h>

#include "DFactory.h"
#include "DArcHit.h"
#include "DMCTrackCandidate.h"

class DQuickFit;


class DFactory_DMCTrackCandidate:public DFactory<DMCTrackCandidate>{

/// Determine track candidates from the MCCheatHits data.
///
/// This is designed to work only on the cheat codes from the
/// Monte Carlo right now. This is because they have x,y,z
/// information with more or less equal error bars. Eventually,
/// this may be evolved into a track finding class for regular
/// detector hits.
///
/// <b>The Algorithm</b>:
/// Tracks are identified using a histogramming method or 
/// Hough transform. Each hit is assumed to fall on a circle
/// in the x/y plane which passes through the origin (beamline).
/// The circle is not completely determined by only two points
/// (the hit and the origin) but the center of the circle is
/// constrained to fall on a line whose parameters can be
/// calculated directly from the hit's x,y coordinates. The
/// lines for all hits are smeared by a density function 
/// (in DArcHit::Density(x,y)) and entered into a 2 dimensional
/// histogram. Peaks in this histogram will emerge where several
/// lines intersect (i.e. the center of a circle). Hits whose
/// lines pass within some distance (DFactory_DMCTrackCandidate::masksize)
/// of the focal point are added to the list of hits for the
/// track candidate. Any hit can be considered part of more than
/// one track.
///
/// The second step is to look at the phi angle of a hit with
/// respect to the circle center as a function of the z position.
/// Hits from the same track will form a line when phi vs. z is
/// plotted. The sign of the slope of the line should match the
/// sign of the charge. The intercept will give the z position
/// of the vertex (assuming the vertex is on the beamline).
/// It can happen that track candidates from step one have hits
/// from more than one track. This would be the case when two
/// particles of equal transverse momentum, but opposite charge
/// go out back-to-back from the vertex. In this case, the two
/// will fall on the same circle in the x,y plane. They will,
/// however, form lines of different slopes on the phi z plane.
/// (note again that phi is in the coordinate system centered
/// on the circle, not on the beamline). By histogramming the
/// the slope for lines formed from all possible pairs of hits
/// in the track candidates, peaks will emerge indicating
/// the slopes of the lines formed by the points. Also, a single
/// peak should appear in a histogram of the intercepts indicating
/// the z-vertex position.

	public:
		DFactory_DMCTrackCandidate();
		~DFactory_DMCTrackCandidate();
		const string toString(void);
	
		derror_t FindCircles(void);
		derror_t FindCirclesHitSub(void);
		derror_t FindCirclesMaskSub(void);
		derror_t FindCirclesInt(void);
		int      FindTracks(float x, float y);
		derror_t ZeroNeighbors(TH2F *hist, int xbin, int ybin);
		int      IntersectionDensity(DArcHit *a, DArcHit *b, float &x, float&y);
		derror_t FillArcDensityHistogram(TH2F *hist);
		derror_t FillSlopeIntDensityHistos(float x0, float y0);
		derror_t DrawPhiZPoints(int which=0);
		
		inline vector<DArcHit*> GetDArcHits(void){return archits;}
		inline vector<TEllipse*> GetCircles(void){return circles;}
		inline vector<DQuickFit*> GetDQuickFits(void){return qfits;}
		derror_t SetMaxDensityHistograms(int N);
		inline int GetNumDensityHistograms(void){return density_histos.size();}
		inline int GetMaxDensityHistograms(void){return max_density_histos;}
		TH2F* GetDensityHistogram(int n);
		TH1F* GetSlopeDensityHistogram(int n);
		TH1F* GetOffsetDensityHistogram(int n);
		DQuickFit* GetQFit(int n);

		float circle_max;		///< max distance of focal point (circle center) from beamline
		float cm_per_bin;		///< cm per bin in one dimension for density histo
		float masksize;		///< max distance a line-of-circle-centers line can
									///< be from a focal point and still be considered
									///< on the circle
		float masksize2;
		int flip_x_axis;		///< change sign of x-coordinate for circles


	private:
		void ClearEvent(void);
		
		vector<DArcHit*> archits;	///< hit info for hits assumed to be on a circle
		vector<TEllipse*> circles;	///< use ellipses to remember circle centers
		vector<TMarker*> markers;	///< For debugging only
		
		vector<TH2F*> density_histos;	///< for debugging
		vector<TH1F*> slope_density_histos;
		vector<TH1F*> offset_density_histos;
		int max_density_histos;			///< maximum number of histos to allocate
		
		vector<DQuickFit*> qfits;
		
		derror_t ThereCanBeOnlyOne(int trk1, int trk2); ///< Aaaah!! the quickening!!!
		derror_t evnt(int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DMCTrackCandidate_

