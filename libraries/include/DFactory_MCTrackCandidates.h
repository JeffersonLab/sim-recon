// $Id$



#ifndef _DMCTRACKCANDIDATES_H_
#define _DMCTRACKCANDIDATES_H_

#include <TH1.h>
#include <TH2.h>
#include <TEllipse.h>
#include <TMarker.h>

#include "DFactory.h"
#include "DArcHit.h"


typedef struct{
	int Nhits;
	int ihit[100];		///< index of hits in MCCheatHits factory
	float x0,y0;		///< center of circle
	float z_vertex;	///< z coordinate of vertex
	float dphidz;		///< dphi/dz in radians per cm
}MCTrackCandidate_t;

class DFactory_MCTrackCandidates:public DFactory{

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
/// lines pass within some distance (DFactory_MCTrackCandidates::masksize)
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
/// (note again that phi means in the coordinate system centered
/// on the circle, not on the beamline). By histogramming the
/// the slope for lines formed from all possible pairs of hits
/// in the track candidates, peaks will emerge indicating
/// the slopes of the lines formed by the points. Also, a single
/// peak should appear in a histogram of the intercepts indicating
/// the z-vertex position.

	public:
		DFactory_MCTrackCandidates(DEvent *event);
		~DFactory_MCTrackCandidates();
		derror_t Print(void);
	
		derror_t FindCircles(void);
		derror_t FindCirclesHitSub(void);
		derror_t FindCirclesMaskSub(void);
		derror_t ZeroNeighbors(TH2F *hist, int xbin, int ybin);
		derror_t FillArcDensityHistogram(TH2F *hist);
		derror_t FillSlopeIntDensityHistos(void);
		derror_t DrawPhiZPoints(void);
		
		inline DArcHit* Getarchits(void){return archit;}
		inline int GetNarchits(void){return Narchits;}
		inline TEllipse* GetCircles(void){return circles;}
		inline int GetNcircles(void){return Ncircles;}
		derror_t SetNumDensityHistograms(int N);
		inline int GetNumDensityHistograms(void){return Ndensity_histos;}
		TH2F* GetDensityHistogram(int n);
		TH1F* GetSlopeDensityHistogram(int n);
		TH1F* GetOffsetDensityHistogram(int n);

		float circle_max;		///< max distance of focal point (circle center) from beamline
		float bins_per_cm;	///< bins per cm in one dimension for density histo
		float masksize;		///< max distance a line-of-circle-centers line can
									///< be from a focal point and still be considered
									///< on the circle
		int flip_x_axis;		///< change sign of x-coordinate for circles

		TH1F *slope_density, *offset_density;

	private:
		DArcHit archit[300];
		int Narchits;
		TEllipse circles[32];	///< use ellipses to remember circle centers
		int Ncircles;
		TMarker *markers;			///< For debugging only
		int Nmarkers;				///< For debugging only
		
		TH2F *density;					///< 2-D density histogram for finding points on a circle
		TH2F *density_histos[8];	///< for debugging
		int Ndensity_histos;			///< number of valid pointers in density_histos[]
		TH1F *slope_density_histos[8];
		TH1F *offset_density_histos[8];

		derror_t evnt(int eventnumber);	///< Called every event.
		
};

#endif //_DTRACKCANDIDATES_H_

