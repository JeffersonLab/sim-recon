// $Id$
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
using namespace std;

#include <math.h>
#include "hddm_s.h"

float RANDOM_MAX = (float)(0x7FFFFFFF);

void SmearCDC(s_HDDM_t *hddm_s);
void AddNoiseHitsCDC(s_HDDM_t *hddm_s);
void SmearFDC(s_HDDM_t *hddm_s);
void AddNoiseHitsFDC(s_HDDM_t *hddm_s);
void SmearFCAL(s_HDDM_t *hddm_s);
void SmearBCAL(s_HDDM_t *hddm_s);
void SmearTOF(s_HDDM_t *hddm_s);
void SmearUPV(s_HDDM_t *hddm_s);
void SmearCherenkov(s_HDDM_t *hddm_s);

double SampleGaussian(double sigma);
double SampleRange(double x1, double x2);

// Do we or do we not add noise hits
bool ADD_NOISE = false;

// The straw diameter is 1.6cm
float CDC_R = 0.8;			// cm

// The stereo pitch is 6 degrees. For hit-based tracks, the
// hit could come anywhere in the 1.6cm diameter straw. Thus,
// the z-range corresponding to the straw diameter is
// 1.6cm * cot(6 degrees) = 15cm. 
float CDC_Z_SIGMA = 7.5;	// cm

// The occupancy of the CDC is used to determine how many noise
// hits to add per event. There are 3240 wires in the CDC. The
// Number of noise hits is sampled from a gaussian with a sigma
// equal to the sqrt(Nwires * CDCOccupancy).
float CDC_AVG_NOISE_HITS = 0.03*3240.0; // 0.03 = 3% occupancy

// The wire spacing in the FDC is 0.5cm. In the actual data,
// multiple planes will be combined into psuedo points since
// each plane only measures in one direction (not counting z).
// At any rate, this should probably be an over estimate.
float FDC_R = 0.25;			// cm

// The z-coordinate of the wire and cathode planes are well
// known, but the hit will have some error on it. This should
// come primarily from the fact that the reconstructed x/y
// coordinates are at the DOCA point which will be off-plane.
// The actual error will come from the error on the track angle
// at the anode plane and the DOCA value measured by drift time.
// For here, we'll just use sigma=1mm.
float FDC_Z = 0.1;			// cm

// The occupancy of the FDC is used to determine how many noise
// hits to add per event. There are 2856 anode wires in the FDC.
// Each anode plane will have two cathode planes that will be
// used to resolve ambiguities and create "hits" that are passed
// on to the tracking code. We base the noise level on the
// number of anode wires.
float FDC_AVG_NOISE_HITS = 0.01*2856.0; // 0.01 = 1% occupancy

//-----------
// Smear
//-----------
void Smear(s_HDDM_t *hddm_s)
{
	SmearCDC(hddm_s);
	if(ADD_NOISE)AddNoiseHitsCDC(hddm_s);
	SmearFDC(hddm_s);
	if(ADD_NOISE)AddNoiseHitsFDC(hddm_s);
	SmearFCAL(hddm_s);
	SmearBCAL(hddm_s);
	SmearTOF(hddm_s);
	SmearUPV(hddm_s);
	SmearCherenkov(hddm_s);
}

//-----------
// SmearCDC
//-----------
void SmearCDC(s_HDDM_t *hddm_s)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->centralDC == HDDM_NULL ||
			hits->centralDC->cdcTruthPoints == HDDM_NULL)continue;
		
		s_CdcTruthPoints_t *cdctruthpoints = hits->centralDC->cdcTruthPoints;
		s_CdcTruthPoint_t *cdctruthpoint = cdctruthpoints->in;
		for(unsigned int j=0; j<cdctruthpoints->mult; j++, cdctruthpoint++){

			// Here we want to move the point to a position
			// randomly sampled from a disc in the x/y
			// plane centered on the actual point.
			float r = cdctruthpoint->r;
			float phi = cdctruthpoint->phi;
			float z = cdctruthpoint->z;
			float x = r*cos(phi);
			float y = r*sin(phi);
								
			// since this needs to be evenly distributed over
			// the disc of radius CDC_R, we may try several times
			float deltaX, deltaY;
			do{
				float sx = (float)random()/RANDOM_MAX;
				deltaX = 2.0*(sx-0.5)*CDC_R;
				float sy = (float)random()/RANDOM_MAX;
				deltaY = 2.0*(sy-0.5)*CDC_R;
			}while(sqrt(deltaX*deltaX + deltaY*deltaY)>CDC_R);
			x+=deltaX;
			y+=deltaY;
			z+= SampleGaussian(CDC_Z_SIGMA);

			cdctruthpoint->r = sqrt(x*x + y*y);
			cdctruthpoint->phi = atan2(y,x);
			if(cdctruthpoint->phi<0.0)cdctruthpoint->phi += 2.0*M_PI;
			cdctruthpoint->z = z;
		}
	}
}

//-----------
// AddNoiseHitsCDC
//-----------
void AddNoiseHitsCDC(s_HDDM_t *hddm_s)
{
	// Since the current pattern finding algorithm only uses values
	// in the s_CdcPoint_t structure, it doesn't really matter what
	// ring,straw we put the noise hits under. Thus, we put them under
	// the first straw of the first ring. If there is no ring/straw, then
	// no hits are added.

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->centralDC == HDDM_NULL ||
			hits->centralDC->cdcTruthPoints == HDDM_NULL)continue;
		
		
		// Get existing hits
		s_CdcTruthPoints_t *old_cdctruthpoints = hits->centralDC->cdcTruthPoints;
		unsigned int Nold = old_cdctruthpoints->mult;
		
		// How many noise hits to add
		unsigned int Nhits = (int)(CDC_AVG_NOISE_HITS + SampleGaussian(sqrt(CDC_AVG_NOISE_HITS)));
		s_CdcTruthPoints_t* cdctruthpoints = make_s_CdcTruthPoints(Nhits + Nold);

		// Add real hits back in first
		cdctruthpoints->mult = 0;
		for(unsigned int j=0; j<Nold; j++){
			cdctruthpoints->in[cdctruthpoints->mult++] = old_cdctruthpoints->in[j];
		}
		
		// Delete memory used for old hits structure and
		// replace pointer in HDDM tree with ours
		free(old_cdctruthpoints);
		hits->centralDC->cdcTruthPoints = cdctruthpoints;
		
		// Add noise hits
		for(unsigned int j=0; j<Nhits; j++){
			s_CdcTruthPoint_t *cdc = &cdctruthpoints->in[cdctruthpoints->mult++];
			cdc->dEdx = 0.0;
			cdc->dradius = 0.0;
			
			// To get a uniform sampling, we have to sample evenly in X/Y
			// until we find a point in the CDC fiducial area
			double x,y,r = 0.0;
			while(r<13.0 || r>59.0){
				x = SampleRange(-60.0, 60.0);
				y = SampleRange(-60.0, 60.0);
				r = sqrt(x*x + y*y);
			}
			
			cdc->phi = M_PI + atan2(y,x);
			cdc->primary = 1;
			cdc->r = r;
			cdc->track = 0;
			cdc->z = SampleRange(17.0, 217.0);
		}
	}
}

//-----------
// SmearFDC
//-----------
void SmearFDC(s_HDDM_t *hddm_s)
{
	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->forwardDC == HDDM_NULL ||
			hits->forwardDC->fdcChambers == HDDM_NULL)continue;

		s_FdcChambers_t* fdcChambers = hits->forwardDC->fdcChambers;
		s_FdcChamber_t *fdcChamber = fdcChambers->in;
		for(unsigned int j=0; j<fdcChambers->mult; j++, fdcChamber++){
			s_FdcTruthPoints_t *fdcTruthPoints = fdcChamber->fdcTruthPoints;
			if(fdcTruthPoints == HDDM_NULL)continue;
			
			s_FdcTruthPoint_t *truth = fdcTruthPoints->in;
			for(unsigned int k=0; k<fdcTruthPoints->mult; k++, truth++){
							
				// Take the simple approach here since
				// we don't know the orientation of the plane
				truth->x += 2.0*((float)random()/RANDOM_MAX-0.5)*FDC_R;
				truth->y += 2.0*((float)random()/RANDOM_MAX-0.5)*FDC_R;
							
				// Z should be well known for the FDC, but smear it slightly anyway
				truth->z += SampleGaussian(FDC_Z);
			}
		}
	}
}

//-----------
// AddNoiseHitsFDC
//-----------
void AddNoiseHitsFDC(s_HDDM_t *hddm_s)
{
	// Since the current pattern finding algorithm only uses values
	// in the s_FdcPoint_t structure, it doesn't really matter what
	// plane,wire we put the noise hits under. Thus, we put them under
	// the first wire of the first plane. If there is no plane/wire, then
	// no hits are added.

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return;
	
	for(unsigned int i=0; i<PE->mult; i++){
		s_HitView_t *hits = PE->in[i].hitView;
		if (hits == HDDM_NULL ||
			hits->forwardDC == HDDM_NULL ||
			hits->forwardDC->fdcChambers == HDDM_NULL)continue;

		s_FdcChambers_t* fdcChambers = hits->forwardDC->fdcChambers;
		s_FdcChamber_t *fdcChamber = fdcChambers->in;
		for(unsigned int j=0; j<fdcChambers->mult; j++, fdcChamber++){
			s_FdcTruthPoints_t *old_fdcTruthPoints = fdcChamber->fdcTruthPoints;
			if(old_fdcTruthPoints == HDDM_NULL)continue;
			
			// Get existing hits
			unsigned int Nold = old_fdcTruthPoints->mult;
			
			// How many noise hits to add
			unsigned int Nhits = (int)(FDC_AVG_NOISE_HITS + SampleGaussian(sqrt(FDC_AVG_NOISE_HITS)));
			s_FdcTruthPoints_t* fdcTruthPoints = make_s_FdcTruthPoints(Nhits + Nold);
			fdcTruthPoints->mult = 0;
	
			// Add real hits back in first
			s_FdcTruthPoint_t *truth = fdcTruthPoints->in;
			for(unsigned int k=0; k<Nold; k++, truth++){
				fdcTruthPoints->in[fdcTruthPoints->mult++] = old_fdcTruthPoints->in[k];
			}
		
			// Delete memory used for old hits structure and
			// replace pointer in HDDM tree with ours
			free(old_fdcTruthPoints);
			fdcChamber->fdcTruthPoints = fdcTruthPoints;

			// Add noise hits
			for(unsigned int k=0; k<Nhits; k++){
				s_FdcTruthPoint_t *fdc = &fdcTruthPoints->in[fdcTruthPoints->mult++];
				fdc->dEdx = 0.0;
				fdc->dradius = 0.0;
				fdc->primary = 1;
				fdc->track = 0;
			
				// To get a uniform sampling, we have to sample evenly in X/Y
				// until we find a point in the FDC fiducial area
				double x,y,r = 0.0;
				while(r<3.5 || r>59.0){
					x = SampleRange(-60.0, 60.0);
					y = SampleRange(-60.0, 60.0);
					r = sqrt(x*x + y*y);
				}
				fdc->x = x;
				fdc->y = y;
				fdc->z = 240.0 + SampleRange(-6.0, +6.0);
				fdc->z += ((float)(random()%4)) * 52.0;
			}
		}
	}
}

//-----------
// SmearFCAL
//-----------
void SmearFCAL(s_HDDM_t *hddm_s)
{



}
//-----------
// SmearBCAL
//-----------
void SmearBCAL(s_HDDM_t *hddm_s)
{



}
//-----------
// SmearTOF
//-----------
void SmearTOF(s_HDDM_t *hddm_s)
{



}
//-----------
// SmearUPV
//-----------
void SmearUPV(s_HDDM_t *hddm_s)
{



}
//-----------
// SmearCherenkov
//-----------
void SmearCherenkov(s_HDDM_t *hddm_s)
{



}
