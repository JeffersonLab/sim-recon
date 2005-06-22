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
void SmearFDC(s_HDDM_t *hddm_s);
void SmearFCAL(s_HDDM_t *hddm_s);
void SmearBCAL(s_HDDM_t *hddm_s);
void SmearTOF(s_HDDM_t *hddm_s);
void SmearUPV(s_HDDM_t *hddm_s);
void SmearCherenkov(s_HDDM_t *hddm_s);

double SampleGaussian(double sigma);

// The straw diameter is 1.6cm
float CDC_R = 0.8;			// cm

// The stereo pitch is 6 degrees. For hit-based tracks, the
// hit could come anywhere in the 1.6cm diameter straw. Thus,
// the z-range corresponding to the straw diameter is
// 1.6cm * cot(6 degrees) = 15cm. 
float CDC_Z_SIGMA = 7.5;	// cm

// The wire spacing in the FDC is 0.5cm. In the actual data,
// multiple planes will be combined into psuedo points since
// each plane only measures in one direction (not counting z).
// At any rate, this should probably be an over estimate.
float FDC_R = 0.25;			// cm

//-----------
// Smear
//-----------
void Smear(s_HDDM_t *hddm_s)
{
	SmearCDC(hddm_s);
	SmearFDC(hddm_s);
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
		// ------------ CdcPoints, Hits --------------
		s_Rings_t *rings=NULL;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->centralDC)
				rings = PE->in[i].hitView->centralDC->rings;
		if(rings){
			for(unsigned int j=0;j<rings->mult;j++){
				//float radius = rings->in[j].radius;
				s_Straws_t *straws = rings->in[j].straws;
				if(straws){
					for(unsigned int k=0;k<straws->mult;k++){
						//float phim = straws->in[k].phim;
						s_CdcPoints_t *cdcpoints = straws->in[k].cdcPoints;
						if(cdcpoints){
							for(unsigned int m=0;m<cdcpoints->mult;m++){
								// Here we want to move the point to a position
								// randomly sampled from a disc in the x/y
								// plane centered on the actual point.
								float r = cdcpoints->in[m].r;
								float phi = cdcpoints->in[m].phi;
								float z = cdcpoints->in[m].z;
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

								cdcpoints->in[m].r = sqrt(x*x + y*y);
								cdcpoints->in[m].phi = atan2(y,x);
								if(cdcpoints->in[m].phi<0.0)cdcpoints->in[m].phi += 2.0*M_PI;
								cdcpoints->in[m].z = z;
							}
						}
					}
				}
			}
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
		s_Chambers_t *chambers = NULL;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->forwardDC)
				chambers = PE->in[i].hitView->forwardDC->chambers;
		if(!chambers)continue;
		
		for(unsigned int j=0;j<chambers->mult;j++){

			s_AnodePlanes_t *anodeplanes = chambers->in[j].anodePlanes;
			if(anodeplanes){
			
				for(unsigned int k=0;k<anodeplanes->mult;k++){
					s_Wires_t *wires = anodeplanes->in[k].wires;
					if(!wires)continue;
				
					for(unsigned int m=0;m<wires->mult;m++){
						s_FdcPoints_t *fdcPoints = wires->in[m].fdcPoints;
						if(!fdcPoints)continue;
						for(unsigned int n=0;n<fdcPoints->mult;n++){
							float x = fdcPoints->in[n].x;
							float y = fdcPoints->in[n].y;
							
							// Take the simple approach here since
							// we don't know the orientation of the plane
							x += 2.0*((float)random()/RANDOM_MAX-0.5)*FDC_R;
							y += 2.0*((float)random()/RANDOM_MAX-0.5)*FDC_R;
							
							// Z should be well known for the FDC
							fdcPoints->in[n].x = x;
							fdcPoints->in[n].y = y;
						}
					}
				}
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
