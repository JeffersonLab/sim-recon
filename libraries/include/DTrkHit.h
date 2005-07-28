// $Id$
//
//    File: DTrkHit.h
// Created: Wed Jul 20 13:43:55 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DTrkHit_
#define _DTrkHit_

#include <math.h>

#include "derror.h"

class DTrkHit{
	public:
		DTrkHit(float x, float y, float z, float r, float phi, int system, int track=0, int ihit=-1);
		virtual ~DTrkHit();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DTrkHit";}
		
		float x,y,z,r,phi;
		int track;
		int system;
		int ihit;
		float phi_circle;
		enum DTrkHit_Flags_t{
			USED				= 0x0001,
			IN_SEED			= 0x0002,
			ON_CIRCLE		= 0x0004,
			ON_TRACK			= 0x0008,
			IGNORE			= 0x0010,
			SYS_CDC			= 0x0100,
			SYS_FDC			= 0x0200,
			SYS_BCAL			= 0x0400,
			SYS_TOF			= 0x0800,
			SYS_CHERENKOV	= 0x1000,
			SYS_FCAL			= 0x2000,
			SYS_UPV			= 0x4000
		};
		unsigned int flags;
		
		inline float DistXY(DTrkHit *hit){return sqrt(DistXY2(hit));}
		inline float DistXY2(DTrkHit *hit){
			float dx = x - hit->x;
			float dy = y - hit->y;
			return dx*dx + dy*dy;
		}
		inline float CosPhiXY(DTrkHit *a, DTrkHit *b){
			float adx = a->x - x;
			float ady = a->y - y;
			float bdx = b->x - x;
			float bdy = b->y - y;
			return (adx*bdx + ady*bdy)/sqrt((adx*adx+ady*ady)*(bdx*bdx+bdy*bdy));
		}
		
	protected:
	
	
	private:

};

#endif // _DTrkHit_

