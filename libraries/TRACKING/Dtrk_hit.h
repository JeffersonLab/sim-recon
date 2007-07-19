// $Id$
//
//    File: Dtrk_hit.h
// Created:Tue May 15 22:56:53 EDT 2007
// Creator: davidl (on Darwin 8.9.1 i386)
//

#ifndef _Dtrk_hit_
#define _Dtrk_hit_

#include <math.h>

#include "JANA/jerror.h"
#include "GlueX.h"
#include "DVector3.h"
#include "DMatrix.h"
#include "DCoordinateSystem.h"

class Dtrk_hit:public DVector3{
	public:
	
		Dtrk_hit(const DVector3 &v):DVector3(v){}
		Dtrk_hit(const DVector3 &v, const Dtrk_hit &hit):DVector3(v){
			covariance = hit.covariance;
			flags = hit.flags;
			phi_circle = hit.phi_circle;
			hitid = hit.hitid;
			wire = hit.wire;
			system = hit.system;
		}
		
		enum Dtrk_hit_Flags_t{
			USED				= 0x0001,
			IN_SEED			= 0x0002,
			ON_CIRCLE		= 0x0004,
			ON_LINE			= 0x0008,
			ON_TRACK			= 0x0010,
			IGNORE			= 0x0020,
		};

		DMatrix covariance;			///< Covariance matrix in XYZ of hit position
		unsigned int flags;			///< see Dtrk_hit_Flags_t enum above
		float phi_circle;				///< phi angle relative to helix center
		oid_t hitid;					///< JANA object id of hit this came from
		const DCoordinateSystem *wire; ///< Wire (if any) this hit came from
		DetectorSystem_t system;	///< detector system this hit came from (see GlueX.h)
		
		inline float DistXY(Dtrk_hit *hit){return sqrt(DistXY2(hit));}
		inline float DistXY2(Dtrk_hit *hit){
			float dx = X() - hit->X();
			float dy = Y() - hit->Y();
			return dx*dx + dy*dy;
		}
		inline float CosPhiXY(Dtrk_hit *a, Dtrk_hit *b){
			float adx = a->X() - X();
			float ady = a->Y() - Y();
			float bdx = b->X() - X();
			float bdy = b->Y() - Y();
			return (adx*bdx + ady*bdy)/sqrt((adx*adx+ady*ady)*(bdx*bdx+bdy*bdy));
		}
		
	protected:
	
	private:

};

#endif // _Dtrk_hit_

