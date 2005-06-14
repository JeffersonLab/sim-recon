// $Id$
//
//    File: DMCReconstructed.cc
// Created: Sat Apr 23 22:25:04 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include "DMCThrown.h"
#include "DMCReconstructed.h"

//------------------
// FindClosestThrown
//------------------
void DMCReconstructed::FindClosestThrown(vector<const DMCThrown*> &mcthrowns)
{
	thrownid = -1;
	thrown_delta_p = 1000.0;
	
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		const DMCThrown *mcthrown = mcthrowns[i];
		
		// Revert to using more common difference in magnitudes
		float delta_p = mcthrown->p - p;
		if(fabs(delta_p) < fabs(thrown_delta_p)){
			thrown_delta_p = delta_p;
			thrownid = i;
		}
		
		#if 0
		// bail early if magnitudes aren't even close
		//if(fabs(mcthrown->p/p - 1.0) > 0.2)continue;
		
		float p_trans_thrown = mcthrown->p*sin(mcthrown->theta);
		float p_trans_recon = p*sin(theta);
		float px = p_trans_thrown*cos(mcthrown->phi) - p_trans_recon*cos(phi);
		float py = p_trans_thrown*sin(mcthrown->phi) - p_trans_recon*sin(phi);
		float pz = mcthrown->p*cos(mcthrown->theta) - p*cos(theta);
		float delta_p2 = px*px + py*py + pz*pz;
		if(delta_p2 < best_delta_p2){
			best_delta_p2 = delta_p2;
			thrownid = i;
		}
		#endif
	}
	
	//thrown_delta_p = sqrt(best_delta_p2);
	
	//if(thrown_delta_p/p > 0.2)thrownid = -1;

}
