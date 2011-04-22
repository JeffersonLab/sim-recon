// $Id$
//
//    File: DBCALHitGroup.h
// Created: Mon Apr 18 15:02:25 EDT 2011
// Creator: davidl (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DBCALHitGroup_
#define _DBCALHitGroup_

#include <vector>

#include <HDDM/hddm_s.h>

class DBCALHitGroup{
	public:
		DBCALHitGroup(s_BcalSiPMUpHit_t *hit):tmin(hit->t),tmax(hit->t){uphits.push_back(hit);}
		DBCALHitGroup(s_BcalSiPMDownHit_t *hit):tmin(hit->t),tmax(hit->t){downhits.push_back(hit);}
		virtual ~DBCALHitGroup(){}
	
		float tmin;
		float tmax;
		std::vector<s_BcalSiPMUpHit_t*> uphits;
		std::vector<s_BcalSiPMDownHit_t*> downhits;
		
		// Add given hit to appropriate container
		void AddHit(s_BcalSiPMUpHit_t *hit){
			uphits.push_back(hit);
			if(hit->t < tmin)tmin = hit->t;
			if(hit->t > tmax)tmax = hit->t;
		}
		void AddHit(s_BcalSiPMDownHit_t *hit){
			downhits.push_back(hit);
			if(hit->t < tmin)tmin = hit->t;
			if(hit->t > tmax)tmax = hit->t;
		}

		// Template so it can handle both up and down hits
		template<typename T>
		bool AddIfInTime(T *hit, float window){
			if(IsInTime(hit, window)){
				AddHit(hit);
				return true;
			}
			return false;
		}
		
		unsigned int GetNhits(void){return uphits.size()+downhits.size();}
		
		// Template so it can handle both up and down hits
		template<typename T>
		bool IsInTime(T *hit, float window){
			// return true if any hit in this group is within the
			// "window" ns of the given hit.
			for(unsigned int i=0; i<uphits.size(); i++){
				if(fabs(uphits[i]->t - hit->t) <= window)return true;
			}
			for(unsigned int i=0; i<downhits.size(); i++){
				if(fabs(downhits[i]->t - hit->t) <= window)return true;
			}
			
			return false;
		}
		
		// MergeHitsIfInRange
		bool MergeHitsIfInRange(DBCALHitGroup *grp){
			if(grp->GetNhits()==0)return false;
			if(GetNhits()==0)return false;
			if(grp->tmax < tmin)return false;
			if(grp->tmin > tmax)return false;
			
			// The other group will get merged into this one
			for(unsigned int i=0; i<grp->uphits.size(); i++)AddHit(grp->uphits[i]);
			for(unsigned int i=0; i<grp->downhits.size(); i++)AddHit(grp->downhits[i]);
			if(grp->tmax > tmax)tmax = grp->tmax;
			if(grp->tmin < tmin)tmin = grp->tmin;

			// Clear the other group
			grp->uphits.clear();
			grp->downhits.clear();
			
			return true;
		}
};

#endif // _DBCALHitGroup_

