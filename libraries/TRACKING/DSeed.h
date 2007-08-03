// $Id$
//
//    File: DSeed.h
// Created: Wed May 23 12:19:35 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp61.jlab.org 8.9.1 i386)
//

#ifndef _DSeed_
#define _DSeed_

#include <TH1.h>

#include "DQuickFit.h"

#include <JANA/jerror.h>

class DSeed{
	public:
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DSeed";}

		vector<Dtrk_hit*> hits_in_seed;
		vector<Dtrk_hit*> hits_on_seed_circle; ///< hits close to circle fit only to seed hits
		vector<Dtrk_hit*> hits_on_circle; ///< hits close to circle fit to all hits_on_seed_circle hits
		vector<Dtrk_hit*> hits_on_circle_with_z; ///< subset of hits_on_circle with z-info
		vector<Dtrk_hit*> hits_on_line; ///< subset of hits_on_circle_with_z that are close to phi-z line
		vector<Dtrk_hit*> hits_on_track;
		
		DQuickFit seed_fit;
		DQuickFit circle_fit1;
		DQuickFit circle_fit2;
		DQuickFit line_fit;
		
		double phizangle;
		double z_vertex;
		
		TH1F *phizangle_hist;
		TH1F *zvertex_hist;
		
		void Clear(void){
			hits_in_seed.clear();
			hits_on_seed_circle.clear();
			hits_on_circle.clear();
			hits_on_circle_with_z.clear();
			hits_on_line.clear();
			hits_on_track.clear();
			
			seed_fit.Clear();
			circle_fit1.Clear();
			circle_fit2.Clear();
			line_fit.Clear();
		}
		
		DSeed(void){phizangle_hist=zvertex_hist=NULL;}
		DSeed(const DSeed& seed){
			if(this==&seed)return;
			Copy(seed);
		}
		void Copy(const DSeed& seed){
			hits_in_seed				= seed.hits_in_seed;
			hits_on_seed_circle		= seed.hits_on_seed_circle;
			hits_on_circle				= seed.hits_on_circle;
			hits_on_circle_with_z	= seed.hits_on_circle_with_z;
			hits_on_line				= seed.hits_on_line;
			hits_on_track				= seed.hits_on_track;
			
			seed_fit						= seed.seed_fit;
			circle_fit1					= seed.circle_fit1;
			circle_fit2					= seed.circle_fit2;
			line_fit						= seed.line_fit;

			phizangle					= seed.phizangle;
			z_vertex						= seed.z_vertex;
			phizangle_hist				= (TH1F*)seed.phizangle_hist->Clone();
			zvertex_hist				= (TH1F*)seed.zvertex_hist->Clone();
		}
		virtual ~DSeed(void){
			if(phizangle_hist)delete phizangle_hist;
			if(zvertex_hist)delete zvertex_hist;
		}
		DSeed& operator=(const DSeed& seed){
			if(this==&seed)return *this;
			Copy(seed);
			return *this;
		}

	protected:	
	private:
};

#endif // _DSeed_

