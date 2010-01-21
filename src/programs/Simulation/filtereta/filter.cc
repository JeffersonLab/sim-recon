// $Id: smear.cc 2432 2007-02-06 04:19:48Z davidl $
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include <math.h>
#include "HDDM/hddm_s.h"


#define rad2deg 57.29577951308

#ifndef _DBG_
#define _DBG_ cerr<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ _DBG_<<endl
#endif


bool DescendedFrom(s_Product_t *start, int pdgtype, s_Product_t *prods, int Nprods);
s_Product_t* GetParent(s_Product_t *start, s_Product_t *prods, int Nprods);

s_HDDM_t *hddm_s_global = NULL;
vector<s_Product_t*> etas;
vector<s_Product_t*> protons;
vector<s_Product_t*> neutrons;
vector<s_Product_t*> pi0s;
vector<s_Product_t*> pips;
vector<s_Product_t*> gammas;
vector<s_Product_t*> gammas_not_from_eta;
vector<s_Product_t*> gammas_not_from_pi0;
vector<s_Product_t*> gammas_not_from_eta_or_pi0;
vector<s_Product_t*> pi0s_not_from_eta;
vector<s_Product_t*> pips_not_from_eta;



//-----------
// Filter
//-----------
bool Filter(s_HDDM_t *hddm_s)
{
	// Return "true" to keep event, "false" to throw it away. We also
	// fill in some global variables with info we find so that we don't
	// have to go through this more than once per event.
	
	// Reset globals
	hddm_s_global = hddm_s; // allows us to verify that global variables hold values for this event
	etas.clear();
	protons.clear();
	neutrons.clear();
	pi0s.clear();
	pips.clear();
	gammas.clear();
	gammas_not_from_eta.clear();
	gammas_not_from_pi0.clear();
	pi0s_not_from_eta.clear();
	pips_not_from_eta.clear();
	gammas_not_from_eta_or_pi0.clear();
	
	// Find pointer to s_Product_t array and number of elements
	s_Product_t *prods = NULL;
	int Nprods = 0;

	s_PhysicsEvents_t* PE = hddm_s->physicsEvents;
	if(!PE) return false;
	int eventnumber = 0;
	
	for(unsigned int i=0; i<PE->mult; i++){
		
		eventnumber = PE->in[i].eventNo;
	
		// ------------ Reactions --------------
		s_Reactions_t *reactions=PE->in[i].reactions;
		if(!reactions)continue;

		for(unsigned int j=0; j<reactions->mult; j++){
			s_Vertices_t *vertices = reactions->in[j].vertices;
			if(vertices){
				for(unsigned int k=0; k<vertices->mult; k++){
					s_Origin_t *origin = vertices->in[k].origin;
					s_Products_t *products = vertices->in[k].products;
					if(products && origin){
						prods = products->in;
						Nprods = products->mult;
					}
				}
			}
		}
	}
	
	// Make lists of each type of particle of interest from all
	// levels of the decay chain
	for(int i=0; i<Nprods; i++){
		// Record eta and nucleons from anywhere in the decay chain
		if(prods[i].pdgtype== 221)etas.push_back(&prods[i]);
		if(prods[i].pdgtype==2212)protons.push_back(&prods[i]);
		if(prods[i].pdgtype==2112)neutrons.push_back(&prods[i]);
		if(prods[i].pdgtype==111 )pi0s.push_back(&prods[i]);
		if(prods[i].pdgtype==211 )pips.push_back(&prods[i]);
		if(prods[i].pdgtype==22  )gammas.push_back(&prods[i]);
	}

	// Make list of gammas that are NOT from eta or pi0 decay chain
	for(unsigned int i=0; i<gammas.size(); i++){
		bool a = DescendedFrom(gammas[i], 221, prods, Nprods);
		bool b = DescendedFrom(gammas[i], 111, prods, Nprods);
		if(!a)gammas_not_from_eta.push_back(gammas[i]);
		if(!b)gammas_not_from_pi0.push_back(gammas[i]);
		if((!a) && (!b))gammas_not_from_eta_or_pi0.push_back(gammas[i]);
	}

	// Make list of pi0s that are NOT from eta decay chain
	for(unsigned int i=0; i<pi0s.size(); i++){
		if(!DescendedFrom(pi0s[i], 221, prods, Nprods))pi0s_not_from_eta.push_back(pi0s[i]);
	}

	// Make list of pips that are NOT from eta decay chain
	for(unsigned int i=0; i<pips.size(); i++){
		if(!DescendedFrom(pips[i], 221, prods, Nprods))pips_not_from_eta.push_back(pips[i]);
	}

#if 0
if(etas.size()==1 && protons.size()==1 && pi0s.size()>0){
_DBG_<<"---------- Event: "<<eventnumber<<" -----------"<<endl;
_DBG_<<"               etas:"<<etas.size()<<endl;
_DBG_<<"            protons:"<<protons.size()<<endl;
_DBG_<<"           neutrons:"<<neutrons.size()<<endl;
_DBG_<<"  pi0s_not_from_eta:"<<pi0s_not_from_eta.size()<<endl;
_DBG_<<"  pips_not_from_eta:"<<pips_not_from_eta.size()<<endl;
_DBG_<<"gammas_not_from_eta:"<<gammas_not_from_eta.size()<<endl;
_DBG_<<"gammas_not_from_pi0:"<<gammas_not_from_pi0.size()<<endl;
_DBG_<<"gammas_not_from_eta_or_pi0:"<<gammas_not_from_eta_or_pi0.size()<<endl;
}
#endif
	return etas.size()>0;
}

//-----------
// DescendedFrom
//-----------
bool DescendedFrom(s_Product_t *start, int pdgtype, s_Product_t *prods, int Nprods)
{
	do{
		if(start->pdgtype == pdgtype)return true;
		start = GetParent(start, prods, Nprods);
	}while(start!=NULL);
	
	return false;
}

//-----------
// GetParent
//-----------
s_Product_t* GetParent(s_Product_t *start, s_Product_t *prods, int Nprods)
{
	// Loop over products until parent is found
	for(int i=0; i<Nprods; i++){
		if(prods[i].id == start->parentid)return &prods[i];
	}
	return NULL;
}

//-----------
// Filter_eta_p
//-----------
bool Filter_eta_p(s_HDDM_t *hddm_s)
{
	return   etas.size()==1 
			&& protons.size()==1
			&& neutrons.size()==0
			&& pi0s_not_from_eta.size()==0
			&& pips_not_from_eta.size()==0
			&& gammas_not_from_eta.size()==0;
}

//-----------
// Filter_eta_p_pi0
//-----------
bool Filter_eta_p_pi0(s_HDDM_t *hddm_s)
{
	return   etas.size()==1 
			&& protons.size()==1
			&& neutrons.size()==0
			&& pi0s_not_from_eta.size()==1
			&& pips_not_from_eta.size()==0
			&& gammas_not_from_eta_or_pi0.size()==0;
}

//-----------
// Filter_eta_n_pip
//-----------
bool Filter_eta_n_pip(s_HDDM_t *hddm_s)
{
	return   etas.size()==1 
			&& protons.size()==0
			&& neutrons.size()==1
			&& pi0s_not_from_eta.size()==0
			&& pips_not_from_eta.size()==1
			&& gammas_not_from_eta_or_pi0.size()==0;
}

//-----------
// Filter_eta_p_gamma
//-----------
bool Filter_eta_p_gamma(s_HDDM_t *hddm_s)
{
	return   etas.size()==1 
			&& protons.size()==1
			&& neutrons.size()==0
			&& pi0s_not_from_eta.size()==0
			&& pips_not_from_eta.size()==0
			&& gammas_not_from_eta_or_pi0.size()==1;
}


