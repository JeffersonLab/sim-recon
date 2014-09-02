// $Id: smear.cc 2432 2007-02-06 04:19:48Z davidl $
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include <math.h>
#include <stdlib.h>
#include "HDDM/hddm_s.hpp"


#define rad2deg 57.29577951308

#ifndef _DBG_
#define _DBG_ std::cerr << __FILE__ << ":" << __LINE__ << " "
#define _DBG__ _DBG_ << std::endl
#endif

typedef hddm_s::ProductList::iterator PlistIter;
bool DescendedFrom(PlistIter start, int pdgtype, hddm_s::ProductList &prods);
PlistIter GetParent(PlistIter start, hddm_s::ProductList &prods);

vector<PlistIter> etas;
vector<PlistIter> protons;
vector<PlistIter> neutrons;
vector<PlistIter> pi0s;
vector<PlistIter> pips;
vector<PlistIter> gammas;
vector<PlistIter> gammas_not_from_eta;
vector<PlistIter> gammas_not_from_pi0;
vector<PlistIter> gammas_not_from_eta_or_pi0;
vector<PlistIter> pi0s_not_from_eta;
vector<PlistIter> pips_not_from_eta;


//-----------
// Filter
//-----------
bool Filter(hddm_s::HDDM &record)
{
   // Return "true" to keep event, "false" to throw it away. We also
   // fill in some global variables with info we find so that we don't
   // have to go through this more than once per event.
   
   // Reset globals
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
   
   // Make lists of each type of particle of interest from all
   // levels of the decay chain
   hddm_s::ProductList prods = record.getProducts();
   hddm_s::ProductList::iterator piter;
   for (piter = prods.begin(); piter != prods.end(); ++piter) {
      // Record eta and nucleons from anywhere in the decay chain
      int pdgtype = piter->getPdgtype();
      if (pdgtype == 221)
         etas.push_back(piter);
      if (pdgtype == 2212)
         protons.push_back(piter);
      if (pdgtype == 2112)
         neutrons.push_back(piter);
      if (pdgtype == 111 )
         pi0s.push_back(piter);
      if (pdgtype == 211 )
         pips.push_back(piter);
      if (pdgtype == 22  )
         gammas.push_back(piter);
   }

   // Make list of gammas that are NOT from eta or pi0 decay chain
   for (unsigned int i=0; i < gammas.size(); i++) {
      bool a = DescendedFrom(gammas[i], 221, prods);
      bool b = DescendedFrom(gammas[i], 111, prods);
      if (!a)
         gammas_not_from_eta.push_back(gammas[i]);
      if (!b)
         gammas_not_from_pi0.push_back(gammas[i]);
      if ((!a) && (!b))
         gammas_not_from_eta_or_pi0.push_back(gammas[i]);
   }

   // Make list of pi0s that are NOT from eta decay chain
   for (unsigned int i=0; i<pi0s.size(); i++) {
      if (! DescendedFrom(pi0s[i], 221, prods))
         pi0s_not_from_eta.push_back(pi0s[i]);
   }

   // Make list of pips that are NOT from eta decay chain
   for (unsigned int i=0; i<pips.size(); i++) {
      if (! DescendedFrom(pips[i], 221, prods))
         pips_not_from_eta.push_back(pips[i]);
   }

#if 0
if (etas.size() == 1 && protons.size() == 1 && pi0s.size()>0) {
_DBG_ << "---------- Event: " << eventnumber << " -----------" << std::endl;
_DBG_ << "               etas:" << etas.size() << std::endl;
_DBG_ << "            protons:" << protons.size() << std::endl;
_DBG_ << "           neutrons:" << neutrons.size() << std::endl;
_DBG_ << "  pi0s_not_from_eta:" << pi0s_not_from_eta.size() << std::endl;
_DBG_ << "  pips_not_from_eta:" << pips_not_from_eta.size() << std::endl;
_DBG_ << "gammas_not_from_eta:" << gammas_not_from_eta.size() << std::endl;
_DBG_ << "gammas_not_from_pi0:" << gammas_not_from_pi0.size() << std::endl;
_DBG_ << "gammas_not_from_eta_or_pi0:" << gammas_not_from_eta_or_pi0.size() << std::endl;
}
#endif
   return (etas.size() > 0);
}

//-----------
// DescendedFrom
//-----------
bool DescendedFrom(PlistIter start, int pdgtype, hddm_s::ProductList &prods)
{
   do {
      if (start->getPdgtype() == pdgtype)
         return true;
      start = GetParent(start, prods);
   } while (start != prods.end());
   
   return false;
}

//-----------
// GetParent
//-----------
PlistIter GetParent(PlistIter start, hddm_s::ProductList &prods)
{
   // Loop over products until parent is found
   PlistIter iter;
   for (iter = prods.begin(); iter != prods.end(); ++iter) {
      if (iter->getId() == start->getParentid())
         return iter;
   }
   return prods.end();
}

//-----------
// Filter_eta_p
//-----------
bool Filter_eta_p(hddm_s::HDDM &record)
{
   return   etas.size() == 1 
         && protons.size() == 1
         && neutrons.size() == 0
         && pi0s_not_from_eta.size() == 0
         && pips_not_from_eta.size() == 0
         && gammas_not_from_eta.size() == 0;
}

//-----------
// Filter_eta_p_pi0
//-----------
bool Filter_eta_p_pi0(hddm_s::HDDM &record)
{
   return   etas.size() == 1 
         && protons.size() == 1
         && neutrons.size() == 0
         && pi0s_not_from_eta.size() == 1
         && pips_not_from_eta.size() == 0
         && gammas_not_from_eta_or_pi0.size() == 0;
}

//-----------
// Filter_eta_n_pip
//-----------
bool Filter_eta_n_pip(hddm_s::HDDM &record)
{
   return   etas.size() == 1 
         && protons.size() == 0
         && neutrons.size() == 1
         && pi0s_not_from_eta.size() == 0
         && pips_not_from_eta.size() == 1
         && gammas_not_from_eta_or_pi0.size() == 0;
}

//-----------
// Filter_eta_p_gamma
//-----------
bool Filter_eta_p_gamma(hddm_s::HDDM &record)
{
   return   etas.size() == 1 
         && protons.size() == 1
         && neutrons.size() == 0
         && pi0s_not_from_eta.size() == 0
         && pips_not_from_eta.size() == 0
         && gammas_not_from_eta_or_pi0.size() == 1;
}
