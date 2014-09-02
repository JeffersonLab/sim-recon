// $Id: smear.cc 2432 2007-02-06 04:19:48Z davidl $
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
using namespace std;

#include <math.h>
#include "HDDM/hddm_s.hpp"


#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "

//-----------
// Filter
//-----------
bool Filter(hddm_s::HDDM &record)
{
   // Return "true" to keep event, "false" to throw it away
   
   // Loop over Physics Events
   hddm_s::PhysicsEventList pes = record.getPhysicsEvents();
   if (pes.size() == 0)
      return false;
   
   //------------- FCAL -------------
   double Efcal = 0.0;
   hddm_s::FcalTruthHitList fcals = record.getFcalTruthHits();
   hddm_s::FcalTruthHitList::iterator fiter;
   for (fiter = fcals.begin(); fiter != fcals.end(); ++fiter) {
      Efcal += fiter->getE();
   }
   //_DBG_ << "Efcal=" << Efcal << std::endl;
   // There must be at least 0.5 GeV in the FCAL to pass the level-1 trigger
   if (Efcal < 0.5)
      return false;
   
   //------------- BCAL -------------
   double Ebcal = 0.0;
   hddm_s::BcalTruthHitList bcals = record.getBcalTruthHits();
   hddm_s::BcalTruthHitList::iterator biter;
   for (biter = bcals.begin(); biter != bcals.end(); ++biter) {
      Ebcal += biter->getE();
   }
   //_DBG_ << "Ebcal=" << Ebcal << endl;
   
   //------------- TOF -------------
   int Ntof_north = 0;
   int Ntof_south = 0;
   hddm_s::FtofTruthHitList ftofs = record.getFtofTruthHits();
   hddm_s::FtofTruthHitList::iterator titer;
   for (titer = ftofs.begin(); titer != ftofs.end(); ++titer) {
      if (titer->getT() < 50.0 && titer->getT() > 0.0) {
         if (titer->getEnd() == 0)
            Ntof_north++;
         else
            Ntof_south++;
      }
   }
      
   // We want the number of TOF coincidences which we'll estimate as the
   // lesser of the north and south hits
   int Ntof = (Ntof_north < Ntof_south)? Ntof_north : Ntof_south;
   //_DBG_ << "Ntof=" << Ntof << std::endl;

   // If there are no hits in the TOF or the BCAL has more energy, then
   // cut the event
   if (Ntof == 0 || Ebcal > Efcal)
      return false;
   
   // Reject events with too many TOF hits
   else if (Ntof > 6)
      return false;
   
   return true;
}
