#ifndef _DBCALTruthHit_
#define _DBCALTruthHit_

/*
 *  DBCALPoint.h
 *
 *  Created by Matthew Shepherd on 3/13/11.
 *
 */

#include "BCAL/DBCALTruthCell.h"

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

using namespace jana;

/** 
   This object gives a z position to BCAL hits using timing information.
   The z position is relative to the target center because higher objects use 
   the position of hits relative to the target in spherical coordinates.

*/

class DBCALTruthHit : public JObject {

public:

  JOBJECT_PUBLIC(DBCALTruthHit);
  
  //DBCALTruthHit(const DBCALTruthCell& cell, double attenutation_length, double c_effective, double track_p0, double track_p1, double track_p2, int end_temp, int chan_indx, double summed_E);
  
DBCALTruthHit(double attenutation_length, double c_effective, double track_p0, double track_p1, double track_p2, int end_temp, int chan_indx, double summed_E);

//  DBCALTruthHitDown(const DBCALTruthCell& cell, double attenutation_length, double c_effective, double track_p0, double track_p1, double track_p2);

  int layer() const{return m_layer;}
  int sector() const{return m_sector;}
  int module() const{return m_module;}
  int end() const{return m_end;}
  int cellId() const{return m_cellId;}

  float E() const{return m_E;}
  float t() const{return m_t;}

  int pulse_peak() const{return m_pulse_peak;}

  void toStrings(vector<pair<string,string> > &items) const {
    AddString(items, "E(GeV)", "%5.3f", m_E);
    AddString(items, "t(ns)", "%5.1f", m_t);
    AddString(items, "module", "%i", m_module);
    AddString(items, "layer", "%i", m_layer);
    AddString(items, "sector", "%i", m_sector);
    AddString(items, "pulse_peak", "%i", m_pulse_peak);
    AddString(items, "end", "%i", m_end);
    AddString(items, "cellId", "%i", m_cellId); 
  }
  
private:
  
  int m_module;
  int m_layer;
  int m_sector;
  int m_cellId;
  double m_t;
  double m_E;
  int m_end;
  int m_pulse_peak;

};

#endif
