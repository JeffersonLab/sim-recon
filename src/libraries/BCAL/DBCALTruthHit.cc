/*
 *  DBCALPoint.cc
 *
 *  Created by Matthew Shepherd on 3/13/11.
 *
 */

#include <iostream>

using namespace std;

#include "BCAL/DBCALTruthHit.h"
#include "BCAL/DBCALTruthCell.h"
#include "BCAL/DBCALGeometry.h"

#include "units.h"

//DBCALTruthHit::DBCALTruthHit(const DBCALTruthCell& cell, double attenuation_length, double c_effective, double track_p0, double track_p1, double track_p2, int end_temp, int chan_indx, double summed_E)
DBCALTruthHit::DBCALTruthHit(double attenuation_length, double c_effective, double track_p0, double track_p1, double track_p2, int end_temp, int chan_indx, double summed_E)
{

  float fibLen = DBCALGeometry::GetBCAL_length();
/*
  int cell_layer = cell.layer;

  if(cell_layer==1) m_layer = 1;
  else if(cell_layer > 1 && cell_layer < 4) m_layer = 2;
  else if(cell_layer > 3 && cell_layer < 7) m_layer = 3;
  else if(cell_layer > 6 && cell_layer < 11) m_layer = 4;

  m_module = cell.module;
  m_sector = cell.sector;

  m_pulse_peak = 1;

  m_cellId = DBCALGeometry::cellId( m_module, m_layer, m_sector );

  if(end_temp==0){

    m_end = 0;

    double dUp = fibLen/2. + cell.zLocal;
  
    m_t = dUp/c_effective + cell.t;
 
    double atten_E = cell.E*exp(-dUp/attenuation_length);

  }

  if(end_temp==1){
 
    m_end = 1;

    double dDown = fibLen/2. - cell.zLocal;

    m_t = dDown/c_effective + cell.t;

    m_E = cell.E*exp(-dDown/attenuation_length);
  }
*/
  m_E = summed_E;
  m_end = end_temp;
  m_t = 1.;

}
/*
DBCALTruthHitDown::DBCALTruthHitUp(const DBCALTruthCell& cell, double attenuation_length, double c_effective, double track_p0, double track_p1, double track_p2)
{

  float fibLen = DBCALGeometry::GetBCAL_length();

  int cell_layer = cell.layer;

  if(cell_layer==1) m_layer = 1;
  else if(cell_layer > 1 && cell_layer < 4) m_layer = 2;
  else if(cell_layer > 3 && cell_layer < 7) m_layer = 3;
  else if(cell_layer > 6 && cell_layer < 11) m_layer = 4;

  m_module = cell.module;
  m_sector = cell.sector;

  m_pulse_peak = 1;

  m_end = 0;

  double dDown = fibLen/2. - cell.zLocal;

  m_t = dDown/c_effective;

  m_E = cell.E*exp(-dDown/attenuation_length);

}*/
