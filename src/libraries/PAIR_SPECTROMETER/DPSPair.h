// $Id$
//
//    File: DPSPair.h
// Created: Fri Mar 20 07:51:31 EDT 2015
// Creator: nsparks (on Linux cua2.jlab.org 2.6.32-431.5.1.el6.x86_64 x86_64)
//

#ifndef _DPSPair_
#define _DPSPair_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include "DPSHit.h"

class DPSPair:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DPSPair);

  ~ DPSPair(){
    //    cout << " DPSPair  destructor is called " << endl;
    delete left;
    delete right;
  }
  
  typedef struct {   
    int    column;       //  tile with the largest energy
    double pulse_peak;   //  pulse peak of the tile with the largest energy
    double integral;     //  total energy in the cluster 
    double t_tile;       //  time of the tile with the largest energy
    int    ntiles;       //  number of tiles in the cluster
    

    double E;            //   Cluster energy
    double t;            //   cluster time
  } PSClust;


  PSClust *left   =  new PSClust;
  PSClust *right  =  new PSClust;

  pair<const PSClust*, const PSClust*> ee;	// first:North(left); second:South(right)	


  void SetPair(int column_l, double pulse_peak_l, double integral_l, double t_tile_l, int ntiles_l, double E_l, double t_l,
	       int column_r, double pulse_peak_r, double integral_r, double t_tile_r, int ntiles_r, double E_r, double t_r){
    
    
    left->column      =  column_l;
    left->pulse_peak  =  pulse_peak_l;
    left->integral    =  integral_l;
    left->t_tile      =  t_tile_l;
    left->ntiles      =  ntiles_l;
    left->E           =  E_l;
    left->t           =  t_l;

    
    right->column      =  column_r;
    right->pulse_peak  =  pulse_peak_r;
    right->integral    =  integral_r;
    right->t_tile      =  t_tile_r;
    right->ntiles      =  ntiles_r;
    right->E           =  E_r;
    right->t           =  t_r;

    ee.first  =  left;
    ee.second =  right;


  }


  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "column_lhit", "%d", ee.first->column);
    AddString(items, "column_rhit", "%d", ee.second->column);
    AddString(items, "E_pair", "%f", ee.first->E+ee.second->E);
    AddString(items, "E_lhit", "%f", ee.first->E);
    AddString(items, "E_rhit", "%f", ee.second->E);
    AddString(items, "t_lhit", "%f", ee.first->t);
    AddString(items, "t_rhit", "%f", ee.second->t);
    AddString(items, "integral_lhit", "%f", ee.first->integral);
    AddString(items, "integral_rhit", "%f", ee.second->integral);
    AddString(items, "pulse_peak_lhit", "%f", ee.first->pulse_peak);
    AddString(items, "pulse_peak_rhit", "%f", ee.second->pulse_peak);
    AddString(items, "t_tile_lhit", "%f", ee.first->t_tile);
    AddString(items, "t_tile_rhit", "%f", ee.second->t_tile);
    AddString(items, "ntiles_lhit", "%d", ee.first->ntiles);
    AddString(items, "ntiles_rhit", "%d", ee.second->ntiles);
  }


		
};

#endif // _DPSPair_

