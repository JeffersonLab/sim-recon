// $Id$
//
//    File: DPSPair_factory.cc
// Created: Fri Mar 20 07:51:31 EDT 2015
// Creator: nsparks (on Linux cua2.jlab.org 2.6.32-431.5.1.el6.x86_64 x86_64)
//
// 
// 10/20/2017  A.S. Significant changes in the PS pair search algorithm, perform hit clusterization
//  Change reported structure to PSClust 


#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

#include "DPSPair_factory.h"
#include "DPSHit.h"
using namespace jana;

inline bool SortByTimeDifference(const DPSPair* pair1, const DPSPair* pair2)
{
  double tdiff1 = fabs(pair1->ee.first->t-pair1->ee.second->t);
  double tdiff2 = fabs(pair2->ee.first->t-pair2->ee.second->t);
  return (tdiff1<tdiff2);
}



//------------------
// init
//------------------
jerror_t DPSPair_factory::init(void)
{
  DELTA_T_CLUST_MAX  =  10.0; // ns
  DELTA_T_PAIR_MAX   =  10.0; // ns

  gPARMS->SetDefaultParameter("PSPair:DELTA_T_CLUST_MAX",DELTA_T_CLUST_MAX,
			      "Maximum difference in ns between hits in a cluster"
			      " in left and right arm of fine PS");

  gPARMS->SetDefaultParameter("PSPair:DELTA_T_PAIR_MAX",DELTA_T_PAIR_MAX,
			      "Maximum difference in ns between a pair of hits"
			      " in left and right arm of fine PS");  
  
  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DPSPair_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DPSPair_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{


  // get fine pair spectrometer hits
  vector<const DPSHit*> hits;
  loop->Get(hits);

  int debug = 0;

  tiles_left.clear();
  tiles_right.clear();

  clust_left.clear();
  clust_right.clear();


  for (unsigned int ii = 0; ii < hits.size(); ii++) {
    if( hits[ii]->arm == 0 ) {
      tile tmp;
      tmp.column      = hits[ii]->column;
      tmp.energy      = hits[ii]->E;
      tmp.integral    = hits[ii]->integral;
      tmp.pulse_peak  = hits[ii]->pulse_peak;
      tmp.time   = hits[ii]->t;
      tmp.used   = 0;
      tiles_left.push_back(tmp);
    }
    if( hits[ii]->arm == 1 ) {
      tile tmp;
      tmp.column      = hits[ii]->column;
      tmp.energy      = hits[ii]->E;
      tmp.integral    = hits[ii]->integral;
      tmp.pulse_peak  = hits[ii]->pulse_peak;
      tmp.time   = hits[ii]->t;
      tmp.used = 0;
      tiles_right.push_back(tmp);
    }
  }
  
  sort(tiles_left.begin(), tiles_left.end(),SortByTile);
  sort(tiles_right.begin(),tiles_right.end(),SortByTile);
  
  int last_tile = -1;
  
  vector<int> my_index;
  

  if(debug)
    cout << " Tiles Left Side = " << tiles_left.size() << endl;
  
  if(tiles_left.size() > 0){
    
    for (unsigned int ii = 0; ii < tiles_left.size(); ii++){
      
      my_index.clear();
      
      if( tiles_left[ii].used == 1) continue;
      
      tiles_left[ii].used = 1;
      
      last_tile = tiles_left[ii].column;
      
      my_index.push_back(ii);
      
      for (unsigned int jj = ii + 1; jj < tiles_left.size(); jj++){
	if(fabs(tiles_left[ii].time - tiles_left[jj].time) < DELTA_T_CLUST_MAX){
	  if( std::abs(tiles_left[jj].column - last_tile) <= 1){
	    
	    tiles_left[jj].used = 1;
	    last_tile = tiles_left[jj].column;
	    
	    my_index.push_back(jj);
	  }
	}
      }
      
      clust tmp;
      tmp.energy = 1;
      tmp.hit_index = my_index;
      
      clust_left.push_back(tmp);  
    }
  }  

  /* Arm B,  South, Right Side */
    
  last_tile = -1;
  
  if(debug){
    cout << endl;
    cout << " Tiles Right Side = " << tiles_right.size() << endl;
  }  

  if(tiles_right.size() > 0){
    
    for (unsigned int ii = 0; ii < tiles_right.size(); ii++){
      
      my_index.clear();
      
      if( tiles_right[ii].used == 1) continue;
	
      tiles_right[ii].used = 1;
      
      last_tile = tiles_right[ii].column;
      
      my_index.push_back(ii);
      
      for (unsigned int jj = ii + 1; jj < tiles_right.size(); jj++){
	if(fabs(tiles_right[ii].time - tiles_right[jj].time) < DELTA_T_CLUST_MAX){
	  if( std::abs(tiles_right[jj].column - last_tile) <= 1){
	    
	    tiles_right[jj].used = 1;
	    last_tile = tiles_right[jj].column;
	    
	    my_index.push_back(jj);
	  }
	}
      }
      
      clust tmp;
      tmp.energy = 1;
      tmp.hit_index = my_index;
      
      clust_right.push_back(tmp);  
    }
  }


  // Cluster energy and time

  double max_integral  =  0.;
  double max_peak      =  0.;
  double max_time      =  0.;
  int    max_column    = -1;

  if( clust_left.size() > 0){
    for(unsigned ii = 0; ii < clust_left.size(); ii++){
      unsigned int nhits = clust_left[ii].hit_index.size();
      if(nhits > 0){
	double norm     = 0.;
	double en_tmp   = 0.;
	double time_tmp = 0.;

	max_integral = 0.;
	
	for(unsigned jj = 0; jj < nhits; jj++){
	  int index = clust_left[ii].hit_index[jj];
	  en_tmp    +=  tiles_left[index].energy*tiles_left[index].integral;
	  time_tmp  +=  tiles_left[index].time*tiles_left[index].integral;
	  norm      +=  tiles_left[index].integral;
	  if(tiles_left[index].integral > max_integral) {
	    max_integral  =  tiles_left[index].integral;
	    max_peak      =  tiles_left[index].pulse_peak;
	    max_time      =  tiles_left[index].time;
	    max_column    =  tiles_left[index].column;
	  }
        }
	if(norm > 0){
	  clust_left[ii].column      =  max_column;
	  clust_left[ii].energy      =  en_tmp/norm;
	  clust_left[ii].time        =  time_tmp/norm;
	  clust_left[ii].integral    =  max_integral;
	  clust_left[ii].pulse_peak  =  max_peak;
	  clust_left[ii].time_tile   =  max_time;
	} else {
	  clust_left[ii].column      =  0;
	  clust_left[ii].energy      =  0;
	  clust_left[ii].time        =  0;
	  clust_left[ii].integral    =  0;
	  clust_left[ii].pulse_peak  =  0;
	  clust_left[ii].time_tile   =  0;
	}	
      }
    }
  }

  max_integral  =  0.;
  max_peak      =  0.;
  max_time      =  0.;
  max_column    = -1;

  if( clust_right.size() > 0){
    for(unsigned ii = 0; ii < clust_right.size(); ii++){
      unsigned int nhits = clust_right[ii].hit_index.size();
      if(nhits > 0){
	double norm     = 0.;
	double en_tmp   = 0.;
	double time_tmp = 0.;


	max_integral = 0.;

	for(unsigned jj = 0; jj < nhits; jj++){
	  int index = clust_right[ii].hit_index[jj];
	  en_tmp    +=  tiles_right[index].energy*tiles_right[index].integral;
	  time_tmp  +=  tiles_right[index].time*tiles_right[index].integral;
	  norm      +=  tiles_right[index].integral;
	  if(tiles_right[index].integral > max_integral) {
	    max_integral  =  tiles_right[index].integral;
	    max_peak      =  tiles_right[index].pulse_peak;
	    max_time      =  tiles_right[index].time;
	    max_column    =  tiles_right[index].column;
	  }	  
	}
	if(norm > 0){
	  clust_right[ii].column      =  max_column;
	  clust_right[ii].energy      =  en_tmp/norm;
	  clust_right[ii].time        =  time_tmp/norm;
	  clust_right[ii].integral    =  max_integral;
	  clust_right[ii].pulse_peak  =  max_peak;
	  clust_right[ii].time_tile   =  max_time;
	} else {
	  clust_right[ii].column      =  0;
	  clust_right[ii].energy      =  0;
	  clust_right[ii].time        =  0;
	  clust_right[ii].integral    =  0;
	  clust_right[ii].pulse_peak  =  0;
	  clust_right[ii].time_tile   =  0;
	}	
      }
    }
  }
  

  if(debug)
    cout << " Number of Left clusters found = " <<  clust_left.size() << endl;

  for(unsigned ii = 0; ii < clust_left.size(); ii++){

    if(debug)
      cout << " Number of tiles inside the cluster = " << 
	clust_left[ii].hit_index.size() << endl;

    for(unsigned jj = 0; jj < clust_left[ii].hit_index.size(); jj++){
      int index_l = clust_left[ii].hit_index[jj];
      if(debug)
	cout << " Index = " << index_l <<  
	  "  Tile     =  " <<  tiles_left[index_l].column  <<  
	  "  Energy   =  " <<  tiles_left[index_l].energy << 
	  "  Integral =  " <<  tiles_left[index_l].integral << endl;    
    }

    if(debug)
      cout << " Cluster ENERGY:    " << ii << "   " << clust_left[ii].energy << 
	"  TIME = "  << clust_left[ii].time << endl;

  }    
  

  if(debug){
    cout << endl;
    cout << endl;
    cout << " Number of Right clusters found = " <<  clust_right.size() << endl;
  }


  for(unsigned ii = 0; ii < clust_right.size(); ii++){
    
    if(debug)
      cout << " Number of tiles inside the cluster = " << 
	clust_right[ii].hit_index.size() << endl;
    
    for(unsigned jj = 0; jj < clust_right[ii].hit_index.size(); jj++){
      int index_r = clust_right[ii].hit_index[jj];
      if(debug)
	cout << " Index = " << index_r <<  
	  "  Tile     =  " <<  tiles_right[index_r].column  <<  
	  "  Energy   =  " <<  tiles_right[index_r].energy << 
	  "  Integral =  " <<  tiles_right[index_r].integral << endl;    
    }
    
    if(debug)
      cout << " Cluster ENERGY:    " << ii << "   " << clust_right[ii].energy << 
	"  TIME = "  << clust_right[ii].time << endl;
    
  }    


  if( (clust_left.size() > 0) && (clust_right.size() > 0)) {
    for(unsigned ii = 0; ii < clust_left.size(); ii++){
      for(unsigned jj = 0; jj < clust_right.size(); jj++){
	
	if(fabs(clust_left[ii].time - clust_right[jj].time) < DELTA_T_PAIR_MAX){
	  
	  if(debug){
	    cout << endl;
	    cout << " PAIR FOUND:  " <<  "  E =  " << 
	      clust_left[ii].energy + clust_right[jj].energy << endl;
	    cout << endl;
	  }	    


	  DPSPair *pair = new DPSPair;

	  pair->SetPair(clust_left[ii].column, clust_left[ii].pulse_peak, clust_left[ii].integral, clust_left[ii].time_tile, 
			clust_left[ii].hit_index.size(), clust_left[ii].energy, clust_left[ii].time, 
			clust_right[jj].column, clust_right[jj].pulse_peak, clust_right[jj].integral, clust_right[jj].time_tile, 
			clust_right[jj].hit_index.size(), clust_right[jj].energy, clust_right[jj].time );
	  
	  _data.push_back(pair);
	  	  
	}
	
	
      }
    }    
  } 

  sort(_data.begin(),_data.end(), SortByTimeDifference);

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DPSPair_factory::erun(void)
{
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DPSPair_factory::fini(void)
{
  return NOERROR;
}

bool  DPSPair_factory::SortByTile(const tile &tile1, const tile &tile2)
{
  if(tile1.column == tile2.column) return tile1.time < tile2.time;
  return (tile1.column < tile2.column);
}
