//-------------------------------------------------------------------------
//
//      Author:     Brad Futch 
//      Date:       Wed 01 Dec 2004 04:37:54 
//
//      Copyright (C) 2004 Brad Futch 
//      
//      This program is free software; you can redistribute it and/or
//      modify it under the terms of the GNU General Public License
//      as published by the Free Software Foundation; either version 2
//      of the License, or (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 59 Temple Place - Suite 330, Boston, 
//      MA  02111-1307, USA.
//
//-------------------------------------------------------------------------

#ifndef __DMAGNETICFIELDMAP_CPP__
#define __DMAGNETICFIELDMAP_CPP__

/*! \file DMagneticFieldMap.cpp
 *  \brief An implementation of a Magnetic Field Map.
 *  @see DMagnetifField.h
 */

#include "DMagneticFieldMap.h"
#include <iostream>
#include <fstream>

using std::cout;
using std::cerr;
using std::ifstream;
using std::endl;


//------------------------------------------------------------------
// NOTE:
//   
//   I make the constructor take input dimensions because 
//   as far as I can tell this information is not in the document.
//   I would think that it would be better to have this information 
//   in the file itself.
//-------------------------------------------------------------------

DMagneticFieldMap::DMagneticFieldMap(int rDim, int zDim)
{
	/// The September 21, 2001 table in dsolenoid.table file has
	/// rDim = 41 and zDim = 251.

  this->rDim = rDim;
  this->zDim = zDim;
  this->phiDim = 0; // place holder for now (1/23/05) D.L.
  this->rMin = 0.0;
  this->rMax = 0.0;
  this->zMin = 0.0;
  this->zMax = 0.0;
  this->constBz =-999.0;
  DEBUG = false;
}

DMagneticFieldMap::DMagneticFieldMap(const float Bz)
{
	// Used if the user just wants a constant magnetic field
	this->constBz = Bz;
}

DMagneticFieldMap::~DMagneticFieldMap()
{
  if (map)
    delete map;
}

void DMagneticFieldMap::setDebug(const bool val)
{
  DEBUG = val;
}

void DMagneticFieldMap::readMap(void)
{
	/// Try opening "dsolenoid.table" in the local directory
	/// If successful, pass the stream pointer to the method below

	ifstream dsolenoid("dsolenoid.table");
	if(!dsolenoid.is_open()){
		cerr<<__FILE__<<":"<<__LINE__<<" Unable to open \"dsolenoid.table!\"";
		cerr<<" Check that it exists in your working directory!"<<endl;
	}else{
		cout<<"Reading Magnetic field map from \"dsolenoid.table\""<<endl;
		readMap(dsolenoid);
		dsolenoid.close();
	}
}

void DMagneticFieldMap::readMap(ifstream &in)
{
  
  if (!in.is_open() )
  {
    cout << "DMagneticFieldMap:readMap: Error reading solenoid map!" 
         << endl;
    exit(0);
  }

  //------------------------------------------------------------------
  // NOTE: Data Format
  //
  // Based on the solenoid table format, read data (comments) until 
  // the first occurance of the characters '=='.  At this point grab 
  // the next 15 lines and throw them away.  The data begins on the
  // next line.
  //
  // The formatted input is of the form:
  //
  // <xcoord> <ycoord> <zcoord>  <Bx> <By> <Bz>
  //
  // Again see the note in the header file about the xy symmetry.
  //
  //------------------------------------------------------------------
  
  int numSkip = 9;
  char c;
  char buf[1000];

  //---------------------------
  // Skip all of the comments
  //---------------------------
  
  while (in.good())
  {
    c = in.get();
    if (c == '=')
    {
      c = in.get();
      if (c == '=')
        break;
    }
  }

  for (int i=0; i < numSkip; i++)
  {
    in.getline(buf, 1000);
  }

  //------------------------------------------------------------
  // At this point we need the dimensions, one way or the other
  //------------------------------------------------------------
  
  if (rDim == 0 || zDim == 0)
  {
    cerr << "DMagneticFieldMap:readMap: Data Dimensions not set, map "
         << "NOT set" << endl;
    return;
  }

  //-------------------------------------------
  // Allocate some space and set up the loop
  //-------------------------------------------

  D3Vector_t pos;
  D3Vector_t mag;
  int        countR;
  int        countZ;

  countR     = 0;
  countZ     = 0;
  map        = new double [rDim*zDim*3];

  while (!in.eof())
  {
    if (countZ >= zDim)
    {
      countZ = 0;
      countR++;
    }
    if (countR >= rDim)
    { 
      //cerr << "DMagneticFieldMap:readMap:Warning  Dimension overrun" 
      //     << endl;
      return;
    }
    
    if (DEBUG)
      printf("%d:%d: ", countR, countZ);
    
    in >> pos.x;
    in >> pos.y;
    in >> pos.z;

    in >> mag.x;
    in >> mag.y;
    in >> mag.z;

    if (DEBUG)
      printf("POS(%.2f, %.2f, %.2f) B(%.2f, %.2f, %.2f)\n", 
              pos.x,pos.y,pos.z,mag.x,mag.y,mag.z);

    int ind; 
    ind = serialize(countR, countZ);
    
    map[ind]     = mag.x;
    map[ind + 1] = mag.y;
    map[ind + 2] = mag.z;

    countZ++;

    //---------------------
    // Set the min and max
    //---------------------

    this->rMin = pos.x < this->rMin ? pos.x : this->rMin;
    this->rMax = pos.x > this->rMax ? pos.x : this->rMax;
    
    this->zMin = pos.z < this->zMin ? pos.z : this->zMin;
    this->zMax = pos.z > this->zMax ? pos.z : this->zMax;
  }
}

//------------------------------------
// Serialize coordinates into the map
//------------------------------------

int DMagneticFieldMap::serialize(int r, int z)
{
  return ((r*zDim) + z)*3;
}

//----------------------------------------------------------------------
// This is going to return the base indices from which to interpolate 
// as well as the percentage of the distance to the next index.
//
// This should make the task of actuall interpolating much easier.
//----------------------------------------------------------------------

void DMagneticFieldMap::getInds(const double &r, const double &z, int ind[2],
                                double &rho, double &zeta)
{

//----------------------------------------------------------------
// First we need to map the the real-space coordinates into [0,1]
//----------------------------------------------------------------

  double rN, zN;
  int baseR, baseZ;

  rN = (r - rMin)/(rMax - rMin);
  zN = (z - zMin)/(zMax - zMin);

  baseR = (int)(rN * (rDim-1));
  baseZ = (int)(zN * (zDim-1));

  rho  = ((rN * (rDim-1)) - baseR)/(rDim-1);
  zeta = ((zN * (zDim-1)) - baseZ)/(zDim-1);
  
  ind[0] = baseR;
  ind[1] = baseZ;

/*
  printf("(rMin, rMax)(zMin, zMax):  (%.2f, %.2f)(%.2f, %.2f)\n", 
          rMin, rMax, zMin, zMax);
  printf("(r,z)->(rN, zN):           (%.2f, %.2f) -> (%f, %f)\n", r, z, rN, zN); 
  printf("(baseR, baseZ)(rho, zeta): (%d, %d)(%f, %f)\n", 
         baseR, baseZ, rho, zeta);
*/

}

//---------------------------------------------------------------
// These are the getQuick methods.  They will fetch the nearest 
// entry in the data structure.
//---------------------------------------------------------------

D3Vector_t DMagneticFieldMap::getQuick(double r, double z)
{
  int ind;
  D3Vector_t vec;
  int inds [2];
  double rho, zeta;
  
  if(constBz>-100.0){
  	// homogeneous magnetic field
	vec.x = vec.y = 0;
	vec.z = constBz;
	return vec;
  }

  getInds(r, z, inds, rho, zeta);

  ind = serialize(inds[0], inds[1]);
  
  int max = rDim*zDim*3;
  if(ind>=0 && ind<max){
  	vec.x = map[ind];
  	vec.y = map[ind + 1];
  	vec.z = map[ind + 2];
  }else{
  	vec.x = 0.0;
	vec.y = 0.0;
	vec.z = 0.0;
  }

  return vec;
}

D3Vector_t DMagneticFieldMap::getQuick(double x, double y, double z)
{
  double r;

  r = sqrt(x*x + y*y);
  return getQuick(r, z);
}

D3Vector_t DMagneticFieldMap::getQuick(const D3Vector_t &vec)
{
  return getQuick(vec.x, vec.y, vec.z);
}

#endif
