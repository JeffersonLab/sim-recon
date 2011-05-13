// $Id$
//
//    File: DBCALGeometry.h
// Created: Thu Nov 17 15:10:51 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#ifndef _DBCALGeometry_
#define _DBCALGeometry_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

// create a single number channel id which is useful in algorithms
// if M L S are module layer sector the bit map looks like:
// MMMM MMMM LLLL SSSS

#define MODULE_SHIFT 8
#define MODULE_MASK 0xFF00
#define LAYER_SHIFT 4
#define LAYER_MASK 0x00F0
#define SECTOR_SHIFT 0
#define SECTOR_MASK 0X000F

#define PI 3.1416

// with this set one will utilize the default summing -- see DBCALGeometry.cc
// for a description of summed geometry.
//
// Author's note:  this is seen as a temporary feature to study effects of
// two different summing schemes.  A preprocessor macro is not the best way
// to change the functionality of code, but mcsmear doesn't use factories to
// provide objects so it is parameters aren't helpful.

//#define BCAL_SUM_CELL

class DBCALGeometry : public JObject {
  
public:
  
  JOBJECT_PUBLIC( DBCALGeometry );
  
  DBCALGeometry();
  
  enum End { kUpstream, kDownstream };
  
  static int NBCALMODS;         ///> number of modules
  static int NBCALLAYS1;        ///> number of layers in first 10 ccm 
  static int NBCALLAYS2;        ///> number of layers in last  15 cm 
  static int NBCALSECS1;        ///> number of sectors in first 10cm of Mod 
  static int NBCALSECS2;        ///> number of sectors in last 15cm of Mod 
  static int NSUMLAYS1;        ///> number of layers summed for digitization in inner region 
  static int NSUMLAYS2;        ///> number of layers summed for digitization in outer region 
  static int NSUMSECS1;        ///> number of sectors summed for digitization in inner region 
  static int NSUMSECS2;        ///> number of sectors summed for digitization in outer region 
  static float BCALINNERRAD;    ///> innner radius of BCAL in cm
  static int BCALMID;         ///> first outer layer (default 7)
  static float m_radius[11];
  float BCALMIDRAD;      ///> mid radius of BCAL in cm
  static float BCALOUTERRAD;    ///> outer radius of BCAL in cm
  static float BCALFIBERLENGTH; ///> BCAL Scintilator fiber lenth in cm
  static float GLOBAL_CENTER;  ///> center of BCAL in gloobal coordinate system
  
  static float ATTEN_LENGTH;   ///> attenuation length
  static float C_EFFECTIVE;    ///> speed of light in fibers 
  
  static bool summingOn() {
    
#ifdef BCAL_SUM_CELL
    return true;
#else
    return false;
#endif
  }
  
  static int module( int cellId );  
  static int layer( int cellId );
  static int sector( int cellId );
  static int fADC_layer( int cellId );
  static int fADC_sector( int cellId );

  static int cellId( int module, int layer, int sector );
  static int fADCId( int module, int layer, int sector );

  static float phi( int cellId );
  static float phiSize( int cellId );
  
  static float r( int cellId );
  static float rSize( int cellId );
  
};

#endif // _DBCALGeometry_
