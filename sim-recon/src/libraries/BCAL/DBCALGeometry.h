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
// provide objects so its parameters aren't helpful.

#define BCAL_SUM_CELL

class DBCALGeometry : public JObject {
  
public:
  
  JOBJECT_PUBLIC( DBCALGeometry );
  
  DBCALGeometry();
  
  enum End { kUpstream, kDownstream };
  
  static const int NBCALMODS=48;         ///> number of modules

  //the distinction between inner layers and outer layers is important, since only the inner layers have TDC readout
#ifdef BCAL_SUM_CELL
  static const int NBCALLAYSIN=3;        ///> number of readout layers in inner BCAL (first 6 SiPM layers)
  static const int NBCALLAYSOUT=1;       ///> number of readout layers in outer BCAL (outer 4 SiPM layers)
#else
  static const int NBCALLAYSIN=6;
  static const int NBCALLAYSOUT=4;
#endif
  static int NSUMLAYSIN[NBCALLAYSIN];        ///> number of radial SiPM layers summed for digitization in each inner readout layer
  static int NSUMLAYSOUT[NBCALLAYSOUT];        ///> number of radial SiPM layers summed for digitization in each outer readout layer
  static int NSUMSECSIN;        ///> for the inner layers, the number of SiPM that will be summed in the azimuthal direction
  static int NSUMSECSOUT;        ///> for the outer layer(s), the number of SiPM that will be summed in the azimuthal direction
  static int NBCALSECSIN;   ///>number of sectors in inner region
  static int NBCALSECSOUT;  ///>number of sectors in outer region
  static float BCALINNERRAD;    ///> innner radius of BCAL in cm

  // Enter the index of the SiPM that designates the first
  // (counting radially outward) of the outer cells (default 7)
  static const int BCALMID=7;         ///> first outer layer (default 7)

  static float m_radius[11];
  static float BCALMIDRAD;     ///> mid radius of BCAL in cm (boundary between inner and outer layers)
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

  //these functions are about encoding/decoding module/layer/sector info in a cellId
  static int cellId( int module, int layer, int sector );  
  static int module( int cellId );  
  static int layer( int cellId );
  static int sector( int cellId );

  //these functions are about finding which readout cell contains a specific SiPM cell
  static int fADC_layer( int SiPM_cellId );
  static int fADC_sector( int SiPM_cellId );
  static int fADCId( int module, int SiPM_layer, int SiPM_sector );
  static int NSiPMs(int fADCId);

  //these functions are about the physical location and dimensions of a readout cell
  static float phi( int fADC_cellId );
  static float phiSize( int fADC_cellId );  
  static float r( int fADC_cellId );
  static float rSize( int fADC_cellId );
  
};

#endif // _DBCALGeometry_
