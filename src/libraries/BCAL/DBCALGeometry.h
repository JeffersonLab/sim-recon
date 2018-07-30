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

class DBCALGeometry : public JObject {
  
public:
  
  JOBJECT_PUBLIC( DBCALGeometry );
  
  DBCALGeometry(int runnumber);
  
  enum End { kUpstream, kDownstream };
  
  // Methods to access and initialize the private variables
  float GetBCAL_inner_rad() const;
  const float* GetBCAL_radii() const;
  float GetBCAL_center() const;
  float GetBCAL_length() const;
  float GetBCAL_phi_shift() const;
  
  float GetBCAL_outer_rad() const { return BCALOUTERRAD; }
  float GetBCAL_middle_rad() const { return BCALMIDRAD; }
  float GetBCAL_middle_cell() const { return BCALMID; }
  float *GetBCAL_cell_radii() { return &(m_radius[0]); }
  
  // as-built geometry
  float GetBCAL_Nmodules() const { return NBCALMODS; }
  float GetBCAL_Nlayers() const { return NBCALLAYERS; }
  float GetBCAL_Nsectors() const { return NBCALSECTORS; }

  float GetBCAL_NInnerLayers() const { return NBCALLAYSIN; }
  float GetBCAL_NOuterLayers() const { return NBCALLAYSOUT; }
  float GetBCAL_NInnerSectors() const { return NBCALSECSIN; }
  float GetBCAL_NOuterSectors() const { return NBCALSECSOUT; }
  // define these for completeness, but they aren't used outside of this class
  // right now, so comment them out
  //vector<float> GetBCAL_NSummedInnerLayers() const { return NSUMLAYSIN; }
  //vector<float> GetBCAL_NSummedOuterLayers() const { return NSUMLAYSOUT; }

  
  // nominal effective velocity
  float GetBCAL_c_effective() const { return C_EFFECTIVE; }

  // nominal attenuation length
  float GetBCAL_attenutation_length() const { return ATTEN_LENGTH; }

  ///these functions are about encoding/decoding module/layer/sector info in a cellId
  int cellId( int module, int layer, int sector ) const;  ///< This object can be used for the SiPM ID or for the fADC ID since they are defined in the same way (4 bits for sector then 4 bits for layer then 8 bits for module.)
  int module( int cellId ) const;  ///< This method can be used for the SiPM ID or for the fADC ID since they are defined in the same way
  int layer( int cellId ) const;   ///< This method can be used for the SiPM ID or for the fADC ID since they are defined in the same way
  int sector( int cellId ) const;  ///< This method can be used for the SiPM ID or for the fADC ID since they are defined in the same way

  ///these functions are about finding which readout cell contains a specific SiPM cell
  int fADC_layer( int SiPM_cellId ) const;
  int fADC_sector( int SiPM_cellId ) const;
  int fADCId( int module, int SiPM_layer, int SiPM_sector ) const;
  int NSiPMs(int fADCId) const;

  ///these functions are about the physical location and dimensions of a readout cell
  float phi( int fADC_cellId ) const;
  float phiSize( int fADC_cellId ) const;  
  float r( int fADC_cellId ) const;
  float rSize( int fADC_cellId ) const;

  ///these are missing functions that fill in some previous gaps.
  int fADCcellId_rphi( float r, float phi ) const;    ///< Method to get the fADC cell ID from an (R, phi) combination.\n  R in cm and phi in radians.
  int getglobalchannelnumber(int module, int layer, int sector, int end) const;  ///< Return a BCAL channel number, in order of significance (module, layer, sector, end).
  int getendchannelnumber(int module, int layer, int sector) const;  ///< Return a channel number for either end, in order of significance (module, layer, sector).
  int getglobalsector(int module, int sector) const;
  int getsector(int globalsector) const;
  int getmodule(int globalsector) const;

private:

  DBCALGeometry();       // forbid the default constructor
  void Initialize(int runnumber);   // this is old, but keep it around for now, make sure no one else can call it

  // as-built geometry
  const int NBCALMODS=48;         ///< number of modules
  const int NBCALLAYERS=4;         ///< number of layers in a module
  const int NBCALSECTORS=4;         ///< number of sectors in a module

  //the distinction between inner layers and outer layers is important, since only the inner layers have TDC readout
  const int NBCALLAYSIN=3;        ///< number of readout layers in inner BCAL (first 6 SiPM layers)
  const int NBCALLAYSOUT=1;       ///< number of readout layers in outer BCAL (outer 4 SiPM layers)

  // On each module there is a 10x4 (r/phi) array of SiPMs
  // 1.2.3.4 summing configuration - This is used in the BCAL as built
  vector<int> NSUMLAYSIN  = {1,2,3};    ///< number of radial SiPM layers summed for digitization in each inner readout layer
  vector<int> NSUMLAYSOUT = {4};   ///< number of radial SiPM layers summed for digitization in each outer readout layer
  const int NBCALSECSIN=4;        ///<number of sectors in inner region
  const int NBCALSECSOUT=4;      ///<number of sectors in outer region
  // the following are completely deprecated
  //const int NSUMSECSIN=1;         ///< for the inner layers, the number of SiPM that will be summed in the azimuthal direction
  //const int NSUMSECSOUT=1;        ///< for the outer layer(s), the number of SiPM that will be summed in the azimuthal direction
  //const int NBCALSECSIN=4/NSUMSECSIN;        ///<number of sectors in inner region
  //const int NBCALSECSOUT=4/NSUMSECSOUT;      ///<number of sectors in outer region

  float BCALINNERRAD=0.;       ///< innner radius of BCAL in cm
  float fADC_radius[5] = {};   ///< BCAL layer radii (4 layers total)
  float GLOBAL_CENTER=0.;      ///< center of BCAL in gloobal coordinate system
  float BCALFIBERLENGTH=0.;    ///< BCAL Scintilator fiber lenth in cm
  float BCAL_PHI_SHIFT=0.;     ///< overall phi roation of BCAL in radians

  // Enter the index of the SiPM that designates the first
  // (counting radially outward) of the outer cells (default 7)
  const int BCALMID=7;         ///< first outer layer (default 7)
  float BCALMIDRAD = m_radius[BCALMID-1];    ///< mid radius of BCAL in cm (boundary between inner and outer layers)
  float BCALOUTERRAD=86.17;     ///< outer radius of BCAL in cm

  float C_EFFECTIVE=16.75;      ///< speed of light in fibers 
  float ATTEN_LENGTH=520.;     ///< attenuation length

  float m_radius[11] = { 64.3, 
				  66.3,
				  68.3,
				  70.3,
				  72.3,
				  74.3,
				  76.3,
				  78.77,
				  81.24,
				  83.70,
				  86.17};
  

};

#endif // _DBCALGeometry_
