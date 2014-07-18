//    File: DBCALGeometry.h
// Created: Fri Nov 26 15:10:51 CST 2010
// Creator: dwbennet

/// MMD: Looking at the cdoe for the DBCALGeometry::phi( int fADC_cell ) method,
/// I have to conclude that sectors are numbered 1 through 4.  
/// Judging by what else I've seen, I have to conclude that layers are numbered 1 through 4 or 1 through 10.  
/// This is now the standard that I will use.


#include <cmath>
#include "DBCALGeometry.h"

#include <HDGEOMETRY/DGeometry.h>

// On each module there is a 10x4 (r/phi) array of SiPMs

#ifdef BCAL_SUM_CELL

// The configuration below has 3x1 (r/phi) summed cells in the inner
// and 2x2 summed cells in the outer region
//int DBCALGeometry::NSUMSECSIN = 1;
//int DBCALGeometry::NSUMSECSOUT = 2;
//int DBCALGeometry::NSUMLAYSIN[] = {3,3};
//int DBCALGeometry::NSUMLAYSOUT[] = {2,2};

// 1.2.3.4 summing configuration
int DBCALGeometry::NSUMSECSIN = 1;
int DBCALGeometry::NSUMSECSOUT = 1;
int DBCALGeometry::NSUMLAYSIN[] = {1,2,3};
int DBCALGeometry::NSUMLAYSOUT[] = {4};

// end of summing configuration

#else

// The configuration below has no summing of SiPMs
int DBCALGeometry::NSUMSECSIN = 1;
int DBCALGeometry::NSUMSECSOUT = 1;
int DBCALGeometry::NSUMLAYSIN[] = {1,1,1,1,1,1};
int DBCALGeometry::NSUMLAYSOUT[] = {1,1,1,1};
// end of no summing configuration

#endif

int DBCALGeometry::NBCALSECSIN = 4/DBCALGeometry::NSUMSECSIN;
int DBCALGeometry::NBCALSECSOUT = 4/DBCALGeometry::NSUMSECSOUT;
float DBCALGeometry::BCAL_PHI_SHIFT = 0.0; // will be overwritten in constructor

float DBCALGeometry::BCALINNERRAD = 64.3;   
float DBCALGeometry::BCALOUTERRAD = 86.17;
float DBCALGeometry::BCALFIBERLENGTH = 390.0;
float DBCALGeometry::GLOBAL_CENTER = 212;
  
float DBCALGeometry::ATTEN_LENGTH = 300.;
float DBCALGeometry::C_EFFECTIVE  = 16.75;

float DBCALGeometry::m_radius[] = { 64.3, 
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

float DBCALGeometry::fADC_radius[] = { 64.3, 
				  66.3,
				  70.3,
				  76.3,
				  86.17};

float DBCALGeometry::BCALMIDRAD = DBCALGeometry::m_radius[DBCALGeometry::BCALMID-1];

DBCALGeometry::DBCALGeometry()  
{
  /// End if groupings do not evenly divide SiPM cells
  bool goodGeometry=true;

  if (NSUMSECSIN <= 0) goodGeometry=false;
  if (4 % NSUMSECSIN != 0) goodGeometry=false;
  if (NSUMSECSOUT <= 0) goodGeometry=false;
  if (4 % NSUMSECSOUT != 0) goodGeometry=false;

  int totalLayersIn=0;
  for (int i=0;i<NBCALLAYSIN;i++) {
    if (NSUMLAYSIN[i] <= 0) goodGeometry=false;
    totalLayersIn += NSUMLAYSIN[i];
  }

  int totalLayersOut=0;
  for (int i=0;i<NBCALLAYSOUT;i++) {
    if (NSUMLAYSOUT[i] <= 0) goodGeometry=false;
    totalLayersOut += NSUMLAYSOUT[i];
  }

  if (totalLayersIn != BCALMID-1) goodGeometry=false;
  if (totalLayersOut != 10-(BCALMID-1)) goodGeometry=false;

  if(!goodGeometry) {
      std::cout<<"ERROR: Bad BCAL fADC groupings, do not evenly divide cells";
      assert (false);
  }
  
  //Get pointer to DGeometry object
  DApplication* dapp=dynamic_cast<DApplication*>(japp);
  const DGeometry *dgeom  = dapp->GetDGeometry(9999);

  // Get overall phi shift of BCAL
  double my_BCAL_PHI_SHIFT;
  dgeom->GetBCALPhiShift(my_BCAL_PHI_SHIFT);
  BCAL_PHI_SHIFT = (float)(my_BCAL_PHI_SHIFT*M_PI/180.0);  // convert to radians

}

//--------------
// module
//--------------
int
DBCALGeometry::module( int cellId ) {
  
  return ( cellId & MODULE_MASK ) >> MODULE_SHIFT;
}

//--------------
// layer
//--------------
int
DBCALGeometry::layer( int cellId ) {
  
  return ( cellId & LAYER_MASK ) >> LAYER_SHIFT;
}

//--------------
// sector
//--------------
int
DBCALGeometry::sector( int cellId ) {
  
  return ( cellId & SECTOR_MASK ) >> SECTOR_SHIFT;
}

//--------------
// cellId
//--------------
int
DBCALGeometry::cellId( int module, int layer, int sector ) {
  
  return ( ( module << MODULE_SHIFT ) | 
           ( layer << LAYER_SHIFT   ) | 
           ( sector << SECTOR_SHIFT ) );
}

//--------------
// fADC_layer
//--------------
int
DBCALGeometry::fADC_layer( int SiPM_cellId ) {
  
  int cell_layer = DBCALGeometry::layer( SiPM_cellId );
  int fADC_layer = 0;
  int tally=0;
  if (cell_layer < BCALMID) {
    for (int i=0;i<NBCALLAYSIN;i++) {
      tally+=NSUMLAYSIN[i];
      if (cell_layer <= tally) {
        fADC_layer=i+1;
        break;
      }
    }
  } else {
    tally=BCALMID-1;
    for (int i=0;i<NBCALLAYSOUT;i++) {
      tally+=NSUMLAYSOUT[i];
      if (cell_layer <= tally) {
        fADC_layer=NBCALLAYSIN+i+1;
        break;
      }
    }
  }
  return fADC_layer;
}

//--------------
// fADC_sector
//--------------
int
DBCALGeometry::fADC_sector( int SiPM_cellId ) {

  int cell_layer = DBCALGeometry::layer( SiPM_cellId );
  int cell_sector = DBCALGeometry::sector( SiPM_cellId );
  int fADC_sector;

  if (cell_layer < BCALMID) {
    fADC_sector = 1 + (cell_sector-1)/NSUMSECSIN;
  } else {
    fADC_sector = 1 + (cell_sector-1)/NSUMSECSOUT;
  }

  return fADC_sector;
}

//--------------
// fADCId
//--------------
int
DBCALGeometry::fADCId( int module, int SiPM_layer, int SiPM_sector ) {
  // This is used to return the readout channel ID which may
  // differ from the cellID if summing is implemented.
  //
  // (5/12/2011 DL)

  int SiPM_cell = cellId(module, SiPM_layer, SiPM_sector);

  int fADC_lay = fADC_layer(SiPM_cell);
  int fADC_sect = fADC_sector(SiPM_cell);

  return cellId(module, fADC_lay, fADC_sect);
}
  
//--------------
// NSiPMs
//--------------
int
DBCALGeometry::NSiPMs(int fADCId)
{
	/// Return the number of SiPMs summed for the given fADCId
	int fadc_lay = layer(fADCId);

	if(fadc_lay<1 || fadc_lay>(NBCALLAYSOUT+NBCALLAYSIN))return 0;

	if(fadc_lay <= NBCALLAYSIN){
		// inner
		return NSUMLAYSIN[fadc_lay-1]*NSUMSECSIN;
	}else{
		// outer
		return NSUMLAYSOUT[fadc_lay-NBCALLAYSIN-1]*NSUMSECSOUT;
	}
}

//--------------
// r
//--------------
float
DBCALGeometry::r( int fADC_cell ) {

  int fADC_lay = layer( fADC_cell );

  float innerRad;
  float outerRad;

  if (fADC_lay <= NBCALLAYSIN) {
    innerRad=m_radius[0];
    outerRad=m_radius[0];
    int innerIndex=0;
    for (int i=0;i<fADC_lay;i++) {
      innerRad=outerRad;
      outerRad=m_radius[innerIndex+NSUMLAYSIN[i]];
      innerIndex=innerIndex+NSUMLAYSIN[i];
    }
  } else {
    innerRad=m_radius[BCALMID-1];
    outerRad=m_radius[BCALMID-1];
    int innerIndex=BCALMID-1;
    for (int i=0;i < (fADC_lay-NBCALLAYSIN);i++) {
      innerRad=outerRad;
      outerRad=m_radius[innerIndex+NSUMLAYSOUT[i]];
      innerIndex=innerIndex+NSUMLAYSOUT[i];
    }
  }
  
  return 0.5 * ( innerRad + outerRad );
}

//--------------
// rSize
//--------------
float
DBCALGeometry::rSize( int fADC_cell ) {

  int fADC_lay = layer( fADC_cell );

  float innerRad;
  float outerRad;

  if (fADC_lay <= NBCALLAYSIN) {
    innerRad=m_radius[0];
    outerRad=m_radius[0];
    int innerIndex=0;
    for (int i=0;i<fADC_lay;i++) {
      innerRad=outerRad;
      outerRad=m_radius[innerIndex+NSUMLAYSIN[i]];
      innerIndex=innerIndex+NSUMLAYSIN[i];
    }
  } else {
    innerRad=m_radius[BCALMID-1];
    outerRad=m_radius[BCALMID-1];
    int innerIndex=BCALMID-1;
    for (int i=0;i < (fADC_lay-NBCALLAYSIN);i++) {
      innerRad=outerRad;
      outerRad=m_radius[innerIndex+NSUMLAYSOUT[i]];
      innerIndex=innerIndex+NSUMLAYSOUT[i];
    }
  }
  
  return ( outerRad - innerRad );
}

//--------------
// phi
//--------------
float
DBCALGeometry::phi( int fADC_cell ) {

  int fADC_lay = layer( fADC_cell );
  int fADC_sect = sector( fADC_cell );
  
  float sectSize;

  int nSects;
  if (fADC_lay <= NBCALLAYSIN) {
    nSects = 4/NSUMSECSIN;
  } else {
    nSects = 4/NSUMSECSOUT;
  }

  fADC_sect += nSects * ( module( fADC_cell ) - 1 );
  sectSize = 2 * PI / (NBCALMODS*nSects);
  
  float my_phi = sectSize * ( (float)fADC_sect - 0.5 );
  my_phi += BCAL_PHI_SHIFT - 2.0*sectSize; // adjust for center of module and overall BCAL shift
  
  return my_phi;
}

//--------------
// phiSize
//--------------
float
DBCALGeometry::phiSize( int fADC_cell ) {

  int fADC_lay = layer( fADC_cell );

  int nSects;
  if (fADC_lay <= NBCALLAYSIN) {
    nSects = 4/NSUMSECSIN;
  } else {
    nSects = 4/NSUMSECSOUT;
  }

  float sectSize = 2 * PI / ( NBCALMODS * nSects );
  
  return sectSize;
}


//--------------
/// fADCcellId_rphi
//--------------
/// Method to get the fADC cell ID from an (R, phi) combination.  R in cm and phi in radians.
int
DBCALGeometry::fADCcellId_rphi( float r, float phi ) {
  int fADC_cellId = 0;
  int SiPM_layer = 0;

  if (r < BCALINNERRAD) return 0;
  else if (r > BCALOUTERRAD) return 0;

  float modulephiSize = (2 * PI) / 48;
  float sectorphiSize = modulephiSize / 4;
  float phi_nooffset = (phi + 2.0*sectorphiSize);
  if (phi_nooffset < 0) return 0;
  else if (phi_nooffset > 2 * PI) return 0;

  if (r < m_radius[1]) SiPM_layer = 1; 
  else if (r < m_radius[2]) SiPM_layer = 2; 
  else if (r < m_radius[3]) SiPM_layer = 3; 
  else if (r < m_radius[4]) SiPM_layer = 4; 
  else if (r < m_radius[5]) SiPM_layer = 5; 
  else if (r < m_radius[6]) SiPM_layer = 6; 
  else if (r < m_radius[7]) SiPM_layer = 7; 
  else if (r < m_radius[8]) SiPM_layer = 8; 
  else if (r < m_radius[9]) SiPM_layer = 9; 
  else  SiPM_layer = 10; 

  float modulefloat = 1 + phi_nooffset / modulephiSize;
  int module = (int)modulefloat;
  float sectorfloat = 1 + ((phi_nooffset - (module-1)*modulephiSize) / sectorphiSize);
  int sector = (int)sectorfloat;
 
  fADC_cellId = DBCALGeometry::fADCId(module, SiPM_layer, sector);  
  return fADC_cellId;
}
 

int DBCALGeometry::getglobalsector(int module, int sector) {
  if (module==0 || sector==0) return 0;
  else return (module-1)*4 + sector;
}
int DBCALGeometry::getsector(int globalsector) {
  if (globalsector<=0) return 0;
  int sector = globalsector%4;
  if (sector==0) sector=4;
  return sector;
}
int DBCALGeometry::getmodule(int globalsector) {
  if (globalsector<=0) return 0;
  else return (globalsector-1)/4+1;
}
