//    File: DBCALGeometry.h
// Created: Fri Nov 26 15:10:51 CST 2010
// Creator: dwbennet

/// MMD: Looking at the code for the DBCALGeometry::phi( int fADC_cell ) method,
/// I have to conclude that sectors are numbered 1 through 4.  
/// Judging by what else I've seen, I have to conclude that layers are numbered 1 through 4 or 1 through 10.  
/// This is now the standard that I will use.


#include <cmath>
#include <TMath.h>
#include "DBCALGeometry.h"

#include <HDGEOMETRY/DGeometry.h>


DBCALGeometry::DBCALGeometry(int runnumber)  
{
  /// End if groupings do not evenly divide SiPM cells
  bool goodGeometry=true;

  //if (NSUMSECSIN <= 0) goodGeometry=false;
  //if (4 % NSUMSECSIN != 0) goodGeometry=false;
  //if (NSUMSECSOUT <= 0) goodGeometry=false;
  //if (4 % NSUMSECSOUT != 0) goodGeometry=false;

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

  // Initialize DBCALGeometry variables
  Initialize(runnumber);

}

void
DBCALGeometry::Initialize(int runnumber) {
  //Get pointer to DGeometry object
  DApplication* dapp=dynamic_cast<DApplication*>(japp);
  const DGeometry *dgeom  = dapp->GetDGeometry(runnumber);

  // Get inner rad of BCAL (including the support plate)
  float my_BCALINNERRAD;
  dgeom->GetBCALRmin(my_BCALINNERRAD);
  BCALINNERRAD = my_BCALINNERRAD;

  // Get layer radii (4 layers)
  vector<float> bcal_fadc_radii;
  dgeom->GetBCALfADCRadii(bcal_fadc_radii);
  if(bcal_fadc_radii.size()==5){
  	for(uint32_t i=0; i<bcal_fadc_radii.size(); i++){
		fADC_radius[i] = bcal_fadc_radii[i];
	}
  }
  else{
  	jerr<<"Did not retrieve 5 values for BCAL fADC radii!!!" << endl;
  	exit(-1);
  }

  // Get BCAL Global Center
  float my_GLOBAL_CENTER;
  dgeom->GetBCALCenterZ(my_GLOBAL_CENTER);
  GLOBAL_CENTER = my_GLOBAL_CENTER;	

  // Get BCAL fiber length
  float my_BCALFIBERLENGTH;
  dgeom->GetBCALLength(my_BCALFIBERLENGTH);
  BCALFIBERLENGTH = my_BCALFIBERLENGTH;

  // Get overall phi shift of BCAL
  float my_BCAL_PHI_SHIFT;
  dgeom->GetBCALPhiShift(my_BCAL_PHI_SHIFT);
  BCAL_PHI_SHIFT = my_BCAL_PHI_SHIFT*TMath::Pi()/180.0;  // convert to radians
}

float
DBCALGeometry::GetBCAL_inner_rad() const {
	return BCALINNERRAD;
}

const float*
DBCALGeometry::GetBCAL_radii() const {
	return fADC_radius;
	//return &(fADC_radius[0]);
}

float
DBCALGeometry::GetBCAL_center() const {
	return GLOBAL_CENTER;
}

float
DBCALGeometry::GetBCAL_length() const {
	return BCALFIBERLENGTH;
}

float
DBCALGeometry::GetBCAL_phi_shift() const {
	return BCAL_PHI_SHIFT;
}

//--------------
// module
//--------------
int
DBCALGeometry::module( int cellId ) const {
  
  return ( cellId & MODULE_MASK ) >> MODULE_SHIFT;
}

//--------------
// layer
//--------------
int
DBCALGeometry::layer( int cellId ) const {
  
  return ( cellId & LAYER_MASK ) >> LAYER_SHIFT;
}

//--------------
// sector
//--------------
int
DBCALGeometry::sector( int cellId ) const {
  
  return ( cellId & SECTOR_MASK ) >> SECTOR_SHIFT;
}

//--------------
// cellId
//--------------
int
DBCALGeometry::cellId( int module, int layer, int sector ) const {
  
  return ( ( module << MODULE_SHIFT ) | 
           ( layer << LAYER_SHIFT   ) | 
           ( sector << SECTOR_SHIFT ) );
}

//--------------
// fADC_layer
//--------------
int
DBCALGeometry::fADC_layer( int SiPM_cellId ) const {
  
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
DBCALGeometry::fADC_sector( int SiPM_cellId ) const {

  int cell_layer = DBCALGeometry::layer( SiPM_cellId );
  int cell_sector = DBCALGeometry::sector( SiPM_cellId );
  int fADC_sector;

  if (cell_layer < BCALMID) {
    //fADC_sector = 1 + (cell_sector-1)/NSUMSECSIN;
    fADC_sector = 1 + (cell_sector-1);
  } else {
    //fADC_sector = 1 + (cell_sector-1)/NSUMSECSOUT;
    fADC_sector = 1 + (cell_sector-1);
  }

  return fADC_sector;
}

//--------------
// fADCId
//--------------
int
DBCALGeometry::fADCId( int module, int SiPM_layer, int SiPM_sector ) const {
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
DBCALGeometry::NSiPMs(int fADCId) const
{
	/// Return the number of SiPMs summed for the given fADCId
	int fadc_lay = layer(fADCId);

	if(fadc_lay<1 || fadc_lay>(NBCALLAYSOUT+NBCALLAYSIN))return 0;

	if(fadc_lay <= NBCALLAYSIN){
		// inner
		//return NSUMLAYSIN[fadc_lay-1]*NSUMSECSIN;
		return NSUMLAYSIN[fadc_lay-1];
	}else{
		// outer
		//return NSUMLAYSOUT[fadc_lay-NBCALLAYSIN-1]*NSUMSECSOUT;
		return NSUMLAYSOUT[fadc_lay-NBCALLAYSIN-1];
	}
}

//--------------
// r
//--------------
float
DBCALGeometry::r( int fADC_cell ) const {

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
DBCALGeometry::rSize( int fADC_cell ) const {

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
DBCALGeometry::phi( int fADC_cell ) const {

  // int fADC_lay = layer( fADC_cell );
  // int fADC_sect = sector( fADC_cell );
  
  // float sectSize;

  // int nSects;
  // if (fADC_lay <= NBCALLAYSIN) {
  //   nSects = 4/NSUMSECSIN;
  // } else {
  //   nSects = 4/NSUMSECSOUT;
  // }

  // fADC_sect += nSects * ( module( fADC_cell ) - 1 );
  // sectSize = 2 * M_PI / (NBCALMODS*nSects);
  
  // float my_phi = sectSize * ( (float)fADC_sect - 0.5 );
  // my_phi += BCAL_PHI_SHIFT - 2.0*sectSize; // adjust for center of module and overall BCAL shift
  
  float globsect = getglobalsector(module(fADC_cell), sector(fADC_cell));
  float sectSize = 2. * M_PI / 48. / 4;
  // The first 2 sectors (half of mudule 1) have negative phi
  float my_phi = sectSize * (globsect-2.5);
  if (my_phi < 0) my_phi += 2 * M_PI;
  return my_phi;

}

//--------------
// phiSize
//--------------
float
DBCALGeometry::phiSize( int fADC_cell ) const {

  // int fADC_lay = layer( fADC_cell );

  // int nSects;
  // if (fADC_lay <= NBCALLAYSIN) {
  //   nSects = 4/NSUMSECSIN;
  // } else {
  //   nSects = 4/NSUMSECSOUT;
  // }

  // float sectSize = 2 * M_PI / ( NBCALMODS * nSects );
  
  float sectSize = 2 * M_PI / 48 / 4;

  return sectSize;
}

//--------------
/// fADCcellId_rphi
//--------------
/// Method to get the fADC cell ID from an (R, phi) combination.  R in cm and phi in radians.
int
DBCALGeometry::fADCcellId_rphi( float r, float phi ) const {
  int fADC_cellId = 0;
  int SiPM_layer = 0;

  if (r < BCALINNERRAD) return 0;
  else if (r > BCALOUTERRAD) return 0;

  float modulephiSize = (2 * M_PI) / 48;
  float sectorphiSize = modulephiSize / 4;
  float phi_nooffset = (phi + 2.0*sectorphiSize);
  if (phi_nooffset < 0) return 0;
  else if (phi_nooffset > 2 * M_PI) return 0;

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

int  DBCALGeometry::getglobalchannelnumber(int module, int layer, int sector, int end) const {
  if (module<=0 || layer<=0 || sector<=0) return 0;
  else return (module-1)*32 + (layer-1)*8 + (sector-1)*2 + end + 1;
}

int DBCALGeometry::getendchannelnumber(int module, int layer, int sector) const {
  if (module<=0 || layer<=0 || sector<=0) return 0;
  else return (module-1)*16 + (layer-1)*4 + sector;
} 

int DBCALGeometry::getglobalsector(int module, int sector) const {
  if (module<=0 || sector<=0) return 0;
  else return (module-1)*4 + sector;
}

int DBCALGeometry::getsector(int globalsector) const {
  if (globalsector<=0) return 0;
  int sector = globalsector%4;
  if (sector==0) sector=4;
  return sector;
}

int DBCALGeometry::getmodule(int globalsector) const {
  if (globalsector<=0) return 0;
  else return (globalsector-1)/4+1;
}
