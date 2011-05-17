//    File: DBCALGeometry.h
// Created: Fri Nov 26 15:10:51 CST 2010
// Creator: dwbennet

#include "DBCALGeometry.h"

int DBCALGeometry::NBCALMODS = 48;

// On each module there is a 10x4 (r/phi) array of SiPMs -- this factory
// allows two different "rectangular" groupings, an inner (1) and an outer (2)

#ifdef BCAL_SUM_CELL

// The configuration below has 3x1 (r/phi) summed cells in the inner
// and 2x2 summed cells in the outer region
int DBCALGeometry::NBCALLAYS1 =  2;
int DBCALGeometry::NBCALLAYS2 =  2; 
int DBCALGeometry::NBCALSECS1 =  4; 
int DBCALGeometry::NBCALSECS2 =  2;
// end of summing configuration

#else

// The configuration below has no summing of SiPMs
int DBCALGeometry::NBCALLAYS1 =  6;
int DBCALGeometry::NBCALLAYS2 =  4; 
int DBCALGeometry::NBCALSECS1 =  4; 
int DBCALGeometry::NBCALSECS2 =  4;
// end of no summing configuration

#endif

int DBCALGeometry::NSUMLAYS1 = 6/DBCALGeometry::NBCALLAYS1;
int DBCALGeometry::NSUMLAYS2 = 4/DBCALGeometry::NBCALLAYS2;
int DBCALGeometry::NSUMSECS1 = 4/DBCALGeometry::NBCALSECS1;
int DBCALGeometry::NSUMSECS2 = 4/DBCALGeometry::NBCALSECS2;

// Enter the index of the SiPM that designates the first
// (counting radially outward) of the outer cells (default 7)
int DBCALGeometry::BCALMID = 7;

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

DBCALGeometry::DBCALGeometry()  
{
  BCALMIDRAD = m_radius[BCALMID-1];
  /// End if groupings do not evenly divide SiPM cells
  if(  (BCALMID-1) % NBCALLAYS1 != 0 
       || (11-BCALMID) % NBCALLAYS2 != 0
       ||  4 % NBCALSECS1 != 0
       ||  4 % NBCALSECS2 != 0)
    {
      std::cout<<"ERROR: Bad BCAL fADC groupings, do not evenly divide cells";
      assert (false);
    }
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
DBCALGeometry::fADC_layer( int cellId ) {
  
  int cell_layer = DBCALGeometry::layer( cellId );
  int fADC_layer = cell_layer;
#ifdef BCAL_SUM_CELL
	if(cell_layer<DBCALGeometry::BCALMID){
		// inner
		fADC_layer = 1+ (cell_layer-1)/NSUMLAYS1;
	}else{
		// outer
		int max_inner_layer = (DBCALGeometry::BCALMID-1)/NSUMLAYS1;
		fADC_layer = 1 + max_inner_layer + (cell_layer-DBCALGeometry::BCALMID)/NSUMLAYS2;
	}
#endif
  return fADC_layer;
}

//--------------
// fADC_sector
//--------------
int
DBCALGeometry::fADC_sector( int cellId ) {

  int cell_layer  = DBCALGeometry::layer( cellId );
  int cell_sector = DBCALGeometry::sector( cellId );
  int fADC_sector = cell_sector;
#ifdef BCAL_SUM_CELL
	if(cell_layer<DBCALGeometry::BCALMID){
		// inner
		fADC_sector = 1 + (cell_sector-1)/NSUMSECS1;
	}else{
		// outer
		fADC_sector = 1 + (cell_sector-1)/NSUMSECS2;
	}
#endif
  return fADC_sector;
}

//--------------
// fADCId
//--------------
int
DBCALGeometry::fADCId( int module, int layer, int sector ) {
	// This is used to return the readout channel ID which may
	// differ from the cellID if summing is implemented.
	//
	// (5/12/2011 DL)
#ifdef BCAL_SUM_CELL
	if(layer<DBCALGeometry::BCALMID){
		// inner
		layer = 1+ (layer-1)/NSUMLAYS1;
		sector = 1 + (sector-1)/NSUMSECS1;
	}else{
		// outer
		int max_inner_layer = (DBCALGeometry::BCALMID-1)/NSUMLAYS1;
		layer = 1 + max_inner_layer + (layer-DBCALGeometry::BCALMID)/NSUMLAYS2;
		sector = 1 + (sector-1)/NSUMSECS2;
	}
#endif
	return cellId(module, layer, sector);
}

//--------------
// r
//--------------
float
DBCALGeometry::r( int cell ) {
 
  int lay = layer( cell );
  
  float innerRad, outerRad;
 
  if( lay <= NBCALLAYS1 ) { // we're in the inner region
   
    // the size of the layer in SiPM units
    int laySize = ( BCALMID - 1 ) / NBCALLAYS1;
    
    innerRad = m_radius[ laySize * ( lay - 1 ) ];
    
    // index below is same as:  laySize * ( lay - 1 ) + laySize
    outerRad = m_radius[ laySize * lay ];
  }
  else{
    
    // find the size of the outer layer in SiPM units
    int laySize = ( 11 - BCALMID ) / NBCALLAYS2;
    
    // subtract the number of inner layers to get
    // layer number within the outer section
    lay -= NBCALLAYS1;
    
    innerRad = m_radius[ BCALMID - 1 + laySize * ( lay - 1 ) ];
    outerRad = m_radius[ BCALMID - 1 + laySize * lay ];
  }
  
  return 0.5 * ( innerRad + outerRad );
}

//--------------
// rSize
//--------------
float
DBCALGeometry::rSize( int cell ) {
  
  int lay = layer( cell );
  
  float innerRad, outerRad;
  
  if( lay <= NBCALLAYS1 ) { // in the inner region
    
    // the size of the layer in SiPM units
    int laySize = ( BCALMID - 1 ) / NBCALLAYS1;
    
    innerRad = m_radius[ laySize * ( lay - 1 ) ];
    
    // index below is same as:  laySize * ( lay - 1 ) + laySize
    outerRad = m_radius[ laySize * lay ];
  }
  else{ // in the outer region
    
    // find the size of the outer layer in SiPM units
    int laySize = ( 11 - BCALMID ) / NBCALLAYS2;
    
    // subtract the number of inner layers to get
    // layer number within the outer section
    lay -= NBCALLAYS1;
    
    innerRad = m_radius[ BCALMID - 1 + laySize * ( lay - 1 ) ];
    outerRad = m_radius[ BCALMID - 1 + laySize * lay ];
  }
  
  return ( outerRad - innerRad );
}

//--------------
// phi
//--------------
float
DBCALGeometry::phi( int cell ) {
 
  int lay = layer( cell );
  int sect = sector( cell );
  
  float sectSize;
  
  if( lay <= NBCALLAYS1 ) { 
    
    sect += NBCALSECS1 * ( module( cell ) - 1 );
    sectSize = 2 * PI / ( NBCALMODS * NBCALSECS1 );
  }
  else {
    
    sect += NBCALSECS2 * ( module( cell ) - 1 );
    sectSize = 2 * PI / ( NBCALMODS * NBCALSECS2 );
  }
  
  return sectSize * ( sect - 0.5 );
}

//--------------
// phiSize
//--------------
float
DBCALGeometry::phiSize( int cell ) {
  
  int lay = layer( cell );
  
  float sectSize;
  
  if( lay <= NBCALLAYS1 ) { 
    
    sectSize = 2 * 3.1416 / ( NBCALMODS * NBCALSECS1 );
  }
  else {
    
    sectSize = 2 * 3.1416 / ( NBCALMODS * NBCALSECS2 );
  }
  
  return sectSize;
}

