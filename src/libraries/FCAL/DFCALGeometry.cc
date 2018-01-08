// $Id: DFCALGeometry.cc 19049 2015-07-16 18:31:31Z staylor $
//
//    File: DFCALGeometry.cc
// Created: Wed Aug 24 10:18:56 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#include <cassert>
#include <math.h>
using namespace std;

#include "DFCALGeometry.h"
#include "DVector2.h"

//---------------------------------
// DFCALGeometry    (Constructor)
//---------------------------------
DFCALGeometry::DFCALGeometry() 
{
  m_numActiveBlocks=0;
  //  double innerRadius = ( kBeamHoleSize - 1 ) / 2. * blockSize() * sqrt(2.);

  // inflate the innner radius by 1% to for "safe" comparison
  // innerRadius *= 1.01;
  //double innerRadius=5.7;
  int calor=0;
  for( int row = 0; row < kBlocksTall; row++ ){
    for( int col = 0; col < kBlocksWide; col++ ){
      // transform to beam axis
      m_positionOnFace[row][col][calor] = 
	DVector2(  ( col - kMidBlock ) * blockSize(calor),
		   ( row - kMidBlock ) * blockSize(calor) );
      float x=m_positionOnFace[row][col][calor].X();
      float y=m_positionOnFace[row][col][calor].Y();
 	
      double thisRadius = m_positionOnFace[row][col][calor].Mod();
      
      if(thisRadius < radius() && (fabs(x)>30.09 || fabs(y)>30.09)){
	m_activeBlock[row][col][calor] = true;
	
	// build the "channel map"
	m_channelNumber[row][col][calor] = m_numActiveBlocks;
	m_row[m_numActiveBlocks] = row;
	m_column[m_numActiveBlocks] =col;
	
	m_numActiveBlocks++;
      }
      else{
	
	m_activeBlock[row][col][calor] = false;
      }
    }
  }
  calor=1;
  for( int row = 0; row < kInnerBlocksTall; row++ ){
    for( int col = 0; col < kInnerBlocksWide; col++ ){
      
      // transform to beam axis
      m_positionOnFace[row][col][calor] = 
	DVector2(  ( col - kInnerMidBlock ) * blockSize(calor),
		   ( row - kInnerMidBlock ) * blockSize(calor) );
      // Carve out a hole for the insert
      float x=m_positionOnFace[row][col][calor].X();
      float y=m_positionOnFace[row][col][calor].Y();
      
      //printf("r %d c %d x %f y %f\n",row,col,x,y);

      if (fabs(x)>5.0 || fabs(y)>5.0){
	m_activeBlock[row][col][calor] = true;
	
	// build the "channel map"
	m_channelNumber[row][col][calor] = m_numActiveBlocks;
	m_row[m_numActiveBlocks] = 100+row;
	m_column[m_numActiveBlocks] =100+col;
	
	m_numActiveBlocks++;
      }
      else{
	
	m_activeBlock[row][col][calor] = false;
      }
    }
  }
}

bool
DFCALGeometry::isBlockActive( int row, int column) const
{
	// I'm inserting these lines to effectively disable the
	// two assert calls below. They are causing all programs
	// (hd_dump, hdview) to exit, even when I'm not interested
	// in the FCAL. This does not fix the underlying problem
	// of why we're getting invalid row/column values.
	// 12/13/05  DL
	//if( row < 0 ||  row >= kBlocksTall )return false;
	//if( column < 0 ||  column >= kBlocksWide )return false;

	// assert(    row >= 0 &&    row < kBlocksTall );
	// assert( column >= 0 && column < kBlocksWide );
	
	int calor=0;
	if (row>=100 && column>=100){
	  row-=100;
	  column-=100;
	  calor=1;
	}
	if (calor==0){
	  if( row < 0 ||  row >= kBlocksTall )return false;
	  if( column < 0 ||  column >= kBlocksWide )return false;
	}
	else{
	  if( row < 0 ||  row >= kInnerBlocksTall )return false;
	  if( column < 0 ||  column >= kInnerBlocksWide )return false;
	}

	return m_activeBlock[row][column][calor];	
}

int
DFCALGeometry::row( float y, int calor ) const 
{
  if (calor==0) return static_cast<int>( y / blockSize(0) + kMidBlock + 0.5);
  return (100+static_cast<int>( y / blockSize(1) + kInnerMidBlock + 0.5));
}

int
DFCALGeometry::column( float x, int calor ) const 
{	
  if (calor==0) return static_cast<int>( x / blockSize(0) + kMidBlock + 0.5);
  return (100+static_cast<int>( x / blockSize(1) + kInnerMidBlock + 0.5));
}

DVector2
DFCALGeometry::positionOnFace( int row, int column) const
{ 
  //	assert(    row >= 0 &&    row < kBlocksTall );
  //	assert( column >= 0 && column < kBlocksWide );
  int calor=0;
  if (row>=100 && column>=100){
    row-=100;
    column-=100;
    calor=1;
  }

	return m_positionOnFace[row][column][calor]; 
}

DVector2 DFCALGeometry::positionOnFace(int channel) const{
  int r=row(channel);
  int c=column(channel);
  int calor=0;
  if (r>=100 && c>=100){
    calor=1;
    r-=100;
    c-=100;
  }
  return m_positionOnFace[r][c][calor];
}


int
DFCALGeometry::channel( int row, int column) const
{
  if( isBlockActive( row, column) ){
    int calor=0;
    if (row>=100 && row>=100){
      calor=1;
      row-=100;
      column-=100;
    }
    return m_channelNumber[row][column][calor]; 
  }
  else{
    
    cerr << "ERROR: request for channel number of inactive block!  row " 
	 << row << " column " <<  column << endl;
    return -1;
  }
}
