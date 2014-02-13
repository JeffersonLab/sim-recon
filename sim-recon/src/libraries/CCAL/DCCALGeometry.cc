// $Id$
//
//    File: DCCALGeometry.cc
// Created: Tue Nov 30 15:42:41 EST 2010
// Creator: davidl (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include <cassert>
#include <math.h>
using namespace std;

#include "DCCALGeometry.h"
#include "DVector2.h"

//---------------------------------
// DCCALGeometry    (Constructor)
//---------------------------------
DCCALGeometry::DCCALGeometry() : 
m_numActiveBlocks( 0 )
{
	double innerRadius = ( kCCALBeamHoleSize - 1 ) / 2. * blockSize() * sqrt(2.);

	// inflate the innner radius by 1% to for "safe" comparison
	innerRadius *= 1.01;
	
	for( int row = 0; row < kCCALBlocksTall; row++ ){
		for( int col = 0; col < kCCALBlocksWide; col++ ){
			
			// transform to beam axis
			m_positionOnFace[row][col] = 
			   DVector2(  ( (double)col - kCCALMidBlock +0.5 ) * blockSize(),
					     ( (double)row - kCCALMidBlock +0.5 ) * blockSize() );
			
			m_activeBlock[row][col] = true;
				
			// build the "channel map"
			m_channelNumber[row][col] = m_numActiveBlocks;
			m_row[m_numActiveBlocks] = row;
			m_column[m_numActiveBlocks] = col;

			m_numActiveBlocks++;
		}
	}
}

bool
DCCALGeometry::isBlockActive( int row, int column ) const
{
	// I'm inserting these lines to effectively disable the
	// two assert calls below. They are causing all programs
	// (hd_dump, hdview) to exit, even when I'm not interested
	// in the FCAL. This does not fix the underlying problem
	// of why we're getting invalid row/column values.
	// 12/13/05  DL
	if( row < 0 ||  row >= kCCALBlocksTall )return false;
	if( column < 0 ||  column >= kCCALBlocksWide )return false;

	assert(    row >= 0 &&    row < kCCALBlocksTall );
	assert( column >= 0 && column < kCCALBlocksWide );
	
	return m_activeBlock[row][column];	
}

int
DCCALGeometry::row( float y ) const 
{	
	return static_cast<int>( y / blockSize() + kCCALMidBlock );
}

int
DCCALGeometry::column( float x ) const 
{	
	return static_cast<int>( x / blockSize() + kCCALMidBlock );
}

DVector2
DCCALGeometry::positionOnFace( int row, int column ) const
{ 
	assert(    row >= 0 &&    row < kCCALBlocksTall );
	assert( column >= 0 && column < kCCALBlocksWide );
	
	return m_positionOnFace[row][column]; 
}

DVector2
DCCALGeometry::positionOnFace( int channel ) const
{
	assert( channel >= 0 && channel < m_numActiveBlocks );
	
	return positionOnFace( m_row[channel], m_column[channel] );
}

int
DCCALGeometry::channel( int row, int column ) const
{
	if( isBlockActive( row, column ) ){
		
		return m_channelNumber[row][column]; 
	}
	else{
		
		cerr << "ERROR: request for channel number of inactive block!" << endl;
		return -1;
	}
}
