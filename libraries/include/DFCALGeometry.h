// $Id$
//
//    File: DFCALGeometry.h
// Created: Wed Aug 24 10:09:27 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#ifndef _DFCALGeometry_
#define _DFCALGeometry_

#include "DFactory.h"
#include "DObject.h"
#include "TVector2.h"
#include "units.h"

/*
class Vector2 {
	
public:
	
	Vector2() : m_x( 0 ), m_y( 0 ) {}
	Vector2( double x, double y ) : m_x( x ), m_y( y ) {}
	
	double Mod( void ) const { return sqrt( m_x * m_x +
											m_y * m_y ); }
	
private:
	
	double m_x;
	double m_y;
};
*/

class DFCALGeometry : public DObject {

public:
	
	HDCLASSDEF(DFCALGeometry);
	
	static const int kBlocksWide      = 53;
	static const int kBlocksTall      = 53;
	static const int kMaxChannels     = kBlocksWide * kBlocksTall;
	static const int kBeamHoleSize    =  3;

	DFCALGeometry();
	~DFCALGeometry();

	static double blockSize() { return 4*k_cm; }
	static double radius() { return 1.08*k_m; }
	
	bool isBlockActive( int row, int column ) const;
	int  numActiveBlocks() const { return m_numActiveBlocks; }
	
	TVector2 positionOnFace( int row, int column ) const;
	TVector2 positionOnFace( int channel ) const;
	
	int channel( int row, int column ) const;

	int row   ( int channel ) const { return m_row[channel];    }
	int column( int channel ) const { return m_column[channel]; }
	
	// get row and column from x and y positions
	int row   ( float y ) const;
	int column( float x ) const;
	
private:

	bool   m_activeBlock[kBlocksTall][kBlocksWide];
	TVector2 m_positionOnFace[kBlocksTall][kBlocksWide];

	int    m_channelNumber[kBlocksTall][kBlocksWide];
	int    m_row[kMaxChannels];
	int    m_column[kMaxChannels];
	
	int    m_numActiveBlocks;
};

#endif // _DFCALGeometry_

