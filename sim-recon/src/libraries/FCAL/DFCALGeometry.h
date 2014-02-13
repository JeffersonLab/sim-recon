// $Id$
//
//    File: DFCALGeometry.h
// Created: Wed Aug 24 10:09:27 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#ifndef _DFCALGeometry_
#define _DFCALGeometry_

#include <JANA/JFactory.h>
#include <JANA/JObject.h>
using namespace jana;

#include "DVector2.h"
#include "units.h"

class DFCALGeometry : public JObject {

#define kBlocksWide 59
#define kBlocksTall 59
#define kMaxChannels kBlocksWide * kBlocksTall
// Do not forget to adjust below formula if number of blocks chage in any direction:
//   this is now used to convert from row/col to coordiantes y/x and back - MK
#define kMidBlock (kBlocksWide-1)/2 			
#define kBeamHoleSize 3

public:
	
	JOBJECT_PUBLIC(DFCALGeometry);

	//static const int kBlocksWide      = 53;
	//static const int kBlocksTall      = 53;
	//static const int kMaxChannels     = kBlocksWide * kBlocksTall;
	//static const int kBeamHoleSize    =  3;

	DFCALGeometry();
	~DFCALGeometry(){}

	static double blockSize()  { return 4*k_cm; }
	static double radius()  { return 1.2*k_m; }
	static double blockLength()  { return 45.0*k_cm; }
	static double fcalFaceZ()  { return 625.3*k_cm; }

        static double fcalMidplane() { return fcalFaceZ() + 0.5 * blockLength() ; } 
	
	bool isBlockActive( int row, int column ) const;
	int  numActiveBlocks() const { return m_numActiveBlocks; }

	
	DVector2 positionOnFace( int row, int column ) const;
	DVector2 positionOnFace( int channel ) const;
	
	int channel( int row, int column ) const;

	int row   ( int channel ) const { return m_row[channel];    }
	int column( int channel ) const { return m_column[channel]; }
	
	// get row and column from x and y positions
	int row   ( float y ) const;
	int column( float x ) const;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "kBlocksWide", "%d", kBlocksWide);
			AddString(items, "kBlocksTall", "%d", kBlocksTall);
			AddString(items, "kMaxChannels", "%d", kMaxChannels);
			AddString(items, "kBeamHoleSize", "%2.3f", kBeamHoleSize);
		}
	
private:

	bool   m_activeBlock[kBlocksTall][kBlocksWide];
	DVector2 m_positionOnFace[kBlocksTall][kBlocksWide];

	int    m_channelNumber[kBlocksTall][kBlocksWide];
	int    m_row[kMaxChannels];
	int    m_column[kMaxChannels];
	
	int    m_numActiveBlocks;
};

#endif // _DFCALGeometry_

