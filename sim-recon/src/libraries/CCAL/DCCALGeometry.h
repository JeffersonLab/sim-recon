// $Id$
//
//    File: DCCALGeometry.h
// Created: Tue Nov 30 15:42:41 EST 2010
// Creator: davidl (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DCCALGeometry_
#define _DCCALGeometry_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

#include "DVector2.h"
#include "units.h"

class DCCALGeometry:public jana::JObject{

#define kCCALBlocksWide 16
#define kCCALBlocksTall 16
#define kCCALMaxChannels kCCALBlocksWide * kCCALBlocksTall
// Do not forget to adjust below formula if number of blocks chage in any direction:
//   this is now used to convert from row/col to coordiantes y/x and back - MK
#define kCCALMidBlock (kCCALBlocksWide)/2 			
#define kCCALBeamHoleSize 2

	public:
		JOBJECT_PUBLIC(DCCALGeometry);
		
		DCCALGeometry();
		~DCCALGeometry(){}

		static double blockSize(){ return 2*k_cm; }
		static double blockLength(){ return 18.0*k_cm; }
		static double ccalFaceZ(){ return 1025.3*k_cm; }
	
		static double ccalMidplane() { return ccalFaceZ() + 0.5 * blockLength() ; } 
	
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
			AddString(items, "kCCALBlocksWide", "%d", kCCALBlocksWide);
			AddString(items, "kCCALBlocksTall", "%d", kCCALBlocksTall);
			AddString(items, "kCCALMaxChannels", "%d", kCCALMaxChannels);
			AddString(items, "kCCALBeamHoleSize", "%2.3f", kCCALBeamHoleSize);
		}
	
	private:

		bool   m_activeBlock[kCCALBlocksTall][kCCALBlocksWide];
		DVector2 m_positionOnFace[kCCALBlocksTall][kCCALBlocksWide];

		int    m_channelNumber[kCCALBlocksTall][kCCALBlocksWide];
		int    m_row[kCCALMaxChannels];
		int    m_column[kCCALMaxChannels];
	
		int    m_numActiveBlocks;
};

#endif // _DCCALGeometry_

