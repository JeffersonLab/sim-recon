// $Id$
//
//    File: DFCALTruthShower.h
// Created: Wed Jan  4 14:43:05 EST 2006
// Creator: davidl (on Linux jlabl1.jlab.org 2.4.21-37.ELsmp i686)
//

#ifndef _DFCALTruthShower_
#define _DFCALTruthShower_

#include "DObject.h"
#include "DFactory.h"

class DFCALTruthShower:public DObject{
	public:
		HDCLASSDEF(DFCALTruthShower);
		DFCALTruthShower(){}
		
		DFCALTruthShower( oid_t id,
						float x, float y, float z,
						float E, float t,
						bool primary, int track) :
			DObject( id ),
			m_x( x ),
			m_y( y ),
			m_z( z ),
			m_E( E ),
			m_t( t ),
			m_primary(primary),
			m_track(track) {}
	
		float x() const { return m_x; }
		float y() const { return m_y; }
		float z() const { return m_z; }
		float t() const { return m_t; }
		float E() const { return m_E; }
		bool primary() const { return m_primary; }
		int track() const { return m_track; }

		void Serialize(JILStream&) const;
		void Deserialize(JILStream&);

	private:
		
		float m_x;
		float m_y;
		float m_z;
		float m_E;
		float m_t;
		bool m_primary;
		int m_track;
};

#endif // _DFCALTruthShower_

