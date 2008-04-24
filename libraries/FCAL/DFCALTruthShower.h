// $Id$
//
//    File: DFCALTruthShower.h
// Created: Wed Jan  4 14:43:05 EST 2006
// Creator: davidl (on Linux jlabl1.jlab.org 2.4.21-37.ELsmp i686)
//

#ifndef _DFCALTruthShower_
#define _DFCALTruthShower_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DFCALTruthShower:public JObject{
	public:
		JOBJECT_PUBLIC(DFCALTruthShower);
		DFCALTruthShower(){}
		
		DFCALTruthShower( oid_t id,
						float x, float y, float z,
						float E, float t,
						bool primary, int track) :
			JObject( id ),
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

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "x(cm)", "%3.1f", x());
			AddString(items, "y(cm)", "%3.1f", y());
			AddString(items, "z(cm)", "%3.1f", z());
			AddString(items, "E(MeV)", "%3.3f", E()*1000.0);
			AddString(items, "t(ns)", "%3.1f", t());
			AddString(items, "primary", "%d", primary());
			AddString(items, "track", "%d", track());
		}

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

