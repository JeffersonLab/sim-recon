// $Id$
//
//    File: DCCALTruthShower.h
// Created: Tue Nov 30 15:02:26 EST 2010
// Creator: davidl (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DCCALTruthShower_
#define _DCCALTruthShower_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DCCALTruthShower:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DCCALTruthShower);
		
		DCCALTruthShower(){}
		
		DCCALTruthShower( oid_t id,
				  float x, float y, float z,
				  float px, float py, float pz,
				  float E, float t,
				  bool primary, int track, int type,
				  int itrack=0) :
		  JObject( id ),
		  m_x( x ),
		  m_y( y ),
		  m_z( z ),
		  m_px( px ),
		  m_py( py ),
		  m_pz( pz ),
		  m_E( E ),
		  m_t( t ),
		  m_primary(primary),
		  m_track(track),
		  m_type(type),
		  m_itrack(itrack){}
	
		float x() const { return m_x; }
		float y() const { return m_y; }
		float z() const { return m_z; }
		float t() const { return m_t; }
		float px() const { return m_px; }
		float py() const { return m_py; }
		float pz() const { return m_pz; }
		float E() const { return m_E; }
		bool primary() const { return m_primary; }
		int track() const { return m_track; }
		int type() const { return m_type; }
		int itrack() const { return m_itrack; }

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "x(cm)", "%3.1f", x());
			AddString(items, "y(cm)", "%3.1f", y());
			AddString(items, "z(cm)", "%3.1f", z());
			AddString(items, "px(MeV/c)", "%3.1f", px()*1000);
			AddString(items, "py(MeV/c)", "%3.1f", py()*1000);
			AddString(items, "pz(MeV/c)", "%3.1f", pz()*1000);
			AddString(items, "E(MeV)", "%3.3f", E()*1000.0);
			AddString(items, "t(ns)", "%3.1f", t());
			AddString(items, "primary", "%d", primary());
			AddString(items, "track", "%d", track());
			AddString(items, "type", "%d", type());
			AddString(items, "itrack", "%d", itrack());
		}

	private:
		
		float m_x;
		float m_y;
		float m_z;
		float m_px;
		float m_py;
		float m_pz;
		float m_E;
		float m_t;
		bool m_primary;
		int m_track;
		int m_type;
		int m_itrack;
};

#endif // _DCCALTruthShower_

