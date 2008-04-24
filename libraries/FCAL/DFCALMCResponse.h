// $Id$
//
//    File: DFCALMCResponse.h
// Created: Wed Aug 24 10:05:31 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#ifndef _DFCALMCResponse_
#define _DFCALMCResponse_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"
#include "DFCALTruthShower.h"

class DFCALMCResponse : public JObject {

public:

	JOBJECT_PUBLIC(DFCALMCResponse);
	DFCALMCResponse(){}
	
	DFCALMCResponse( oid_t id, int channel, 
					 double E, double t ) :
		JObject( id ),
		m_channel( channel ),
		m_E( E ),
		m_t( t ) {}

	int channel() const { return m_channel; }
	
	double E() const { return m_E; }
	double t() const { return m_t; }
	
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "channel", "%d", channel());
			AddString(items, "E(MeV)", "%3.1f", E());
			AddString(items, "t(ns)", "%2.3f", t());
		}

private:
		
	int m_channel;

	double m_E;
	double m_t;
};

#endif // _DFCALMCResponse_

