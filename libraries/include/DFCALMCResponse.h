// $Id$
//
//    File: DFCALMCResponse.h
// Created: Wed Aug 24 10:05:31 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#ifndef _DFCALMCResponse_
#define _DFCALMCResponse_

#include "DFactory.h"
#include "DObject.h"
#include "DHDDMForwardShower.h"

class DFCALMCResponse : public DObject {

public:

	HDCLASSDEF(DFCALMCResponse);
	DFCALMCResponse(){}
	
	DFCALMCResponse( identifier_t id, int channel, 
					 double E, double t ) :
		DObject( id ),
		m_channel( channel ),
		m_E( E ),
		m_t( t ) {}

	int channel() const { return m_channel; }
	
	double E() const { return m_E; }
	double t() const { return m_t; }
	
private:
		
	int m_channel;

	double m_E;
	double m_t;
};

#endif // _DFCALMCResponse_

