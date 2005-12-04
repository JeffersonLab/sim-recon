// $Id$
//
//    File: DHDDMForwardShower.h
// Created: Mon Aug 29 15:14:08 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#ifndef _DHDDMForwardShower_
#define _DHDDMForwardShower_

#include "DObject.h"
#include "DFactory.h"

class DHDDMForwardShower : public DObject{

public:

	HDCLASSDEF(DHDDMForwardShower);
	DHDDMForwardShower(){}
	
	DHDDMForwardShower( identifier_t id,
						float x, float y,
						float E, float t ) :
		DObject( id ),
		m_x( x ),
		m_y( y ),
		m_E( E ),
		m_t( t ) {}
	
	float x() const { return m_x; }
	float y() const { return m_y; }
	float t() const { return m_t; }
	float E() const { return m_E; }
	
private:
		
	float m_x;
	float m_y;
	float m_E;
	float m_t;
	
};

#endif // _DHDDMForwardShower_

