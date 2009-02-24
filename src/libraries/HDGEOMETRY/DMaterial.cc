// $Id$
//
//    File: DMaterial.cc
// Created: Fri Apr 25 15:44:12 EDT 2008
// Creator: davidl (on Darwin swire-d95.jlab.org 8.11.1 i386)
//

#include <sstream>
using namespace std;

#include "DMaterial.h"

//---------------------------------
// DMaterial    (Constructor)
//---------------------------------
DMaterial::DMaterial(string &name, double A, double Z, double density, double radlen)
{
	this->kname = name;
	this->kA = A;
	this->kZ = Z;
	this->kdensity = density;
	this->kXo = radlen;
	
	last_ptype = Unknown;
	last_p = -1.0;
	last_dEdx = -1.0;
}

//---------------------------------
// toString
//---------------------------------
string DMaterial::toString(void)
{
	stringstream ss;
	ss<<"DMaterial: "<<kname<<endl;
	ss<<"============================"<<endl;
	ss<<"name: "<<kname<<endl;
	ss<<"A: "<<kA<<" g/mol"<<endl;
	ss<<"Z: "<<kZ<<endl;
	ss<<"density: "<<kdensity<<" g/cm^3"<<endl;
	ss<<"radlen: "<<kXo<<" g/cm^2"<<endl;
	ss<<endl;
	
	return ss.str();
}

//---------------------------------
// GetdEdx
//---------------------------------
double DMaterial::GetdEdx(Particle_t ptype, double p)
{

	return 0.0;
}
