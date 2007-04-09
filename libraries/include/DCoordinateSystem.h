// $Id$
//
//    File: DCoordinateSystem.h
// Created: Thu Nov 16 10:15:12 EST 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.8.0 powerpc)
//

#ifndef _DCoordinateSystem_
#define _DCoordinateSystem_

#include <DVector3.h>

class DCoordinateSystem{
	public:
		DCoordinateSystem(){}
		virtual ~DCoordinateSystem(){}
		
		DVector3 origin;	//< x,y,z-coordinate of origin in lab coordinates
		DVector3 sdir;		//< direction of s-axis in lab coordinates
		DVector3 tdir;		//< direction of t-axis in lab coordinates
		DVector3 udir;		//< direction of u-axis in lab coordinates
		double L;			//< length of wire (if this represents a wire)
		
		inline void ToLab(double &x, double &y, double &z);
		inline void FromLab(double &x, double &y, double &z);
		inline void ToLab(DVector3 &pos);
		inline void FromLab(DVector3 &pos);
};

//---------------------------------
// ToLab
//---------------------------------
inline void DCoordinateSystem::ToLab(double &x, double &y, double &z)
{
	/// Transform the given vector from this coordinate
	/// system into the lab coordinate system.
	DVector3 pos(x,y,z);
	ToLab(pos);
	x= pos.x();
	y= pos.y();
	z= pos.z();
}

//---------------------------------
// FromLab
//---------------------------------
inline void DCoordinateSystem::FromLab(double &x, double &y, double &z)
{
	/// Transform the given vector from the lab coordinate
	/// system into this coordinate system.
	DVector3 pos(x,y,z);
	FromLab(pos);
	x= pos.x();
	y= pos.y();
	z= pos.z();
}

//---------------------------------
// ToLab
//---------------------------------
inline void DCoordinateSystem::ToLab(DVector3 &pos)
{
	/// Transform the given vector from this coordinate
	/// system into the lab coordinate system.
	pos = origin + pos.x()*sdir + pos.y()*tdir + pos.z()*udir;
}

//---------------------------------
// FromLab
//---------------------------------
inline void DCoordinateSystem::FromLab(DVector3 &pos)
{
	/// Transform the given vector from the lab coordinate
	/// system into this coordinate system.
	pos -= origin;
	pos.SetXYZ(pos.Dot(sdir), pos.Dot(tdir), pos.Dot(udir));
}


#endif // _DCoordinateSystem_

