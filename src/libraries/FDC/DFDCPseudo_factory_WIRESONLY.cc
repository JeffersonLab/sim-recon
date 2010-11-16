// $Id$
//
//    File: DFDCPseudo_factory_WIRESONLY.cc
// Created: Fri Nov  9 09:57:12 EST 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//

#include <cmath>
using namespace std;

#include "DVector2.h"
#include "FDC/DFDCIntersection.h"

#include "DFDCPseudo_factory_WIRESONLY.h"

//------------------
// evnt
//------------------
jerror_t DFDCPseudo_factory_WIRESONLY::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DFDCIntersection*> fdcintersections;
	loop->Get(fdcintersections);

	for(unsigned int i=0; i<fdcintersections.size(); i++){
		const DFDCIntersection *fdcintersection = fdcintersections[i];
		
		// Here we add 2 pseudo points for each intersection point
		MakePseudo(fdcintersection->hit1, fdcintersection->wire1, fdcintersection->pos);
		MakePseudo(fdcintersection->hit2, fdcintersection->wire2, fdcintersection->pos);
		
	}

	return NOERROR;
}

//------------------
// MakePseudo
//------------------
void DFDCPseudo_factory_WIRESONLY::MakePseudo(const DFDCHit *hit, const DFDCWire *wire, const DVector3 &pos)
{
	DFDCPseudo *pseudo = new DFDCPseudo;
	
	DVector2 R(pos.X(), pos.Y());
	DVector2 udir(wire->udir.X(), wire->udir.Y());
	DVector2 a(wire->origin.X(), wire->origin.Y());
	
	pseudo->w = R*a;
	pseudo->dw = 1.0/sqrt(12.0); // cm
	pseudo->s = (R-a)*udir;
	pseudo->ds = 1.0/sqrt(12.0); // cm
	pseudo->wire = wire;
	pseudo->time = hit->t;
	pseudo->status = 1;  // 1 external hit used to find intersection
	pseudo->x = pos.X();
	pseudo->y = pos.Y();
	
	// Intialize covariance matrix in w/s coordinate system
	DMatrix cov(2,2);
	cov(0,0) = pseudo->dw*pseudo->dw;	cov(1,0) = 0.0;
	cov(0,1) = 0.0;							cov(1,1) = pseudo->ds*pseudo->ds;
	
	// Build rotation matrix for rotating into x/y coordinate system
	DMatrix Rot(2,2);
	double cos_angle = cos(wire->angle);
	double sin_angle = sin(wire->angle);
	Rot(0,0) = cos_angle;		Rot(1,0) = sin_angle;
	Rot(0,1) = -sin_angle;		Rot(1,1) = cos_angle;
	DMatrix RotT(DMatrix::kTransposed, Rot);
	DMatrix RotCov=RotT*cov*Rot;
	
	pseudo->covxx=RotCov(0,0);
	pseudo->covxy=RotCov(1,0);
	pseudo->covyy=RotCov(1,1);

	_data.push_back(pseudo);
}

