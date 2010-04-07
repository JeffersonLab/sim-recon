// $Id$
//
//    File: DVertexCalculator.h
// Created: Wed Apr  7 10:56:12 EDT 2010
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//

#ifndef _DVertexCalculator_
#define _DVertexCalculator_

#include <JANA/jerror.h>

#include <PID/DVertex.h>
class DKinematicData;

class DVertexCalculator{
	public:
		DVertexCalculator();
		virtual ~DVertexCalculator();
		
		void InitializeGeometry(DGeometry *geom);
		void FindVertex(vector<const DKinematicData*> &trks, DVertex *vertex);
		
	private:

		double target_length;
		double target_z;

};

#endif // _DVertexCalculator_

