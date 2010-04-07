// $Id$
//
//    File: DVertexCalculator.cc
// Created: Wed Apr  7 10:56:12 EDT 2010
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//

#include "DVertexCalculator.h"

//---------------------------------
// DVertexCalculator    (Constructor)
//---------------------------------
DVertexCalculator::DVertexCalculator()
{
	// Initialize to reasonable defaults
	target_length = 0.0;
	target_z = 0.0;
}

//---------------------------------
// InitializeGeometry
//---------------------------------
void DVertexCalculator::InitializeGeometry(DGeometry *geom)
{
	// Get Target parameters from XML
	if(geom){
		geom->GetTargetLength(target_length);
		geom->GetTargetZ(target_z);
	}
}

//------------------
// FindVertex
//------------------
void DVertexCalculator::FindVertex(vector<const DKinematicData*> &trks, DVertex *vertex)
{
	// Initialize vertex parameters to center of target as fall-back
	vertex->x.SetXYZ(0.0, 0.0, target_z);
	vertex->cov.ResizeTo(3,3);
	vertex->beamline_used = true;
	
	// If there are no tracks assume center of target which is already set so just return
	if(trks.size()==0)return;
	
	// Simply average POCA to beamline for all tracks
	vertex->x.SetXYZ(0.0, 0.0, 0.0);
	for(unsigned int i=0; i<trks.size(); i++){
		vertex->x += trks[i]->position();
		vertex->AddAssociatedObject(trks[i]);
	}
	vertex->x *= (1.0/(double)trks.size());
	
	// Make sure the vertex is finite
	if(!finite(vertex->x.Mag2())){
		vertex->x.SetXYZ(0.0, 0.0, target_z);
		// Remove any associated objects we added
		for(unsigned int i=0; i<trks.size(); i++){
			vertex->RemoveAssociatedObject(trks[i]);
		}
	}

	return;
}

//---------------------------------
// ~DVertexCalculator    (Destructor)
//---------------------------------
DVertexCalculator::~DVertexCalculator()
{

}
