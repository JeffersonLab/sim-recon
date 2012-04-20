// $Id$
//
//    File: DFDCIntersection_factory.cc
// Created: Tue Oct 30 11:24:53 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//


#include "DFDCIntersection_factory.h"
#include "FDC/DFDCGeometry.h"
#include "HDGEOMETRY/DGeometry.h"

//------------------
// init
//------------------
jerror_t DFDCIntersection_factory::init(void)
{
	MAX_DIST2 = 2.0*2.0;

	// Create skeleton of data structure to hold hits
	vector<const DFDCHit*> mt_trkhits;
	vector<vector<const DFDCHit*> > mt_trkhits_by_layer;
	for(int i=0; i<6; i++)mt_trkhits_by_layer.push_back(mt_trkhits);
	for(int i=0; i<4; i++)fdchits_by_package.push_back(mt_trkhits_by_layer);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DFDCIntersection_factory::brun(JEventLoop *loop, int runnumber)
{
  // Get pointer to DGeometry object
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  DGeometry *dgeom  = dapp->GetDGeometry(runnumber);
  if (!dgeom->GetFDCWires(fdcwires)){
    _DBG_<< "FDC geometry not available!" <<endl;
    USE_FDC=false;
  }

  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DFDCIntersection_factory::evnt(JEventLoop *loop, int eventnumber)
{
  if (!USE_FDC) return NOERROR;

	// Clear fdchits_by_package structure
	for(int package=1; package<=4; package++){
		for(int layer=1; layer<=6; layer++){
			fdchits_by_package[package-1][layer-1].clear();
		}
	}

	// Get raw hits
	vector<const DFDCHit*> fdchits;
	loop->Get(fdchits);

	// For events with a very large number of hits, assume
	// we can't reconstruct them so bail early
	// Feb. 8, 2008  D.L.
	if(fdchits.size()>(5.0+5.0+1.0)*25.0*6.0){
		_DBG_<<"Too many hits in FDC! Intersection point reconstruction in FDC bypassed for event "<<loop->GetJEvent().GetEventNumber()<<endl;
		return NOERROR;
	}
	
	// Sort wire hits by package and layer
	for(unsigned int i=0; i<fdchits.size(); i++){
		const DFDCHit *fdchit = fdchits[i];
		if(fdchit->type != 0)continue; // filter out cathode hits
		if(fdchit->gLayer<1 || fdchit->gLayer>24){
			_DBG_<<"FDC gLayer out of range! ("<<fdchit->gLayer<<" is not between 1 and 24)"<<endl;
			continue;
		}

		int package = (fdchit->gLayer-1)/6 + 1;
		fdchits_by_package[package-1][(fdchit->gLayer-1)%6].push_back(fdchit);
	}

	// Loop over packages, creating intersection points
	for(unsigned int package=1; package<=4; package++){
		//MakeIntersectionPoints(fdchits_by_package[package-1]);
		MakeRestrictedIntersectionPoints(fdchits_by_package[package-1]);
	}

	return NOERROR;
}

//------------------
// MakeIntersectionPoints
//------------------
void DFDCIntersection_factory::MakeIntersectionPoints(vector<vector<const DFDCHit*> >&hits_by_layer)
{
	// hits_by_layer should contain all of the wire hits in a single
	// FDC package. They are stored in set of nested vectors with the
	// outer vector being the list of layers (always 6 in length)
	// and the inner vector the hits within the layer.

	// Loop over layers
	for(unsigned int i=0; i<hits_by_layer.size()-1; i++){
		vector<const DFDCHit*> &layer1 = hits_by_layer[i];
		vector<const DFDCHit*> &layer2 = hits_by_layer[i+1];
		
		FindIntersections(layer1, layer2, _data);
	}
}

//------------------
// MakeRestrictedIntersectionPoints
//------------------
void DFDCIntersection_factory::MakeRestrictedIntersectionPoints(vector<vector<const DFDCHit*> >&hits_by_layer)
{
	// hits_by_layer should contain all of the wire hits in a single
	// FDC package. They are stored in a set of nested vectors with the
	// outer vector being the list of layers (always 6 in length)
	// and the inner vector the hits within the layer.
	
	// This method will first make a list of all the intersection points
	// between adjacent layers in the package. It will then make a list
	// which intersection points to keep by finding those that meet
	// one of the following criteria:
	//
	// 1. For one of the wires used, this is the only intersection
	//    point at this z location.
	//
	// 2. The intersection point has a corresponding intersection
	//    point in an adjacent plane (i.e. a triple wire coincidence)
	// 
	// Note that there is a class of hits this method fails to accept
	// at the moment. This is when 2 wire planes have 2 hits each leading
	// to 4 intersections, only 2 of which are real. If the adjacent
	// plane has only one hit that is in coincidence with only
	// 1 of the points, the other "true" point could be inferred. At this
	// point, the second true point will not be kept in these cases
	// (or similar ones with more tracks).

	// First, get all of the intersection points
	vector<vector<DFDCIntersection *> > intersections(hits_by_layer.size()-1);
	for(unsigned int i=0; i<hits_by_layer.size()-1; i++){
		
		// Find intersection points between the current layer and next one downstream
		FindIntersections(hits_by_layer[i], hits_by_layer[i+1], intersections[i]);
	}
	
	// Here we need to make 2 lists for every wire hit. The first is of
	// intersection points made from intersections with its upstream
	// neighbor and the second from intersections with its downstream
	// neighbor.
	map<const DFDCWire*, vector<DFDCIntersection*> > upstr;
	map<const DFDCWire*, vector<DFDCIntersection*> > dnstr;
	for(unsigned int i=0; i<intersections.size(); i++){
		vector<DFDCIntersection *> &layer = intersections[i];
		
		for(unsigned int j=0; j<layer.size(); j++){
			DFDCIntersection *fdcinter = layer[j];
			// Note here that wire2 is in the downstream layer for the hit
			// and wire1 is the upstream layer. Therefore, this intersection
			// should be upstream for wire2 and downstream for wire 1.
			upstr[fdcinter->wire2].push_back(fdcinter);
			dnstr[fdcinter->wire1].push_back(fdcinter);
		}
	}
	
	// Keep track of the "good" intersection objects in a map.
	map<DFDCIntersection *, bool> intersects_to_keep;
	
	// Loop over both the upstr and dnstr lists to find wires
	// with only one hit in either its upstream or downstream
	// plane.
	map<const DFDCWire*, vector<DFDCIntersection*> >::iterator iter;
	for(iter=upstr.begin(); iter!=upstr.end(); iter++){
		vector<DFDCIntersection*> &hits = iter->second;
		if(hits.size()==1)intersects_to_keep[hits[0]] = true;
	}
	for(iter=dnstr.begin(); iter!=dnstr.end(); iter++){
		vector<DFDCIntersection*> &hits = iter->second;
		if(hits.size()==1)intersects_to_keep[hits[0]] = true;
	}

	// Find triple wire plane coincidences
	// Loop over intersection layers
	for(unsigned int i=0; i<intersections.size()-1; i++){
		
		// Loop over intersections in this intersection layer
		vector<DFDCIntersection *> &layer1 = intersections[i];
		for(unsigned int j=0; j<layer1.size(); j++){
			DFDCIntersection *int1 = layer1[j];			

			// Loop over hits in downstream layer
			vector<DFDCIntersection *> &layer2 = intersections[i+1];
			for(unsigned int k=0; k<layer2.size(); k++){
				DFDCIntersection *int2 = layer2[k];
				// Only interested in hits that share a wire
				if(int1->wire2 != int2->wire1)continue;
			
				// Find distance between hits
				double dx = int1->pos.X() - int2->pos.X();
				double dy = int1->pos.Y() - int2->pos.Y();
				double dist2 = dx*dx + dy*dy;
				
				if(dist2<MAX_DIST2){
					intersects_to_keep[int1] = true;
					intersects_to_keep[int2] = true;
				}
			}
		}
	}
	
	// Loop over all intersection points. Any that are in "intersects_to_keep"
	// copy to _data. All others, delete.
	for(unsigned int i=0; i<intersections.size(); i++){
		for(unsigned int j=0; j<intersections[i].size(); j++){
			DFDCIntersection *intersection = intersections[i][j];
			bool keep = intersects_to_keep[intersection];
			if(keep){
				_data.push_back(intersection);
			}else{
				delete intersection;
			}
		}
	}
}

//------------------
// FindIntersections
//------------------
void DFDCIntersection_factory::FindIntersections(vector<const DFDCHit*> &layer1, vector<const DFDCHit*> &layer2, vector<DFDCIntersection*> &intersections)
{
	// Loop over hits in layer1
	for(unsigned int j=0; j<layer1.size(); j++){
		const DFDCHit* hit1 = layer1[j];
		const DFDCWire *wire1 =fdcwires[hit1->gLayer-1][hit1->element-1];
		if(!wire1){
			_DBG_<<"No wire for layer="<<hit1->gLayer<<" wire="<<hit1->element<<" !!"<<endl;
			continue;
		}
		DVector2 a(wire1->origin.X(), wire1->origin.Y());
		DVector2 b(wire1->udir.X(), wire1->udir.Y());
		b /= b.Mod();
		
		// Loop over hits in layer2
		for(unsigned int k=0; k<layer2.size(); k++){
			const DFDCHit* hit2 = layer2[k];
			const DFDCWire *wire2 =fdcwires[hit2->gLayer-1][hit2->element-1];
			if(!wire2){
				_DBG_<<"No wire for layer="<<hit2->gLayer<<" wire="<<hit2->element<<" !!"<<endl;
				continue;
			}
			DVector2 c(wire2->origin.X(), wire2->origin.Y());
			DVector2 d(wire2->udir.X(), wire2->udir.Y());
			d /= d.Mod();
			
			// Find intersection of the projections of wires 1 and 2 on the X/Y plane
			DVector2 e = a - c;
			double gamma = d*b;
			double beta = (d*e - gamma*(b*e))/(1.0-gamma*gamma);
			if(!finite(beta))continue; // wires must be parallel
			if(fabs(beta) > wire2->L/2.0)continue; // intersection is past end of wire 
			double alpha = beta*gamma - b*e;
			if(fabs(alpha) > wire1->L/2.0)continue; // intersection is past end of wire 
			
			DVector2 x = c + beta*d; // intersection point in X/Y
			
			// Add the intersection point
			DFDCIntersection *fdcint = new DFDCIntersection;
			fdcint->hit1 = hit1;
			fdcint->hit2 = hit2;
			fdcint->wire1 = wire1;
			fdcint->wire2 = wire2;
			fdcint->pos.SetXYZ(x.X(), x.Y(), (wire1->origin.Z()+wire2->origin.Z())/2.0);
			
			intersections.push_back(fdcint);
		}
	}
}

