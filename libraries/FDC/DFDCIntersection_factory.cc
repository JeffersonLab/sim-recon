// $Id$
//
//    File: DFDCIntersection_factory.cc
// Created: Tue Oct 30 11:24:53 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//


#include "DFDCIntersection_factory.h"
#include "FDC/DFDCGeometry.h"
#include "DVector2.h"

//------------------
// evnt
//------------------
jerror_t DFDCIntersection_factory::evnt(JEventLoop *loop, int eventnumber)
{
	// Get raw hits
	vector<const DFDCHit*> fdchits;
	loop->Get(fdchits);

	// Create skeleton of data structure to hold hits
	vector<vector<vector<const DFDCHit*> > > fdchits_by_package; ///< fdchits_by_package[package][layer][hit]
	vector<const DFDCHit*> mt_trkhits;
	vector<vector<const DFDCHit*> > mt_trkhits_by_layer;
	for(int i=0; i<6; i++)mt_trkhits_by_layer.push_back(mt_trkhits);
	for(int i=0; i<4; i++)fdchits_by_package.push_back(mt_trkhits_by_layer);
	
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
		MakeIntersectionPoints(fdchits_by_package[package-1]);
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
		for(unsigned int j=0; j<layer1.size(); j++){
			const DFDCHit* hit1 = layer1[j];
			const DFDCWire *wire1 = DFDCGeometry::GetDFDCWire(hit1->gLayer, hit1->element);
			if(!wire1){
				_DBG_<<"No wire for layer="<<hit1->gLayer<<" wire="<<hit1->element<<" !!"<<endl;
				continue;
			}
			DVector2 a(wire1->origin.X(), wire1->origin.Y());
			DVector2 b(wire1->udir.X(), wire1->udir.Y());
			b /= b.Mod();
			
			for(unsigned int k=0; k<layer2.size(); k++){
				const DFDCHit* hit2 = layer2[k];
				const DFDCWire *wire2 = DFDCGeometry::GetDFDCWire(hit2->gLayer, hit2->element);
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

				_data.push_back(fdcint);
			}
		}
	}
}

//------------------
// toString
//------------------
const string DFDCIntersection_factory::toString(void)
{
#if 0
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	// Put the class specific code to produce nicely formatted ASCII here.
	// The DFactory_base class has several methods defined to help. They
	// rely on positions of colons (:) in the header. Here's an example:
	//
	printheader("row:    x:     y:");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DFDCIntersection *myDFDCIntersection = _data[i];
	
		printnewrow();
		printcol("%d",	i);
//		printcol("%1.3f",	myDFDCIntersection->x);
//		printcol("%3.2f",	myDFDCIntersection->y);
		printrow();
	}
#endif

	return _table;
}
