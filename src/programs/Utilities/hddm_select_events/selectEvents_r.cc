// $Id$
//
// Created Oct 10, 2013  David Lawrence

#include "hddm_select_events.h"
#include <HDDM/hddm_r.hpp>
#include "particleType.h"

using namespace std;
using std::distance;

//-----------------------------------------------------------------
// selectEvent_s
//-----------------------------------------------------------------
bool selectEvent_r(int select_type, hddm_r::HDDM *record, int nevents, bool debug){

  // Select Lambda -> p pi- events
  if(select_type==1){

    bool foundLambda = false;
    int  idLambda    = -999;
    int  idProton    = -999;
    bool foundProton = false;
    if(debug) cout << "event: " << nevents << "----------------------------------------------------------------" << endl;

    hddm_r::ReconstructedPhysicsEvent &re = record->getReconstructedPhysicsEvent();
    hddm_r::ReactionList reactions = re.getReactions();
    // Loop over ReactionList
    hddm_r::ReactionList::iterator iter;
    for (iter = reactions.begin(); iter != reactions.end(); ++iter) {
      int ireaction = std::distance(reactions.begin(), iter);
      
      hddm_r::VertexList vertices = iter->getVertices();
      hddm_r::VertexList::iterator iter_vertex;
      // Loop over VertexList
      for (iter_vertex = vertices.begin(); iter_vertex != vertices.end(); ++iter_vertex){
	int ivertex = std::distance(vertices.begin(), iter_vertex);

	hddm_r::ProductList products = iter_vertex->getProducts();
	hddm_r::ProductList::iterator iter_product;
	// Loop over ProductList
	for (iter_product = products.begin(); iter_product != products.end(); ++iter_product){
	  int iparticle = std::distance(products.begin(), iter_product);

	  if(debug) cout << "reaction: " << ireaction << " vertex: " << ivertex << " particle: " << iparticle
			 << " id: " << setw(3) << iter_product->getId() << " parent: " << iter_product->getParentId()
			 << " PDG type: " << setw(7) << iter_product->getPdgtype()
			 << " type: " << ParticleType((Particle_t)PDGtoPType((Particle_t)iter_product->getPdgtype()))
			 << " vertex: (" << iter_vertex->getOrigin().getVx() << ", " << iter_vertex->getOrigin().getVy() << ", " << iter_vertex->getOrigin().getVz() << ")"
			 << endl;
	  
	  // Find Lambda
	  if(iter_product->getPdgtype() == PDGtype(Lambda)){
	    foundLambda = true;
	    idLambda = iter_product->getId();
	  }

	  // Find decay products of Lambda
	  if(foundLambda){
	    // if parent was Lambda
	    if(iter_product->getParentId() == idLambda){
	      // find what type it is
	      if(iter_product->getPdgtype() == PDGtype(Proton)){
		foundProton = true;
		idProton = iter_product->getId();
	      }
	    }
	  }

	  // If we found a proton decaying from Lambda,
	  // save this event
	  if(foundProton){
	    return true;
	  }

	} // end of loop over ProductList
      }	// end of loop over VertexList
    } // end of loop over ReactiionList
  }
  // select Lambda -> p pi-, eta -> gamma gamma events
  else if(select_type==2){

    bool foundLambda = false;
    int  idLambda    = -999;
    int  idProton    = -999;
    bool foundProton = false;

    bool foundEta = false;
    int  idEta    = -999;
    int  idPhoton[2];
    int  nPhoton = 0;

    if(debug) cout << "event: " << nevents << "----------------------------------------------------------------" << endl;

    hddm_r::ReconstructedPhysicsEvent &re = record->getReconstructedPhysicsEvent();
    hddm_r::ReactionList reactions = re.getReactions();
    // Loop over ReactionList
    hddm_r::ReactionList::iterator iter;
    for (iter = reactions.begin(); iter != reactions.end(); ++iter) {
      int ireaction = std::distance(reactions.begin(), iter);
      
      hddm_r::VertexList vertices = iter->getVertices();
      hddm_r::VertexList::iterator iter_vertex;
      // Loop over VertexList
      for (iter_vertex = vertices.begin(); iter_vertex != vertices.end(); ++iter_vertex){
	int ivertex = std::distance(vertices.begin(), iter_vertex);

	hddm_r::ProductList products = iter_vertex->getProducts();
	hddm_r::ProductList::iterator iter_product;
	// Loop over ProductList
	for (iter_product = products.begin(); iter_product != products.end(); ++iter_product){
	  int iparticle = std::distance(products.begin(), iter_product);

	  if(debug) cout << "reaction: " << ireaction << " vertex: " << ivertex << " particle: " << iparticle
			 << " id: " << setw(3) << iter_product->getId() << " parent: " << iter_product->getParentId()
			 << " PDG type: " << setw(7) << iter_product->getPdgtype()
			 << " type: " << ParticleType((Particle_t)PDGtoPType((Particle_t)iter_product->getPdgtype()))
			 << " vertex: (" << iter_vertex->getOrigin().getVx() << ", " << iter_vertex->getOrigin().getVy() << ", " << iter_vertex->getOrigin().getVz() << ")"
			 << endl;
	  
	  // Find Lambda
	  if(iter_product->getPdgtype() == PDGtype(Lambda)){
	    foundLambda = true;
	    idLambda = iter_product->getId();
	  }

	  // Find decay products of Lambda
	  if(foundLambda){
	    // if parent was Lambda
	    if(iter_product->getParentId() == idLambda){
	      // find what type it is
	      if(iter_product->getPdgtype() == PDGtype(Proton)){
		foundProton = true;
		idProton = iter_product->getId();
	      }
	    }
	  }

	  // Find Eta
	  if(iter_product->getPdgtype() == PDGtype(Eta)){
	    foundEta = true;
	    idEta = iter_product->getId();
	  }

	  // Find decay products of Eta
	  if(foundEta){
	    // if parent was Eta
	    if(iter_product->getParentId() == idEta){
	      // find what type it is
	      if(iter_product->getPdgtype() == PDGtype(Gamma)){
		idPhoton[nPhoton] = iter_product->getId();
		nPhoton++;
	      }
	    }
	  }

	  // If we found a proton decaying from Lambda,
	  // save this event
	  if(foundProton && nPhoton==2){
	    return true;
	  }

	} // end of loop over ProductList
      }	// end of loop over VertexList
    } // end of loop over ReactiionList
  }
  //__________________________________________________________________________________________________

  // Select Lambda -> p pi- events, decay length is less than VERTEXDIFFMAX
  else if(select_type==3){
    const double VERTEXDIFFMAX = 10.0;

    bool foundLambda = false;
    int  idLambda    = -999;
    int  idProton    = -999;
    bool foundProton = false;
    double vertexLambda[3];
    double vertexProton[3];
    for(int v=0;v<3;v++){
      vertexLambda[v] = -999;
      vertexProton[v] = 999;
    }
    double vertexdiff = 999;

    hddm_r::ReconstructedPhysicsEvent &re = record->getReconstructedPhysicsEvent();
    hddm_r::ReactionList reactions = re.getReactions();
    // Loop over ReactionList
    hddm_r::ReactionList::iterator iter;
    for (iter = reactions.begin(); iter != reactions.end(); ++iter) {
      int ireaction = std::distance(reactions.begin(), iter);
      
      hddm_r::VertexList vertices = iter->getVertices();
      hddm_r::VertexList::iterator iter_vertex;
      // Loop over VertexList
      for (iter_vertex = vertices.begin(); iter_vertex != vertices.end(); ++iter_vertex){
	int ivertex = std::distance(vertices.begin(), iter_vertex);

	hddm_r::ProductList products = iter_vertex->getProducts();
	hddm_r::ProductList::iterator iter_product;
	// Loop over ProductList
	for (iter_product = products.begin(); iter_product != products.end(); ++iter_product){
	  int iparticle = std::distance(products.begin(), iter_product);

	  if(debug) cout << "reaction: " << ireaction << " vertex: " << ivertex << " particle: " << iparticle
			 << " id: " << setw(3) << iter_product->getId() << " parent: " << iter_product->getParentId()
			 << " PDG type: " << setw(7) << iter_product->getPdgtype()
			 << " type: " << ParticleType((Particle_t)PDGtoPType((Particle_t)iter_product->getPdgtype()))
			 << " vertex: (" << iter_vertex->getOrigin().getVx() << ", " << iter_vertex->getOrigin().getVy() << ", " << iter_vertex->getOrigin().getVz() << ")"
			 << endl;
	  
	  // Find Lambda
	  if(iter_product->getPdgtype() == PDGtype(Lambda)){
	    foundLambda = true;
	    idLambda = iter_product->getId();
	    vertexLambda[0] = iter_vertex->getOrigin().getVx();
	    vertexLambda[1] = iter_vertex->getOrigin().getVy();
	    vertexLambda[2] = iter_vertex->getOrigin().getVz();
	  }

	  // Find decay products of Lambda
	  if(foundLambda){
	    // if parent was Lambda
	    if(iter_product->getParentId() == idLambda){
	      // find what type it is
	      if(iter_product->getPdgtype() == PDGtype(Proton)){
		foundProton = true;
		idProton = iter_product->getId();
		vertexProton[0] = iter_vertex->getOrigin().getVx();
		vertexProton[1] = iter_vertex->getOrigin().getVy();
		vertexProton[2] = iter_vertex->getOrigin().getVz();
		vertexdiff = sqrt(pow(vertexLambda[0] - vertexProton[0],2.) + pow(vertexLambda[1] - vertexProton[1],2.) + pow(vertexLambda[2] - vertexProton[2],2.));
	      }
	    }
	  }

	  // If we found a proton decaying from Lambda,
	  // and the vertex difference is less than VERTEXDIFFMAX,
	  // save this event
	  if(foundProton && vertexdiff < VERTEXDIFFMAX){
	    return true;
	  }

	} // end of loop over ProductList
      }	// end of loop over VertexList
    } // end of loop over ReactiionList
  } // end of select_type 3
  //__________________________________________________________________________________________________

  // default is to return false
  return false;
}
