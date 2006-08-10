//********************************************************
// DFDCPseudo_factory.cc - factory producing first-order 
// reconstructed points for the FDC
// Author: Craig Bookwalter (craigb at jlab.org)
// Date:   March 2006
//********************************************************

#include "DFDCPseudo_factory.h"

#include <TMatrixD.h>

#define X0 0
#define QA 1
#define K2 2
#define ITER_MAX 100
#define TOLX 1e-4
#define TOLF 1e-4
#define A_OVER_H 0.475


///
/// DFDCCathodeCluster_gLayer_cmp(): 
/// non-member function passed to std::sort() to sort DFDCCathodeCluster pointers 
/// by their gLayer attributes.
///
bool DFDCCathodeCluster_gLayer_cmp(const DFDCCathodeCluster* a, 
								   const DFDCCathodeCluster* b) {
	return a->gLayer < b->gLayer;
}

///
/// DFDCPseudo_factory::DFDCPseudo_factory():
/// default constructor -- initializes log file
///
DFDCPseudo_factory::DFDCPseudo_factory() {
	logFile = new ofstream("DFDCPseudo_factory.log");
	_log = new JStreamLog(*logFile, "PSEUDO");
	*_log << "File initialized." << endMsg;
}

///
/// DFDCPseudo_factory::~DFDCPseudo_factory():
/// default destructor -- closes log file
///
DFDCPseudo_factory::~DFDCPseudo_factory() {
	logFile->close();
	delete logFile;
	delete _log;
}

///
/// DFDCPseudo_factory::evnt():
/// this is the place that anode hits and DFDCCathodeClusters are organized into pseudopoints.
/// For now, this is done purely by geometry, with no drift or peak-finding. See also
/// DFDCPseudo_factory::makePseudo().
///
jerror_t DFDCPseudo_factory::evnt(JEventLoop* eventLoop, int eventNo) {
	vector<const DFDCHit*> fdcHits;
	vector<const DFDCHit*> xHits;
	vector<const DFDCCathodeCluster*> cathClus;
	vector<const DFDCCathodeCluster*> uClus;
	vector<const DFDCCathodeCluster*> oneLayerU;
	vector<const DFDCCathodeCluster*> vClus;
	vector<const DFDCCathodeCluster*> oneLayerV;
	std::map<int, const DFDCHit*> oneLayerX;
	DFDCGeometry geo;
	float angle = 0.0;
	float pi	= 3.1415926;
	
	eventLoop->Get(fdcHits);
	eventLoop->Get(cathClus);
	
	// Ensure clusters are in order of ascending Z position.
	std::sort(cathClus.begin(), cathClus.end(), DFDCCathodeCluster_gLayer_cmp);
	
	// Sift through hits and select out anode hits.
	for (unsigned int i=0; i < fdcHits.size(); i++)
		if (fdcHits[i]->type == 0)
			xHits.push_back(fdcHits[i]);
			
	// Sift through clusters and put U and V clusters into respective vectors.
	for (unsigned int i=0; i < cathClus.size(); i++) {
		if (cathClus[i]->plane == 1)
			vClus.push_back(cathClus[i]);
		else
			uClus.push_back(cathClus[i]);
	}

	vector<const DFDCCathodeCluster*>::iterator uIt = uClus.begin();
	vector<const DFDCCathodeCluster*>::iterator vIt = vClus.begin();
	vector<const DFDCHit*>::iterator xIt = xHits.begin();
	
	// For each layer, get its sets of V, X, and U hits, and then pass them to the geometrical
	// organization routine, DFDCPseudo_factory::makePseudo()
	for (int iLayer=1; iLayer <= 24; iLayer++) {
		for (; ((uIt != uClus.end() && (*uIt)->gLayer == iLayer)); uIt++)
			oneLayerU.push_back(*uIt);
		for (; ((vIt != vClus.end() && (*vIt)->gLayer == iLayer)); vIt++)
			oneLayerV.push_back(*vIt);
		for (; ((xIt != xHits.end() && (*xIt)->gLayer == iLayer)); xIt++)
			oneLayerX.insert(std::pair<const int, const DFDCHit*>((*xIt)->element,
																  *xIt));
		makePseudo(oneLayerX, oneLayerU, oneLayerV, geo.getLayerRotation(iLayer), 	
					iLayer);
		angle += pi/3.0;
		oneLayerU.clear();
		oneLayerV.clear();
		oneLayerX.clear();
	}
	
	return NOERROR;
}

/// 
/// DFDCPseudo_factory::makePseudo():
/// performs UV+X matching to create pseudopoints
///
void DFDCPseudo_factory::makePseudo(	std::map<int, const DFDCHit*>& x,
										vector<const DFDCCathodeCluster*>& u,
										vector<const DFDCCathodeCluster*>& v,
										float angle,
										int layer)
{
	if ((u.size() == 0) || (v.size() == 0) || (x.size() == 0))
		return;
	std::map<int, const DFDCHit*>::iterator xIt;
	vector<const DFDCCathodeCluster*>::iterator uIt;
	vector<const DFDCCathodeCluster*>::iterator vIt;

	// Loop over all U and V clusters looking for peaks
	upeaks.clear();
	for (uIt = u.begin(); uIt != u.end(); uIt++){
	  vector<const DFDCHit*>::const_iterator strip;
	  for (strip=(*uIt)->members.begin();
	       strip!=(*uIt)->members.end();strip++){
	    FindCentroid((*uIt)->members,strip,upeaks);
	  }
	}	
	vpeaks.clear();
	for (vIt = v.begin(); vIt != v.end(); vIt++){
	  vector<const DFDCHit*>::const_iterator strip;
	  for (strip=(*vIt)->members.begin();
	       strip!=(*vIt)->members.end();strip++){
	    FindCentroid((*vIt)->members,strip,vpeaks);
	  }
	}


	// Loop over all U and V clusters
	for (uIt = u.begin(); uIt != u.end(); uIt++) {
		for (vIt = v.begin(); vIt != v.end(); vIt++) {
		 
	  

			// Find the left and right edges of the cluster
			float leftX 	= intersectX((*uIt)->beginStrip, (*vIt)->beginStrip);
			float rightX 	= intersectX((*uIt)->endStrip, (*vIt)->endStrip);

			// Find the Y coordinate of the U and V strips with maxmimum dE
			float y			= intersectY((*uIt)->maxStrip, (*vIt)->maxStrip);
			float res 		= ((*uIt)->width + (*vIt)->width / 2) / 2;

			// Since wire numbers are integers, cast left and right bounds to integers.
			int left 		= static_cast<int>(floor(leftX + 60.0));
			int right 		= static_cast<int>(ceil(rightX + 60.0));
			
			// Ask the map if it has a hit for wire numbers between left and right;
			// if it does, create a new pseudopoint.			
			for (int i = left; i <= right; i++) 
			{
				xIt = x.find(i);
				if (xIt == x.end())
					continue;
				float xc = (*xIt).first - 60;
				xc = xc*cos(angle) - y*sin(angle);
				y = xc*sin(angle) + y*sin(angle);
				DFDCPseudo* newPseu = new DFDCPseudo(xc, y, layer, res);
				_data.push_back(newPseu);
			}
		}
	}
}			

///
/// DFDCPseudo_factory::intersectX():
/// finds the X coordinate of a U-V intersection
///
float DFDCPseudo_factory::intersectX(int u, int v) {	
	float uYint = (119 - u)*sqrt(2.0) + (1/sqrt(2.0));
	float vYint = (v - 119)*sqrt(2.0) - (1/sqrt(2.0)); 	

	return (uYint - vYint) / 2;
}			

///
/// DFDCPseudo_factory::intersectY():
/// finds the Y coordinate of a U-V intersection
///
float DFDCPseudo_factory::intersectY(int u, int v) {
	float uYint = (119 - u)*sqrt(2.0) + (1/sqrt(2.0));
	float vYint = (v - 119)*sqrt(2.0) - (1/sqrt(2.0));
	
	return (uYint + vYint) / 2;
}



//
/// DFDCPseudo_factory::FindCentroid()
//   Uses the Newton-Raphson method for solving the set of non-linear
// equations describing the charge distribution over 3 strips for the peak 
// position x0, the anode charge qa, and the "width" parameter k2.  
// See Numerical Recipes in C p.379-383.
// Updates list of centroids. 
///

jerror_t DFDCPseudo_factory::FindCentroid(const vector<const DFDCHit*>& H, 
		       vector<const DFDCHit *>::const_iterator peak,
                       vector<centroid_t>&centroids){

  // Check for a peak with charge on adjacent strips 
  if (peak>H.begin() && peak+1!=H.end() && (*peak)->dE > (*(peak-1))->dE 
      && (*peak)->dE>(*(peak+1))->dE){
    // Define some matrices for use in the Newton-Raphson iteration
    TMatrixD J(3,3);  //Jacobean matrix
    TMatrixD F(3,1),N(3,1),X(3,1),par(3,1),dpar(3,1);
    int i=0;
    double sum=0.;
 
    // Initialize the matrices to some suitable starting values
    par(X0,0)=double((*peak)->element);
    par(K2,0)=1.;
    for (vector<const DFDCHit*>::const_iterator j=peak-1;j<=peak+1;j++){
      X(i,0)=double((*j)->element);
      N(i++,0)=double((*j)->dE);
      sum+=double((*j)->dE);
    }
    par(QA,0)=2.*sum;

    // Newton-Raphson procedure
    double errf=0.,errx=0.;
    for (int iter=1;iter<=ITER_MAX;iter++){
      errf=0.;
      errx=0.;
      for (i=0;i<3;i++){
	 double argp=par(K2,0)*(par(X0,0)-X(i,0)+A_OVER_H);
	 double argm=par(K2,0)*(par(X0,0)-X(i,0)-A_OVER_H);

	 //Find the Jacobian matrix:  J_ij = dF_i/dx_j. 
	 J(i,QA)=-(tanh(argp)-tanh(argm))/4.;
	 J(i,K2)=-par(QA,0)/4.*(argp/par(K2,0)*(1.-tanh(argp)*tanh(argp))
				-argm/par(K2,0)*(1.-tanh(argm)*tanh(argm)));
	 J(i,X0)=-par(QA,0)*par(K2,0)/4.
	   *(tanh(argm)*tanh(argm)-tanh(argp)*tanh(argp));
	 
	 // update F_i
	 F(i,0)=N(i,0)-par(QA,0)/4.*(tanh(argp)-tanh(argm));

	 errf+=fabs(F(i,0));
      }
      // Check for convergence
      if (errf<TOLF){
	centroid_t temp;
	temp.pos=par(X0,0);
	temp.q=par(QA,0);
	centroids.push_back(temp);
	
	return NOERROR;
      }
      // Find the corrections to the vector par:
      TMatrixD InvJ(TMatrixD::kInverted,J);      
      dpar=(-1)*InvJ*F;
      // calculate the improved values of the parameters
      par+=dpar;
      
      //Check for convergence
      for (i=0;i<3;i++){
	errx+=fabs(dpar(i,0));
      }
      if (errx<TOLX){	
	centroid_t temp;
	temp.pos=par(X0,0);
	temp.q=par(QA,0);
	centroids.push_back(temp);
	
	return NOERROR;
      }
    }
  }
  return NOERROR;
}
