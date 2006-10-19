//********************************************************
// DFDCPseudo_factory.cc - factory producing first-order 
// reconstructed points for the FDC
// Author: Craig Bookwalter (craigb at jlab.org)
// Date:   March 2006
// UVX cathode-wire matching revised by Simon Taylor, Aug 2006
//********************************************************

#include "DFDCPseudo_factory.h"
#include "DFDCGeometry.h"

#include <TMatrixD.h>
#include <TDecompLU.h>

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
///
jerror_t DFDCPseudo_factory::evnt(JEventLoop* eventLoop, int eventNo) {
	vector<const DFDCHit*> fdcHits;
	vector<const DFDCHit*> xHits;
	vector<const DFDCCathodeCluster*> cathClus;
	vector<const DFDCCathodeCluster*> uClus;
	vector<const DFDCCathodeCluster*> oneLayerU;
	vector<const DFDCCathodeCluster*> vClus;
	vector<const DFDCCathodeCluster*> oneLayerV;
	vector<const DFDCHit*> oneLayerX;
	
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
			oneLayerX.push_back(*xIt);
		makePseudo(oneLayerX, oneLayerU, oneLayerV,iLayer);
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
void DFDCPseudo_factory::makePseudo(vector<const DFDCHit*>& x,
				    vector<const DFDCCathodeCluster*>& u,
				    vector<const DFDCCathodeCluster*>& v,
				    int layer)
{
  if ((u.size() == 0) || (v.size() == 0) || (x.size() == 0))
    return;
  vector<const DFDCHit*>::iterator xIt;
  vector<const DFDCCathodeCluster*>::iterator uIt;
  vector<const DFDCCathodeCluster*>::iterator vIt;
  centroid_t temp;
  float E1,E2,pos1,pos2;
  
  // Loop over all U and V clusters looking for peaks
  upeaks.clear();
  for (uIt = u.begin(); uIt != u.end(); uIt++){
    vector<const DFDCHit*>::const_iterator strip=(*uIt)->members.begin();
    switch((*uIt)->members.size()){
    case 0: // Make sure we have data!!
      break;
    case 1: // One isolated hit in the cluster:  use element number itself
      temp.pos=(*strip)->element;
      temp.q=2.*((*strip)->dE);  // Each cathode view should see half the 
                                 // anode charge
      upeaks.push_back(temp);
      break;
    case 2: //Two adjacent hits: use average for the centroid
      pos1=(*strip)->element;
      pos2=(*(strip+1))->element;
      E1=(*strip)->dE;
      E2=(*(strip+1))->dE;      
      temp.pos=(pos1*E1+pos2*E2)/(E1+E2);
      temp.q=2.*(E1+E2);
      upeaks.push_back(temp);
      break;
    default:      
      for (strip=(*uIt)->members.begin();
	   strip!=(*uIt)->members.end();strip++){
	FindCentroid((*uIt)->members,strip,upeaks);
      }
      break;
    }
  }	
  vpeaks.clear();
  for (vIt = v.begin(); vIt != v.end(); vIt++){
    vector<const DFDCHit*>::const_iterator strip=(*vIt)->members.begin();
    switch((*vIt)->members.size()){
    case 0: // Make sure we have data!!
      break;
    case 1: // One isolated hit in the cluster:  use element number itself
      temp.pos=(*strip)->element;
      temp.q=2.*((*strip)->dE);
      temp.numstrips=1;
      vpeaks.push_back(temp);
      break;
    case 2: //Two adjacent hits: use average for the centroid
      pos1=(*strip)->element;
      pos2=(*(strip+1))->element;
      E1=(*strip)->dE;
      E2=(*(strip+1))->dE;      
      temp.pos=(pos1*E1+pos2*E2)/(E1+E2);
      temp.q=2.*(E1+E2);
      temp.numstrips=2;
      vpeaks.push_back(temp);
      break;
    default:      
      for (strip=(*vIt)->members.begin();
	   strip!=(*vIt)->members.end();strip++){	
	FindCentroid((*vIt)->members,strip,vpeaks);
      }
      break;
    }
  }
  if (upeaks.size()*vpeaks.size()>0){
    //Loop over all u and v centroids looking for matches with wires
    for (unsigned int i=0;i<upeaks.size();i++){
      for (unsigned int j=0;j<vpeaks.size();j++){
	// In the layer local coordinate system, wires are quantized 
	// in the x-direction and y is along the wire.
	float x_from_strips=DFDCGeometry::getXLocalStrips(upeaks[i].pos,
							  vpeaks[j].pos);
	float y_from_strips=DFDCGeometry::getYLocalStrips(upeaks[i].pos,
							  vpeaks[j].pos);
	for(xIt=x.begin();xIt!=x.end();xIt++){
	  float x_from_wire=DFDCGeometry::getWireR(*xIt);
	  if (fabs(x_from_wire-x_from_strips)<WIRE_SPACING/2.){
	    int status=upeaks[i].numstrips+vpeaks[j].numstrips;
	    float xres=WIRE_SPACING/2./sqrt(12.);
	    float yres=fabs(x_from_wire-x_from_strips);
	    
	    yres=x_from_wire-x_from_strips;
	   
	    DFDCPseudo* newPseu = new DFDCPseudo(x_from_wire,xres,
						 y_from_strips,
						 yres,layer,(*xIt)->element,
						 (*xIt)->t,status);
	    _data.push_back(newPseu);
	  } // match in x
	} // xIt loop
      } // vpeaks loop
    } // upeaks loop
  } // if we have peaks in both u and v views
}			

//
/// DFDCPseudo_factory::FindCentroid()
///   Uses the Newton-Raphson method for solving the set of non-linear
/// equations describing the charge distribution over 3 strips for the peak 
/// position x0, the anode charge qa, and the "width" parameter k2.  
/// See Numerical Recipes in C p.379-383.
/// Updates list of centroids. 
///

jerror_t DFDCPseudo_factory::FindCentroid(const vector<const DFDCHit*>& H, 
		       vector<const DFDCHit *>::const_iterator peak,
                       vector<centroid_t>&centroids){
  centroid_t temp; 
  // Make sure we do not exceed the range of the vector
  if (peak>H.begin() && peak+1!=H.end()){
    float err_diff1=0.,err_diff2=0.;

    // Some code for checking for significance of fluctuations.
    // Currently disabled.
    //float dq1=(*(peak-1))->dq;
    //float dq2=(*peak)->dq;
    //float dq3=(*(peak+1))->dq;
    //err_diff1=sqrt(dq1*dq1+dq2*dq2);
    //err_diff2=sqrt(dq2*dq2+dq3*dq3);
  
    // Check for a peak in three adjacent strips
    if ((*peak)->dE-(*(peak-1))->dE > err_diff1
		&& (*peak)->dE-(*(peak+1))->dE > err_diff2){
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
	  temp.pos=par(X0,0);
	  temp.q=par(QA,0);
	  temp.numstrips=3;
	  centroids.push_back(temp);
	  
	  return NOERROR;
	}
	// Check that J is invertible
	TDecompLU lu(J);
	if (lu.Decompose()==false) return UNRECOVERABLE_ERROR; // error placeholder

	// Invert the J matrix
	TMatrixD InvJ(TMatrixD::kInverted,J);
	
	// Find the corrections to the vector par:
	dpar=InvJ*F;
	// calculate the improved values of the parameters
	par-=dpar;
	
	//Check for convergence
	for (i=0;i<3;i++){
	  errx+=fabs(dpar(i,0));
	}
	if (errx<TOLX){	
	  temp.pos=par(X0,0);
	  temp.q=par(QA,0);
	  temp.numstrips=3;
	  centroids.push_back(temp);
	  
	  return NOERROR;
	}
      }
    }
  }
  return INFINITE_RECURSION; // error placeholder
}
