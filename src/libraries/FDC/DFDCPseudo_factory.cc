//********************************************************
// DFDCPseudo_factory.cc - factory producing first-order 
// reconstructed points for the FDC
// Author: Craig Bookwalter (craigb at jlab.org)
// Date:   March 2006
// UVX cathode-wire matching revised by Simon Taylor, Aug 2006
//********************************************************

#include "DFDCPseudo_factory.h"
#include "HDGEOMETRY/DGeometry.h"
#include "DFDCGeometry.h"
#include <TRACKING/DTrackHitSelectorTHROWN.h>
#include <TROOT.h>

#define HALF_CELL 0.5
#define MAX_DEFLECTION 0.15
#define X0 0
#define QA 1
#define K2 2
#define ITER_MAX 100
#define TOLX 1e-4
#define TOLF 1e-4
#define A_OVER_H 0.4
#define ALPHA 1e-4 // rate parameter for Newton step backtracking algorithm
#define W_EFF 0.0302 // keV
#define GAS_GAIN 8e4
#define ELECTRON_CHARGE 1.6022e-4 // fC


///
/// DFDCAnode_gLayer_cmp(): 
/// non-member function passed to std::sort() to sort DFDCHit pointers 
/// for the anode wires by their gLayer attributes.
///
bool DFDCAnode_gLayer_cmp(const DFDCHit* a, const DFDCHit* b) {
	return a->gLayer < b->gLayer;
}

///
/// DFDCPseudo_factory::DFDCPseudo_factory():
/// default constructor -- initializes log file
///
DFDCPseudo_factory::DFDCPseudo_factory() {
	_log = new JStreamLog(std::cout, "FDC PSEUDO >>");
	*_log << "File initialized." << endMsg;
}

///
/// DFDCPseudo_factory::~DFDCPseudo_factory():
/// default destructor -- closes log file
///
DFDCPseudo_factory::~DFDCPseudo_factory() {
	delete _log;
}

//------------------
// init
//------------------
jerror_t DFDCPseudo_factory::init(void)
{
  RIN_FIDUCIAL = 1.5;
  ROUT_FIDUCIAL=48.0;
  MAX_ALLOWED_FDC_HITS=(5+5+1)*24*10;
  STRIP_ANODE_TIME_CUT=10.;

  gPARMS->SetDefaultParameter("FDC:ROUT_FIDUCIAL",ROUT_FIDUCIAL, "Outer fiducial radius of FDC in cm"); 
  gPARMS->SetDefaultParameter("FDC:RIN_FIDUCIAL",RIN_FIDUCIAL, "Inner fiducial radius of FDC in cm");
  gPARMS->SetDefaultParameter("FDC:MAX_ALLOWED_FDC_HITS",MAX_ALLOWED_FDC_HITS, "Max. number of FDC hits (includes both cathode strips and wires hits) to allow before considering event too busy to attempt FDC tracking");
  gPARMS->SetDefaultParameter("FDC:STRIP_ANODE_TIME_CUT",STRIP_ANODE_TIME_CUT,
			      "maximum time difference between strips and wires (in ns)"); 

  DEBUG_HISTS = false;
  gPARMS->SetDefaultParameter("FDC:DEBUG_HISTS",DEBUG_HISTS);


  return NOERROR;
}


//------------------
// brun
//------------------
jerror_t DFDCPseudo_factory::brun(JEventLoop *loop, int runnumber)
{
  // Get pointer to DGeometry object
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  const DGeometry *dgeom  = dapp->GetDGeometry(runnumber);
  dgeom->GetFDCWires(fdcwires);

  if(DEBUG_HISTS){
    dapp->Lock();

    // Histograms may already exist. (Another thread may have created them)
    // Try and get pointers to the existing ones.

    qa_vs_qc= (TH2F*) gROOT->FindObject("qa_vs_qc");
    if (!qa_vs_qc) qa_vs_qc=new TH2F("qa_vs_qc","Anode charge from wire versus charge derived from cathodes",200,0,10,200,0,10);

    qa_qc_diff=(TH2F*)gROOT->FindObject("qa_qc_diff");
    if (!qa_qc_diff) qa_qc_diff=new TH2F("qa_qc_diff","(qa-qcv) vs (qa-qcu)",
					 200,-25,25,200,-25,25);

     dtv_vs_dtu= (TH2F*) gROOT->FindObject("dtv_vs_dtu");
    if (!dtv_vs_dtu) dtv_vs_dtu=new TH2F("dtv_vs_dtu","t(wire)-t(v) vs t(wire)-t(u)",100,-50,50,100,-50,50);

    u_wire_dt_vs_wire=(TH2F *) gROOT->FindObject("u_wire_dt_vs_wire");
    if (!u_wire_dt_vs_wire) u_wire_dt_vs_wire=new TH2F("u_wire_dt_vs_wire","wire/u cathode time difference vs wire number",
			   96,0.5,96.5,100,-50,50);
    v_wire_dt_vs_wire=(TH2F *) gROOT->FindObject("v_wire_dt_vs_wire");
    if (!v_wire_dt_vs_wire) v_wire_dt_vs_wire=new TH2F("v_wire_dt_vs_wire","wire/v cathode time difference vs wire number",
						       96,0.5,96.5,100,-50,50);
 uv_dt_vs_u=(TH2F *) gROOT->FindObject("uv_dt_vs_u");
    if (!uv_dt_vs_u) uv_dt_vs_u=new TH2F("uv_dt_vs_u","uv time difference vs u",
			   192,0.5,192.5,100,-50,50); 
    uv_dt_vs_v=(TH2F *) gROOT->FindObject("uv_dt_vs_v");
    if (!uv_dt_vs_v) uv_dt_vs_v=new TH2F("uv_dt_vs_v","uv time difference vs v",
			   192,0.5,192.5,100,-50,50);
    dapp->Unlock();
  }

  return NOERROR;
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

	// Get all FDC hits (anode and cathode)
	eventLoop->Get(fdcHits);

	// For events with a very large number of hits, assume
	// we can't reconstruct them so bail early
	// Feb. 8, 2008  D.L. (updated to config param. Nov. 18, 2010 D.L.)
	if(fdcHits.size()>MAX_ALLOWED_FDC_HITS){
		_DBG_<<"Too many hits in FDC ("<<fdcHits.size()<<", max="<<MAX_ALLOWED_FDC_HITS<<")! Pseudopoint reconstruction in FDC bypassed for event "<<eventLoop->GetJEvent().GetEventNumber()<<endl;
		return NOERROR;
	}

	// Get cathode clusters
	eventLoop->Get(cathClus);

	// Sift through hits and select out anode hits.
	for (unsigned int i=0; i < fdcHits.size(); i++)
		if (fdcHits[i]->type == 0)
			xHits.push_back(fdcHits[i]);
	// Make sure the wires are also in order of ascending z position
	std::sort(xHits.begin(), xHits.end(), DFDCAnode_gLayer_cmp);
			
	// Sift through clusters and put U and V clusters into respective vectors.
	for (unsigned int i=0; i < cathClus.size(); i++) {
		if (cathClus[i]->plane == 1)
			vClus.push_back(cathClus[i]);
		else
			uClus.push_back(cathClus[i]);
	}
	
	// If this is simulated data then we want to match up the truth hit
	// with this "real" hit. Ideally, this would be done at the
	// DFDCHit object level, but the organization of the data in HDDM
	// makes that difficult. Here we have the full wire definition so
	// we make the connection here.
	vector<const DMCTrackHit*> mctrackhits;
	eventLoop->Get(mctrackhits);

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
	  if (oneLayerU.size()>0 && oneLayerV.size()>0 && oneLayerX.size()>0)
	    makePseudo(oneLayerX, oneLayerU, oneLayerV,iLayer, mctrackhits);
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
				    int layer,
				    vector<const DMCTrackHit*> &mctrackhits)
{
  vector<const DFDCHit*>::iterator xIt;
  vector<centroid_t>upeaks;
  vector<centroid_t>vpeaks;


  //printf("---------u clusters --------\n");
  // Loop over all U and V clusters looking for peaks
  for (unsigned int i=0;i<u.size();i++){
    //printf("Cluster %d\n",i);
    if (u[i]->members.size()>2)
      {
      for (vector<const DFDCHit*>::const_iterator strip=u[i]->members.begin();
	   strip!=u[i]->members.end();strip++){  
	//printf("  %d %f %f\n",(*strip)->element,(*strip)->q,(*strip)->t);
	if (FindCentroid(u[i]->members,strip,upeaks)==NOERROR){
	  upeaks[upeaks.size()-1].cluster=i;
	}
      }
    }
  }  
  //printf("---------v cluster --------\n");	
  for (unsigned int i=0;i<v.size();i++){
    //printf("Cluster %d\n",i);
    if (v[i]->members.size()>2)
      {
      for (vector<const DFDCHit*>::const_iterator strip=v[i]->members.begin();
	   strip!=v[i]->members.end();strip++){		
	//printf("  %d %f %f\n",(*strip)->element,(*strip)->q,(*strip)->t);
	if (FindCentroid(v[i]->members,strip,vpeaks)==NOERROR){
	  vpeaks[vpeaks.size()-1].cluster=i;
	}
      }
    }
  }
  if (upeaks.size()*vpeaks.size()>0){
    //Loop over all u and v centroids looking for matches with wires
    for (unsigned int i=0;i<upeaks.size();i++){
      for (unsigned int j=0;j<vpeaks.size();j++){
	// In the layer local coordinate system, wires are quantized 
	// in the x-direction and y is along the wire.
	double x_from_strips=DFDCGeometry::getXLocalStrips(upeaks[i].pos,
						 vpeaks[j].pos);
	double y_from_strips=DFDCGeometry::getYLocalStrips(upeaks[i].pos,
							  vpeaks[j].pos);
	for(xIt=x.begin();xIt!=x.end();xIt++){
	  if ((*xIt)->element<=WIRES_PER_PLANE && (*xIt)->element>0){
	    double x_from_wire=DFDCGeometry::getWireR(*xIt);
	    // Test radial value for checking whether or not the hit is within
	    // the fiducial region of the detector
	    double rtest=sqrt(x_from_wire*x_from_wire
			     +y_from_strips*y_from_strips);
	    double delta_x=x_from_wire-x_from_strips;
	 
	    if (fabs(delta_x)<0.5*WIRE_SPACING && rtest<ROUT_FIDUCIAL
		&& rtest>RIN_FIDUCIAL){
	      double dt1 = (*xIt)->t - upeaks[i].t;
	      double dt2 = (*xIt)->t - vpeaks[j].t;
	      double q_cathodes=0.5*(upeaks[i].q+vpeaks[j].q);
	      double charge_to_energy=W_EFF/(GAS_GAIN*ELECTRON_CHARGE);
	      double dE=charge_to_energy*q_cathodes;

	      if (DEBUG_HISTS){
		qa_vs_qc->Fill(dE,((*xIt)->q));
		
		qa_qc_diff->Fill(((*xIt)->q-charge_to_energy*upeaks[i].q),
				 ((*xIt)->q-charge_to_energy*vpeaks[i].q));
		dtv_vs_dtu->Fill(dt1,dt2);
		u_wire_dt_vs_wire->Fill((*xIt)->element,(*xIt)->t-upeaks[i].t);
		v_wire_dt_vs_wire->Fill((*xIt)->element,(*xIt)->t-vpeaks[j].t);
		uv_dt_vs_u->Fill(upeaks[i].pos,upeaks[i].t-vpeaks[j].t);
		uv_dt_vs_v->Fill(vpeaks[j].pos,upeaks[i].t-vpeaks[j].t);
	      }
	      if (sqrt(dt1*dt1+dt2*dt2)>STRIP_ANODE_TIME_CUT) continue;
	      
	      int status=upeaks[i].numstrips+vpeaks[j].numstrips;
	      double xres=WIRE_SPACING/2./sqrt(12.);
	      double cosangle,sinangle;	   

	      DFDCPseudo* newPseu = new DFDCPseudo;
	      newPseu->u = upeaks[i].pos;
	      newPseu->v = vpeaks[j].pos;
	      newPseu->w      = x_from_wire;
	      newPseu->dw     = xres;
	      newPseu->s      = y_from_strips;
	      newPseu->ds     = delta_x;
	      newPseu->wire   = fdcwires[layer-1][(*xIt)->element-1];
	      newPseu->time   = (*xIt)->t;
	      newPseu->status = status;
	      newPseu->itrack = (*xIt)->itrack;

	      newPseu->AddAssociatedObject(v[vpeaks[j].cluster]);
	      newPseu->AddAssociatedObject(u[upeaks[i].cluster]);
	      
	      newPseu->dE = dE;
	      
	    
	      // It can occur (although rarely) that newPseu->wire is NULL
	      // which causes us to crash below. In these cases, we can't really
	      // make a psuedo point so we delete the current object
	      // and just go on to the next one.
	      if(newPseu->wire==NULL){
		_DBG_<<"newPseu->wire=NULL! This shouldn't happen. Complain to staylor@jlab.org"<<endl;
		delete newPseu;
		continue;
	      }
	      sinangle=newPseu->wire->udir(0);
	      cosangle=newPseu->wire->udir(1);

	      newPseu->xy.Set((newPseu->w)*cosangle+(newPseu->s)*sinangle,
			       -(newPseu->w)*sinangle+(newPseu->s)*cosangle);
	      
	      double sigx2=HALF_CELL*HALF_CELL/3.;
	      double sigy2=MAX_DEFLECTION*MAX_DEFLECTION/3.;
	      newPseu->covxx=sigx2*cosangle*cosangle+sigy2*sinangle*sinangle;
	      newPseu->covyy=sigx2*sinangle*sinangle+sigy2*cosangle*cosangle;
	      newPseu->covxy=(sigy2-sigx2)*sinangle*cosangle;

	      // Try matching truth hit with this "real" hit.
	      const DMCTrackHit *mctrackhit = DTrackHitSelectorTHROWN::GetMCTrackHit(newPseu->wire, DRIFT_SPEED*newPseu->time, mctrackhits);
	      if(mctrackhit)newPseu->AddAssociatedObject(mctrackhit);
	      
	      _data.push_back(newPseu);
	    } // match in x
	  } else _DBG_ << "Bad wire " << (*xIt)->element <<endl;
	} // xIt loop
      } // vpeaks loop
    } // upeaks loop
  } // if we have peaks in both u and v views
}			

//
/// DFDCPseudo_factory::CalcMeanTime()
/// Calculate the mean and rms of the times of the hits passed in "H".
/// The contents of H should be pointers to a single cluster in a
/// cathode plane. 
//  1/2/2008 D.L.
void DFDCPseudo_factory::CalcMeanTime(const vector<const DFDCHit*>& H, double &t, double &t_rms)
{
	// Calculate mean
	t=0.0;
	for(unsigned int i=0; i<H.size(); i++)t+=H[i]->t;
	if(H.size()>0)t/=(double)H.size();
	
	// Calculate RMS
	t_rms=0.0;
	for(unsigned int i=0; i<H.size(); i++)t_rms+=pow((double)(H[i]->t-t),2.0);
	if(H.size()>0)t_rms = sqrt(t_rms/(double)H.size());
}

// Find the mean time and rms for a group of 3 hits with a maximum in the 
// center hit
void DFDCPseudo_factory::CalcMeanTime(vector<const DFDCHit *>::const_iterator peak, double &t, double &t_rms)
{
  // Calculate mean
  t=0.0;
  for (vector<const DFDCHit*>::const_iterator j=peak-1;j<=peak+1;j++){
    t+=(*j)->t;
  }
  t/=3.;
  
  // Calculate RMS
  t_rms=0.0;
  for (vector<const DFDCHit*>::const_iterator j=peak-1;j<=peak+1;j++){
    t_rms+=((*j)->t-t)*((*j)->t-t);
  }
  t_rms=sqrt(t_rms/3.);
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
    double err_diff1=0.,err_diff2=0.;
	 
    // Fill in time info in temp  1/2/2008 D.L.
    //CalcMeanTime(H, temp.t, temp.t_rms);

    // Some code for checking for significance of fluctuations.
    // Currently disabled.
    //double dq1=(*(peak-1))->dq;
    //double dq2=(*peak)->dq;
    //double dq3=(*(peak+1))->dq;
    //err_diff1=sqrt(dq1*dq1+dq2*dq2);
    //err_diff2=sqrt(dq2*dq2+dq3*dq3);
   
    // Check for a peak in three adjacent strips
    if ((*peak)->q-(*(peak-1))->q > err_diff1
                && (*peak)->q-(*(peak+1))->q > err_diff2){
      // Define some matrices for use in the Newton-Raphson iteration
      DMatrix J(3,3);  //Jacobean matrix
      DMatrix F(3,1),N(3,1),X(3,1),par(3,1);
      DMatrix newpar(3,1);
      int i=0;
      double sum=0.;
 
      // Initialize the matrices to some suitable starting values
      par(X0,0)=double((*peak)->element);
      par(K2,0)=1.;
      for (vector<const DFDCHit*>::const_iterator j=peak-1;j<=peak+1;j++){
        X(i,0)=double((*j)->element);
        N(i++,0)=double((*j)->q);
        sum+=double((*j)->q);
      }
      par(QA,0)=2.*sum;
      newpar=par;
    
      // Newton-Raphson procedure
      double errf=0.,errx=0;
      for (int iter=1;iter<=ITER_MAX;iter++){
        errf=0.;
        errx=0.;

        // Compute Jacobian matrix: J_ij = dF_i/dx_j.
        for (i=0;i<3;i++){
          double argp=par(K2,0)*(par(X0,0)-X(i,0)+A_OVER_H);
          double argm=par(K2,0)*(par(X0,0)-X(i,0)-A_OVER_H);
           
          J(i,QA)=-(tanh(argp)-tanh(argm))/4.;
          J(i,K2)=-par(QA,0)/4.*(argp/par(K2,0)*(1.-tanh(argp)*tanh(argp))
                                 -argm/par(K2,0)*(1.-tanh(argm)*tanh(argm)));
          J(i,X0)=-par(QA,0)*par(K2,0)/4.
            *(tanh(argm)*tanh(argm)-tanh(argp)*tanh(argp)); 
	  F(i,0)=N(i,0)-par(QA,0)/4.*(tanh(argp)-tanh(argm));
	  if (fabs(F(i,0))>errf) errf=fabs(F(i,0));
        }
	// Check for convergence
	if (errf<TOLF){	
	  temp.pos=par(X0,0);
	  temp.q=par(QA,0);
	  temp.numstrips=3;  
	  temp.t=(*peak)->t;
	  temp.t_rms=0.;
	  //CalcMeanTime(peak,temp.t,temp.t_rms);
	  centroids.push_back(temp);
  
	  return NOERROR;
	}

        // Check that J is invertible
        TDecompLU lu(J);
        if (lu.Decompose()==false){
          *_log << ":FindCentroid: Singular matrix"<< endMsg;
          return UNRECOVERABLE_ERROR; // error placeholder
        }

	// Find the new set of parameters
        jerror_t error;
	if ((error=FindNewParmVec(N,X,F,J,par,newpar))!=NOERROR){
	  *_log << ":FindNewParmVec:  error=" << error << endMsg;

	  return error;
	}

        //Check for convergence
        for (i=0;i<3;i++){
          if (fabs(par(i,0)-newpar(i,0))>errx) errx=fabs(par(i,0)-newpar(i,0));
	}
        if (errx<TOLX){
          temp.pos=par(X0,0);
          temp.q=par(QA,0);
          temp.numstrips=3;
	  temp.t=(*peak)->t;
	  temp.t_rms=0.;
	  //CalcMeanTime(peak,temp.t,temp.t_rms);
          centroids.push_back(temp);
        
          return NOERROR;
        }
	par=newpar;
      } // iterations
    }
  }
  else{
    return VALUE_OUT_OF_RANGE;
  }
  return INFINITE_RECURSION; // error placeholder
}
     


///
/// DFDCPseudo_factory::FindNewParmVec()
///   Routine used by FindCentroid to backtrack along the direction
/// of the Newton step if the step is too large in order to avoid conditions
/// where the iteration procedure starts to diverge.
/// The procedure uses a scalar quantity f= 1/2 F.F for this purpose.
/// Algorithm described in Numerical Recipes in C, pp. 383-389.
///
 
jerror_t DFDCPseudo_factory::FindNewParmVec(DMatrix N,DMatrix X,DMatrix F,
				       DMatrix J,DMatrix par,
				       DMatrix &newpar){
  // Invert the J matrix
  DMatrix InvJ(DMatrix::kInverted,J);
  
  // Find the full Newton step
  DMatrix fullstep(3,1);
  fullstep=InvJ*F; fullstep*=-1;

  // find the rate of decrease for the Newton-Raphson step
  DMatrix slope(1,1);
  slope=(-1.0)*DMatrix(DMatrix::kTransposed,F)*F; //dot product

  // This should be a negative number...
  if (slope(0,0)>=0){
    return VALUE_OUT_OF_RANGE;
  }
  
  double lambda=1.0;  // Start out with full Newton step
  double lambda_temp,lambda2=lambda;
  DMatrix f(1,1),f2(1,1),newF(3,1),newf(1,1);

  // Compute starting values for f=1/2 F.F 
  f=-0.5*slope;

  for (;;){
    newpar=par+lambda*fullstep;

    // Compute the value of the vector F and f=1/2 F.F with the current step
    for (int i=0;i<3;i++){
      double argp=newpar(K2,0)*(newpar(X0,0)-X(i,0)+A_OVER_H);
      double argm=newpar(K2,0)*(newpar(X0,0)-X(i,0)-A_OVER_H);
      newF(i,0)=N(i,0)-newpar(QA,0)/4.*(tanh(argp)-tanh(argm));
    }
    newf=0.5*(DMatrix(DMatrix::kTransposed,newF)*newF); // dot product
 
    if (lambda<0.1) {  // make sure the step is not too small
      newpar=par;
      return NOERROR;
    } // Check if we have sufficient function decrease
    else if (newf(0,0)<=f(0,0)+ALPHA*lambda*slope(0,0)){
      return NOERROR;
    }
    else{
      // g(lambda)=f(par+lambda*fullstep)
      if (lambda==1.0){//first attempt: quadratic approximation for g(lambda)
        lambda_temp=-slope(0,0)/2./(newf(0,0)-f(0,0)-slope(0,0));
      }
      else{ // cubic approximation for g(lambda)
        double temp1=newf(0,0)-f(0,0)-lambda*slope(0,0);
        double temp2=f2(0,0)-f(0,0)-lambda2*slope(0,0);
        double a=(temp1/lambda/lambda-temp2/lambda2/lambda2)/(lambda-lambda2);
        double b=(-lambda2*temp1/lambda/lambda+lambda*temp2/lambda2/lambda2)
          /(lambda-lambda2);
        if (a==0.0) lambda_temp=-slope(0,0)/2./b;
        else{
          double disc=b*b-3.0*a*slope(0,0);
          if (disc<0.0) lambda_temp=0.5*lambda;
          else if (b<=0.0) lambda_temp=(-b+sqrt(disc))/3./a;
          else lambda_temp=-slope(0,0)/(b+sqrt(disc));
        }
        // ensure that we are headed in the right direction...
        if (lambda_temp>0.5*lambda) lambda_temp=0.5*lambda;
      }
    }
    lambda2=lambda;
    f2=newf;
    // Make sure that new version of lambda is not too small
    lambda=(lambda_temp>0.1*lambda ? lambda_temp : 0.1*lambda);
  } 
}

