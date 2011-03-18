//************************************************************************
// DFDCSegment_factory.cc - factory producing track segments from pseudopoints
//************************************************************************

#include "DFDCSegment_factory.h"
#include "DANA/DApplication.h"
//#include "HDGEOMETRY/DLorentzMapCalibDB.h"
#include <math.h>

#define HALF_CELL 0.5
#define MAX_DEFLECTION 0.15
#define EPS 1e-8
#define KILL_RADIUS 5.0 
#define Z_TARGET 65.0
#define Z_VERTEX_CUT 25.0
#define MATCH_RADIUS 5.0
#define ADJACENT_MATCH_RADIUS 1.0
#define SIGN_CHANGE_CHISQ_CUT 10.0
#define BEAM_VARIANCE 0.01 // cm^2
#define USED_IN_SEGMENT 0x8
#define CORRECTED 0x10
#define MAX_ITER 10
#define TARGET_LENGTH 30.0 //cm
#define MIN_TANL 2.0
#define ONE_THIRD  0.33333333333333333
#define SQRT3      1.73205080756887719

// Variance for position along wire using PHENIX angle dependence, transverse
// diffusion, and an intrinsic resolution of 127 microns.
inline double fdc_y_variance(double alpha,double x){
  return 0.00060+0.0064*tan(alpha)*tan(alpha)+0.0004*fabs(x);
}

bool DFDCSegment_package_cmp(const DFDCPseudo* a, const DFDCPseudo* b) {
  return a->wire->layer>b->wire->layer;
}

DFDCSegment_factory::DFDCSegment_factory() {
        _log = new JStreamLog(std::cout, "FDCSEGMENT >>");
        *_log << "File initialized." << endMsg;
}


///
/// DFDCSegment_factory::~DFDCSegment_factory():
/// default destructor -- closes log file
///
DFDCSegment_factory::~DFDCSegment_factory() {
        delete _log;	
}
///
/// DFDCSegment_factory::brun():
/// Initialization: read in deflection map, get magnetic field map
///
jerror_t DFDCSegment_factory::brun(JEventLoop* eventLoop, int eventNo) { 
  DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  bfield = dapp->GetBfield();
  lorentz_def=dapp->GetLorentzDeflections();

  *_log << "Table of Lorentz deflections initialized." << endMsg;
  return NOERROR;
}

///
/// DFDCSegment_factory::evnt():
/// Routine where pseudopoints are combined into track segments
///
jerror_t DFDCSegment_factory::evnt(JEventLoop* eventLoop, int eventNo) {
  myeventno=eventNo;

  vector<const DFDCPseudo*>pseudopoints;
  eventLoop->Get(pseudopoints);  
 
  // Skip segment finding if there aren't enough points to form a sensible 
  // segment 
  if (pseudopoints.size()>=3){
    // Group pseudopoints by package
    vector<const DFDCPseudo*>package[4];
    for (vector<const DFDCPseudo*>::const_iterator i=pseudopoints.begin();
	 i!=pseudopoints.end();i++){
      package[((*i)->wire->layer-1)/6].push_back(*i);
    }
    
    // Find the segments in each package
    for (int j=0;j<4;j++){
      std::sort(package[j].begin(),package[j].end(),DFDCSegment_package_cmp);
      // We need at least 3 points to define a circle, so skip if we don't 
      // have enough points.
      if (package[j].size()>2) FindSegments(package[j]);
    } 
  } // pseudopoints>2
  
  return NOERROR;
}

// Riemann Line fit: linear regression of s on z to determine the tangent of 
// the dip angle and the z position of the closest approach to the beam line.
// Also returns predicted positions along the helical path.
//
jerror_t DFDCSegment_factory::RiemannLineFit(vector<const DFDCPseudo *>points,
					     DMatrix &CR,vector<xyz_t>&XYZ){
  unsigned int n=points.size()+1;
  vector<int>bad(n);  // Keep track of "bad" intersection points
  // Fill matrix of intersection points
  for (unsigned int m=0;m<n-1;m++){
    double r2=points[m]->xy.Mod2();
    double denom= N[0]*N[0]+N[1]*N[1];
    double numer=dist_to_origin+r2*N[2];
    double ratio=numer/denom;

    DVector2 xy_int0(-N[0]*ratio,-N[1]*ratio);
    double temp=denom*r2-numer*numer;
    if (temp<0){    
      bad[m]=1;
      XYZ[m].xy=xy_int0;
      continue; 
    }
    temp=sqrt(temp)/denom;
    
    // Choose sign of square root based on proximity to actual measurements
    DVector2 delta(N[1]*temp,-N[0]*temp);
    DVector2 xy1=xy_int0+delta;
    DVector2 xy2=xy_int0-delta;
    double diff1=(xy1-points[m]->xy).Mod2();
    double diff2=(xy2-points[m]->xy).Mod2();
    if (diff1>diff2){
      XYZ[m].xy=xy2;
    }
    else{
      XYZ[m].xy=xy1;
    }
  }
  // Fake target point
  XYZ[n-1].xy.Set(0.,0.);

  // All arc lengths are measured relative to some reference plane with a hit.
  // Don't use a "bad" hit for the reference...
  unsigned int start=0;
  for (unsigned int i=0;i<bad.size();i++){
    if (!bad[i]){
      start=i;
      break;
    }
  }
  if ((start!=0 && ref_plane==0) || (start!=2&&ref_plane==2)) ref_plane=start;

  // Linear regression to find z0, tanl   
  double sumv=0.,sumx=0.;
  double sumy=0.,sumxx=0.,sumxy=0.;
  double sperp=0.,sperp_old=0.,ratio,Delta;
  double z=0,zlast=0;
  DVector2 oldxy=XYZ[start].xy;
  for (unsigned int k=start;k<n;k++){
    zlast=z;
    sperp_old=sperp;
    if (!bad[k]){
      DVector2 diffxy=XYZ[k].xy-oldxy;
      ratio=diffxy.Mod()/(2.*rc);
      // Make sure the argument for the arcsin does not go out of range...
      sperp=sperp_old+(ratio>1?2.*rc*(M_PI/2.):2.*rc*asin(ratio));
      //z=XYZ(k,2);
      z=XYZ[k].z;

      // Assume errors in s dominated by errors in R 
      double inv_var=1./CR(k,k);
      sumv+=inv_var;
      sumy+=sperp*inv_var;
      sumx+=z*inv_var;
      sumxx+=z*z*inv_var;
      sumxy+=sperp*z*inv_var;
      
      // Save the current x and y coordinates
      //oldx=XYZ(k,0);
      //oldy=XYZ(k,1);
      oldxy=XYZ[k].xy;
    }
  }
  Delta=sumv*sumxx-sumx*sumx;
  // Track parameters z0 and tan(lambda)
  tanl=-Delta/(sumv*sumxy-sumy*sumx); 
  //z0=(sumxx*sumy-sumx*sumxy)/Delta*tanl;
  // Error in tanl 
  var_tanl=sumv/Delta*(tanl*tanl*tanl*tanl);

  // Vertex position
  sperp-=sperp_old; 
  if (tanl<0){ 
    // a negative tanl makes no sense for FDC segments if we 
    // assume that the particle came from the target
    zvertex=Z_TARGET;
    tanl=(zlast-zvertex)/sperp;
  }
  else{
    zvertex=zlast-tanl*sperp;
  }

  return NOERROR;
}

// Use predicted positions (propagating from plane 1 using a helical model) to
// update the R and RPhi covariance matrices.
//
jerror_t DFDCSegment_factory::UpdatePositionsAndCovariance(unsigned int n,
     double r1sq,vector<xyz_t>&XYZ,DMatrix &CRPhi,DMatrix &CR){
  double delta_x=XYZ[ref_plane].xy.X()-xc; 
  double delta_y=XYZ[ref_plane].xy.Y()-yc;
  double r1=sqrt(r1sq);
  double denom=delta_x*delta_x+delta_y*delta_y;

  // Auxiliary matrices for correcting for non-normal track incidence to FDC
  // The correction is 
  //    CRPhi'= C*CRPhi*C+S*CR*S, where S(i,i)=R_i*kappa/2
  //                                and C(i,i)=sqrt(1-S(i,i)^2)  
  DMatrix C(n,n);
  DMatrix S(n,n);

  // Predicted positions
  Phi1=atan2(delta_y,delta_x);
  double z1=XYZ[ref_plane].z;
  double y1=XYZ[ref_plane].xy.X();
  double x1=XYZ[ref_plane].xy.Y();
  double var_R1=CR(ref_plane,ref_plane);
  for (unsigned int k=0;k<n;k++){       
    double sperp=charge*(XYZ[k].z-z1)/tanl;
    double sinp=sin(Phi1+sperp/rc);
    double cosp=cos(Phi1+sperp/rc);
    XYZ[k].xy.Set(xc+rc*cosp,yc+rc*sinp);
 
    // Error analysis.  We ignore errors in N because there doesn't seem to 
    // be any obvious established way to estimate errors in eigenvalues for 
    // small eigenvectors.  We assume that most of the error comes from 
    // the measurement for the reference plane radius and the dip angle anyway.
    double Phi=XYZ[k].xy.Phi();
    double sinPhi=sin(Phi);
    double cosPhi=cos(Phi);
    double dRPhi_dx=Phi*cosPhi-sinPhi;
    double dRPhi_dy=Phi*sinPhi+cosPhi;
    
    double dx_dx1=rc*sinp*delta_y/denom;
    double dx_dy1=-rc*sinp*delta_x/denom;
    double dx_dtanl=sinp*sperp/tanl;
    
    double dy_dx1=-rc*cosp*delta_y/denom;
    double dy_dy1=rc*cosp*delta_x/denom;
    double dy_dtanl=-cosp*sperp/tanl;
 
    double dRPhi_dx1=dRPhi_dx*dx_dx1+dRPhi_dy*dy_dx1;
    double dRPhi_dy1=dRPhi_dx*dx_dy1+dRPhi_dy*dy_dy1;
    double dRPhi_dtanl=dRPhi_dx*dx_dtanl+dRPhi_dy*dy_dtanl;

    double dR_dx1=cosPhi*dx_dx1+sinPhi*dy_dx1;
    double dR_dy1=cosPhi*dx_dy1+sinPhi*dy_dy1;
    double dR_dtanl=cosPhi*dx_dtanl+sinPhi*dy_dtanl;
      	
    double cdist=dist_to_origin+r1sq*N[2];
    double n2=N[0]*N[0]+N[1]*N[1];

    double ydenom=y1*n2+N[1]*cdist;
    double dy1_dr1=-r1*(2.*N[1]*N[2]*y1+2.*N[2]*cdist-N[0]*N[0])/ydenom;
    double var_y1=dy1_dr1*dy1_dr1*var_R1;

    double xdenom=x1*n2+N[0]*cdist;
    double dx1_dr1=-r1*(2.*N[0]*N[2]*x1+2.*N[2]*cdist-N[1]*N[1])/xdenom;
    double var_x1=dx1_dr1*dx1_dr1*var_R1;

    CRPhi(k,k)=dRPhi_dx1*dRPhi_dx1*var_x1+dRPhi_dy1*dRPhi_dy1*var_y1
      +dRPhi_dtanl*dRPhi_dtanl*var_tanl;
    CR(k,k)=dR_dx1*dR_dx1*var_x1+dR_dy1*dR_dy1*var_y1
      +dR_dtanl*dR_dtanl*var_tanl;
    
    // double stemp=sqrt(XYZ(k,0)*XYZ(k,0)+XYZ(k,1)*XYZ(k,1))/(4.*rc);
    double stemp=XYZ[k].xy.Mod()/(4.*rc);
    double ctemp=1.-stemp*stemp;
    if (ctemp>0){
      S(k,k)=stemp;
      C(k,k)=sqrt(ctemp);
    }
    else{
      S(k,k)=0.;
      C(k,k)=1.;
    }
  }
  
  // Correction for non-normal incidence of track on FDC 
  CRPhi=C*CRPhi*C+S*CR*S;
 
  return NOERROR;
}

// Riemann Circle fit:  points on a circle in x,y project onto a plane cutting
// the circular paraboloid surface described by (x,y,x^2+y^2).  Therefore the
// task of fitting points in (x,y) to a circle is transormed to the taks of
// fitting planes in (x,y, w=x^2+y^2) space
//
jerror_t DFDCSegment_factory::RiemannCircleFit(vector<const DFDCPseudo *>points,
					       DMatrix &CRPhi){
  unsigned int n=points.size()+1;
  DMatrix X(n,3);
  DMatrix Xavg(1,3);
  DMatrix A(3,3);
  // vector of ones
  DMatrix OnesT(1,n);
  double W_sum=0.;
  DMatrix W(n,n);

  // The goal is to find the eigenvector corresponding to the smallest 
  // eigenvalue of the equation
  //            lambda=n^T (X^T W X - W_sum Xavg^T Xavg)n
  // where n is the normal vector to the plane slicing the cylindrical 
  // paraboloid described by the parameterization (x,y,w=x^2+y^2),
  // and W is the weight matrix, assumed for now to be diagonal.
  // In the absence of multiple scattering, W_sum is the sum of all the 
  // diagonal elements in W.
  // At this stage we ignore the multiple scattering.
  for (unsigned int i=0;i<n-1;i++){
    X(i,0)=points[i]->xy.X();
    X(i,1)=points[i]->xy.Y();
    X(i,2)=points[i]->xy.Mod2();
    OnesT(0,i)=1.;
    W(i,i)=1./CRPhi(i,i);
    W_sum+=W(i,i);
  }
  unsigned int last_index=n-1;
  OnesT(0,last_index)=1.;
  W(last_index,last_index)=1./CRPhi(last_index,last_index);
  W_sum+=W(last_index,last_index);
  var_avg=1./W_sum;
  Xavg=var_avg*(OnesT*(W*X));
  // Store in private array for use in other routines
  xavg[0]=Xavg(0,0);
  xavg[1]=Xavg(0,1);
  xavg[2]=Xavg(0,2);
  
  A=DMatrix(DMatrix::kTransposed,X)*(W*X)
    -W_sum*(DMatrix(DMatrix::kTransposed,Xavg)*Xavg);
  if(!A.IsValid())return UNRECOVERABLE_ERROR;

  // The characteristic equation is 
  //   lambda^3+B2*lambda^2+lambda*B1+B0=0 
  //
  double B2=-(A(0,0)+A(1,1)+A(2,2));
  double B1=A(0,0)*A(1,1)-A(1,0)*A(0,1)+A(0,0)*A(2,2)-A(2,0)*A(0,2)
    +A(1,1)*A(2,2)-A(2,1)*A(1,2);
  double B0=-A.Determinant();
  if(B0==0 || !finite(B0))return UNRECOVERABLE_ERROR;

  // The roots of the cubic equation are given by 
  //        lambda1= -B2/3 + S+T
  //        lambda2= -B2/3 - (S+T)/2 + i sqrt(3)/2. (S-T)
  //        lambda3= -B2/3 - (S+T)/2 - i sqrt(3)/2. (S-T)
  // where we define some temporary variables:
  //        S= (R+sqrt(Q^3+R^2))^(1/3)
  //        T= (R-sqrt(Q^3+R^2))^(1/3)
  //        Q=(3*B1-B2^2)/9
  //        R=(9*B2*B1-27*B0-2*B2^3)/54
  //        sum=S+T;
  //        diff=i*(S-T)
  // We divide Q and R by a safety factor to prevent multiplying together 
  // enormous numbers that cause unreliable results.

  double Q=(3.*B1-B2*B2)/9.e4; 
  double R=(9.*B2*B1-27.*B0-2.*B2*B2*B2)/54.e6;
  double Q1=Q*Q*Q+R*R;
  if (Q1<0) Q1=sqrt(-Q1);
  else{
    return VALUE_OUT_OF_RANGE;
  }

  // DeMoivre's theorem for fractional powers of complex numbers:  
  //      (r*(cos(theta)+i sin(theta)))^(p/q)
  //                  = r^(p/q)*(cos(p*theta/q)+i sin(p*theta/q))
  //
  //double temp=100.*pow(R*R+Q1*Q1,0.16666666666666666667);
  double temp=100.*sqrt(cbrt(R*R+Q1*Q1));
  double theta1=ONE_THIRD*atan2(Q1,R);
  double sum_over_2=temp*cos(theta1);
  double diff_over_2=-temp*sin(theta1);
  // Third root
  double lambda_min=-ONE_THIRD*B2-sum_over_2+SQRT3*diff_over_2;

  // Normal vector to plane
  N[0]=1.;
  N[1]=(A(1,0)*A(0,2)-(A(0,0)-lambda_min)*A(1,2))
    /(A(0,1)*A(2,1)-(A(1,1)-lambda_min)*A(0,2));
  N[2]=(A(2,0)*(A(1,1)-lambda_min)-A(1,0)*A(2,1))
    /(A(1,2)*A(2,1)-(A(2,2)-lambda_min)*(A(1,1)-lambda_min));
  
  // Normalize: n1^2+n2^2+n3^2=1
  double denom=sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]);
  for (int i=0;i<3;i++){
    N[i]/=denom;
  }
      
  // Distance to origin
  dist_to_origin=-(N[0]*Xavg(0,0)+N[1]*Xavg(0,1)+N[2]*Xavg(0,2));

  // Center and radius of the circle
  xc=-N[0]/2./N[2];
  yc=-N[1]/2./N[2];
  rc=sqrt(1.-N[2]*N[2]-4.*dist_to_origin*N[2])/2./fabs(N[2]);

  return NOERROR;
}

// Riemann Helical Fit based on transforming points on projection to x-y plane 
// to a circular paraboloid surface combined with a linear fit of the arc 
// length versus z. Uses RiemannCircleFit and RiemannLineFit.
//
jerror_t DFDCSegment_factory::RiemannHelicalFit(vector<const DFDCPseudo*>points,
						DMatrix &CR,vector<xyz_t>&XYZ){
  double Phi;
  unsigned int num_points=points.size()+1;
  DMatrix CRPhi(num_points,num_points); 
 
  // Fill initial matrices for R and RPhi measurements
  XYZ[num_points-1].z=Z_TARGET;
  for (unsigned int m=0;m<points.size();m++){
    XYZ[m].z=points[m]->wire->origin.z();

    //Phi=atan2(points[m]->y,points[m]->x);
    Phi=points[m]->xy.Phi();
    double cosPhi=cos(Phi);
    double sinPhi=sin(Phi);
    double Phi_cosPhi=Phi*cosPhi;
    double Phi_sinPhi=Phi*sinPhi;
    CRPhi(m,m)
      =(Phi_cosPhi-sinPhi)*(Phi_cosPhi-sinPhi)*points[m]->covxx
      +(Phi_sinPhi+cosPhi)*(Phi_sinPhi+cosPhi)*points[m]->covyy
      +2.*(Phi_sinPhi+cosPhi)*(Phi_cosPhi-sinPhi)*points[m]->covxy;

    CR(m,m)=cosPhi*cosPhi*points[m]->covxx
      +sinPhi*sinPhi*points[m]->covyy
      +2.*sinPhi*cosPhi*points[m]->covxy;
  }
  CR(points.size(),points.size())=BEAM_VARIANCE;
  CRPhi(points.size(),points.size())=BEAM_VARIANCE;

  // Reference track:
  jerror_t error=NOERROR;  
  // First find the center and radius of the projected circle
  error=RiemannCircleFit(points,CRPhi); 
  if (error!=NOERROR) return error;
  
  // Get reference track estimates for z0 and tanl and intersection points
  // (stored in XYZ)
  error=RiemannLineFit(points,CR,XYZ);
  if (error!=NOERROR) return error;

  // Guess particle charge (+/-1);
  charge=GetCharge(points.size(),XYZ,CR,CRPhi);

  double r1sq=XYZ[ref_plane].xy.Mod2();
  UpdatePositionsAndCovariance(num_points,r1sq,XYZ,CRPhi,CR);
  
  // Preliminary circle fit 
  error=RiemannCircleFit(points,CRPhi); 
  if (error!=NOERROR) return error;

  // Preliminary line fit
  error=RiemannLineFit(points,CR,XYZ);
  if (error!=NOERROR) return error;

  // Guess particle charge (+/-1);
  charge=GetCharge(points.size(),XYZ,CR,CRPhi);
  
  r1sq=XYZ[ref_plane].xy.Mod2();
  UpdatePositionsAndCovariance(num_points,r1sq,XYZ,CRPhi,CR);
  
  // Final circle fit 
  error=RiemannCircleFit(points,CRPhi); 
  if (error!=NOERROR) return error;

  // Final line fit
  error=RiemannLineFit(points,CR,XYZ);
  if (error!=NOERROR) return error;

  // Guess particle charge (+/-1)
  charge=GetCharge(points.size(),XYZ,CR,CRPhi);
  
  // Final update to covariance matrices
  ref_plane=0;
  r1sq=XYZ[ref_plane].xy.Mod2();
  UpdatePositionsAndCovariance(num_points,r1sq,XYZ,CRPhi,CR);

  // Store residuals and path length for each measurement
  chisq=0.;
  for (unsigned int m=0;m<points.size();m++){
    double sperp=charge*(XYZ[m].z-XYZ[ref_plane].z)/tanl;
    double phi_s=Phi1+sperp/rc;
    DVector2 XY(xc+rc*cos(phi_s),yc+rc*sin(phi_s));
    chisq+=(XY-points[m]->xy).Mod2()/CR(m,m);
  }
  return NOERROR;
}


// DFDCSegment_factory::FindSegments
//  Associate nearest neighbors within a package with track segment candidates.
// Provide guess for (seed) track parameters
//
jerror_t DFDCSegment_factory::FindSegments(vector<const DFDCPseudo*>points){
  // Clear all the "used" flags;
  used.clear();

  // Put indices for the first point in each plane before the most downstream
  // plane in the vector x_list.
  double old_z=points[0]->wire->origin.z();
  vector<unsigned int>x_list;
  x_list.push_back(0);
  for (unsigned int i=0;i<points.size();i++){
    used.push_back(false);
    if (points[i]->wire->origin.z()!=old_z){
      x_list.push_back(i);
    }
    old_z=points[i]->wire->origin.z();
  }
  x_list.push_back(points.size()); 

  unsigned int start=0;
  // loop over the start indices, starting with the first plane
  while (start<x_list.size()-1){
    int num_used=0;
    // For each iteration, count how many hits we have used already in segments
    for (unsigned int i=0;i<points.size();i++){
      if(used[i]==true) num_used++;
    }
    // Break out of loop if we don't have enough hits left to fit a circle
    if (points.size()-num_used<3) break;

    // Now loop over the list of track segment start points
    for (unsigned int i=x_list[start];i<x_list[start+1];i++){
      if (used[i]==false){
	used[i]=true;
	chisq=1.e8;

	// Clear track parameters
	tanl=D=z0=phi0=Phi1=xc=yc=rc=0.;
	charge=1.;
	
	// Point in the current plane in the package 
	DVector2 XY=points[i]->xy;
	
	// Create list of nearest neighbors
	vector<const DFDCPseudo*>neighbors;
	neighbors.push_back(points[i]);
	unsigned int match=0;
	double delta,delta_min=1000.;
	for (unsigned int k=0;k<x_list.size()-1;k++){
	  delta_min=1000.;
	  match=0;
	  for (unsigned int m=x_list[k];m<x_list[k+1];m++){
	    delta=(XY-points[m]->xy).Mod();
	    if (delta<delta_min && delta<MATCH_RADIUS){
	      delta_min=delta;
	      match=m;
	    }
	  }	
	  if (match!=0 
	      && used[match]==false
	      ){
	    XY=points[match]->xy;
	    used[match]=true;
	    neighbors.push_back(points[match]);	  
	  }
	}
	unsigned int num_neighbors=neighbors.size();
	
	// Skip to next segment seed if we don't have enough points to fit a 
	// circle
	if (num_neighbors<3) continue;    

	bool do_sort=false;
	// Look for hits adjacent to the ones we have in our segment candidate
	for (unsigned int k=0;k<points.size();k++){
	  for (unsigned int j=0;j<num_neighbors;j++){
	    delta=(points[k]->xy-neighbors[j]->xy).Mod();
	    if (delta<ADJACENT_MATCH_RADIUS && 
		abs(neighbors[j]->wire->wire-points[k]->wire->wire)==1
		&& neighbors[j]->wire->origin.z()==points[k]->wire->origin.z()){
	      used[k]=true;
	      neighbors.push_back(points[k]);
	      do_sort=true;
	    }      
	  }
	}
	// Sort if we added another point
	if (do_sort)
	  std::sort(neighbors.begin(),neighbors.end(),DFDCSegment_package_cmp);
    
	// list of points on track and the corresponding covariances
	vector<xyz_t>XYZ(neighbors.size()+1);
	DMatrix CR(neighbors.size()+1,neighbors.size()+1);
   	
	// Arc lengths in helical model are referenced relative to the plane
	// ref_plane within a segment.  For a 6 hit segment, ref_plane=2 is 
	// roughly in the center of the segment.
	ref_plane=2; 
	
	// Perform the Riemann Helical Fit on the track segment
	jerror_t error=RiemannHelicalFit(neighbors,CR,XYZ);   /// initial hit-based fit
		
	if (error==NOERROR){  
	  // Estimate for azimuthal angle
	  phi0=atan2(-xc,yc); 
	  if (charge<0) phi0+=M_PI;
	  
	  // Look for distance of closest approach nearest to target
	  D=-charge*rc-xc/sin(phi0);
	  
	  // Creat a new segment
	  DFDCSegment *segment = new DFDCSegment;	

	  // Initialize seed track parameters
	  segment->q=charge; //charge 
	  segment->phi0=phi0;      // Phi
	  segment->D=D;       // D=distance of closest approach to origin   
	  segment->tanl=tanl;     // tan(lambda), lambda=dip angle
	  segment->z_vertex=zvertex;// z-position at closest approach to origin

	  segment->hits=neighbors;
	  segment->xc=xc;
	  segment->yc=yc;
	  segment->rc=rc;
	  segment->Phi1=Phi1;
	  segment->chisq=chisq;
	  
	  _data.push_back(segment);
	}
      }
    }
    // Look for a new plane to start looking for a segment
    while (start<x_list.size()-1){
      if (used[x_list[start]]==false) break;
      start++;
    }
  } //while loop over x_list planes

  return NOERROR;
}

// Linear regression to find charge
double DFDCSegment_factory::GetCharge(unsigned int n,vector<xyz_t>&XYZ, 
				      DMatrix &CR, 
				      DMatrix &CRPhi){
  double inv_var=0.; 
  double sumv=0.;
  double sumy=0.;
  double sumx=0.;
  double sumxx=0.,sumxy=0,Delta;
  double slope,r2;
  double phi_old=XYZ[0].xy.Phi();
  for (unsigned int k=0;k<n;k++){   
    double tempz=XYZ[k].z;
    double phi_z=XYZ[k].xy.Phi();
    // Check for problem regions near +pi and -pi
    if (fabs(phi_z-phi_old)>M_PI){  
      if (phi_old<0) phi_z-=2.*M_PI;
      else phi_z+=2.*M_PI;
    }
    r2=XYZ[k].xy.Mod2();
    inv_var=r2/(CRPhi(k,k)+phi_z*phi_z*CR(k,k));
    sumv+=inv_var;
    sumy+=phi_z*inv_var;
    sumx+=tempz*inv_var;
    sumxx+=tempz*tempz*inv_var;
    sumxy+=phi_z*tempz*inv_var;
    phi_old=phi_z;
  }
  Delta=sumv*sumxx-sumx*sumx;
  slope=(sumv*sumxy-sumy*sumx)/Delta; 
  
  // Guess particle charge (+/-1);
  if (slope<0.) return -1.;
  return 1.;
}

//----------------------------------------------------------------------------
// The following routine is no longer in use: 
// Correct avalanche position along wire and incorporate drift data for 
// coordinate away from the wire using results of preliminary hit-based fit
//
//#define R_START 7.6
//#define Z_TOF 617.4
//#include "HDGEOMETRY/DLorentzMapCalibDB.h
//#define SC_V_EFF 15.
//#define SC_LIGHT_GUIDE     140.
//#define SC_CENTER          38.75   
//#define TOF_BAR_LENGTH      252.0 
//#define TOF_V_EFF 15.
//#define FDC_X_RESOLUTION 0.028
//#define FDC_Y_RESOLUTION 0.02 //cm
/*
jerror_t DFDCSegment_factory::CorrectPoints(vector<DFDCPseudo*>points,
					    DMatrix XYZ){
  // dip angle
  double lambda=atan(tanl);  
  double alpha=M_PI/2.-lambda;
  
  if (alpha == 0. || rc==0.) return VALUE_OUT_OF_RANGE;

  // Get Bfield, needed to guess particle momentum
  double Bx,By,Bz,B;
  double x=XYZ(ref_plane,0);
  double y=XYZ(ref_plane,1);
  double z=XYZ(ref_plane,2);
  
  bfield->GetField(x,y,z,Bx,By,Bz);
  B=sqrt(Bx*Bx+By*By+Bz*Bz);

  // Momentum and beta
  double p=0.002998*B*rc/cos(lambda);
  double beta=p/sqrt(p*p+0.140*0.140);

  // Correct start time for propagation from (0,0)
  double my_start_time=0.;
  if (use_tof){
    //my_start_time=ref_time-(Z_TOF-Z_TARGET)/sin(lambda)/beta/29.98;
    // If we need to use the tof, the angle relative to the beam line is 
    // small enough that sin(lambda) ~ 1.
    my_start_time=ref_time-(Z_TOF-Z_TARGET)/beta/29.98;
    //my_start_time=0;
  }
  else{
    double ratio=R_START/2./rc;
    if (ratio<=1.)
      my_start_time=ref_time
	-2.*rc*asin(R_START/2./rc)*(1./cos(lambda)/beta/29.98);
    else
      my_start_time=ref_time
	-rc*M_PI*(1./cos(lambda)/beta/29.98);
    
  }
  //my_start_time=0.;

  for (unsigned int m=0;m<points.size();m++){
    DFDCPseudo *point=points[m];

    double cosangle=point->wire->udir(1);
    double sinangle=point->wire->udir(0);
    x=XYZ(m,0);
    y=XYZ(m,1);
    z=point->wire->origin.z();
    double delta_x=0,delta_y=0;   
    // Variances based on expected resolution
    double sigx2=FDC_X_RESOLUTION*FDC_X_RESOLUTION;
   
    // Find difference between point on helical path and wire
    double w=x*cosangle-y*sinangle-point->w;
     // .. and use it to determine which sign to use for the drift time data
    double sign=(w>0?1.:-1.);

    // Correct the drift time for the flight path and convert to distance units
    // assuming the particle is a pion
    delta_x=sign*(point->time-fdc_track[m].s/beta/29.98-my_start_time)*55E-4;

    // Variance for position along wire. Includes angle dependence from PHENIX
    // and transverse diffusion
    double sigy2=fdc_y_variance(alpha,delta_x);
    
    // Next find correction to y from table of deflections
    delta_y=lorentz_def->GetLorentzCorrection(x,y,z,alpha,delta_x);

    // Fill x and y elements with corrected values
    point->ds =-delta_y;     
    point->dw =delta_x;
    point->x=(point->w+point->dw)*cosangle+(point->s+point->ds)*sinangle;
    point->y=-(point->w+point->dw)*sinangle+(point->s+point->ds)*cosangle;
    point->covxx=sigx2*cosangle*cosangle+sigy2*sinangle*sinangle;
    point->covyy=sigx2*sinangle*sinangle+sigy2*cosangle*cosangle;
    point->covxy=(sigy2-sigx2)*sinangle*cosangle;
    point->status|=CORRECTED;
  }
  return NOERROR;
}
*/
