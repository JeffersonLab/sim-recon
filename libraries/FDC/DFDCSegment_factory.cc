//************************************************************************
// DFDCSegment_factory.cc - factory producing track segments from pseudopoints
//************************************************************************

#include "DFDCSegment_factory.h"
#include "TRACKING/DQuickFit.h"
#include <math.h>
#define HALF_CELL 0.5
#define MAX_DEFLECTION 0.15
#define EPS 1e-3
#define KILL_RADIUS 5.0 
#define NOT_CORRECTED 0x8
#define Z_TARGET 65.0
#define MATCH_RADIUS 5.0

static bool got_deflection_file=true;
static bool warn=true;

bool DFDCSegment_package_cmp(const DFDCPseudo* a, const DFDCPseudo* b) {
  return a->wire->layer>b->wire->layer;
}


// Locate a position in array xx given x
void locate(double *xx,int n,double x,int *j){
  int ju,jm,jl;
  int ascnd;
  
  jl=-1;
  ju=n;
  ascnd=(xx[n-1]>=xx[0]);
  while(ju-jl>1){
    jm=(ju+jl)>>1;
    if (x>=xx[jm]==ascnd)
      jl=jm;
    else
      ju=jm;
  }
  if (x==xx[0]) *j=0;
  else if (x==xx[n-1]) *j=n-2;
  else *j=jl; 
}


// Polynomial interpolation on a grid.
// Adapted from Numerical Recipes in C (2nd Edition), pp. 121-122.
void polint(double *xa, double *ya,int n,double x, double *y,double *dy){
  int i,m,ns=0;
  double den,dif,dift,ho,hp,w;
  
  double *c= new double[n];
  double *d= new double[n];

  dif=fabs(x-xa[0]);
  for (i=0;i<n;i++){
    if ((dift=fabs(x-xa[i]))<dif){
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];

  for (m=1;m<n;m++){
    for (i=1;i<=n-m;i++){
      ho=xa[i-1]-x;
      hp=xa[i+m-1]-x;
      w=c[i+1-1]-d[i-1];
      if ((den=ho-hp)==0.0) return;
      
      den=w/den;
      d[i-1]=hp*den;
      c[i-1]=ho*den;      
    }    
    *y+=(*dy=(2*ns<(n-m) ?c[ns+1]:d[ns--]));
  }
  delete []c;
  delete []d;
}

DFDCSegment_factory::DFDCSegment_factory() {
        logFile = new ofstream("DFDCSegment_factory.log");
        _log = new JStreamLog(*logFile, "SEGMENT");
        *_log << "File initialized." << endMsg;
}


///
/// DFDCSegment_factory::~DFDCSegment_factory():
/// default destructor -- closes log file
///
DFDCSegment_factory::~DFDCSegment_factory() {
        logFile->close();
        delete logFile;
        delete _log;
}
///
/// DFDCSegment_factory::brun():
/// Initialization: read in deflection map
///
jerror_t DFDCSegment_factory::brun(JEventLoop* eventLoop, int eventNo) { 
  ifstream fdc_deflection_file;

  fdc_deflection_file.open("fdc_deflections.dat");
  if (!fdc_deflection_file.is_open()){ 
    *_log << "Failed to open fdc_deflection file." << endMsg;
    got_deflection_file=false;
    return RESOURCE_UNAVAILABLE;
  }

  double bx;
  double bz;
  BField[0]=BField[1]=BField[2]=BField[3]=0.;

  char buf[80];
  fdc_deflection_file.getline(buf,80,'\n');
  fdc_deflection_file.setf(ios::skipws); 
  for (int i=0;i<LORENTZ_X_POINTS;i++){
    for (int j=0;j<LORENTZ_Z_POINTS;j++){
      fdc_deflection_file >> lorentz_x[i];
      fdc_deflection_file >> lorentz_z[j];
      fdc_deflection_file >> bx;
      fdc_deflection_file >> bz;
      fdc_deflection_file >> lorentz_nx[i][j];
      fdc_deflection_file >> lorentz_nz[i][j];  

    BField[j/PACKAGE_Z_POINTS]+=sqrt(bx*bx+bz*bz)
	/double(PACKAGE_Z_POINTS*LORENTZ_X_POINTS);  
    }
  }
  fdc_deflection_file.close();

  *_log << "Table of Lorentz deflections initialized." << endMsg;
  return NOERROR;
}

///
/// DFDCSegment_factory::evnt():
/// Routine where pseudopoints are converted into space points
///
jerror_t DFDCSegment_factory::evnt(JEventLoop* eventLoop, int eventNo) {
  vector<const DFDCPseudo*>pseudopoints;
  eventLoop->Get(pseudopoints);  

  // Group pseudopoints by package
  vector<DFDCPseudo*>package[4];
  for (vector<const DFDCPseudo*>::iterator i=pseudopoints.begin();
       i!=pseudopoints.end();i++){
    package[((*i)->wire->layer-1)/6].push_back((DFDCPseudo*)(*i));

  }
  for (int j=0;j<4;j++){
    std::sort(package[j].begin(), package[j].end(), DFDCSegment_package_cmp);
    KalmanFilter(package[j]);
  
  }

  // Correct for the Lorentz effect given DOCAs
  if (got_deflection_file){
    for (unsigned int k=0;k<_data.size();k++){
      DFDCSegment *segment=_data[k];
      for (unsigned int m=0;m<segment->hits.size();m++)
	CorrectPointY(segment->S,segment->hits[m]);
    }
  }
  else if (warn){
    fprintf(stderr,
	    "DFDCSegment_factory: No Lorentz deflection file found!\n");
      fprintf(stderr,"DFDCSegment_factory: --> Pseudo-points will not be corrected for Lorentz effect.");
      warn=false;      
  }


  return NOERROR;
}

//
// Kalman Filter
//
jerror_t DFDCSegment_factory::KalmanFilter(vector<DFDCPseudo*>points){
  // Check that we have data to fit
  if (points.size()<2) return OBJECT_NOT_AVAILABLE;
  DMatrix Seed(5,1); //  Seed track state vector
  DMatrix S(5,1),Old_S(5,1);  // State vector 
  DMatrix M(2,1);  // measurement vector
  DMatrix M_pred(2,1); // prediction for hit position
  DMatrix H(2,5);  // Track projection matrix
  DMatrix H_T(5,2); // Transpose of track projection matrix
  DMatrix F(5,5);  // State vector propagation matrix
  DMatrix F_T(5,5); // Transpose of state vector propagation matrix
  DMatrix C(5,5);  // State vector covariance matrix
  DMatrix Q(5,5);  // Process noise covariance matrix
  DMatrix K(5,2);  // Kalman gain matrix
  DMatrix V(2,2);  // Measurement covariance matrix
  DMatrix X(2,1);  // Position on helical track
 
  // Put indices for the first point in each plane before the most downstream
  // plane in the vector x_list.
  double old_z=points[0]->wire->origin(2);
  vector<unsigned int>x_list;
  for (unsigned int i=0;i<points.size();i++){
    if (points[i]->wire->origin(2)!=old_z){
      x_list.push_back(i);
    }
    old_z=points[i]->wire->origin(2);
  }
  x_list.push_back(points.size()); 

  // Now loop over the list of track segment end points
  for (unsigned int i=0;i<x_list[0];i++){
    // Creat a new segment
    DFDCSegment *segment = new DFDCSegment;
 
    // Put first point in list of hits belonging to this segment
    segment->hits.push_back(points[i]);

    // Point in the last plane in the package 
    double x=points[i]->x;
    double y=points[i]->y; 
    double z=points[i]->wire->origin(2);
    double r=sqrt(x*x+y*y);
    double Phi0=atan2(y,x),Phi1;
  
    // Create list of nearest neighbors
    vector<DFDCPseudo*>neighbors;
    neighbors.push_back(points[i]);
    unsigned int match=0;
    double delta,delta_min=1000.,y1=0,x1=0,xtemp,ytemp;
    for (unsigned int k=0;k<x_list.size()-1;k++){
      delta_min=1000.;
      match=0;
      for (unsigned int m=x_list[k];m<x_list[k+1];m++){
	xtemp=points[m]->x;
	ytemp=points[m]->y;
	delta=sqrt((x-xtemp)*(x-xtemp)+(y-ytemp)*(y-ytemp));
        if (delta<delta_min && delta<MATCH_RADIUS){
          delta_min=delta;
          match=m;
        }
      }	
      if (match!=0){
        x=points[match]->x;
        y=points[match]->y;
        neighbors.push_back(points[match]);
      }
    }
    x1=neighbors[neighbors.size()-1]->x;
    y1=neighbors[neighbors.size()-1]->y;
    Phi1=atan2(y1,x1);

    // Guess particle charge (+/-1);
    double q=1.;     
    if (Phi0-Phi1<0) q=-1.;
    
    // Guess for the center and radius of the circle from the projection of 
    // the helical path
    double xc,yc,rc;
    DQuickFit *fit= new DQuickFit();
    for (unsigned int k=0;k<neighbors.size();k++){
      fit->AddHitXYZ(neighbors[k]->x,neighbors[k]->y,
		     neighbors[k]->wire->origin(2));
    }
    fit->FitCircle();
    xc=fit->x0;
    yc=fit->y0;
    rc=sqrt(xc*xc+yc*yc);
  
    // Initialization of covariance matrix 
    for (int j=0;j<5;j++){
      C(j,j)=1.;
    }
    
    // Initial guess for state vector
    double kappa,phi,tanl=0.,sperp,z0;
    kappa=q/2./rc;
    phi=atan2(-xc,yc);
    if (q<0) phi+=M_PI;
    
    // Linear regression to find z0, tanl   
    double sigx2=HALF_CELL*HALF_CELL/3.;
    double sigy2=MAX_DEFLECTION*MAX_DEFLECTION/3.;
    double var_s=sigx2+sigy2; 
    double sum_var=1./var_s;
    double sum_s=0;
    double sum_z=Z_TARGET/var_s;
    double sum_zz=Z_TARGET*Z_TARGET/var_s,sum_sz=0,Delta;
    double cosangle,sinangle;
    for (unsigned int k=0;k<neighbors.size();k++){   
      double tempz=neighbors[k]->wire->origin(2);
      r=sqrt(neighbors[k]->x*neighbors[k]->x+neighbors[k]->y*neighbors[k]->y);
      double ratio=r/2./rc; 
      // Make sure the argument for the arcsin does not go out of range...
      if (ratio>1.) 
	sperp=2.*rc*(M_PI/2.);
      else
	sperp=2.*rc*asin(ratio);
      // Assume errors in sperp domininated by errors in x and y at the current
      // z-position.
      var_s=(neighbors[k]->w*neighbors[k]->w*sigx2
	     +neighbors[k]->s*neighbors[k]->s*sigy2)/r/r/(1.-r*r/4./rc/rc);
      sum_var+=1./var_s;
      sum_s+=sperp/var_s;
      sum_z+=tempz/var_s;
      sum_zz+=tempz*tempz/var_s;
      sum_sz+=sperp*tempz/var_s;
    }
    Delta=sum_var*sum_zz-sum_z*sum_z;
    tanl=Delta/(sum_var*sum_sz-sum_s*sum_z); 
    z0=-(sum_zz*sum_s-sum_z*sum_sz)/Delta*tanl;

    Seed(0,0)=kappa;    // Curvature (kappa)
    Seed(1,0)=phi;      // Phi
    Seed(2,0)=0.;       // D   
    Seed(3,0)=tanl;     // tan(lambda)
    Seed(4,0)=z0; // z0
    S=Seed;
    
    // Loop over list of neighboring points, updating the state vector at each
    // step using the Kalman algorithm
    for (unsigned int k=0;k<neighbors.size();k++){      
      // The next measurement 
      M(0,0)=neighbors[k]->x;
      M(1,0)=neighbors[k]->y; 
     
      // ... and its covariance matrix taking into account rotation    
      cosangle=neighbors[k]->wire->udir(1);
      sinangle=neighbors[k]->wire->udir(0);
      V(0,0)=sigx2*cosangle*cosangle+sigy2*sinangle*sinangle;
      V(1,1)=sigx2*sinangle*sinangle+sigy2*cosangle*cosangle;
      V(0,1)=V(1,0)=(sigy2-sigx2)*sinangle*cosangle;
           
      // Put the current point in the list of pseudopoints belonging to this
      // segment.
      //segment->hits.push_back(neighbors[k]);
      old_z=z;
      z=neighbors[k]->wire->origin(2);

      // Compute Process noise covariance matrix
      GetProcessNoiseCovariance(S,segment->hits,Q);	

      // position on the seed trajectory
      double xold,yold;
      GetHelicalTrackPosition(old_z,Seed,xold,yold);  
      X(0,0)=xold;
      X(1,0)=yold;

      // position on current trajectory
      GetHelicalTrackPosition(z,S,x,y); 
       
      // Get the state vector from the non-linear state propagation function
      DMatrix S1(5,1);
      GetStateVector(xold,yold,old_z,x,y,z,S,S1);
 
      // Get state transport matrix describing how S evolves
      GetStateTransportMatrix(xold,yold,x,y,S,F);

      // Transport the covariance matrix
      F_T=DMatrix(DMatrix::kTransposed,F);
      C=F*(C*F_T)+Q;
	
      // Get predicted estimate for S
      S=Seed+F*(S1-S);
    
      // Get the projection matrix to propagate the track to the next z-plane
      GetTrackProjectionMatrix(z,S1,H);
      // Get transpose of projection matrix
      H_T = DMatrix(DMatrix::kTransposed,H);
    
      // Get prediction for next measurement
      M_pred=X+H*(S-Seed);
      
      // Compute new V and check that it is invertible
      DMatrix V_new(2,2);
      V_new=V+H*(C*H_T);
      TDecompLU lu(V_new);
      if (lu.Decompose()==false){
	*_log << "Singular matrix"<< endMsg;
	continue;
      }
      
      // Compute Kalman gain matrix
      K=C*(H_T*DMatrix(DMatrix::kInverted,V_new));

      // Update the state vector 
      S=S+K*(M-M_pred); 
      
      // Update state vector covariance matrix
      C=C-K*(H*C);
    }
      
    segment->S.ResizeTo(S);
    segment->S=S;
    segment->cov.ResizeTo(C);
    segment->cov=C;
    segment->hits=neighbors;
    
    _data.push_back(segment);
  }
  return NOERROR;
}

// Get state vector in translated coordinate system
jerror_t DFDCSegment_factory::GetStateVector(double xold, double yold,
    double old_z,double x, double y,double z,DMatrix S,DMatrix &S1){
  double phi=S(1,0);
  double D=S(2,0);
  double dx=x-xold+D*sin(phi);
  double dy=y-yold-D*cos(phi);
  double dz=z-old_z;
  double kappa=S(0,0);
  double xc=dx-(D+1./2./kappa)*sin(phi);
  double yc=dy+(D+1./2./kappa)*cos(phi);
  double q=kappa/fabs(kappa);
  double rc=sqrt(xc*xc+yc*yc);
  double tanl=S(3,0);
  double z0=S(4,0);
  double phi1=atan2(-xc,yc);

  if (q<0) phi1+=M_PI;

  S1(0,0)=S(0,0);
  S1(1,0)=phi1;
  S1(2,0)=q*rc-1./2./kappa;
  S1(3,0)=tanl;
  S1(4,0)=z0+dz+tanl*(phi1-phi)/2./kappa;
  
  return NOERROR;
}


// Computation of transport matrix for S
jerror_t DFDCSegment_factory:: GetStateTransportMatrix(double xold, 
	double yold, double x, double y,DMatrix S, DMatrix &F){
  // Old version of state vector 
  double tanl=S(3,0);      
  double kappa=S(0,0);
  double phi=S(1,0);
  double D=S(2,0);

  double dx=x-xold+D*sin(phi);
  double dy=y-yold-D*cos(phi);
  double xc=dx-(D+1./2./kappa)*sin(phi);
  double yc=dy+(D+1./2./kappa)*cos(phi);
  double phi1=atan2(-xc,yc);
  double q=kappa/fabs(kappa);
  double rc2=xc*xc+yc*yc;
  double rc=sqrt(rc2);

  if (q<0) phi1+=M_PI;

  // Some temperorary variables for computing F
  double a=yc*sin(phi)+xc*cos(phi);
  double c=yc*cos(phi)-xc*sin(phi);
  
  // Transport matrix for S
  F(0,0)=F(3,3)=F(4,4)=1.;
  F(1,0)=-a/2./kappa/kappa/rc2;
  F(1,1)= (D+1./2./kappa)*c/rc2;
  F(1,2)=a/rc2;
  F(2,0)=(1.-q*c/rc)/2./kappa/kappa;
  F(2,1)=-q*a/rc*(D+1./2./kappa);
  F(2,2)=q*c/rc;
  F(4,0)=tanl/2./kappa*((phi-phi1)/kappa+F(0,1));
  F(4,1)=tanl/2./kappa*(F(1,1)-1.);
  F(4,2)=tanl/2./kappa*F(1,2);
  F(4,3)=(phi1-phi)/2./kappa;

  return NOERROR;

}


jerror_t DFDCSegment_factory::GetProcessNoiseCovariance(DMatrix S,
			    vector<DFDCPseudo*>points,DMatrix &Q){
  if (!got_deflection_file) return RESOURCE_UNAVAILABLE;

  // Track parameters
  double tanl=S(3,0);
  double cosl=cos(atan(tanl));
  double kappa=S(0,0);
  double D=S(2,0);
  double phi=S(1,0);
  double sigma2_ms=0;  	

  // Momentum
  double p=0.003*BField[points[0]->wire->layer/PACKAGE_Z_POINTS]
    /2./fabs(kappa)/cosl;
  double beta=p/sqrt(p*p+0.140*0.140); // assume pion 

  //Materials: copper, Kapton, Mylar
  double thickness[3]={4e-4,50e-4,13e-4};
  double density[3]={8.96,1.42,1.39};
  double X0[3]={12.86,40.56,39.95};
  double material_sum=0.;
  for (int i=0;i<3;i++){
    material_sum+=thickness[i]*density[i]/X0[i];
  }	
  // RMS from multiple scattering
  sigma2_ms=0.0136*0.0136/p/p/beta/beta*material_sum;

  Q(0,0)=kappa*kappa*tanl*tanl;
  Q(0,3)=Q(3,0)=kappa*tanl*(1.+tanl*tanl);
  Q(1,1)=1./cosl/cosl;
  Q(2,2)=D*D/tan(phi)/tan(phi)/cosl/cosl;
  Q(3,3)=(1.+tanl*tanl)*(1.+tanl*tanl);
  
  Q*=sigma2_ms;
  
  return NOERROR;
}

// Compute track projection matrix
jerror_t DFDCSegment_factory::GetTrackProjectionMatrix(double z,
						       DMatrix S,DMatrix &H){
  double tanl=S(3,0);      
  double kappa=S(0,0);
  double phi=S(1,0);
  double D=S(2,0);
  double z0=S(4,0);
  double sperp=(z-z0)/tanl;

  H(0,0)=1./2./kappa/kappa*(cos(phi)*(2.*kappa*sperp*cos(2.*kappa*sperp)
				      -sin(2.*kappa*sperp))
			    +sin(phi)*(1.-cos(2.*kappa*sperp)
				       -2.*kappa*sperp*sin(2.*kappa*sperp)));
  H(1,0)=1./2./kappa/kappa*(sin(phi)*(2.*kappa*sperp*cos(2.*kappa*sperp)
				      -sin(2.*kappa*sperp))
			    -cos(phi)*(1.-cos(2.*kappa*sperp)
				       -2.*kappa*sperp*sin(2.*kappa*sperp)));
  H(0,1)=-D*cos(phi)-1./2./kappa*sin(phi)*sin(2.*kappa*sperp)
    -1./2./kappa*cos(phi)*(1.-cos(2.*kappa*sperp));  
  H(1,1)=-D*sin(phi)+1./2./kappa*cos(phi)*sin(2.*kappa*sperp)
    -1./2./kappa*sin(phi)*(1.-cos(2.*kappa*sperp));
  H(0,2)=-sin(phi);
  H(1,2)=cos(phi);
  H(0,3)=(sin(phi)*sin(2.*kappa*sperp)-cos(phi)*cos(2.*kappa*sperp))*sperp
    /tanl;
  H(1,3)=-(sin(phi)*sin(2.*kappa*sperp)+cos(phi)*cos(2.*kappa*sperp))*sperp
    /tanl;
  H(0,4)=(sin(phi)*sin(2.*kappa*sperp)-cos(phi)*cos(2.*kappa*sperp))/tanl;
  H(1,4)=-(sin(phi)*sin(2.*kappa*sperp)+cos(phi)*cos(2.*kappa*sperp))/tanl;
  
  return NOERROR;
}


jerror_t DFDCSegment_factory::GetHelicalTrackPosition(double z,DMatrix S,
			  double  &xpos, double &ypos){
  // Track parameters
  double kappa=S(0,0);
  double phi=S(1,0);
  double D=S(2,0);
  double tanl=S(3,0);
  double z0=S(4,0);
  double sperp=(z-z0)/tanl;
  
  xpos=-D*sin(phi)+(1./2./kappa)*cos(phi)*sin(2.*kappa*sperp)
    -(1./2./kappa)*sin(phi)*(1.-cos(2.*kappa*sperp)); 
  ypos=D*cos(phi)+(1./2./kappa)*sin(phi)*sin(2.*kappa*sperp)
    +(1./2./kappa)*cos(phi)*(1.-cos(2.*kappa*sperp));

  return NOERROR;
}

// Correct avalanche position along wire from doca data
jerror_t DFDCSegment_factory::CorrectPointY(DMatrix S,
					  DFDCPseudo *point){
  
  DVector3 P1,P2,P3,R;
  double z=point->wire->origin(2);
  double cosangle=point->wire->udir(1);
  double sinangle=point->wire->udir(0);
  double x=point->x;
  double y=point->y;
  double delta_x=0,delta_y=0;
  double xhelix,yhelix;
 
  // Find the distance of closest approach between the point on the track
  // segment and the wire
  P1=point->wire->origin;
  P2.SetXYZ(x,y,z);
  GetHelicalTrackPosition(z,S,xhelix,yhelix);
  P3.SetXYZ(xhelix,yhelix,z);
  R=P1-P3+((P3-P1).Dot(P2-P1))/(P2-P1).Dot(P2-P1)*(P2-P1);
  delta_x=R(0)*cosangle-R(1)*sinangle;

  // Next find correction to y from table of deflections
  double tanr,tanz,r=sqrt(x*x+y*y);
  double phi=atan2(y,x);
  double ytemp[LORENTZ_X_POINTS],ytemp2[LORENTZ_X_POINTS],dummy;
  int imin,imax, ind,ind2;

  // Locate positions in x and z arrays given r and z
  locate(lorentz_x,LORENTZ_X_POINTS,r,&ind);
  locate(lorentz_z,LORENTZ_Z_POINTS,z,&ind2);

  // First do interpolation in z direction 
  imin=PACKAGE_Z_POINTS*(ind2/PACKAGE_Z_POINTS); // Integer division...
  for (int j=0;j<LORENTZ_X_POINTS;j++){
    polint(&lorentz_z[imin],&lorentz_nx[j][imin],PACKAGE_Z_POINTS,z,
	   &ytemp[j],&dummy);
    polint(&lorentz_z[imin],&lorentz_nz[j][imin],PACKAGE_Z_POINTS,z,
	   &ytemp2[j],&dummy);
  }
  // Then do final interpolation in x direction 
  imin=(ind>0)?(ind-1):0;
  imax=(ind<LORENTZ_X_POINTS-2)?(ind+2):(LORENTZ_X_POINTS-1);
  polint(&lorentz_x[imin],ytemp,imax-imin+1,r,&tanr,&dummy);
  polint(&lorentz_x[imin],ytemp2,imax-imin+1,r,&tanz,&dummy);

  // Compute correction
  double alpha=atan(S(3,0));  // Angle of incidence in x relative to z
  delta_y=tanr*delta_x*cos(alpha)-tanz*delta_x*sin(alpha)*cos(phi);

  // Check to see if distance to the wire from the track segment is consistent
  // with the cell size.
  if (fabs(delta_x)>=0.5){
    delta_y=0.;
    point->status|=NOT_CORRECTED;
  }
  point->s -=delta_y; 

  return NOERROR;
}




//------------------
// toString
//------------------
const string DFDCSegment_factory::toString(void)
{
        // Ensure our Get method has been called so _data is up to date
        Get();
        if(_data.size()<=0)return string(); // don't print anything if we have no data!

	return string();

}
