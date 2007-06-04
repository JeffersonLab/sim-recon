//************************************************************************
// DFDCSegment_factory.cc - factory producing track segments from pseudopoints
//************************************************************************

#include "DFDCSegment_factory.h"
#define HALF_CELL 0.5
#define MAX_DEFLECTION 0.15
#define EPS 1e-3
#define KILL_RADIUS 5.0 
#define NOT_CORRECTED 0x8
#define Z_TARGET 65.0

static bool got_deflection_file=true;
static bool warn=true;

bool DFDCSegment_package_cmp(const DFDCPseudo* a, const DFDCPseudo* b) {
  return a->wire->layer>b->wire->layer;
}


// Locate a position in array xx given x
void locate(float *xx,int n,float x,int *j){
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
void polint(float *xa, float *ya,int n,float x, float *y,float *dy){
  int i,m,ns=0;
  float den,dif,dift,ho,hp,w;
  
  float *c= new float[n];
  float *d= new float[n];

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

jerror_t FindCircle(float x, float y,vector<DFDCPseudo *>points,
		    unsigned int &start,
		    float &xc, float &yc, 
		    float &rc){
   // Find nearest neighbor in second-most downsteam plane with hits
  float x1,y1,xkeep=0,ykeep=0,diff,diff_min=1000.;
  for (unsigned int k=start;k<points.size();k++){
    float sinangle=points[k]->wire->udir(0);
    float cosangle=points[k]->wire->udir(1);
    x1=(points[k]->w)*cosangle+(points[k]->s)*sinangle;
    y1=-(points[k]->w)*sinangle+(points[k]->s)*cosangle;
    diff = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1));
    if (diff<diff_min){
      xkeep=x1;
      ykeep=y1;
      diff_min=diff;
      start=k;
    }
  }
  
  //Compute the center and radius of the circle through (x,y),(y1,y1) and (0,0)
  DMatrix A(3,3),D(3,3),B(3,3),F(3,3);
  A(0,0)=B(0,1)=F(0,1)=x;
  A(1,0)=B(1,1)=F(1,1)=xkeep;
  A(0,1)=D(0,1)=F(0,2)=y;
  A(1,1)=D(1,1)=F(1,2)=ykeep;
  D(0,0)=B(0,0)=F(0,0)=x*x+y*y; 
  D(1,0)=B(1,0)=F(1,0)=xkeep*xkeep+ykeep*ykeep;
  A(0,2)=A(1,2)=A(2,2)=B(0,2)=B(1,2)=B(2,2)=D(0,2)=D(1,2)=D(2,2)=1.0;
  
  float d,a,b,f;
  b=B.Determinant();
  d=-D.Determinant();
  a=A.Determinant();
  f=-F.Determinant(); 

  xc=-d/2/a;
  yc=-b/2/a;
  rc=sqrt((d*d+b*b)/4/a/a-f/a);

	return NOERROR;
 }




DFDCSegment_factory::DFDCSegment_factory() {
        logFile = new ofstream("DFDCSegment_factory.log");
        _log = new JStreamLog(*logFile, "POINT");
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

  float bx;
  float bz;
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
	/float(PACKAGE_Z_POINTS*LORENTZ_X_POINTS);  
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

  return NOERROR;
}

//
// Kalman Filter
//
jerror_t DFDCSegment_factory::KalmanFilter(vector<DFDCPseudo*>points){
  // Check that we have data to fit
  if (points.size()<2) return OBJECT_NOT_AVAILABLE;

  DMatrix S(5,1),Old_S(5,1);  // State vector 
  DMatrix M(2,1);  // measurement vector
  DMatrix H(2,5);  // Track projection matrix
  DMatrix F(5,5);  // State vector propagation matrix
  DMatrix C(5,5);  // State vector covariance matrix
  DMatrix Q(5,5);  // Process noise covariance matrix
  DMatrix K(5,2);  // Kalman gain matrix
  DMatrix V(2,2);  // Measurement covariance matrix
  DMatrix X(2,1);  // position on helical path

  // Each point in the most downstream plane of a package represents the 
  // end of track segment.  Put indices for these in the vector x_list.
  float old_z=points[0]->wire->origin(2);
  unsigned int start=1;
  vector<unsigned int>x_list;
  for (unsigned int i=0;i<points.size();i++){
    if (points[i]->wire->origin(2)!=old_z){
      start=i;
      break;
    }
    x_list.push_back(i);
    old_z=points[i]->wire->origin(2);
  }
  
  // Now loop over the list of track segment end points
  for (unsigned int i=0;i<x_list.size();i++){
    unsigned int new_start=start;
    
    DFDCSegment *segment = new DFDCSegment;
 
    // Put first point in list of hits
    segment->hits.push_back(points[i]);

    // Point in the last plane in the package 
    float sinangle=points[i]->wire->udir(0);
    float cosangle=points[i]->wire->udir(1);
    float x=(points[i]->w)*cosangle+(points[i]->s)*sinangle;
    float y=-(points[i]->w)*sinangle+(points[i]->s)*cosangle; 
    float z=points[i]->wire->origin(2);
    float Phi1=atan2(y,x),Phi2;

    // Find closest point in previous plane
    float delta_min=1000.,ybest=y,xbest=x;
    old_z=points[start]->wire->origin(2);
    for (unsigned int k=start;k<points.size();k++){     
      sinangle=points[k]->wire->udir(0);
      cosangle=points[k]->wire->udir(1);
      float x1=(points[k]->w)*cosangle+(points[k]->s)*sinangle;
      float y1=-(points[k]->w)*sinangle+(points[k]->s)*cosangle;
      float delta=sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1));
      if (points[k]->wire->origin(2)!=old_z){
	break;
      }
      if (delta<delta_min){
	delta_min=delta;
	xbest=x1;
	ybest=y1;
	new_start=k;
      }
      old_z=points[k]->wire->origin(2);
    }
    Phi2=atan2(ybest,xbest);

    // Guess particle charge (+/-1);
    float q=1.;      
    if (Phi1-Phi2<0) q=-1.;
    
    // Guess for the center and radius of the circle from the projection of 
    // the helical path
    float xc,yc,rc,rc2;
    FindCircle(x,y,points,new_start,xc,yc,rc);
   
    // Initialization of covariance matrix 
    for (int j=0;j<5;j++){
      C(j,j)=1.;
    }
 
    // Some useful tracking parameter variables
    float sperp,tanl;
    float z0=Z_TARGET;
    float kappa=q/2./rc;
    float phi=atan2(-xc,yc);
    float D=0.;
    sperp=asin(2.*kappa*(x*cos(phi)+y*sin(phi)))/2./kappa;

    // Initial guess for state vector
    S(0,0)=kappa;
    S(1,0)=phi;
    S(2,0)=D;
    S(3,0)=tanl=(z-z0)/sperp;
    S(4,0)=z0;
    Old_S=S; 

    // Now loop over the remaining points in the other planes
    for (unsigned int k=new_start;k<points.size();k++){
      // Old version of state vector before adding current measurement
      tanl=Old_S(3,0);      
      kappa=Old_S(0,0);
      phi=Old_S(1,0);
      D=Old_S(2,0);
      z0=Old_S(4,0);

      // Track position
      float xpos,ypos;
      GetHelicalTrackPosition(z,Old_S,xpos,ypos);
      X(0,0)=xpos;
      X(1,0)=ypos;
 
      // The next measurement 
      sinangle=points[k]->wire->udir(0);
      cosangle=points[k]->wire->udir(1);
      M(0,0)=(points[k]->w)*cosangle+(points[k]->s)*sinangle;
      M(1,0)=-(points[k]->w)*sinangle+(points[k]->s)*cosangle; 

      // ... and its covariance matrix taking into account rotation
      float sigx2=HALF_CELL*HALF_CELL/12.;
      float sigy2=MAX_DEFLECTION*MAX_DEFLECTION/12.;
      V(0,0)=sigx2*cosangle*cosangle+sigy2*sinangle*sinangle;
      V(1,1)=sigx2*sinangle*sinangle+sigy2*cosangle*cosangle;
      V(0,1)=V(1,0)=(sigy2-sigx2)*sinangle*cosangle;
     
      // Compute track projection matrix
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
      
      // Make sure the next hit is consistent with the track segment we are
      // working on...
      DMatrix M_pred=X+H*(S-Old_S);
      float delta_mx=M_pred(0,0)-M(0,0);
      float delta_my=M_pred(1,0)-M(1,0);
      if (fabs(delta_mx)<KILL_RADIUS && fabs(delta_my)<KILL_RADIUS){ 
	// Put the current point in the list of pseudopoints belonging to this
	// segment.
	segment->hits.push_back(points[k]);

        // Arc length
	z=points[k]->wire->origin(2);
	sperp=(z-z0)/tanl;	
	
        // Center of the circle
        xc=-(D+1./2./kappa)*sin(phi);
	yc=(D+1./2./kappa)*cos(phi);
	rc2=xc*xc+yc*yc;

	// Charge (+/-1)
	q=kappa/fabs(kappa);

	// Some temperorary variables for computing F
	float a=yc*sin(phi)+xc*cos(phi);
	float c=yc*cos(phi)-xc*sin(phi);

	// Transport matrix for S
	F(0,0)=F(3,3)=F(4,4)=1.;
	F(0,1)=-a/2./kappa/kappa/rc2;
	F(1,1)=	(D-1./2./kappa)*c/rc2;
	F(2,1)=a/rc2;
	F(0,2)=(1.-q*c/sqrt(rc2))/2./kappa/kappa;
	F(1,2)=-q*a/sqrt(rc)*(D+1./2./kappa);
	F(2,2)=q*c/sqrt(rc2);
	F(0,4)=tanl/2./kappa*((Old_S(1,0)-S(1,0))/kappa+F(1,0));
	F(1,4)=tanl/2./kappa*(F(1,1)-1.);
	F(2,4)=tanl/2./kappa*F(1,2);
	F(3,4)=(S(1,0)-Old_S(1,0))/kappa;

	// Compute Process noise covariance matrix
	GetProcessNoiseCovariance(S,segment->hits,Q);

	// Transport the covariance matrix
	DMatrix F_T(5,5);
	F_T=DMatrix(DMatrix::kTransposed,F);
	C=F*(C*F_T)+Q;

	// Compute Kalman gain matrix
	DMatrix H_T(5,2);
	H_T = DMatrix(DMatrix::kTransposed,H);
	K=C*(H_T*DMatrix(DMatrix::kInverted,V+H*(C*H_T)));

	// Update the state vector 
	Old_S=S; 
	S=Old_S+K*(M-M_pred); 

	// Update state vector covariance matrix
	C=C-K*(H*C);
      }
    }
  
    // Correct for the Lorentz effect given DOCAs
    if (got_deflection_file)
      for (unsigned int k=0;k<segment->hits.size();k++){
	CorrectPointY(z0,S,segment->hits[k]);
      }
    else if (warn){
      fprintf(stderr,
	      "DFDCSegment_factory: No Lorentz deflection file found!\n");
      fprintf(stderr,"DFDCSegment_factory: --> Pseudo-points will not be corrected for Lorentz effect.");
      warn=false;      
    }
    
    segment->S.ResizeTo(S);
    segment->S=S;
    segment->cov.ResizeTo(C);
    segment->cov=C;
    
    _data.push_back(segment);
  }
  return NOERROR;
}


jerror_t DFDCSegment_factory::GetProcessNoiseCovariance(DMatrix S,
			    vector<DFDCPseudo*>points,DMatrix &Q){
  if (!got_deflection_file) return RESOURCE_UNAVAILABLE;

  // Track parameters
  float tanl=S(3,0);
  float cosl=cos(atan(tanl));
  float kappa=S(0,0);
  float D=S(2,0);
  float phi=S(1,0);
  float sigma2_ms=0;  	

  // Momentum
  float p=0.003*BField[points[0]->wire->layer/PACKAGE_Z_POINTS]
    /2./fabs(kappa)/cosl;
  float beta=p/sqrt(p*p+0.140*0.140); // assume pion 

  //Materials: copper, Kapton, Mylar
  float thickness[3]={4e-4,50e-4,13e-4};
  float density[3]={8.96,1.42,1.39};
  float X0[3]={12.86,40.56,39.95};
  float material_sum=0.;
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

jerror_t DFDCSegment_factory::GetHelicalTrackPosition(float z,DMatrix S,
			  float  &xpos, float &ypos){
  // Track parameters
  float kappa=S(0,0);
  float phi=S(1,0);
  float D=S(2,0);
  float tanl=S(3,0);
  float z0=S(4,0);
  float sperp=(z-z0)/tanl;
  
  xpos=-D*sin(phi)+1./2./kappa*cos(phi)*sin(2.*kappa*sperp)
    -1./2./kappa*sin(phi)*(1.-cos(2.*kappa*sperp)); 
  ypos=D*cos(phi)+1./2./kappa*sin(phi)*sin(2.*kappa*sperp)
    +1./2./kappa*cos(phi)*(1.-cos(2.*kappa*sperp));
  return NOERROR;
}

// Correct avalanche position along wire from doca data
jerror_t DFDCSegment_factory::CorrectPointY(float z0, DMatrix S,
					  DFDCPseudo *point){
  
  DVector3 P1,P2,P3,R;
  float z=point->wire->origin(2);
  float cosangle=point->wire->udir(1);
  float sinangle=point->wire->udir(0);
  float x=(point->w)*cosangle+(point->s)*sinangle;
  float y=-(point->w)*sinangle+(point->s)*cosangle; 
  float delta_x=0,delta_y=0;
  float xhelix,yhelix;
 
  // Find the distance of closest approach between the point on the track
  // segment and the wire
  P1=point->wire->origin;
  P2.SetXYZ(x,y,z);
  GetHelicalTrackPosition(z,S,xhelix,yhelix);
  P3.SetXYZ(xhelix,yhelix,z);
  R=P1-P3+((P3-P1).Dot(P2-P1))/(P2-P1).Dot(P2-P1)*(P2-P1);
  delta_x=R(0)*cosangle-R(1)*sinangle;

  // Next find correction to y from table of deflections
  float tanr,tanz,r=sqrt(x*x+y*y);
  float phi=atan2(y,x);
  float ytemp[LORENTZ_X_POINTS],ytemp2[LORENTZ_X_POINTS],dummy;
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
  float alpha=atan(S(3,0));  // Angle of incidence in x relative to z
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
