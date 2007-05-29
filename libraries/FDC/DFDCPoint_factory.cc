//************************************************************************
// DFDCPoint_factory.cc - factory producing track segments from pseudopoints
//************************************************************************

#include "DFDCPoint_factory.h"
#define HALF_CELL 0.5
#define MAX_DEFLECTION 0.15
#define EPS 1e-3
#define KILL_RADIUS 1
#define NOT_CORRECTED 0x8

bool DFDCPoint_package_cmp(const DFDCPseudo* a, const DFDCPseudo* b) {
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



DFDCPoint_factory::DFDCPoint_factory() {
        logFile = new ofstream("DFDCPoint_factory.log");
        _log = new JStreamLog(*logFile, "POINT");
        *_log << "File initialized." << endMsg;
}


///
/// DFDCPoint_factory::~DFDCPoint_factory():
/// default destructor -- closes log file
///
DFDCPoint_factory::~DFDCPoint_factory() {
        logFile->close();
        delete logFile;
        delete _log;
}
///
/// DFDCPoint_factory::brun():
/// Initialization: read in deflection map
///
jerror_t DFDCPoint_factory::brun(JEventLoop* eventLoop, int eventNo) { 
  ifstream fdc_deflection_file;

  fdc_deflection_file.open("fdc_deflections.dat");
  if (!fdc_deflection_file.is_open()){ 
    *_log << "Failed to open fdc_deflection file." << endMsg;
    return RESOURCE_UNAVAILABLE;
  }

  float dummy;
  char buf[80];
  fdc_deflection_file.getline(buf,80,'\n');
  fdc_deflection_file.setf(ios::skipws); 
  for (int i=0;i<LORENTZ_X_POINTS;i++){
    for (int j=0;j<LORENTZ_Z_POINTS;j++){
      fdc_deflection_file >> lorentz_x[i];
      fdc_deflection_file >> lorentz_z[j];
      fdc_deflection_file >> dummy;
      fdc_deflection_file >> dummy;
      fdc_deflection_file >> lorentz_nx[i][j];
      fdc_deflection_file >> lorentz_nz[i][j];     
    }
  }
  fdc_deflection_file.close();

  *_log << "Table of Lorentz deflections initialized." << endMsg;
  return NOERROR;
}

///
/// DFDCPoint_factory::evnt():
/// Routine where pseudopoints are converted into space points
///
jerror_t DFDCPoint_factory::evnt(JEventLoop* eventLoop, int eventNo) {
  vector<const DFDCPseudo*>pseudopoints;
  eventLoop->Get(pseudopoints);  

  // Group pseudopoints by package
  vector<DFDCPseudo*>package[4];
  for (vector<const DFDCPseudo*>::iterator i=pseudopoints.begin();
       i!=pseudopoints.end();i++){
    package[((*i)->wire->layer-1)/6].push_back((DFDCPseudo*)(*i));

  }
  for (int j=0;j<4;j++){
    std::sort(package[j].begin(), package[j].end(), DFDCPoint_package_cmp);
    KalmanFilter(package[j]);
  
  }

  return NOERROR;
}

//
// Kalman Filter
//
jerror_t DFDCPoint_factory::KalmanFilter(vector<DFDCPseudo*>points){
  // Check that we have data to fit
  if (points.size()<2) return OBJECT_NOT_AVAILABLE;

  DMatrix S(4,1);  // State vector = {x0, y0, cos(theta_x), cos(theta_y)}
  DMatrix M(2,1);  // measurement vector
  DMatrix H(2,4);  // Track projection matrix
  DMatrix F(4,4);  // State vector propagation matrix
  DMatrix C(4,4);  // State vector covariance matrix
  DMatrix Q(4,4);  // Process noise covariance matrix
  DMatrix K(4,2);  // Kalman gain matrix
  DMatrix V(2,2);  // Measurement covariance matrix

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
    DFDCPoint *segment = new DFDCPoint;
    // Put first point in list of hits
    segment->hits.push_back(points[i]); 

    // Compute various geometrical information 
    // (including plane-to-plane rotations)
    float z0=DFDCGeometry::GetZpackage(points[i]->wire->layer);
    float sinangle=points[i]->wire->udir(0);
    float cosangle=points[i]->wire->udir(1);
    float x=(points[i]->w)*cosangle+(points[i]->s)*sinangle;
    float y=-(points[i]->w)*sinangle+(points[i]->s)*cosangle; 
    float r=sqrt(x*x+y*y);
    float z=points[i]->wire->origin(2);
    float delta_z=z-z0;
    float phi=atan2(y,x);
    float gamma=(z-65.0)/sqrt(r*r+(z-65.0)*(z-65.0));
    float sintheta=r/sqrt(r*r+(z-65.0)*(z-65.0));
    float alpha=cos(phi)*sintheta;
    float beta=sin(phi)*sintheta;
    
    // Initial guess for state vector
    S(0,0)=x-alpha*delta_z/gamma;
    S(1,0)=y-beta*delta_z/gamma;
    S(2,0)=alpha;
    S(3,0)=beta;
  
    // Initialization 
    for (int j=0;j<4;j++){
      C(j,j)=1.;
      F(j,j)=1.;
    }
    H(0,0)=1;
    H(1,1)=1;

    // Now loop over the remaining points in the other planes
    for (unsigned int k=start;k<points.size();k++){ 
      // The measurement
      sinangle=points[k]->wire->udir(0);
      cosangle=points[k]->wire->udir(1);
      M(0,0)=(points[k]->w)*cosangle+(points[k]->s)*sinangle;
      M(1,0)=-(points[k]->w)*sinangle+(points[k]->s)*cosangle; 
          
      // Compute track projection matrix
      alpha=S(2,0);
      beta=S(3,0);
      gamma=sqrt(1-alpha*alpha-beta*beta);
      delta_z=points[k]->wire->origin(2)-z0;
      H(0,2)=H(1,3)=delta_z/gamma;
      
      // Compute measurement covariance matrix taking into account rotation
      float sigx2=HALF_CELL*HALF_CELL/12.;
      float sigy2=MAX_DEFLECTION*MAX_DEFLECTION/12.;
      V(0,0)=sigx2*cosangle*cosangle+sigy2*sinangle*sinangle;
      V(1,1)=sigx2*sinangle*sinangle+sigy2*cosangle*cosangle;
      V(0,1)=V(1,0)=(sigy2-sigx2)*sinangle*cosangle;
      
      // Make sure the next hit is consistent with the track segment we are
      // working on...
      DMatrix M_pred=H*S;
      if (fabs(M_pred(0,0)-M(0,0))<KILL_RADIUS 
	  && fabs(M_pred(1,0)-M(1,0))<KILL_RADIUS){ 
	// Put the current point in the list of pseudopoints belonging to this
	// segment.
	segment->hits.push_back(points[k]);

	// Propagate the covariance matrix
	DMatrix F_T(4,4);
	F_T=DMatrix(DMatrix::kTransposed,F);
	C=F*(C*F_T);
	
	// Compute Kalman gain matrix
	DMatrix H_T(4,2);
	H_T = DMatrix(DMatrix::kTransposed,H);
	K=C*(H_T*DMatrix(DMatrix::kInverted,V+H*(C*H_T)));
	
	// Update the state vector 
	S=F*S+K*(M-H*(F*S));
	
	// Update state vector covariance matrix
	C=C-K*(H*C);
      }
    }

    // Correct for the Lorentz effect given DOCAs
    for (unsigned int k=0;k<segment->hits.size();k++){
      CorrectPointY(z0,S,segment->hits[k]);
    }
    
    segment->pos.SetXYZ(S(0,0),S(1,0),z0);
    segment->dir.SetXYZ(S(2,0),S(3,0),0);
    
    _data.push_back(segment);
  }
  return NOERROR;
}

// Correct avalanche position along wire from doca data
jerror_t DFDCPoint_factory::CorrectPointY(float z0, DMatrix S,
					  DFDCPseudo *point){
  DVector3 P1,P2,P3,R;
  float z=point->wire->origin(2);
  float cosangle=point->wire->udir(1);
  float sinangle=point->wire->udir(0);
  float x=(point->w)*cosangle+(point->s)*sinangle;
  float y=-(point->w)*sinangle+(point->s)*cosangle; 
  float delta_z=z-z0;
  float delta_x=0,delta_y=0;

  // Find the distance of closest approach between the point on the track
  // segment and the wire
  P1=point->wire->origin;
  P2.SetXYZ(x,y,z);
  P3.SetXYZ(S(0,0)+S(2,0)*delta_z/sqrt(1-S(2,0)*S(2,0)-S(3,0)*S(3,0)),
	    S(1,0)+S(3,0)*delta_z/sqrt(1-S(2,0)*S(2,0)-S(3,0)*S(3,0)),z);
  R=P1-P3+((P3-P1).Dot(P2-P1))/(P2-P1).Dot(P2-P1)*(P2-P1);
  delta_x=R(0)*cosangle-R(1)*sinangle;

  //
  // Next find correction to y from table of deflections
  //
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
  float alpha=asin(S(2,0));  // Angle of incidence in x relative to z
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
const string DFDCPoint_factory::toString(void)
{
        // Ensure our Get method has been called so _data is up to date
        Get();
        if(_data.size()<=0)return string(); // don't print anything if we have no data!

	return string();

}
