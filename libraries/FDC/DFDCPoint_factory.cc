//************************************************************************
// DFDCPoint_factory.cc - factory producing space points from pseudopoints
//************************************************************************

#include "DFDCPoint_factory.h"
#define HALF_CELL 0.5
#define MAX_DEFLECTION 0.15

bool DFDCPoint_package_cmp(const DFDCPoint* a, const DFDCPoint* b) {
  return a->wire->layer>b->wire->layer;
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
/// DFDCPoint_factory::evnt():
/// Routine where pseudopoints are converted into space points
///
jerror_t DFDCPoint_factory::evnt(JEventLoop* eventLoop, int eventNo) {
  vector<const DFDCPseudo*>pseudopoints;
  eventLoop->Get(pseudopoints);  

  vector<const DFDCPoint*>package[4];
  for (vector<const DFDCPseudo*>::iterator i=pseudopoints.begin();
       i!=pseudopoints.end();i++){
    DFDCPoint *newPoint = new DFDCPoint();
    
    float sinalpha=(*i)->wire->udir(0);
    float cosalpha=(*i)->wire->udir(1);
    float x=(*i)->w;
    float y=(*i)->s;
     
    newPoint->wire=(*i)->wire;
    newPoint->time=(*i)->time;
    newPoint->dist=(*i)->dist;
    newPoint->x=x*cosalpha+y*sinalpha;
    newPoint->y=-x*sinalpha+y*cosalpha; 
    
    package[((*i)->wire->layer-1)/6].push_back(newPoint);

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
jerror_t DFDCPoint_factory::KalmanFilter(vector<const DFDCPoint *>points){
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

  float z0=DFDCGeometry::GetZpackage(points[0]->wire->layer);
  float r=sqrt(points[0]->x*points[0]->x+points[0]->y*points[0]->y);
  float z=points[0]->wire->origin(2);
  float delta_z=z-z0;
  float phi=atan2(points[0]->y,points[0]->x);
  float gamma=(z-65.0)/sqrt(r*r+(z-65.0)*(z-65.0));
  float sintheta=r/sqrt(r*r+(z-65.0)*(z-65.0));
  float alpha=sin(phi)*sintheta;
  float beta=cos(phi)*sintheta;

  // Initial guess for state vector
  S(0,0)=points[0]->x-alpha*delta_z/gamma;
  S(1,0)=points[0]->y-beta*delta_z/gamma;
  S(2,0)=alpha;
  S(3,0)=beta;

  // Initialization 
  for (int i=0;i<4;i++){
    C(i,i)=1.;
    F(i,i)=1.;
  }
  H(0,0)=1;
  H(1,1)=1;
 
  for (unsigned int k=1;k<points.size();k++){ 
    // The measurement
    M(0,0)=points[k]->x;
    M(1,0)=points[k]->y;
 
    // Compute track projection matrix
    alpha=S(2,0);
    beta=S(3,0);
    gamma=sqrt(1-alpha*alpha-beta*beta);
    delta_z=points[k]->wire->origin(2)-z0;
    H(0,2)=H(1,3)=delta_z/gamma;

    // Compute measurement covariance matrix taking into account rotation
    float cosangle=points[k]->wire->udir(1);
    float sinangle=points[k]->wire->udir(0);
    float sigx2=HALF_CELL*HALF_CELL/12.;
    float sigy2=MAX_DEFLECTION*MAX_DEFLECTION/12.;
    V(0,0)=sigx2*cosangle*cosangle+sigy2*sinangle*sinangle;
    V(1,1)=sigx2*sinangle*sinangle+sigy2*cosangle*cosangle;
    V(0,1)=V(1,0)=(sigy2-sigx2)*sinangle*cosangle;

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
  for (unsigned int i=0;i<points.size();i++){
    CorrectPointY(z0,S,points[i]);
  }

  return NOERROR;
}

// Correct avalanche position along wire from doca data
jerror_t DFDCPoint_factory::CorrectPointY(float z0, DMatrix S,
				       const DFDCPoint *point){
  DVector3 P1,P2,P3,R;
  float z=point->wire->origin(2);
  float cosangle=point->wire->udir(1);
  float sinangle=point->wire->udir(0);
  float delta_z=z-z0;
  float delta_x=0;

  P1=point->wire->origin;
  P2.SetXYZ(point->x,point->y,z);
  P3.SetXYZ(S(0,0)+S(2,0)*delta_z/sqrt(1-S(2,0)*S(2,0)-S(3,0)*S(3,0)),
	    S(1,0)+S(3,0)*delta_z/sqrt(1-S(2,0)*S(2,0)-S(3,0)*S(3,0)),z);
  R=P1-P3+((P3-P1).Dot(P2-P1))/(P2-P1).Dot(P2-P1)*(P2-P1);

  delta_x=R(0)*cosangle-R(1)*sinangle;

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

        printheader("layer: wire: time(ns):     x(cm):     y(cm):   status:");

        for(unsigned int i=0; i<_data.size(); i++){
                DFDCPoint *hit = _data[i];
        
                printnewrow();
                printcol("%d",          hit->wire->layer);
                printcol("%d",          hit->wire->wire);
                printcol("%3.1f",       hit->time);
                printcol("%3.1f",       hit->x);
                printcol("%1.4f",       hit->y);
                printcol("%d",          hit->status);
                printrow();
        }

        return _table;

}
