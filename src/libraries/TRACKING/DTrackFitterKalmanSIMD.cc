//************************************************************************
// DTrackFitterKalmanSIMD.cc
//************************************************************************

#include "DTrackFitterKalmanSIMD.h"
#include "CDC/DCDCTrackHit.h"
#include "HDGEOMETRY/DLorentzDeflections.h"
#include "HDGEOMETRY/DMaterialMap.h"
#include "HDGEOMETRY/DRootGeom.h"
#include "DANA/DApplication.h"
#include <JANA/JCalibration.h>
#include "PID/DParticleID.h"

#include <TH2F.h>
#include <TH1I.h>
#include <TROOT.h>
#include <TMath.h>
#include <DMatrix.h>

#include <iomanip>
#include <math.h>

#define MAX_TB_PASSES 20
#define MAX_WB_PASSES 20
#define MAX_P 12.0
#define ALPHA 1./137.
#define CHISQ_DELTA 0.01


// Local boolean routines for sorting
//bool static DKalmanSIMDHit_cmp(DKalmanSIMDHit_t *a, DKalmanSIMDHit_t *b){
//  return a->z<b->z;
//}

inline bool static DKalmanSIMDFDCHit_cmp(DKalmanSIMDFDCHit_t *a, DKalmanSIMDFDCHit_t *b){
  if (fabs(a->z-b->z)<EPS){
    if (fabs(a->t-b->t)<EPS){
      double tsum_1=a->hit->t_u+a->hit->t_v;
      double tsum_2=b->hit->t_u+b->hit->t_v;
      if (fabs(tsum_1-tsum_2)<EPS){
	return (a->dE>b->dE);
      }
      return (tsum_1<tsum_2);
    }
    return(a->t<b->t);
  }
   return a->z<b->z;
}
inline bool static DKalmanSIMDCDCHit_cmp(DKalmanSIMDCDCHit_t *a, DKalmanSIMDCDCHit_t *b){
   if (a==NULL || b==NULL){
      cout << "Null pointer in CDC hit list??" << endl;
      return false;
   }
   const DCDCWire *wire_a= a->hit->wire;
   const DCDCWire *wire_b= b->hit->wire;
   if(wire_b->ring == wire_a->ring){
      return wire_b->straw < wire_a->straw;
   }

   return (wire_b->ring>wire_a->ring);
}


// Locate a position in array xx given x
void DTrackFitterKalmanSIMD::locate(const double *xx,int n,double x,int *j){
   int jl=-1;
   int ju=n;
   int ascnd=(xx[n-1]>=xx[0]);
   while(ju-jl>1){
      int jm=(ju+jl)>>1;
      if ( (x>=xx[jm])==ascnd)
         jl=jm;
      else
         ju=jm;
   }
   if (x==xx[0]) *j=0;
   else if (x==xx[n-1]) *j=n-2;
   else *j=jl; 
}



// Locate a position in vector xx given x
unsigned int DTrackFitterKalmanSIMD::locate(vector<double>&xx,double x){
   int n=xx.size();
   if (x==xx[0]) return 0;
   else if (x==xx[n-1]) return n-2;

   int jl=-1;
   int ju=n;
   int ascnd=(xx[n-1]>=xx[0]);
   while(ju-jl>1){
      int jm=(ju+jl)>>1;
      if ( (x>=xx[jm])==ascnd)
         jl=jm;
      else
         ju=jm;
   }
   return jl;
}

// Crude approximation for the variance in drift distance due to smearing
double DTrackFitterKalmanSIMD::fdc_drift_variance(double t){
   //return FDC_ANODE_VARIANCE;
   if (t<5.) t=5.;
   double sigma=DRIFT_RES_PARMS[0]/(t+1.)+DRIFT_RES_PARMS[1]+DRIFT_RES_PARMS[2]*t*t;

   return sigma*sigma;
}

// Convert time to distance for the cdc and compute variance
void DTrackFitterKalmanSIMD::ComputeCDCDrift(double dphi,double delta,double t,
      double B,
      double &d, double &V, double &tcorr){
   //d=0.39; // initialize at half-cell
   //V=0.0507; // initialize with (cell size)/12.
   tcorr = t; // Need this value even when t is negative for calibration
   if (t>0){
      //double dtc =(CDC_DRIFT_BSCALE_PAR1 + CDC_DRIFT_BSCALE_PAR2 * B)* t;
      //tcorr=t-dtc;

      //      CDC_RES_PAR2=0.005;
      double sigma=CDC_RES_PAR1/(tcorr+1.) + CDC_RES_PAR2 + CDC_RES_PAR3*tcorr;

      // Variables to store values for time-to-distance functions for delta=0
      // and delta!=0
      double f_0=0.;
      double f_delta=0.;
      // Derivative of d with respect to t, needed to add t0 variance 
      // dependence to sigma
      double dd_dt=0;
      // Scale factor to account for affect of B-field on maximum drift time
      double Bscale=long_drift_Bscale_par1+long_drift_Bscale_par2*B;
      tcorr=t*Bscale;

      //	if (delta>0)
      if (delta>-EPS2){
         double a1=long_drift_func[0][0];
         double a2=long_drift_func[0][1];
         double b1=long_drift_func[1][0];
         double b2=long_drift_func[1][1];
         double c1=long_drift_func[2][0];
         double c2=long_drift_func[2][1];
         double c3=long_drift_func[2][2];

         // use "long side" functional form
         double my_t=0.001*tcorr;
         double sqrt_t=sqrt(my_t);
         double t3=my_t*my_t*my_t;
         double delta_mag=fabs(delta);
         double a=a1+a2*delta_mag;
         double b=b1+b2*delta_mag;
         double c=c1+c2*delta_mag+c3*delta*delta;
         f_delta=a*sqrt_t+b*my_t+c*t3;
         f_0=a1*sqrt_t+b1*my_t+c1*t3;

         dd_dt=0.001*(0.5*a/sqrt_t+b+3.*c*my_t*my_t);
      }
      else{
         double my_t=0.001*tcorr;
         double sqrt_t=sqrt(my_t);
         double delta_mag=fabs(delta);

         // use "short side" functional form
         double a1=short_drift_func[0][0];
         double a2=short_drift_func[0][1];
         double a3=short_drift_func[0][2];
         double b1=short_drift_func[1][0];
         double b2=short_drift_func[1][1];
         double b3=short_drift_func[1][2];

         double delta_sq=delta*delta;
         double a=a1+a2*delta_mag+a3*delta_sq;
         double b=b1+b2*delta_mag+b3*delta_sq;
         f_delta=a*sqrt_t+b*my_t;
         f_0=a1*sqrt_t+b1*my_t;

         dd_dt=0.001*(0.5*a/sqrt_t+b);
      }

      unsigned int max_index=cdc_drift_table.size()-1;
      if (tcorr>cdc_drift_table[max_index]){
         //_DBG_ << "t: " << tcorr <<" d " << f_delta <<endl;
         d=f_delta;
         V=sigma*sigma+mVarT0*dd_dt*dd_dt;

         return;
      }

      // Drift time is within range of table -- interpolate...
      unsigned int index=0;
      index=locate(cdc_drift_table,tcorr);
      double dt=cdc_drift_table[index+1]-cdc_drift_table[index];
      double frac=(tcorr-cdc_drift_table[index])/dt;
      double d_0=0.01*(double(index)+frac); 

      if (fabs(delta) < EPS2){
         d=d_0;
         V=sigma*sigma+mVarT0*dd_dt*dd_dt;
      }
      else{
         double P=0.;
         double tcut=250.0; // ns
         if (tcorr<tcut) {
            P=(tcut-tcorr)/tcut;
         }
         d=f_delta*(d_0/f_0*P+1.-P);
         V=sigma*sigma+mVarT0*dd_dt*dd_dt;
      }
   }
   else { // Time is negative, or exactly zero, choose position at wire, with error of t=0 hit
      d=0.;
      double sigma = CDC_RES_PAR1+CDC_RES_PAR2;
      double dt=cdc_drift_table[1]-cdc_drift_table[0];
      V=sigma*sigma+mVarT0*0.0001/(dt*dt);
      //V=0.0507; // straw radius^2 / 12
   }

}

#define FDC_T0_OFFSET 17.6
// Interpolate on a table to convert time to distance for the fdc
/*
   double DTrackFitterKalmanSIMD::fdc_drift_distance(double t,double Bz){
   double a=93.31,b=4.614,Bref=2.143;
   t*=(a+b*Bref)/(a+b*Bz);
   int id=int((t+FDC_T0_OFFSET)/2.);
   if (id<0) id=0;
   if (id>138) id=138;
   double d=fdc_drift_table[id];  
   if (id!=138){
   double frac=0.5*(t+FDC_T0_OFFSET-2.*double(id));
   double dd=fdc_drift_table[id+1]-fdc_drift_table[id];
   d+=frac*dd;
   }

   return d;
   }
   */

// parametrization of time-to-distance for FDC
double DTrackFitterKalmanSIMD::fdc_drift_distance(double time,double Bz){
  if (time<0.) return 0.;
  double d=0.; 
  time/=1.+FDC_DRIFT_BSCALE_PAR1+FDC_DRIFT_BSCALE_PAR2*Bz*Bz;
  double tsq=time*time;
  double t_high=DRIFT_FUNC_PARMS[4];
  
  if (time<t_high){
    d=DRIFT_FUNC_PARMS[0]*sqrt(time)+DRIFT_FUNC_PARMS[1]*time
      +DRIFT_FUNC_PARMS[2]*tsq+DRIFT_FUNC_PARMS[3]*tsq*time;
  }
  else{
    double t_high_sq=t_high*t_high;
    d=DRIFT_FUNC_PARMS[0]*sqrt(t_high)+DRIFT_FUNC_PARMS[1]*t_high
      +DRIFT_FUNC_PARMS[2]*t_high_sq+DRIFT_FUNC_PARMS[3]*t_high_sq*t_high;
    d+=DRIFT_FUNC_PARMS[5]*(time-t_high);
  }
    
  return d;
}


DTrackFitterKalmanSIMD::DTrackFitterKalmanSIMD(JEventLoop *loop):DTrackFitter(loop){
   FactorForSenseOfRotation=(bfield->GetBz(0.,0.,65.)>0.)?-1.:1.;

   // Some useful values
   two_m_e=2.*ELECTRON_MASS;
   m_e_sq=ELECTRON_MASS*ELECTRON_MASS;

   // Get the position of the CDC downstream endplate from DGeometry
   double endplate_rmin,endplate_rmax;
   geom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);
   endplate_z-=0.5*endplate_dz;
   endplate_z_downstream=endplate_z+endplate_dz;
   endplate_rmin+=0.1;  // put just inside CDC
   endplate_r2min=endplate_rmin*endplate_rmin;
   endplate_r2max=endplate_rmax*endplate_rmax;

   // Beginning of the cdc
   vector<double>cdc_center;
   vector<double>cdc_upstream_endplate_pos; 
   vector<double>cdc_endplate_dim;
   geom->Get("//posXYZ[@volume='CentralDC'/@X_Y_Z",cdc_origin);
   geom->Get("//posXYZ[@volume='centralDC']/@X_Y_Z",cdc_center);
   geom->Get("//posXYZ[@volume='CDPU']/@X_Y_Z",cdc_upstream_endplate_pos);
   geom->Get("//tubs[@name='CDPU']/@Rio_Z",cdc_endplate_dim);
   cdc_origin[2]+=cdc_center[2]+cdc_upstream_endplate_pos[2]
      +0.5*cdc_endplate_dim[2];
   
   geom->GetCDCWires(cdcwires);
   //   geom->GetCDCRmid(cdc_rmid); // THIS ISN'T IMPLEMENTED!!
   // extract the "mean" radius of each ring from the wire data
   for(int ring=0; ring<cdcwires.size(); ring++)
  		cdc_rmid.push_back( cdcwires[ring][0]->origin.Perp() );
      
   // Outer detector geometry parameters
   geom->GetFCALZ(dFCALz); 
   if (geom->GetDIRCZ(dDIRCz)==false) dDIRCz=1000.;
   vector<double>tof_face;
   geom->Get("//section/composition/posXYZ[@volume='ForwardTOF']/@X_Y_Z",
	      tof_face);
   vector<double>tof_plane;  
   geom->Get("//composition[@name='ForwardTOF']/posXYZ[@volume='forwardTOF']/@X_Y_Z/plane[@value='0']", tof_plane);
   dTOFz=tof_face[2]+tof_plane[2]; 
   geom->Get("//composition[@name='ForwardTOF']/posXYZ[@volume='forwardTOF']/@X_Y_Z/plane[@value='1']", tof_plane);
   dTOFz+=tof_face[2]+tof_plane[2];
   dTOFz*=0.5;  // mid plane between tof planes

   
   
   // Get start counter geometry;
   if (geom->GetStartCounterGeom(sc_pos, sc_norm)){
     // Create vector of direction vectors in scintillator planes
     for (int i=0;i<30;i++){
       vector<DVector3>temp;
       for (unsigned int j=0;j<sc_pos[i].size()-1;j++){
	 double dx=sc_pos[i][j+1].x()-sc_pos[i][j].x();
	 double dy=sc_pos[i][j+1].y()-sc_pos[i][j].y();
	 double dz=sc_pos[i][j+1].z()-sc_pos[i][j].z();
	 temp.push_back(DVector3(dx/dz,dy/dz,1.));
       }
       sc_dir.push_back(temp);
     }
     SC_END_NOSE_Z=sc_pos[0][12].z();
     SC_BARREL_R2=sc_pos[0][0].Perp2();
     SC_PHI_SECTOR1=sc_pos[0][0].Phi();
   }

   // Get z positions of fdc wire planes
   geom->GetFDCZ(fdc_z_wires);
   // for now, assume the z extent of a package is the difference between the positions
   // of the two wire planes.  save half of this distance
   fdc_package_size = (fdc_z_wires[1]-fdc_z_wires[0]) / 2.;
   geom->GetFDCRmin(fdc_rmin_packages);
   geom->GetFDCRmax(fdc_rmax);

   ADD_VERTEX_POINT=false; 
   gPARMS->SetDefaultParameter("KALMAN:ADD_VERTEX_POINT", ADD_VERTEX_POINT);
  
   THETA_CUT=60.0; 
   gPARMS->SetDefaultParameter("KALMAN:THETA_CUT", THETA_CUT);

   RING_TO_SKIP=0;
   gPARMS->SetDefaultParameter("KALMAN:RING_TO_SKIP",RING_TO_SKIP);

   PLANE_TO_SKIP=0;
   gPARMS->SetDefaultParameter("KALMAN:PLANE_TO_SKIP",PLANE_TO_SKIP);

   MIN_HITS_FOR_REFIT=8; 
   gPARMS->SetDefaultParameter("KALMAN:MIN_HITS_FOR_REFIT", MIN_HITS_FOR_REFIT);
   PHOTON_ENERGY_CUTOFF=0.125; 
   gPARMS->SetDefaultParameter("KALMAN:PHOTON_ENERGY_CUTOFF",
         PHOTON_ENERGY_CUTOFF); 

   USE_FDC_HITS=true;
   gPARMS->SetDefaultParameter("TRKFIT:USE_FDC_HITS",USE_FDC_HITS);
   USE_CDC_HITS=true;
   gPARMS->SetDefaultParameter("TRKFIT:USE_CDC_HITS",USE_CDC_HITS);

   // Flag to enable calculation of alignment derivatives
   ALIGNMENT=false;
   gPARMS->SetDefaultParameter("TRKFIT:ALIGNMENT",ALIGNMENT);

   ALIGNMENT_FORWARD=false;
   gPARMS->SetDefaultParameter("TRKFIT:ALIGNMENT_FORWARD",ALIGNMENT_FORWARD);

   ALIGNMENT_CENTRAL=false;
   gPARMS->SetDefaultParameter("TRKFIT:ALIGNMENT_CENTRAL",ALIGNMENT_CENTRAL);

   if(ALIGNMENT){ALIGNMENT_FORWARD=true;ALIGNMENT_CENTRAL=true;}

   DEBUG_HISTS=false; 
   gPARMS->SetDefaultParameter("KALMAN:DEBUG_HISTS", DEBUG_HISTS);

   DEBUG_LEVEL=0;
   gPARMS->SetDefaultParameter("KALMAN:DEBUG_LEVEL", DEBUG_LEVEL); 

   USE_T0_FROM_WIRES=0;
   gPARMS->SetDefaultParameter("KALMAN:USE_T0_FROM_WIRES", USE_T0_FROM_WIRES);

   ESTIMATE_T0_TB=false;
   gPARMS->SetDefaultParameter("KALMAN:ESTIMATE_T0_TB",ESTIMATE_T0_TB);

   ENABLE_BOUNDARY_CHECK=true;
   gPARMS->SetDefaultParameter("GEOM:ENABLE_BOUNDARY_CHECK",
         ENABLE_BOUNDARY_CHECK);

   USE_MULS_COVARIANCE=true;
   gPARMS->SetDefaultParameter("TRKFIT:USE_MULS_COVARIANCE",
         USE_MULS_COVARIANCE);  

   USE_PASS1_TIME_MODE=false;
   gPARMS->SetDefaultParameter("KALMAN:USE_PASS1_TIME_MODE",USE_PASS1_TIME_MODE); 

   USE_FDC_DRIFT_TIMES=true;
   gPARMS->SetDefaultParameter("TRKFIT:USE_FDC_DRIFT_TIMES",
         USE_FDC_DRIFT_TIMES);

   RECOVER_BROKEN_TRACKS=true;
   gPARMS->SetDefaultParameter("KALMAN:RECOVER_BROKEN_TRACKS",RECOVER_BROKEN_TRACKS);

   NUM_CDC_SIGMA_CUT=3.5;
   NUM_FDC_SIGMA_CUT=3.5;
   gPARMS->SetDefaultParameter("KALMAN:NUM_CDC_SIGMA_CUT",NUM_CDC_SIGMA_CUT,
         "maximum distance in number of sigmas away from projection to accept cdc hit");
   gPARMS->SetDefaultParameter("KALMAN:NUM_FDC_SIGMA_CUT",NUM_FDC_SIGMA_CUT,
         "maximum distance in number of sigmas away from projection to accept fdc hit"); 

   ANNEAL_SCALE=1.5;
   ANNEAL_POW_CONST=20.0;
   gPARMS->SetDefaultParameter("KALMAN:ANNEAL_SCALE",ANNEAL_SCALE,
         "Scale factor for annealing");
   gPARMS->SetDefaultParameter("KALMAN:ANNEAL_POW_CONST",ANNEAL_POW_CONST,
         "Annealing parameter"); 
   FORWARD_ANNEAL_SCALE=1.5;
   FORWARD_ANNEAL_POW_CONST=20.0;
   gPARMS->SetDefaultParameter("KALMAN:FORWARD_ANNEAL_SCALE",
         FORWARD_ANNEAL_SCALE,
         "Scale factor for annealing");
   gPARMS->SetDefaultParameter("KALMAN:FORWARD_ANNEAL_POW_CONST",
         FORWARD_ANNEAL_POW_CONST,
         "Annealing parameter");

   FORWARD_PARMS_COV=false;
   gPARMS->SetDefaultParameter("KALMAN:FORWARD_PARMS_COV",FORWARD_PARMS_COV); 

   CDC_VAR_SCALE_FACTOR=1.;
   gPARMS->SetDefaultParameter("KALMAN:CDC_VAR_SCALE_FACTOR",CDC_VAR_SCALE_FACTOR); 
   CDC_T_DRIFT_MIN=-8.; // One f125 clock
   gPARMS->SetDefaultParameter("KALMAN:CDC_T_DRIFT_MIN",CDC_T_DRIFT_MIN);

   MOLIERE_FRACTION=0.9;
   gPARMS->SetDefaultParameter("KALMAN:MOLIERE_FRACTION",MOLIERE_FRACTION);    
   MS_SCALE_FACTOR=2.0;
   gPARMS->SetDefaultParameter("KALMAN:MS_SCALE_FACTOR",MS_SCALE_FACTOR);
   MOLIERE_RATIO1=0.5/(1.-MOLIERE_FRACTION);
   MOLIERE_RATIO2=MS_SCALE_FACTOR*1.e-6/(1.+MOLIERE_FRACTION*MOLIERE_FRACTION); //scale_factor/(1+F*F)

   COVARIANCE_SCALE_FACTOR_CENTRAL=20.0;
   gPARMS->SetDefaultParameter("KALMAN:COVARIANCE_SCALE_FACTOR_CENTRAL",
			       COVARIANCE_SCALE_FACTOR_CENTRAL);
 
   COVARIANCE_SCALE_FACTOR_FORWARD=2.0;
   gPARMS->SetDefaultParameter("KALMAN:COVARIANCE_SCALE_FACTOR_FORWARD",
                               COVARIANCE_SCALE_FACTOR_FORWARD);


   DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
   JCalibration *jcalib = dapp->GetJCalibration((loop->GetJEvent()).GetRunNumber());
   vector< map<string, double> > tvals;
   cdc_drift_table.clear();
   if (jcalib->Get("CDC/cdc_drift_table", tvals)==false){    
      for(unsigned int i=0; i<tvals.size(); i++){
         map<string, double> &row = tvals[i];
         cdc_drift_table.push_back(1000.*row["t"]);
      }
   }
   else{
      jerr << " CDC time-to-distance table not available... bailing..." << endl;
      exit(0);
   }

   int straw_number[28]={42,42,54,54,66,66,80,80,93,93,106,106,
      123,123,135,135,146,146,158,158,170,170,
      182,182,197,197,209,209};
   max_sag.clear();
   sag_phi_offset.clear();
   int straw_count=0,ring_count=0;
   if (jcalib->Get("CDC/sag_parameters", tvals)==false){
      vector<double>temp,temp2;
      for(unsigned int i=0; i<tvals.size(); i++){
         map<string, double> &row = tvals[i];

         temp.push_back(row["offset"]);
         temp2.push_back(row["phi"]);

         straw_count++;
         if (straw_count==straw_number[ring_count]){
            max_sag.push_back(temp);
            sag_phi_offset.push_back(temp2);
            temp.clear();
            temp2.clear();
            straw_count=0;
            ring_count++;
         }
      }
   }

   if (jcalib->Get("CDC/drift_parameters", tvals)==false){
      map<string, double> &row = tvals[0]; // long drift side
      long_drift_func[0][0]=row["a1"];
      long_drift_func[0][1]=row["a2"];
      long_drift_func[0][2]=row["a3"];  
      long_drift_func[1][0]=row["b1"];
      long_drift_func[1][1]=row["b2"];
      long_drift_func[1][2]=row["b3"];
      long_drift_func[2][0]=row["c1"];
      long_drift_func[2][1]=row["c2"];
      long_drift_func[2][2]=row["c3"];
      long_drift_Bscale_par1=row["B1"];
      long_drift_Bscale_par2=row["B2"];

      row = tvals[1]; // short drift side
      short_drift_func[0][0]=row["a1"];
      short_drift_func[0][1]=row["a2"];
      short_drift_func[0][2]=row["a3"];  
      short_drift_func[1][0]=row["b1"];
      short_drift_func[1][1]=row["b2"];
      short_drift_func[1][2]=row["b3"];
      short_drift_func[2][0]=row["c1"];
      short_drift_func[2][1]=row["c2"];
      short_drift_func[2][2]=row["c3"];
      short_drift_Bscale_par1=row["B1"];
      short_drift_Bscale_par2=row["B2"];
   }

   map<string, double> cdc_drift_parms;
   jcalib->Get("CDC/cdc_drift_parms", cdc_drift_parms);
   CDC_DRIFT_BSCALE_PAR1 = cdc_drift_parms["bscale_par1"];
   CDC_DRIFT_BSCALE_PAR2 = cdc_drift_parms["bscale_par2"];

   map<string, double> cdc_res_parms;
   jcalib->Get("CDC/cdc_resolution_parms_v2", cdc_res_parms);
   CDC_RES_PAR1 = cdc_res_parms["res_par1"];
   CDC_RES_PAR2 = cdc_res_parms["res_par2"];
   CDC_RES_PAR3 = cdc_res_parms["res_par3"];

   // Parameters for correcting for deflection due to Lorentz force
   map<string,double>lorentz_parms;
   jcalib->Get("FDC/lorentz_deflection_parms",lorentz_parms);
   LORENTZ_NR_PAR1=lorentz_parms["nr_par1"];
   LORENTZ_NR_PAR2=lorentz_parms["nr_par2"];
   LORENTZ_NZ_PAR1=lorentz_parms["nz_par1"];
   LORENTZ_NZ_PAR2=lorentz_parms["nz_par2"];

   // Parameters for accounting for variation in drift distance from FDC
   map<string,double>drift_res_parms;
   jcalib->Get("FDC/drift_resolution_parms",drift_res_parms); 
   DRIFT_RES_PARMS[0]=drift_res_parms["p0"];   
   DRIFT_RES_PARMS[1]=drift_res_parms["p1"];
   DRIFT_RES_PARMS[2]=drift_res_parms["p2"]; 

   // Time-to-distance function parameters for FDC
   map<string,double>drift_func_parms;
   jcalib->Get("FDC/drift_function_parms",drift_func_parms); 
   DRIFT_FUNC_PARMS[0]=drift_func_parms["p0"];   
   DRIFT_FUNC_PARMS[1]=drift_func_parms["p1"];
   DRIFT_FUNC_PARMS[2]=drift_func_parms["p2"]; 
   DRIFT_FUNC_PARMS[3]=drift_func_parms["p3"];
   DRIFT_FUNC_PARMS[4]=1000.;
   DRIFT_FUNC_PARMS[5]=0.;
   map<string,double>drift_func_ext;
   if (jcalib->Get("FDC/drift_function_ext",drift_func_ext)==false){
     DRIFT_FUNC_PARMS[4]=drift_func_ext["p4"]; 
     DRIFT_FUNC_PARMS[5]=drift_func_ext["p5"]; 
   }
   // Factors for taking care of B-dependence of drift time for FDC
   map<string, double> fdc_drift_parms;
   jcalib->Get("FDC/fdc_drift_parms", fdc_drift_parms);
   FDC_DRIFT_BSCALE_PAR1 = fdc_drift_parms["bscale_par1"];
   FDC_DRIFT_BSCALE_PAR2 = fdc_drift_parms["bscale_par2"];


   /*
      if (jcalib->Get("FDC/fdc_drift2", tvals)==false){
      for(unsigned int i=0; i<tvals.size(); i++){
      map<string, float> &row = tvals[i];
      iter_float iter = row.begin();
      fdc_drift_table[i] = iter->second;
      }
      }
      else{
      jerr << " FDC time-to-distance table not available... bailing..." << endl;
      exit(0);
      }
      */

   for (unsigned int i=0;i<5;i++)I5x5(i,i)=1.;


   // center of the target
   map<string, double> targetparms;
   if (jcalib->Get("TARGET/target_parms",targetparms)==false){
      TARGET_Z = targetparms["TARGET_Z_POSITION"];
   }
   else{
      geom->GetTargetZ(TARGET_Z);
   }
   if (ADD_VERTEX_POINT){
     gPARMS->SetDefaultParameter("KALMAN:VERTEX_POSITION",TARGET_Z);
   }

   // Beam position and direction
   map<string, double> beam_vals;
   jcalib->Get("PHOTON_BEAM/beam_spot",beam_vals);
   beam_center.Set(beam_vals["x"],beam_vals["y"]); 
   beam_dir.Set(beam_vals["dxdz"],beam_vals["dydz"]);
   jout << " Beam spot: x=" << beam_center.X() << " y=" << beam_center.Y()
	<< " z=" << beam_vals["z"]
	<< " dx/dz=" << beam_dir.X() << " dy/dz=" << beam_dir.Y() << endl;
   beam_z0=beam_vals["z"];
   
   // Inform user of some configuration settings
   static bool config_printed = false;
   if(!config_printed){
      config_printed = true;
      stringstream ss;
      ss << "vertex constraint: " ;
      if(ADD_VERTEX_POINT){
         ss << "z = " << TARGET_Z << "cm" << endl;
      }else{
         ss << "<none>" << endl;
      }
      jout << ss.str();
   } // config_printed

  if(DEBUG_HISTS){ 
   for (auto i=0; i < 46; i++){
      double min = -10., max=10.;
      if(i%23<12) {min=-5; max=5;}
      if(i<23)alignDerivHists[i]=new TH1I(Form("CentralDeriv%i",i),Form("CentralDeriv%i",i),200, min, max);
      else alignDerivHists[i]=new TH1I(Form("ForwardDeriv%i",i%23),Form("ForwardDeriv%i",i%23),200, min, max);
   }
   brentCheckHists[0]=new TH2I("ForwardBrentCheck","DOCA vs ds", 100, -5., 5., 100, 0.0, 1.5);
   brentCheckHists[1]=new TH2I("CentralBrentCheck","DOCA vs ds", 100, -5., 5., 100, 0.0, 1.5);
  }
   
	dResourcePool_TMatrixFSym = std::make_shared<DResourcePool<TMatrixFSym>>();
	dResourcePool_TMatrixFSym->Set_ControlParams(20, 20, 50);

	my_fdchits.reserve(24);
	my_cdchits.reserve(28);
	fdc_updates.reserve(24);
	cdc_updates.reserve(28);
	cdc_used_in_fit.reserve(28);
	fdc_used_in_fit.reserve(24);
}

//-----------------
// ResetKalmanSIMD
//-----------------
void DTrackFitterKalmanSIMD::ResetKalmanSIMD(void)
{
   last_material_map=0;

   for (unsigned int i=0;i<my_cdchits.size();i++){
      delete my_cdchits[i];
   } 
   for (unsigned int i=0;i<my_fdchits.size();i++){
      delete my_fdchits[i];
   }
   central_traj.clear();
   forward_traj.clear();
   my_fdchits.clear();
   my_cdchits.clear();
   fdc_updates.clear();
   cdc_updates.clear();
   fdc_used_in_fit.clear();
   cdc_used_in_fit.clear();

   cov.clear();
   fcov.clear();

   len = 0.0;
   ftime=0.0;
   var_ftime=0.0;
   x_=0.,y_=0.,tx_=0.,ty_=0.,q_over_p_ = 0.0;
   z_=0.,phi_=0.,tanl_=0.,q_over_pt_ =0, D_= 0.0;
   chisq_ = 0.0;
   ndf_ = 0;
   MASS=0.13957;
   mass2=MASS*MASS;
   Bx=By=0.;
   Bz=2.;
   dBxdx=0.,dBxdy=0.,dBxdz=0.,dBydx=0.,dBydy=0.,dBydy=0.,dBzdx=0.,dBzdy=0.,dBzdz=0.;
   // Step sizes
   mStepSizeS=2.0;
   mStepSizeZ=2.0;
   //mStepSizeZ=0.5;

   //if (fit_type==kTimeBased){
   //   mStepSizeS=0.5;
   //   mStepSizeZ=0.5;
   // }


   mT0=0.,mT0MinimumDriftTime=1e6;
   mVarT0=25.;

   mCDCInternalStepSize=0.5;
   //mCDCInternalStepSize=1.0;
   //mCentralStepSize=0.75;
   mCentralStepSize=0.75;

   mT0Detector=SYS_CDC;

   IsHadron=true;
   IsElectron=false;
   IsPositron=false;
   
   PT_MIN=0.01;
   Q_OVER_P_MAX=100.;

   // These variables provide the approximate location along the trajectory
   // where there is an indication of a kink in the track      
   break_point_fdc_index=0;
   break_point_cdc_index=0;
   break_point_step_index=0;

}

//-----------------
// FitTrack
//-----------------
DTrackFitter::fit_status_t DTrackFitterKalmanSIMD::FitTrack(void)
{
   // Reset member data and free an memory associated with the last fit,
   // but some of which only for wire-based fits 
   ResetKalmanSIMD();

   // Check that we have enough FDC and CDC hits to proceed
   if (cdchits.size()+fdchits.size()<6) return kFitNotDone;

   // Copy hits from base class into structures specific to DTrackFitterKalmanSIMD  
   if (USE_CDC_HITS && cdchits.size()>=MIN_CDC_HITS) 
      for(unsigned int i=0; i<cdchits.size(); i++)AddCDCHit(cdchits[i]);
   if (USE_FDC_HITS && fdchits.size()>=MIN_FDC_HITS)
      for(unsigned int i=0; i<fdchits.size(); i++)AddFDCHit(fdchits[i]);

   unsigned int num_good_cdchits=my_cdchits.size();
   unsigned int num_good_fdchits=my_fdchits.size(); 

   // keep track of the range of detector elements that could be hit
   // for calculating the number of expected hits later on
   //int min_cdc_ring=-1, max_cdc_ring=-1;

   // Order the cdc hits by ring number
   if (num_good_cdchits>0){
      stable_sort(my_cdchits.begin(),my_cdchits.end(),DKalmanSIMDCDCHit_cmp);

	  //min_cdc_ring = my_cdchits[0]->hit->wire->ring;
	  //max_cdc_ring = my_cdchits[my_cdchits.size()-1]->hit->wire->ring;

      // Look for multiple hits on the same wire
      for (unsigned int i=0;i<my_cdchits.size()-1;i++){
         if (my_cdchits[i]->hit->wire->ring==my_cdchits[i+1]->hit->wire->ring && 
               my_cdchits[i]->hit->wire->straw==my_cdchits[i+1]->hit->wire->straw){
            num_good_cdchits--;	
            if (my_cdchits[i]->tdrift<my_cdchits[i+1]->tdrift){
               my_cdchits[i+1]->status=late_hit;
            }
            else{
               my_cdchits[i]->status=late_hit;
            }
         }
      }

   }
   // Order the fdc hits by z
   if (num_good_fdchits>0){
      stable_sort(my_fdchits.begin(),my_fdchits.end(),DKalmanSIMDFDCHit_cmp);

      // Look for multiple hits on the same wire 
      for (unsigned int i=0;i<my_fdchits.size()-1;i++){
         if (my_fdchits[i]->hit->wire->layer==my_fdchits[i+1]->hit->wire->layer &&
               my_fdchits[i]->hit->wire->wire==my_fdchits[i+1]->hit->wire->wire){
            num_good_fdchits--;
	    if (fabs(my_fdchits[i]->t-my_fdchits[i+1]->t)<EPS){
	      double tsum_1=my_fdchits[i]->hit->t_u+my_fdchits[i]->hit->t_v;
	      double tsum_2=my_fdchits[i+1]->hit->t_u+my_fdchits[i+1]->hit->t_v;
	      if (tsum_1<tsum_2){
		my_fdchits[i+1]->status=late_hit;
	      }
	      else{
		my_fdchits[i]->status=late_hit;
	      }
	    }
            else if (my_fdchits[i]->t<my_fdchits[i+1]->t){
               my_fdchits[i+1]->status=late_hit;
            }
            else{
               my_fdchits[i]->status=late_hit;
            }
         }
      }
   }
   if(num_good_cdchits+num_good_fdchits<6) return kFitNotDone;

   // Create vectors of updates (from hits) to S and C
   if (my_cdchits.size()>0){
      cdc_updates=vector<DKalmanUpdate_t>(my_cdchits.size());
      // Initialize vector to keep track of whether or not a hit is used in 
      // the fit
      cdc_used_in_fit=vector<bool>(my_cdchits.size());
   }
   if (my_fdchits.size()>0){
      fdc_updates=vector<DKalmanUpdate_t>(my_fdchits.size());
      // Initialize vector to keep track of whether or not a hit is used in 
      // the fit
      fdc_used_in_fit=vector<bool>(my_fdchits.size());
   }

   // start time and variance
   if (fit_type==kTimeBased){
      mT0=input_params.t0();
      switch(input_params.t0_detector()){
      case SYS_TOF:
	mVarT0=0.01;
	break;
      case SYS_CDC:
	mVarT0=7.5;
	break;
      case SYS_FDC:
	mVarT0=7.5;
	break;
      case SYS_BCAL:
	mVarT0=0.25;
	break;
      default:
	mVarT0=0.09;
	break;
      }
   }
   
   //_DBG_ << SystemName(input_params.t0_detector()) << " " << mT0 <<endl;

   //Set the mass
   MASS=input_params.mass();
   mass2=MASS*MASS;
   m_ratio=ELECTRON_MASS/MASS;
   m_ratio_sq=m_ratio*m_ratio;

   // Is this particle an electron or positron?
   if (MASS<0.001){
      IsHadron=false;
      if (input_params.charge()<0.) IsElectron=true;
      else IsPositron=true;
   }
   if (DEBUG_LEVEL>0)
     {
       _DBG_ << "------Starting " 
         <<(fit_type==kTimeBased?"Time-based":"Wire-based") 
         << " Fit with " << my_fdchits.size() << " FDC hits and " 
         << my_cdchits.size() << " CDC hits.-------" <<endl;
      if (fit_type==kTimeBased){
	_DBG_ << " Using t0=" << mT0 << " from DET=" 
            << input_params.t0_detector() <<endl;
      }
   }
   // Do the fit
   jerror_t error = KalmanLoop();
   if (error!=NOERROR){
      if (DEBUG_LEVEL>0)
         _DBG_ << "Fit failed with Error = " << error <<endl;
      return kFitFailed;
   }

   // Copy fit results into DTrackFitter base-class data members
   DVector3 mom,pos;
   GetPosition(pos);
   GetMomentum(mom);
   double charge = GetCharge();
   fit_params.setPosition(pos);
   fit_params.setMomentum(mom);
   fit_params.setTime(mT0MinimumDriftTime);
   fit_params.setPID(IDTrack(charge, MASS));
   fit_params.setT0(mT0MinimumDriftTime,4.,mT0Detector);

   if (DEBUG_LEVEL>0){
      _DBG_ << "----- Pass: " 
         << (fit_type==kTimeBased?"Time-based ---":"Wire-based ---") 
         << " Mass: " << MASS 
         << " p=" 	<< mom.Mag()
         << " theta="  << 90.0-180./M_PI*atan(tanl_)
         << " vertex=(" << x_ << "," << y_ << "," << z_<<")"
         << " chi2=" << chisq_
         <<endl;
      if(DEBUG_LEVEL>1){
         //Dump pulls
         for (unsigned int iPull = 0; iPull < pulls.size(); iPull++){
            if (pulls[iPull].cdc_hit != NULL){
	      _DBG_ << " ring: " <<  pulls[iPull].cdc_hit->wire->ring
		    << " straw: " << pulls[iPull].cdc_hit->wire->straw  
		    << " Residual: " << pulls[iPull].resi
		    << " Err: " << pulls[iPull].err
		    << " tdrift: " << pulls[iPull].tdrift
		    << " doca: " << pulls[iPull].d
		    << " docaphi: " << pulls[iPull].docaphi
		    << " z: " << pulls[iPull].z
		    << " cos(theta_rel): " << pulls[iPull].cosThetaRel
		    << " tcorr: " << pulls[iPull].tcorr 
		    << endl;
            }
         }
      }
   }

   DMatrixDSym errMatrix(5);
   // Fill the tracking error matrix and the one needed for kinematic fitting
   if (fcov.size()!=0){      
      // We MUST fill the entire matrix (not just upper right) even though 
      // this is a DMatrixDSym
      for (unsigned int i=0;i<5;i++){
         for (unsigned int j=0;j<5;j++){
            errMatrix(i,j)=fcov[i][j];
         }
      }
      if (FORWARD_PARMS_COV){
         fit_params.setForwardParmFlag(true);    
         fit_params.setTrackingStateVector(x_,y_,tx_,ty_,q_over_p_);

         // Compute and fill the error matrix needed for kinematic fitting
         fit_params.setErrorMatrix(Get7x7ErrorMatrixForward(errMatrix));
      }
      else {
         fit_params.setForwardParmFlag(false); 
         fit_params.setTrackingStateVector(q_over_pt_,phi_,tanl_,D_,z_);

         // Compute and fill the error matrix needed for kinematic fitting
         fit_params.setErrorMatrix(Get7x7ErrorMatrix(errMatrix));
      }
   }
   else if (cov.size()!=0){
      fit_params.setForwardParmFlag(false);

      // We MUST fill the entire matrix (not just upper right) even though 
      // this is a DMatrixDSym
      for (unsigned int i=0;i<5;i++){
         for (unsigned int j=0;j<5;j++){
            errMatrix(i,j)=cov[i][j];
         }
      }
      fit_params.setTrackingStateVector(q_over_pt_,phi_,tanl_,D_,z_);

      // Compute and fill the error matrix needed for kinematic fitting
      fit_params.setErrorMatrix(Get7x7ErrorMatrix(errMatrix));
   }
   auto locTrackingCovarianceMatrix = dResourcePool_TMatrixFSym->Get_SharedResource();
   locTrackingCovarianceMatrix->ResizeTo(5, 5);
   for(unsigned int loc_i = 0; loc_i < 5; ++loc_i)
   {
      for(unsigned int loc_j = 0; loc_j < 5; ++loc_j)
         (*locTrackingCovarianceMatrix)(loc_i, loc_j) = errMatrix(loc_i, loc_j);

   }
   fit_params.setTrackingErrorMatrix(locTrackingCovarianceMatrix);
   this->chisq = GetChiSq();
   this->Ndof = GetNDF();
   fit_status = kFitSuccess;

   // figure out the number of expected hits for this track based on the final fit
	set<const DCDCWire *> expected_hit_straws;
	set<int> expected_hit_fdc_planes;

	for(int i=0; i<extrapolations[SYS_CDC].size(); i++) {
		// figure out the radial position of the point to see which ring it's in
		double r = extrapolations[SYS_CDC][i].position.Perp();
		int ring=0;
		for(; ring<cdc_rmid.size(); ring++) {
			//_DBG_ << "Rs = " << r << " " << cdc_rmid[ring] << endl;
			if( (r<cdc_rmid[ring]-0.78) || (fabs(r-cdc_rmid[ring])<0.78) )
				break;
		}
		if(ring == cdc_rmid.size()) ring--;
		//_DBG_ << "ring = " << ring << endl;
		//_DBG_ << "ring = " << ring << "  stereo = " << cdcwires[ring][0]->stereo << endl;
		int best_straw=0;
		double best_dist_diff=fabs((extrapolations[SYS_CDC][i].position 
			- cdcwires[ring][0]->origin).Mag());		
	    // match based on straw center
	    for(int straw=1; straw<cdcwires[ring].size(); straw++) {
	    	DVector3 wire_position = cdcwires[ring][straw]->origin;  // start with the nominal wire center
	    	// now take into account the z dependence due to the stereo angle
	    	double dz = extrapolations[SYS_CDC][i].position.Z() - cdcwires[ring][straw]->origin.Z();
	    	double ds = dz*tan(cdcwires[ring][straw]->stereo);
	    	wire_position += DVector3(-ds*sin(cdcwires[ring][straw]->origin.Phi()), ds*cos(cdcwires[ring][straw]->origin.Phi()), dz);
	    	double diff = fabs((extrapolations[SYS_CDC][i].position
				- wire_position).Mag());
			if( diff < best_dist_diff )
				best_straw = straw;
	    }
	    
	    expected_hit_straws.insert(cdcwires[ring][best_straw]);
	}
	
	for(int i=0; i<extrapolations[SYS_FDC].size(); i++) {
		// check to make sure that the track goes through the sensitive region of the FDC
		// assume one hit per plane
		double z = extrapolations[SYS_FDC][i].position.Z();
		double r = extrapolations[SYS_FDC][i].position.Perp();

		// see if we're in the "sensitive area" of a package
		for(int plane=0; plane<fdc_z_wires.size(); plane++) {
			int package = plane/6;
			if(fabs(z-fdc_z_wires[plane]) < fdc_package_size) {
				if( r<fdc_rmax && r>fdc_rmin_packages[package]) {
					expected_hit_fdc_planes.insert(plane);
				}
				break; // found the right plane
			}
 		}
	}
	
	potential_cdc_hits_on_track = expected_hit_straws.size();
	potential_fdc_hits_on_track = expected_hit_fdc_planes.size();

    if(DEBUG_LEVEL>0) {
   		_DBG_ << " CDC hits/potential hits " << my_cdchits.size() << "/" << potential_cdc_hits_on_track 
        	 << "  FDC hits/potential hits " << my_fdchits.size() << "/" << potential_fdc_hits_on_track  << endl;
	}
	
   //_DBG_  << "========= done!" << endl;

   return fit_status;
}

//-----------------
// ChiSq
//-----------------
double DTrackFitterKalmanSIMD::ChiSq(fit_type_t fit_type, DReferenceTrajectory *rt, double *chisq_ptr, int *dof_ptr, vector<pull_t> *pulls_ptr)
{
   // This simply returns whatever was left in for the chisq/NDF from the last fit.
   // Using a DReferenceTrajectory is not really appropriate here so the base class'
   // requirement of it should be reviewed.
   double chisq = GetChiSq();
   unsigned int ndf = GetNDF();

   if(chisq_ptr)*chisq_ptr = chisq;
   if(dof_ptr)*dof_ptr = int(ndf);
   if(pulls_ptr)*pulls_ptr = pulls;

   return chisq/double(ndf);
}

// Return the momentum at the distance of closest approach to the origin.
inline void DTrackFitterKalmanSIMD::GetMomentum(DVector3 &mom){
   double pt=1./fabs(q_over_pt_);
   mom.SetXYZ(pt*cos(phi_),pt*sin(phi_),pt*tanl_);
}

// Return the "vertex" position (position at which track crosses beam line)
inline void DTrackFitterKalmanSIMD::GetPosition(DVector3 &pos){
   pos.SetXYZ(x_,y_,z_);
}

// Add FDC hits
jerror_t DTrackFitterKalmanSIMD::AddFDCHit(const DFDCPseudo *fdchit){
   DKalmanSIMDFDCHit_t *hit= new DKalmanSIMDFDCHit_t;

   hit->t=fdchit->time;
   hit->uwire=fdchit->w;
   hit->vstrip=fdchit->s;
   hit->vvar=fdchit->ds*fdchit->ds;
   hit->z=fdchit->wire->origin.z();
   hit->cosa=cos(fdchit->wire->angle);
   hit->sina=sin(fdchit->wire->angle);
   hit->phiX=fdchit->wire->angles.X();
   hit->phiY=fdchit->wire->angles.Y();
   hit->phiZ=fdchit->wire->angles.Z();

   hit->nr=0.;
   hit->nz=0.;
   hit->dE=1e6*fdchit->dE;
   hit->hit=fdchit;
   hit->status=good_hit;

   my_fdchits.push_back(hit);

   return NOERROR;
}

//  Add CDC hits
jerror_t DTrackFitterKalmanSIMD::AddCDCHit (const DCDCTrackHit *cdchit){
   DKalmanSIMDCDCHit_t *hit= new DKalmanSIMDCDCHit_t;

   hit->hit=cdchit;
   hit->status=good_hit;
   hit->origin.Set(cdchit->wire->origin.x(),cdchit->wire->origin.y());
   double one_over_uz=1./cdchit->wire->udir.z();
   hit->dir.Set(one_over_uz*cdchit->wire->udir.x(),
         one_over_uz*cdchit->wire->udir.y());
   hit->z0wire=cdchit->wire->origin.z();
   hit->cosstereo=cos(cdchit->wire->stereo);
   hit->tdrift=cdchit->tdrift;
   my_cdchits.push_back(hit);

   return NOERROR;
}

// Calculate the derivative of the state vector with respect to z
jerror_t DTrackFitterKalmanSIMD::CalcDeriv(double z,
      const DMatrix5x1 &S, 
      double dEdx, 
      DMatrix5x1 &D){
   double tx=S(state_tx),ty=S(state_ty);
   double q_over_p=S(state_q_over_p);

   // Turn off dEdx if the magnitude of the momentum drops below some cutoff
   if (fabs(q_over_p)>Q_OVER_P_MAX){
      dEdx=0.;
   }
   // Try to keep the direction tangents from heading towards 90 degrees
   if (fabs(tx)>TAN_MAX) tx=TAN_MAX*(tx>0.0?1.:-1.); 
   if (fabs(ty)>TAN_MAX) ty=TAN_MAX*(ty>0.0?1.:-1.);

   // useful combinations of terms
   double kq_over_p=qBr2p*q_over_p;
   double tx2=tx*tx;
   double ty2=ty*ty;
   double txty=tx*ty;
   double one_plus_tx2=1.+tx2;
   double dsdz=sqrt(one_plus_tx2+ty2);
   double dtx_Bfac=ty*Bz+txty*Bx-one_plus_tx2*By;
   double dty_Bfac=Bx*(1.+ty2)-txty*By-tx*Bz;
   double kq_over_p_dsdz=kq_over_p*dsdz;

   // Derivative of S with respect to z
   D(state_x)=tx;
   D(state_y)=ty;
   D(state_tx)=kq_over_p_dsdz*dtx_Bfac;
   D(state_ty)=kq_over_p_dsdz*dty_Bfac;

   D(state_q_over_p)=0.;
   if (CORRECT_FOR_ELOSS && fabs(dEdx)>EPS){
      double q_over_p_sq=q_over_p*q_over_p;
      double E=sqrt(1./q_over_p_sq+mass2); 
      D(state_q_over_p)=-q_over_p_sq*q_over_p*E*dEdx*dsdz;
   }
   return NOERROR;
}

// Reference trajectory for forward tracks in CDC region
// At each point we store the state vector and the Jacobian needed to get to 
//this state along z from the previous state.
jerror_t DTrackFitterKalmanSIMD::SetCDCForwardReferenceTrajectory(DMatrix5x1 &S){
   int i=0,forward_traj_length=forward_traj.size();
   double z=z_;
   double r2=0.;
   bool stepped_to_boundary=false;

   // Magnetic field and gradient at beginning of trajectory
   //bfield->GetField(x_,y_,z_,Bx,By,Bz); 
   bfield->GetFieldAndGradient(x_,y_,z_,Bx,By,Bz,
         dBxdx,dBxdy,dBxdz,dBydx,
         dBydy,dBydz,dBzdx,dBzdy,dBzdz);

   // Reset cdc status flags
   for (unsigned int i=0;i<my_cdchits.size();i++){
      if (my_cdchits[i]->status!=late_hit) my_cdchits[i]->status=good_hit;
   }

   // Continue adding to the trajectory until we have reached the endplate
   // or the maximum radius
   while(z<endplate_z && z>cdc_origin[2] &&
         r2<endplate_r2max && fabs(S(state_q_over_p))<Q_OVER_P_MAX
         && fabs(S(state_tx))<TAN_MAX
         && fabs(S(state_ty))<TAN_MAX
        ){
      if (PropagateForwardCDC(forward_traj_length,i,z,r2,S,stepped_to_boundary)
            !=NOERROR) return UNRECOVERABLE_ERROR;
   }

   // Only use hits that would fall within the range of the reference trajectory
   /*
   for (unsigned int i=0;i<my_cdchits.size();i++){
      DKalmanSIMDCDCHit_t *hit=my_cdchits[i];
      double my_r2=(hit->origin+(z-hit->z0wire)*hit->dir).Mod2();
      if (my_r2>r2) hit->status=bad_hit;
   }
   */

   // If the current length of the trajectory deque is less than the previous 
   // trajectory deque, remove the extra elements and shrink the deque
   if (i<(int)forward_traj.size()){
      forward_traj_length=forward_traj.size();
      for (int j=0;j<forward_traj_length-i;j++){
         forward_traj.pop_front();
      }
   }

   // return an error if there are not enough entries in the trajectory
   if (forward_traj.size()<2) return RESOURCE_UNAVAILABLE;

   if (DEBUG_LEVEL>20)
   {
      cout << "--- Forward cdc trajectory ---" <<endl;
      for (unsigned int m=0;m<forward_traj.size();m++){
         //      DMatrix5x1 S=*(forward_traj[m].S); 
         DMatrix5x1 S=(forward_traj[m].S);
         double tx=S(state_tx),ty=S(state_ty);
         double phi=atan2(ty,tx);
         double cosphi=cos(phi);
         double sinphi=sin(phi);
         double p=fabs(1./S(state_q_over_p));
         double tanl=1./sqrt(tx*tx+ty*ty);
         double sinl=sin(atan(tanl));
         double cosl=cos(atan(tanl));
         cout
            << setiosflags(ios::fixed)<< "pos: " << setprecision(4) 
            << forward_traj[m].S(state_x) << ", "
            << forward_traj[m].S(state_y) << ", "
            << forward_traj[m].z << "  mom: "
            << p*cosphi*cosl<< ", " << p*sinphi*cosl << ", " 
            << p*sinl << " -> " << p
            <<"  s: " << setprecision(3) 	   
            << forward_traj[m].s 
            <<"  t: " << setprecision(3) 	   
            << forward_traj[m].t/SPEED_OF_LIGHT 
            <<"  B: " << forward_traj[m].B 
            << endl;
      }
   }

   // Current state vector
   S=forward_traj[0].S;

   // position at the end of the swim
   x_=forward_traj[0].S(state_x);
   y_=forward_traj[0].S(state_y);
   z_=forward_traj[0].z;

   return NOERROR;
}

// Routine that extracts the state vector propagation part out of the reference
// trajectory loop
jerror_t DTrackFitterKalmanSIMD::PropagateForwardCDC(int length,int &index,
      double &z,double &r2,
      DMatrix5x1 &S,
      bool &stepped_to_boundary){
   DMatrix5x5 J,Q;
   DKalmanForwardTrajectory_t temp;
   int my_i=0;
   temp.h_id=0;
   temp.num_hits=0;
   double dEdx=0.;
   double s_to_boundary=1e6;
   double dz_ds=1./sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));

   // current position 
   DVector3 pos(S(state_x),S(state_y),z);
   temp.z=z;
   // squared radius 
   r2=pos.Perp2();

   temp.s=len;  
   temp.t=ftime;
   temp.S=S;

   // Kinematic variables
   double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
   double one_over_beta2=1.+mass2*q_over_p_sq;
   if (one_over_beta2>BIG) one_over_beta2=BIG;

   // get material properties from the Root Geometry
   if (ENABLE_BOUNDARY_CHECK && fit_type==kTimeBased){
     DVector3 mom(S(state_tx),S(state_ty),1.);
     if(geom->FindMatKalman(pos,mom,temp.K_rho_Z_over_A,
			    temp.rho_Z_over_A,temp.LnI,temp.Z,
			    temp.chi2c_factor,temp.chi2a_factor,
			    temp.chi2a_corr,
			    last_material_map,
			    &s_to_boundary)!=NOERROR){
       return UNRECOVERABLE_ERROR;
     }
   }
   else
     {
       if(geom->FindMatKalman(pos,temp.K_rho_Z_over_A,
			      temp.rho_Z_over_A,temp.LnI,temp.Z,
			      temp.chi2c_factor,temp.chi2a_factor,
			      temp.chi2a_corr,
			      last_material_map)!=NOERROR){
	 return UNRECOVERABLE_ERROR;
       }
     }

   // Get dEdx for the upcoming step
   if (CORRECT_FOR_ELOSS){
      dEdx=GetdEdx(S(state_q_over_p),temp.K_rho_Z_over_A,temp.rho_Z_over_A,
            temp.LnI,temp.Z); 
   }

   index++; 
   if (index<=length){
      my_i=length-index;
      forward_traj[my_i].s=temp.s;
      forward_traj[my_i].t=temp.t;
      forward_traj[my_i].h_id=temp.h_id;
      forward_traj[my_i].num_hits=0;
      forward_traj[my_i].z=temp.z;
      forward_traj[my_i].rho_Z_over_A=temp.rho_Z_over_A;
      forward_traj[my_i].K_rho_Z_over_A=temp.K_rho_Z_over_A;
      forward_traj[my_i].LnI=temp.LnI;
      forward_traj[my_i].Z=temp.Z;
      forward_traj[my_i].S=S;
   } 

   // Determine the step size based on energy loss 
   //double step=mStepSizeS*dz_ds; 
   double max_step_size
      =(z<endplate_z&& r2>endplate_r2min)?mCDCInternalStepSize:mStepSizeS;
   double ds=mStepSizeS;
   if (z<endplate_z && r2<endplate_r2max && z>cdc_origin[2]){
      if (!stepped_to_boundary){
         stepped_to_boundary=false;
         if (fabs(dEdx)>EPS){
            ds=DE_PER_STEP/fabs(dEdx);
         }  
         if (ds>max_step_size) ds=max_step_size;  
         if (s_to_boundary<ds){
            ds=s_to_boundary+EPS3;
            stepped_to_boundary=true;
         }
         if(ds<MIN_STEP_SIZE)ds=MIN_STEP_SIZE; 	   
      }
      else{
         ds=MIN_STEP_SIZE;
         stepped_to_boundary=false;
      }
   }
   double newz=z+ds*dz_ds; // new z position  

   // Store magnetic field
   temp.B=sqrt(Bx*Bx+By*By+Bz*Bz);

   // Step through field
   ds=FasterStep(z,newz,dEdx,S);

   // update path length
   len+=fabs(ds);

   // Update flight time
   ftime+=ds*sqrt(one_over_beta2);// in units where c=1

   // Get the contribution to the covariance matrix due to multiple 
   // scattering
   GetProcessNoise(z,ds,temp.chi2c_factor,temp.chi2a_factor,temp.chi2a_corr,
         temp.S,Q);

   // Energy loss straggling
   if (CORRECT_FOR_ELOSS){
      double varE=GetEnergyVariance(ds,one_over_beta2,temp.K_rho_Z_over_A);
      Q(state_q_over_p,state_q_over_p)=varE*q_over_p_sq*q_over_p_sq*one_over_beta2;   
   }

   // Compute the Jacobian matrix and its transpose
   StepJacobian(newz,z,S,dEdx,J);

   // update the trajectory
   if (index<=length){
      forward_traj[my_i].B=temp.B;
      forward_traj[my_i].Q=Q;
      forward_traj[my_i].J=J;
   }
   else{	
      temp.Q=Q;
      temp.J=J;
      temp.Ckk=Zero5x5;
      temp.Skk=Zero5x1;
      forward_traj.push_front(temp);    
   }

   //update z
   z=newz;

   return NOERROR;
}

// Routine that extracts the state vector propagation part out of the reference
// trajectory loop
jerror_t DTrackFitterKalmanSIMD::PropagateCentral(int length, int &index,
      DVector2 &my_xy,
      double &var_t_factor,
      DMatrix5x1 &Sc,
      bool &stepped_to_boundary){
   DKalmanCentralTrajectory_t temp;
   DMatrix5x5 J;  // State vector Jacobian matrix 
   DMatrix5x5 Q;  // Process noise covariance matrix

   //Initialize some variables needed later
   double dEdx=0.;
   double s_to_boundary=1e6; 
   int my_i=0;
   // Kinematic variables
   double q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
   double q_over_p_sq=q_over_p*q_over_p;
   double one_over_beta2=1.+mass2*q_over_p_sq;
   if (one_over_beta2>BIG) one_over_beta2=BIG;

   // Reset D to zero
   Sc(state_D)=0.;

   temp.xy=my_xy;	
   temp.s=len;
   temp.t=ftime;
   temp.h_id=0;
   temp.S=Sc;

   // Store magnitude of magnetic field
   temp.B=sqrt(Bx*Bx+By*By+Bz*Bz);

   // get material properties from the Root Geometry
   DVector3 pos3d(my_xy.X(),my_xy.Y(),Sc(state_z));
   if (ENABLE_BOUNDARY_CHECK && fit_type==kTimeBased){
     DVector3 mom(cos(Sc(state_phi)),sin(Sc(state_phi)),Sc(state_tanl));
     if(geom->FindMatKalman(pos3d,mom,temp.K_rho_Z_over_A,
			    temp.rho_Z_over_A,temp.LnI,temp.Z,
			    temp.chi2c_factor,temp.chi2a_factor,
			    temp.chi2a_corr,
			    last_material_map,
			    &s_to_boundary)
	!=NOERROR){
       return UNRECOVERABLE_ERROR;
     }
   }
   else if(geom->FindMatKalman(pos3d,temp.K_rho_Z_over_A,
			       temp.rho_Z_over_A,temp.LnI,temp.Z,
			       temp.chi2c_factor,temp.chi2a_factor,
			       temp.chi2a_corr,
			       last_material_map)!=NOERROR){
     return UNRECOVERABLE_ERROR;
   }

   if (CORRECT_FOR_ELOSS){
      dEdx=GetdEdx(q_over_p,temp.K_rho_Z_over_A,temp.rho_Z_over_A,temp.LnI,
            temp.Z);
   }

   // If the deque already exists, update it
   index++; 
   if (index<=length){
      my_i=length-index;
      central_traj[my_i].B=temp.B;
      central_traj[my_i].s=temp.s;
      central_traj[my_i].t=temp.t;
      central_traj[my_i].h_id=0;
      central_traj[my_i].xy=temp.xy;
      central_traj[my_i].rho_Z_over_A=temp.rho_Z_over_A;
      central_traj[my_i].K_rho_Z_over_A=temp.K_rho_Z_over_A;
      central_traj[my_i].LnI=temp.LnI;
      central_traj[my_i].Z=temp.Z;
      central_traj[my_i].S=Sc;
   } 

   // Adjust the step size
   double step_size=mStepSizeS;    
   if (stepped_to_boundary){
      step_size=MIN_STEP_SIZE;
      stepped_to_boundary=false;
   }
   else{
      if (fabs(dEdx)>EPS){
         step_size=DE_PER_STEP/fabs(dEdx);
      } 
      if(step_size>mStepSizeS) step_size=mStepSizeS; 
      if (s_to_boundary<step_size){
         step_size=s_to_boundary+EPS3;
         stepped_to_boundary=true;
      }
      if(step_size<MIN_STEP_SIZE)step_size=MIN_STEP_SIZE;
   } 
   double r2=my_xy.Mod2();
   if (r2>endplate_r2min 
         && step_size>mCDCInternalStepSize) step_size=mCDCInternalStepSize;
   // Propagate the state through the field
   FasterStep(my_xy,step_size,Sc,dEdx);

   // update path length
   len+=step_size;

   // Update flight time
   double dt=step_size*sqrt(one_over_beta2); // in units of c=1
   double one_minus_beta2=1.-1./one_over_beta2;
   ftime+=dt;
   var_ftime+=dt*dt*one_minus_beta2*one_minus_beta2*0.0004;
   var_t_factor=dt*dt*one_minus_beta2*one_minus_beta2;

   //printf("t %f sigt %f\n",TIME_UNIT_CONVERSION*ftime,TIME_UNIT_CONVERSION*sqrt(var_ftime));

   // Multiple scattering    
   GetProcessNoiseCentral(step_size,temp.chi2c_factor,temp.chi2a_factor,
         temp.chi2a_corr,temp.S,Q);

   // Energy loss straggling in the approximation of thick absorbers    
   if (CORRECT_FOR_ELOSS){
      double varE
         =GetEnergyVariance(step_size,one_over_beta2,temp.K_rho_Z_over_A);    
      Q(state_q_over_pt,state_q_over_pt)
         +=varE*temp.S(state_q_over_pt)*temp.S(state_q_over_pt)*one_over_beta2
         *q_over_p_sq;
   }

   // B-field and gradient at current (x,y,z)
   bfield->GetFieldAndGradient(my_xy.X(),my_xy.Y(),Sc(state_z),Bx,By,Bz,
         dBxdx,dBxdy,dBxdz,dBydx,
         dBydy,dBydz,dBzdx,dBzdy,dBzdz);

   // Compute the Jacobian matrix and its transpose
   StepJacobian(my_xy,temp.xy-my_xy,-step_size,Sc,dEdx,J);

   // Update the trajectory info
   if (index<=length){
      central_traj[my_i].Q=Q;
      central_traj[my_i].J=J;
   }
   else{
      temp.Q=Q;
      temp.J=J;
      temp.Ckk=Zero5x5;
      temp.Skk=Zero5x1;
      central_traj.push_front(temp);    
   }

   return NOERROR;
}



// Reference trajectory for central tracks
// At each point we store the state vector and the Jacobian needed to get to this state 
// along s from the previous state.
// The tricky part is that we swim out from the target to find Sc and pos along the trajectory 
// but we need the Jacobians for the opposite direction, because the filter proceeds from 
// the outer hits toward the target.
jerror_t DTrackFitterKalmanSIMD::SetCDCReferenceTrajectory(const DVector2 &xy,
      DMatrix5x1 &Sc){
   bool stepped_to_boundary=false;
   int i=0,central_traj_length=central_traj.size();
   // factor for scaling momentum resolution to propagate variance in flight 
   // time
   double var_t_factor=0; 

   // Magnetic field and gradient at beginning of trajectory
   //bfield->GetField(x_,y_,z_,Bx,By,Bz);
   bfield->GetFieldAndGradient(x_,y_,z_,Bx,By,Bz,
         dBxdx,dBxdy,dBxdz,dBydx,
         dBydy,dBydz,dBzdx,dBzdy,dBzdz);

   // Copy of initial position in xy
   DVector2 my_xy=xy;
   double r2=xy.Mod2(),z=z_;

   // Reset cdc status flags
   for (unsigned int j=0;j<my_cdchits.size();j++){
      if (my_cdchits[j]->status!=late_hit)my_cdchits[j]->status=good_hit;
   }

   // Continue adding to the trajectory until we have reached the endplate
   // or the maximum radius
   while(z<endplate_z && z>=Z_MIN && r2<endplate_r2max
         && fabs(Sc(state_q_over_pt))<Q_OVER_PT_MAX
         && fabs(Sc(state_tanl))<TAN_MAX
        ){
      if (PropagateCentral(central_traj_length,i,my_xy,var_t_factor,Sc,stepped_to_boundary)
            !=NOERROR) return UNRECOVERABLE_ERROR;    
      z=Sc(state_z);
      r2=my_xy.Mod2();
   }

   // If the current length of the trajectory deque is less than the previous 
   // trajectory deque, remove the extra elements and shrink the deque
   if (i<(int)central_traj.size()){
      int central_traj_length=central_traj.size();
      for (int j=0;j<central_traj_length-i;j++){
         central_traj.pop_front();
      }
   }

   // Only use hits that would fall within the range of the reference trajectory
   /*for (unsigned int j=0;j<my_cdchits.size();j++){
      DKalmanSIMDCDCHit_t *hit=my_cdchits[j];
      double my_r2=(hit->origin+(z-hit->z0wire)*hit->dir).Mod2();
      if (my_r2>r2) hit->status=bad_hit;
   }
   */

   // return an error if there are not enough entries in the trajectory
   if (central_traj.size()<2) return RESOURCE_UNAVAILABLE;

   if (DEBUG_LEVEL>20)
   {
      cout << "---------" << central_traj.size() <<" entries------" <<endl;
      for (unsigned int m=0;m<central_traj.size();m++){
         DMatrix5x1 S=central_traj[m].S;
         double cosphi=cos(S(state_phi));
         double sinphi=sin(S(state_phi));
         double pt=fabs(1./S(state_q_over_pt));
         double tanl=S(state_tanl);

         cout
            << m << "::"
            << setiosflags(ios::fixed)<< "pos: " << setprecision(4) 
            << central_traj[m].xy.X() << ", "
            << central_traj[m].xy.Y() << ", "
            << central_traj[m].S(state_z) << "  mom: "
            << pt*cosphi << ", " << pt*sinphi << ", " 
            << pt*tanl << " -> " << pt/cos(atan(tanl))
            <<"  s: " << setprecision(3) 	   
            << central_traj[m].s 
            <<"  t: " << setprecision(3) 	   
            << central_traj[m].t/SPEED_OF_LIGHT 
            <<"  B: " << central_traj[m].B 
            << endl;
      }
   }

   // State at end of swim
   Sc=central_traj[0].S;

   return NOERROR;
}

// Routine that extracts the state vector propagation part out of the reference
// trajectory loop
jerror_t DTrackFitterKalmanSIMD::PropagateForward(int length,int &i,
      double &z,double zhit,
      DMatrix5x1 &S, bool &done,
      bool &stepped_to_boundary,
      bool &stepped_to_endplate){
   DMatrix5x5 J,Q;    
   DKalmanForwardTrajectory_t temp;

   // Initialize some variables
   temp.h_id=0;
   temp.num_hits=0;
   int my_i=0;
   double s_to_boundary=1e6;
   double dz_ds=1./sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));

   // current position
   DVector3 pos(S(state_x),S(state_y),z);
   double r2=pos.Perp2();

   temp.s=len;
   temp.t=ftime;
   temp.z=z;
   temp.S=S;

   // Kinematic variables  
   double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
   double one_over_beta2=1.+mass2*q_over_p_sq;
   if (one_over_beta2>BIG) one_over_beta2=BIG;

   // get material properties from the Root Geometry
   if (ENABLE_BOUNDARY_CHECK && fit_type==kTimeBased){
     DVector3 mom(S(state_tx),S(state_ty),1.);
     if (geom->FindMatKalman(pos,mom,temp.K_rho_Z_over_A,
			     temp.rho_Z_over_A,temp.LnI,temp.Z,
			     temp.chi2c_factor,temp.chi2a_factor,
			     temp.chi2a_corr,
			     last_material_map,
			     &s_to_boundary)
	 !=NOERROR){
       return UNRECOVERABLE_ERROR;      
     }  
   }
   else
     {
       if (geom->FindMatKalman(pos,temp.K_rho_Z_over_A,
			       temp.rho_Z_over_A,temp.LnI,temp.Z,
			       temp.chi2c_factor,temp.chi2a_factor,
			       temp.chi2a_corr,
			       last_material_map)!=NOERROR){
	 return UNRECOVERABLE_ERROR;      
       }       
     }

   // Get dEdx for the upcoming step
   double dEdx=0.;
   if (CORRECT_FOR_ELOSS){
      dEdx=GetdEdx(S(state_q_over_p),temp.K_rho_Z_over_A,
            temp.rho_Z_over_A,temp.LnI,temp.Z);
   }
   i++;
   my_i=length-i;
   if (i<=length){
      forward_traj[my_i].s=temp.s;
      forward_traj[my_i].t=temp.t;
      forward_traj[my_i].h_id=temp.h_id;
      forward_traj[my_i].num_hits=temp.num_hits;
      forward_traj[my_i].z=temp.z;
      forward_traj[my_i].rho_Z_over_A=temp.rho_Z_over_A;
      forward_traj[my_i].K_rho_Z_over_A=temp.K_rho_Z_over_A;
      forward_traj[my_i].LnI=temp.LnI;
      forward_traj[my_i].S=S;
   } 
   else{
      temp.S=S;
   }

   // Determine the step size based on energy loss 
   //  step=mStepSizeS*dz_ds;
   double max_step_size
      =(z<endplate_z&& r2>endplate_r2min)?mCentralStepSize:mStepSizeS;
   double ds=mStepSizeS;
   if (z>cdc_origin[2]){
      if (!stepped_to_boundary){
         stepped_to_boundary=false;
         if (fabs(dEdx)>EPS){
            ds=DE_PER_STEP/fabs(dEdx);
         } 
         if (ds>max_step_size) ds=max_step_size;
         if (s_to_boundary<ds){
            ds=s_to_boundary+EPS3;
            stepped_to_boundary=true;
         }
         if(ds<MIN_STEP_SIZE)ds=MIN_STEP_SIZE;

      }
      else{
         ds=MIN_STEP_SIZE;
         stepped_to_boundary=false;
      }
   }

   double dz=stepped_to_endplate ? endplate_dz : ds*dz_ds;
   double newz=z+dz; // new z position  
   // Check if we are stepping through the CDC endplate
   if (newz>endplate_z && z<endplate_z){
      //  _DBG_ << endl;
      newz=endplate_z+EPS3;
      stepped_to_endplate=true;
   }

   // Check if we are about to step to one of the wire planes
   done=false;
   if (newz>zhit){ 
      newz=zhit;
      done=true;
   }

   // Store magnitude of magnetic field
   temp.B=sqrt(Bx*Bx+By*By+Bz*Bz);

   // Step through field
   ds=FasterStep(z,newz,dEdx,S);

   // update path length
   len+=ds;

   // update flight time
   ftime+=ds*sqrt(one_over_beta2); // in units where c=1

   // Get the contribution to the covariance matrix due to multiple 
   // scattering
   GetProcessNoise(z,ds,temp.chi2c_factor,temp.chi2a_factor,temp.chi2a_corr,
         temp.S,Q);

   // Energy loss straggling  
   if (CORRECT_FOR_ELOSS){
      double varE=GetEnergyVariance(ds,one_over_beta2,temp.K_rho_Z_over_A);	
      Q(state_q_over_p,state_q_over_p)=varE*q_over_p_sq*q_over_p_sq*one_over_beta2;
   }

   // Compute the Jacobian matrix and its transpose
   StepJacobian(newz,z,S,dEdx,J);

   // update the trajectory data
   if (i<=length){
      forward_traj[my_i].B=temp.B;
      forward_traj[my_i].Q=Q;
      forward_traj[my_i].J=J;
   }
   else{
      temp.Q=Q;
      temp.J=J;
      temp.Ckk=Zero5x5;
      temp.Skk=Zero5x1;
      forward_traj.push_front(temp);
   }

   // update z
   z=newz;

   return NOERROR;
}

// Reference trajectory for trajectories with hits in the forward direction
// At each point we store the state vector and the Jacobian needed to get to this state 
// along z from the previous state.
jerror_t DTrackFitterKalmanSIMD::SetReferenceTrajectory(DMatrix5x1 &S){   

   // Magnetic field and gradient at beginning of trajectory
   //bfield->GetField(x_,y_,z_,Bx,By,Bz);
   bfield->GetFieldAndGradient(x_,y_,z_,Bx,By,Bz,
         dBxdx,dBxdy,dBxdz,dBydx,
         dBydy,dBydz,dBzdx,dBzdy,dBzdz);

   // progress in z from hit to hit
   double z=z_;
   int i=0;

   int forward_traj_length=forward_traj.size();
   // loop over the fdc hits   
   double zhit=0.,old_zhit=0.;
   bool stepped_to_boundary=false;
   bool stepped_to_endplate=false;
   unsigned int m=0;
   for (m=0;m<my_fdchits.size();m++){
      if (fabs(S(state_q_over_p))>Q_OVER_P_MAX
	  || fabs(S(state_tx))>TAN_MAX
	  || fabs(S(state_ty))>TAN_MAX
	  || S(state_x)*S(state_x)+S(state_y)*S(state_y)>50.*50.  
	  || z>400. || z<Z_MIN
         ){
         break;
      }

      zhit=my_fdchits[m]->z;
      if (fabs(old_zhit-zhit)>EPS){
         bool done=false;
         while (!done){
	   if (fabs(S(state_q_over_p))>=Q_OVER_P_MAX
	       || fabs(S(state_tx))>TAN_MAX
	       || fabs(S(state_ty))>TAN_MAX  
	       || S(state_x)*S(state_x)+S(state_y)*S(state_y)>50.*50.
	       || z>400. || z< Z_MIN
               ){
               break;
            }

            if (PropagateForward(forward_traj_length,i,z,zhit,S,done,
                     stepped_to_boundary,stepped_to_endplate)
                  !=NOERROR)
               return UNRECOVERABLE_ERROR;
         } 
      }
      old_zhit=zhit;
   }

   // If m<2 then no useable FDC hits survived the check on the magnitude on the 
   // momentum
   if (m<2) return UNRECOVERABLE_ERROR;

   // Make sure the reference trajectory goes one step beyond the most 
   // downstream hit plane
   if (m==my_fdchits.size()){
      bool done=false;  
      if (PropagateForward(forward_traj_length,i,z,400.,S,done,
               stepped_to_boundary,stepped_to_endplate)
            !=NOERROR)
         return UNRECOVERABLE_ERROR;  
      if (PropagateForward(forward_traj_length,i,z,400.,S,done,
               stepped_to_boundary,stepped_to_endplate)
            !=NOERROR)
         return UNRECOVERABLE_ERROR;   
   }

   // Shrink the deque if the new trajectory has less points in it than the 
   // old trajectory
   if (i<(int)forward_traj.size()){
      int mylen=forward_traj.size();
      //_DBG_ << "Shrinking: " << mylen << " to " << i << endl;
      for (int j=0;j<mylen-i;j++){
         forward_traj.pop_front();
      }
      //    _DBG_ << " Now " << forward_traj.size() << endl;
   }

   // If we lopped off some hits on the downstream end, truncate the trajectory to 
   // the point in z just beyond the last valid hit
   unsigned int my_id=my_fdchits.size();
   if (m<my_id){
      if (zhit<z) my_id=m;
      else my_id=m-1;
      zhit=my_fdchits[my_id-1]->z;
      //_DBG_ << "Shrinking: " << forward_traj.size()<< endl;
      for (;;){
         z=forward_traj[0].z;
         if (z<zhit+EPS2) break;
         forward_traj.pop_front();
      }
      //_DBG_ << " Now " << forward_traj.size() << endl;

      // Temporory structure keeping state and trajectory information
      DKalmanForwardTrajectory_t temp;
      temp.h_id=0;
      temp.num_hits=0;
      temp.B=0.;
      temp.K_rho_Z_over_A=temp.rho_Z_over_A=temp.LnI=0.;
      temp.Q=Zero5x5; 

      // last S vector
      S=forward_traj[0].S;

      // Step just beyond the last hit 
      double newz=z+0.01;
      double ds=Step(z,newz,0.,S); // ignore energy loss for this small step
      temp.s=forward_traj[0].s+ds;
      temp.z=newz;
      temp.S=S;

      // Flight time
      double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
      double one_over_beta2=1.+mass2*q_over_p_sq;
      if (one_over_beta2>BIG) one_over_beta2=BIG;
      temp.t=forward_traj[0].t+ds*sqrt(one_over_beta2); // in units where c=1

      // Covariance and state vector needed for smoothing code
      temp.Ckk=Zero5x5;
      temp.Skk=Zero5x1;

      // Jacobian matrices 
      temp.J=I5x5;

      forward_traj.push_front(temp);
   }

   // return an error if there are not enough entries in the trajectory
   if (forward_traj.size()<2) return RESOURCE_UNAVAILABLE;

   // Fill in Lorentz deflection parameters
   for (unsigned int m=0;m<forward_traj.size();m++){
      if (my_id>0){
         unsigned int hit_id=my_id-1;
         double z=forward_traj[m].z;
         if (fabs(z-my_fdchits[hit_id]->z)<EPS2){
            forward_traj[m].h_id=my_id;

            // Get the magnetic field at this position along the trajectory
            bfield->GetField(forward_traj[m].S(state_x),forward_traj[m].S(state_y),
                  z,Bx,By,Bz);
            double Br=sqrt(Bx*Bx+By*By);

            // Angle between B and wire
            double my_phi=0.;
            if (Br>0.) my_phi=acos((Bx*my_fdchits[hit_id]->sina 
                     +By*my_fdchits[hit_id]->cosa)/Br);
            /*
               lorentz_def->GetLorentzCorrectionParameters(forward_traj[m].pos.x(),
               forward_traj[m].pos.y(),
               forward_traj[m].pos.z(),
               tanz,tanr);
               my_fdchits[hit_id]->nr=tanr;
               my_fdchits[hit_id]->nz=tanz;
               */

            my_fdchits[hit_id]->nr=LORENTZ_NR_PAR1*Bz*(1.+LORENTZ_NR_PAR2*Br);
            my_fdchits[hit_id]->nz=(LORENTZ_NZ_PAR1+LORENTZ_NZ_PAR2*Bz)*(Br*cos(my_phi));


            my_id--;

            unsigned int num=1;
            while (hit_id>0 
                  && fabs(my_fdchits[hit_id]->z-my_fdchits[hit_id-1]->z)<EPS){
               hit_id=my_id-1;
               num++;
               my_id--;
            }
            forward_traj[m].num_hits=num;
         }

      }
   }

   if (DEBUG_LEVEL>20)
   {
      cout << "--- Forward fdc trajectory ---" <<endl;
      for (unsigned int m=0;m<forward_traj.size();m++){
         DMatrix5x1 S=(forward_traj[m].S);
         double tx=S(state_tx),ty=S(state_ty);
         double phi=atan2(ty,tx);
         double cosphi=cos(phi);
         double sinphi=sin(phi);
         double p=fabs(1./S(state_q_over_p));
         double tanl=1./sqrt(tx*tx+ty*ty);
         double sinl=sin(atan(tanl));
         double cosl=cos(atan(tanl));
         cout
            << setiosflags(ios::fixed)<< "pos: " << setprecision(4) 
            << forward_traj[m].S(state_x) << ", "
            << forward_traj[m].S(state_y) << ", "
            << forward_traj[m].z << "  mom: "
            << p*cosphi*cosl<< ", " << p*sinphi*cosl << ", " 
            << p*sinl << " -> " << p
            <<"  s: " << setprecision(3) 	   
            << forward_traj[m].s 
            <<"  t: " << setprecision(3) 	   
            << forward_traj[m].t/SPEED_OF_LIGHT 
            <<"  id: " << forward_traj[m].h_id
            << endl;
      }
   }


   // position at the end of the swim
   z_=z;
   x_=S(state_x);
   y_=S(state_y);

   return NOERROR;
}

// Step the state vector through the field from oldz to newz.
// Uses the 4th-order Runga-Kutte algorithm.
double DTrackFitterKalmanSIMD::Step(double oldz,double newz, double dEdx,
      DMatrix5x1 &S){
   double delta_z=newz-oldz;
   if (fabs(delta_z)<EPS) return 0.; // skip if the step is too small

   // Direction tangents
   double tx=S(state_tx);
   double ty=S(state_ty);
   double ds=sqrt(1.+tx*tx+ty*ty)*delta_z;

   double delta_z_over_2=0.5*delta_z;
   double midz=oldz+delta_z_over_2;
   DMatrix5x1 D1,D2,D3,D4;

   //B-field and gradient at  at (x,y,z)
   bfield->GetFieldAndGradient(S(state_x),S(state_y),oldz,Bx,By,Bz, 
         dBxdx,dBxdy,dBxdz,dBydx,
         dBydy,dBydz,dBzdx,dBzdy,dBzdz);
   double Bx0=Bx,By0=By,Bz0=Bz;

   // Calculate the derivative and propagate the state to the next point
   CalcDeriv(oldz,S,dEdx,D1);
   DMatrix5x1 S1=S+delta_z_over_2*D1;

   // Calculate the field at the first intermediate point
   double dx=S1(state_x)-S(state_x);
   double dy=S1(state_y)-S(state_y);
   Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*delta_z_over_2;
   By=By0+dBydx*dx+dBydy*dy+dBydz*delta_z_over_2;
   Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*delta_z_over_2;

   // Calculate the derivative and propagate the state to the next point
   CalcDeriv(midz,S1,dEdx,D2);
   S1=S+delta_z_over_2*D2;

   // Calculate the field at the second intermediate point
   dx=S1(state_x)-S(state_x);
   dy=S1(state_y)-S(state_y);
   Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*delta_z_over_2;
   By=By0+dBydx*dx+dBydy*dy+dBydz*delta_z_over_2;
   Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*delta_z_over_2;

   // Calculate the derivative and propagate the state to the next point
   CalcDeriv(midz,S1,dEdx,D3);
   S1=S+delta_z*D3;

   // Calculate the field at the final point
   dx=S1(state_x)-S(state_x);
   dy=S1(state_y)-S(state_y);
   Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*delta_z;
   By=By0+dBydx*dx+dBydy*dy+dBydz*delta_z;
   Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*delta_z;

   // Final derivative
   CalcDeriv(newz,S1,dEdx,D4);

   //  S+=delta_z*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);
   double dz_over_6=delta_z*ONE_SIXTH;
   double dz_over_3=delta_z*ONE_THIRD;
   S+=dz_over_6*D1;
   S+=dz_over_3*D2;
   S+=dz_over_3*D3;
   S+=dz_over_6*D4;

   // Don't let the magnitude of the momentum drop below some cutoff
   //if (fabs(S(state_q_over_p))>Q_OVER_P_MAX) 
   //  S(state_q_over_p)=Q_OVER_P_MAX*(S(state_q_over_p)>0.0?1.:-1.);
   // Try to keep the direction tangents from heading towards 90 degrees
   //if (fabs(S(state_tx))>TAN_MAX) 
   //  S(state_tx)=TAN_MAX*(S(state_tx)>0.0?1.:-1.); 
   //if (fabs(S(state_ty))>TAN_MAX) 
   //  S(state_ty)=TAN_MAX*(S(state_ty)>0.0?1.:-1.);

   return ds;
}
// Step the state vector through the field from oldz to newz.
// Uses the 4th-order Runga-Kutte algorithm.
// Uses the gradient to compute the field at the intermediate and last 
// points.
double DTrackFitterKalmanSIMD::FasterStep(double oldz,double newz, double dEdx,
      DMatrix5x1 &S){
   double delta_z=newz-oldz;
   if (fabs(delta_z)<EPS) return 0.; // skip if the step is too small

   // Direction tangents
   double tx=S(state_tx);
   double ty=S(state_ty);
   double ds=sqrt(1.+tx*tx+ty*ty)*delta_z;

   double delta_z_over_2=0.5*delta_z;
   double midz=oldz+delta_z_over_2;
   DMatrix5x1 D1,D2,D3,D4;
   double Bx0=Bx,By0=By,Bz0=Bz;

   // The magnetic field at the beginning of the step is assumed to be 
   // obtained at the end of the previous step through StepJacobian

   // Calculate the derivative and propagate the state to the next point
   CalcDeriv(oldz,S,dEdx,D1);
   DMatrix5x1 S1=S+delta_z_over_2*D1;

   // Calculate the field at the first intermediate point
   double dx=S1(state_x)-S(state_x);
   double dy=S1(state_y)-S(state_y);
   Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*delta_z_over_2;
   By=By0+dBydx*dx+dBydy*dy+dBydz*delta_z_over_2;
   Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*delta_z_over_2;

   // Calculate the derivative and propagate the state to the next point
   CalcDeriv(midz,S1,dEdx,D2);
   S1=S+delta_z_over_2*D2;

   // Calculate the field at the second intermediate point
   dx=S1(state_x)-S(state_x);
   dy=S1(state_y)-S(state_y);
   Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*delta_z_over_2;
   By=By0+dBydx*dx+dBydy*dy+dBydz*delta_z_over_2;
   Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*delta_z_over_2;

   // Calculate the derivative and propagate the state to the next point
   CalcDeriv(midz,S1,dEdx,D3);
   S1=S+delta_z*D3;

   // Calculate the field at the final point
   dx=S1(state_x)-S(state_x);
   dy=S1(state_y)-S(state_y);
   Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*delta_z;
   By=By0+dBydx*dx+dBydy*dy+dBydz*delta_z;
   Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*delta_z;

   // Final derivative
   CalcDeriv(newz,S1,dEdx,D4);

   //  S+=delta_z*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);
   double dz_over_6=delta_z*ONE_SIXTH;
   double dz_over_3=delta_z*ONE_THIRD;
   S+=dz_over_6*D1;
   S+=dz_over_3*D2;
   S+=dz_over_3*D3;
   S+=dz_over_6*D4;

   // Don't let the magnitude of the momentum drop below some cutoff
   //if (fabs(S(state_q_over_p))>Q_OVER_P_MAX) 
   //  S(state_q_over_p)=Q_OVER_P_MAX*(S(state_q_over_p)>0.0?1.:-1.);
   // Try to keep the direction tangents from heading towards 90 degrees
   //if (fabs(S(state_tx))>TAN_MAX) 
   //  S(state_tx)=TAN_MAX*(S(state_tx)>0.0?1.:-1.); 
   //if (fabs(S(state_ty))>TAN_MAX) 
   //  S(state_ty)=TAN_MAX*(S(state_ty)>0.0?1.:-1.);

   return ds;
}



// Compute the Jacobian matrix for the forward parametrization.
jerror_t DTrackFitterKalmanSIMD::StepJacobian(double oldz,double newz,
      const DMatrix5x1 &S,
      double dEdx,DMatrix5x5 &J){
   // Initialize the Jacobian matrix
   //J.Zero();
   //for (int i=0;i<5;i++) J(i,i)=1.;
   J=I5x5;

   // Step in z
   double delta_z=newz-oldz;
   if (fabs(delta_z)<EPS) return NOERROR; //skip if the step is too small 

   // Current values of state vector variables
   double x=S(state_x), y=S(state_y),tx=S(state_tx),ty=S(state_ty);
   double q_over_p=S(state_q_over_p);

   //B-field and field gradient at (x,y,z)
   //if (get_field) 
   bfield->GetFieldAndGradient(x,y,oldz,Bx,By,Bz,dBxdx,dBxdy,
         dBxdz,dBydx,dBydy,
         dBydz,dBzdx,dBzdy,dBzdz);

   // Don't let the magnitude of the momentum drop below some cutoff
   if (fabs(q_over_p)>Q_OVER_P_MAX){
      q_over_p=Q_OVER_P_MAX*(q_over_p>0.0?1.:-1.);
      dEdx=0.;
   }
   // Try to keep the direction tangents from heading towards 90 degrees
   if (fabs(tx)>TAN_MAX) tx=TAN_MAX*(tx>0.0?1.:-1.); 
   if (fabs(ty)>TAN_MAX) ty=TAN_MAX*(ty>0.0?1.:-1.);
   // useful combinations of terms
   double kq_over_p=qBr2p*q_over_p;
   double tx2=tx*tx;
   double twotx2=2.*tx2;
   double ty2=ty*ty;
   double twoty2=2.*ty2;
   double txty=tx*ty;
   double one_plus_tx2=1.+tx2;
   double one_plus_ty2=1.+ty2;
   double one_plus_twotx2_plus_ty2=one_plus_ty2+twotx2;
   double one_plus_twoty2_plus_tx2=one_plus_tx2+twoty2;
   double dsdz=sqrt(1.+tx2+ty2);
   double ds=dsdz*delta_z;
   double kds=qBr2p*ds;
   double kqdz_over_p_over_dsdz=kq_over_p*delta_z/dsdz;
   double kq_over_p_ds=kq_over_p*ds;
   double dtx_Bdep=ty*Bz+txty*Bx-one_plus_tx2*By;
   double dty_Bdep=Bx*one_plus_ty2-txty*By-tx*Bz;
   double Bxty=Bx*ty;
   double Bytx=By*tx;
   double Bztxty=Bz*txty;
   double Byty=By*ty;
   double Bxtx=Bx*tx;

   // Jacobian
   J(state_x,state_tx)=J(state_y,state_ty)=delta_z;
   J(state_tx,state_q_over_p)=kds*dtx_Bdep;
   J(state_ty,state_q_over_p)=kds*dty_Bdep;
   J(state_tx,state_tx)+=kqdz_over_p_over_dsdz*(Bxty*(one_plus_twotx2_plus_ty2)
         -Bytx*(3.*one_plus_tx2+twoty2)
         +Bztxty);
   J(state_tx,state_x)=kq_over_p_ds*(ty*dBzdx+txty*dBxdx-one_plus_tx2*dBydx);
   J(state_ty,state_ty)+=kqdz_over_p_over_dsdz*(Bxty*(3.*one_plus_ty2+twotx2)
         -Bytx*(one_plus_twoty2_plus_tx2)
         -Bztxty);
   J(state_ty,state_y)= kq_over_p_ds*(one_plus_ty2*dBxdy-txty*dBydy-tx*dBzdy);
   J(state_tx,state_ty)=kqdz_over_p_over_dsdz
      *((Bxtx+Bz)*(one_plus_twoty2_plus_tx2)-Byty*one_plus_tx2);
   J(state_tx,state_y)= kq_over_p_ds*(tx*dBzdy+txty*dBxdy-one_plus_tx2*dBydy);
   J(state_ty,state_tx)=-kqdz_over_p_over_dsdz*((Byty+Bz)*(one_plus_twotx2_plus_ty2)
         -Bxtx*one_plus_ty2);
   J(state_ty,state_x)=kq_over_p_ds*(one_plus_ty2*dBxdx-txty*dBydx-tx*dBzdx);
   if (CORRECT_FOR_ELOSS && fabs(dEdx)>EPS){
      double one_over_p_sq=q_over_p*q_over_p;
      double E=sqrt(1./one_over_p_sq+mass2); 
      J(state_q_over_p,state_q_over_p)-=dEdx*ds/E*(2.+3.*mass2*one_over_p_sq);
      double temp=-(q_over_p*one_over_p_sq/dsdz)*E*dEdx*delta_z;
      J(state_q_over_p,state_tx)=tx*temp;
      J(state_q_over_p,state_ty)=ty*temp;
   }

   return NOERROR;
}

// Calculate the derivative for the alternate set of parameters {q/pT, phi, 
// tan(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::CalcDeriv(DVector2 &dpos,const DMatrix5x1 &S,
      double dEdx,DMatrix5x1 &D1){
   //Direction at current point
   double tanl=S(state_tanl);
   // Don't let tanl exceed some maximum
   if (fabs(tanl)>TAN_MAX) tanl=TAN_MAX*(tanl>0.0?1.:-1.);  

   double phi=S(state_phi);
   double cosphi=cos(phi);
   double sinphi=sin(phi);
   double lambda=atan(tanl);
   double sinl=sin(lambda);
   double cosl=cos(lambda);
   // Other parameters
   double q_over_pt=S(state_q_over_pt);
   double pt=fabs(1./q_over_pt);

   // Turn off dEdx if the pt drops below some minimum
   if (pt<PT_MIN) {
      dEdx=0.;
   }
   double kq_over_pt=qBr2p*q_over_pt;

   // Derivative of S with respect to s
   double By_cosphi_minus_Bx_sinphi=By*cosphi-Bx*sinphi;
   D1(state_q_over_pt)
      =kq_over_pt*q_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
   double one_over_cosl=1./cosl;
   if (CORRECT_FOR_ELOSS && fabs(dEdx)>EPS){    
      double p=pt*one_over_cosl;
      double p_sq=p*p;
      double E=sqrt(p_sq+mass2);
      D1(state_q_over_pt)-=q_over_pt*E*dEdx/p_sq;
   }
   //  D1(state_phi)
   //  =kq_over_pt*(Bx*cosphi*sinl+By*sinphi*sinl-Bz*cosl);
   D1(state_phi)=kq_over_pt*((Bx*cosphi+By*sinphi)*sinl-Bz*cosl);
   D1(state_tanl)=kq_over_pt*By_cosphi_minus_Bx_sinphi*one_over_cosl;  
   D1(state_z)=sinl;

   // New direction
   dpos.Set(cosl*cosphi,cosl*sinphi);

   return NOERROR;
}

// Calculate the derivative and Jacobian matrices for the alternate set of 
// parameters {q/pT, phi, tan(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::CalcDerivAndJacobian(const DVector2 &xy,
      DVector2 &dxy,
      const DMatrix5x1 &S,
      double dEdx,
      DMatrix5x5 &J1,
      DMatrix5x1 &D1){  
   //Direction at current point
   double tanl=S(state_tanl);
   // Don't let tanl exceed some maximum
   if (fabs(tanl)>TAN_MAX) tanl=TAN_MAX*(tanl>0.0?1.:-1.);  

   double phi=S(state_phi);
   double cosphi=cos(phi);
   double sinphi=sin(phi);
   double lambda=atan(tanl);
   double sinl=sin(lambda);
   double cosl=cos(lambda);
   double cosl2=cosl*cosl;
   double cosl3=cosl*cosl2;
   double one_over_cosl=1./cosl;
   // Other parameters
   double q_over_pt=S(state_q_over_pt);
   double pt=fabs(1./q_over_pt);
   double q=pt*q_over_pt;

   // Turn off dEdx if pt drops below some minimum
   if (pt<PT_MIN) {
      dEdx=0.;
   }
   double kq_over_pt=qBr2p*q_over_pt;

   // B-field and gradient at (x,y,z)
   bfield->GetFieldAndGradient(xy.X(),xy.Y(),S(state_z),Bx,By,Bz,
         dBxdx,dBxdy,dBxdz,dBydx,
         dBydy,dBydz,dBzdx,dBzdy,dBzdz);

   // Derivative of S with respect to s
   double By_cosphi_minus_Bx_sinphi=By*cosphi-Bx*sinphi;
   double By_sinphi_plus_Bx_cosphi=By*sinphi+Bx*cosphi;
   D1(state_q_over_pt)=kq_over_pt*q_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
   D1(state_phi)=kq_over_pt*(By_sinphi_plus_Bx_cosphi*sinl-Bz*cosl);
   D1(state_tanl)=kq_over_pt*By_cosphi_minus_Bx_sinphi*one_over_cosl;
   D1(state_z)=sinl;

   // New direction
   dxy.Set(cosl*cosphi,cosl*sinphi);

   // Jacobian matrix elements
   J1(state_phi,state_phi)=kq_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
   J1(state_phi,state_q_over_pt)
      =qBr2p*(By_sinphi_plus_Bx_cosphi*sinl-Bz*cosl);
   J1(state_phi,state_tanl)=kq_over_pt*(By_sinphi_plus_Bx_cosphi*cosl
         +Bz*sinl)*cosl2;
   J1(state_phi,state_z)
      =kq_over_pt*(dBxdz*cosphi*sinl+dBydz*sinphi*sinl-dBzdz*cosl);

   J1(state_tanl,state_phi)=-kq_over_pt*By_sinphi_plus_Bx_cosphi*one_over_cosl;
   J1(state_tanl,state_q_over_pt)=qBr2p*By_cosphi_minus_Bx_sinphi*one_over_cosl;
   J1(state_tanl,state_tanl)=kq_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
   J1(state_tanl,state_z)=kq_over_pt*(dBydz*cosphi-dBxdz*sinphi)*one_over_cosl;  
   J1(state_q_over_pt,state_phi)
      =-kq_over_pt*q_over_pt*sinl*By_sinphi_plus_Bx_cosphi;  
   J1(state_q_over_pt,state_q_over_pt)
      =2.*kq_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
   J1(state_q_over_pt,state_tanl)
      =kq_over_pt*q_over_pt*cosl3*By_cosphi_minus_Bx_sinphi;
   if (CORRECT_FOR_ELOSS && fabs(dEdx)>EPS){  
      double p=pt*one_over_cosl;
      double p_sq=p*p;
      double m2_over_p2=mass2/p_sq;
      double E=sqrt(p_sq+mass2);

      D1(state_q_over_pt)-=q_over_pt*E/p_sq*dEdx;
      J1(state_q_over_pt,state_q_over_pt)-=dEdx*(2.+3.*m2_over_p2)/E;
      J1(state_q_over_pt,state_tanl)+=q*dEdx*sinl*(1.+2.*m2_over_p2)/(p*E);
   }
   J1(state_q_over_pt,state_z)
      =kq_over_pt*q_over_pt*sinl*(dBydz*cosphi-dBxdz*sinphi);
   J1(state_z,state_tanl)=cosl3;

   return NOERROR;
}

// Step the state and the covariance matrix through the field
jerror_t DTrackFitterKalmanSIMD::StepStateAndCovariance(DVector2 &xy,
      double ds,
      double dEdx,
      DMatrix5x1 &S,
      DMatrix5x5 &J,
      DMatrix5x5 &C){
   //Initialize the Jacobian matrix
   J=I5x5;
   if (fabs(ds)<EPS) return NOERROR; // break out if ds is too small

   // B-field and gradient at current (x,y,z)
   bfield->GetFieldAndGradient(xy.X(),xy.Y(),S(state_z),Bx,By,Bz,
         dBxdx,dBxdy,dBxdz,dBydx,
         dBydy,dBydz,dBzdx,dBzdy,dBzdz);
   double Bx0=Bx,By0=By,Bz0=Bz;

   // Matrices for intermediate steps
   DMatrix5x1 D1,D2,D3,D4;
   DMatrix5x1 S1;
   DMatrix5x5 J1;
   DVector2 dxy1,dxy2,dxy3,dxy4;
   double ds_2=0.5*ds;

   // Find the derivative at the first point, propagate the state to the 
   // first intermediate point and start filling the Jacobian matrix
   CalcDerivAndJacobian(xy,dxy1,S,dEdx,J1,D1);
   S1=S+ds_2*D1; 

   // Calculate the field at the first intermediate point
   double dz=S1(state_z)-S(state_z);
   double dx=ds_2*dxy1.X();
   double dy=ds_2*dxy1.Y();  
   Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
   By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
   Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

   // Calculate the derivative and propagate the state to the next point
   CalcDeriv(dxy2,S1,dEdx,D2);
   S1=S+ds_2*D2; 

   // Calculate the field at the second intermediate point
   dz=S1(state_z)-S(state_z);
   dx=ds_2*dxy2.X();
   dy=ds_2*dxy2.Y();  
   Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
   By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
   Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

   // Calculate the derivative and propagate the state to the next point
   CalcDeriv(dxy3,S1,dEdx,D3);
   S1=S+ds*D3;

   // Calculate the field at the final point
   dz=S1(state_z)-S(state_z);
   dx=ds*dxy3.X();
   dy=ds*dxy3.Y();  
   Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
   By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
   Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

   // Final derivative
   CalcDeriv(dxy4,S1,dEdx,D4);

   // Position vector increment
   //DVector3 dpos
   //  =ds*(ONE_SIXTH*dpos1+ONE_THIRD*dpos2+ONE_THIRD*dpos3+ONE_SIXTH*dpos4);
   double ds_over_6=ds*ONE_SIXTH;
   double ds_over_3=ds*ONE_THIRD;
   DVector2 dxy=ds_over_6*dxy1;
   dxy+=ds_over_3*dxy2;
   dxy+=ds_over_3*dxy3;
   dxy+=ds_over_6*dxy4;

   // New Jacobian matrix
   J+=ds*J1;

   // Deal with changes in D
   double B=sqrt(Bx0*Bx0+By0*By0+Bz0*Bz0);
   //double qrc_old=qpt/(qBr2p*Bz_);
   double qpt=1./S(state_q_over_pt);
   double q=(qpt>0.)?1.:-1.;
   double qrc_old=qpt/(qBr2p*B);
   double sinphi=sin(S(state_phi));
   double cosphi=cos(S(state_phi));
   double qrc_plus_D=S(state_D)+qrc_old;
   dx=dxy.X();
   dy=dxy.Y();
   double dx_sinphi_minus_dy_cosphi=dx*sinphi-dy*cosphi;
   double rc=sqrt(dxy.Mod2()
         +2.*qrc_plus_D*(dx_sinphi_minus_dy_cosphi)
         +qrc_plus_D*qrc_plus_D);
   double q_over_rc=q/rc;

   J(state_D,state_D)=q_over_rc*(dx_sinphi_minus_dy_cosphi+qrc_plus_D);
   J(state_D,state_q_over_pt)=qpt*qrc_old*(J(state_D,state_D)-1.);
   J(state_D,state_phi)=q_over_rc*qrc_plus_D*(dx*cosphi+dy*sinphi);

   // New xy vector
   xy+=dxy;

   // New state vector
   //S+=ds*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);
   S+=ds_over_6*D1;
   S+=ds_over_3*D2;
   S+=ds_over_3*D3;
   S+=ds_over_6*D4;

   // Don't let the pt drop below some minimum
   //if (fabs(1./S(state_q_over_pt))<PT_MIN) {
   //  S(state_q_over_pt)=(1./PT_MIN)*(S(state_q_over_pt)>0.0?1.:-1.);
   // }
   // Don't let tanl exceed some maximum
   if (fabs(S(state_tanl))>TAN_MAX){
      S(state_tanl)=TAN_MAX*(S(state_tanl)>0.0?1.:-1.);
   }

   // New covariance matrix
   // C=J C J^T
   //C=C.SandwichMultiply(J);
   C=J*C*J.Transpose();

   return NOERROR;
}

// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
// Uses the gradient to compute the field at the intermediate and last 
// points.
jerror_t DTrackFitterKalmanSIMD::FasterStep(DVector2 &xy,double ds,
      DMatrix5x1 &S,
      double dEdx){  
   if (fabs(ds)<EPS) return NOERROR; // break out if ds is too small

   // Matrices for intermediate steps
   DMatrix5x1 D1,D2,D3,D4;
   DMatrix5x1 S1;
   DVector2 dxy1,dxy2,dxy3,dxy4;
   double ds_2=0.5*ds;
   double Bx0=Bx,By0=By,Bz0=Bz;

   // The magnetic field at the beginning of the step is assumed to be 
   // obtained at the end of the previous step through StepJacobian

   // Calculate the derivative and propagate the state to the next point
   CalcDeriv(dxy1,S,dEdx,D1);
   S1=S+ds_2*D1; 

   // Calculate the field at the first intermediate point
   double dz=S1(state_z)-S(state_z);
   double dx=ds_2*dxy1.X();
   double dy=ds_2*dxy1.Y();  
   Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
   By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
   Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

   // Calculate the derivative and propagate the state to the next point
   CalcDeriv(dxy2,S1,dEdx,D2);
   S1=S+ds_2*D2; 

   // Calculate the field at the second intermediate point
   dz=S1(state_z)-S(state_z);
   dx=ds_2*dxy2.X();
   dy=ds_2*dxy2.Y();  
   Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
   By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
   Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

   // Calculate the derivative and propagate the state to the next point
   CalcDeriv(dxy3,S1,dEdx,D3);
   S1=S+ds*D3;

   // Calculate the field at the final point
   dz=S1(state_z)-S(state_z);
   dx=ds*dxy3.X();
   dy=ds*dxy3.Y();  
   Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
   By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
   Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

   // Final derivative
   CalcDeriv(dxy4,S1,dEdx,D4);

   // New state vector
   //  S+=ds*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);
   double ds_over_6=ds*ONE_SIXTH;
   double ds_over_3=ds*ONE_THIRD;
   S+=ds_over_6*D1;
   S+=ds_over_3*D2;
   S+=ds_over_3*D3;
   S+=ds_over_6*D4;

   // New position
   //pos+=ds*(ONE_SIXTH*dpos1+ONE_THIRD*dpos2+ONE_THIRD*dpos3+ONE_SIXTH*dpos4);
   xy+=ds_over_6*dxy1;
   xy+=ds_over_3*dxy2;
   xy+=ds_over_3*dxy3;
   xy+=ds_over_6*dxy4;

   // Don't let the pt drop below some minimum
   //if (fabs(1./S(state_q_over_pt))<PT_MIN) {
   //  S(state_q_over_pt)=(1./PT_MIN)*(S(state_q_over_pt)>0.0?1.:-1.);
   //}
   // Don't let tanl exceed some maximum
   if (fabs(S(state_tanl))>TAN_MAX){
      S(state_tanl)=TAN_MAX*(S(state_tanl)>0.0?1.:-1.);
   }

   return NOERROR;
}

// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::Step(DVector2 &xy,double ds,
      DMatrix5x1 &S,
      double dEdx){  
   if (fabs(ds)<EPS) return NOERROR; // break out if ds is too small

   // B-field and gradient at current (x,y,z)
   bfield->GetFieldAndGradient(xy.X(),xy.Y(),S(state_z),Bx,By,Bz,
         dBxdx,dBxdy,dBxdz,dBydx,
         dBydy,dBydz,dBzdx,dBzdy,dBzdz);
   double Bx0=Bx,By0=By,Bz0=Bz;

   // Matrices for intermediate steps
   DMatrix5x1 D1,D2,D3,D4;
   DMatrix5x1 S1;
   DVector2 dxy1,dxy2,dxy3,dxy4;
   double ds_2=0.5*ds;

   // Find the derivative at the first point, propagate the state to the 
   // first intermediate point
   CalcDeriv(dxy1,S,dEdx,D1);
   S1=S+ds_2*D1; 

   // Calculate the field at the first intermediate point
   double dz=S1(state_z)-S(state_z);
   double dx=ds_2*dxy1.X();
   double dy=ds_2*dxy1.Y();  
   Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
   By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
   Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

   // Calculate the derivative and propagate the state to the next point
   CalcDeriv(dxy2,S1,dEdx,D2);
   S1=S+ds_2*D2; 

   // Calculate the field at the second intermediate point
   dz=S1(state_z)-S(state_z);
   dx=ds_2*dxy2.X();
   dy=ds_2*dxy2.Y();  
   Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
   By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
   Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

   // Calculate the derivative and propagate the state to the next point
   CalcDeriv(dxy3,S1,dEdx,D3);
   S1=S+ds*D3;

   // Calculate the field at the final point
   dz=S1(state_z)-S(state_z);
   dx=ds*dxy3.X();
   dy=ds*dxy3.Y();  
   Bx=Bx0+dBxdx*dx+dBxdy*dy+dBxdz*dz;
   By=By0+dBydx*dx+dBydy*dy+dBydz*dz;
   Bz=Bz0+dBzdx*dx+dBzdy*dy+dBzdz*dz;

   // Final derivative
   CalcDeriv(dxy4,S1,dEdx,D4);

   // New state vector
   //  S+=ds*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);
   double ds_over_6=ds*ONE_SIXTH;
   double ds_over_3=ds*ONE_THIRD;
   S+=ds_over_6*D1;
   S+=ds_over_3*D2;
   S+=ds_over_3*D3;
   S+=ds_over_6*D4;

   // New position
   //pos+=ds*(ONE_SIXTH*dpos1+ONE_THIRD*dpos2+ONE_THIRD*dpos3+ONE_SIXTH*dpos4);
   xy+=ds_over_6*dxy1;
   xy+=ds_over_3*dxy2;
   xy+=ds_over_3*dxy3;
   xy+=ds_over_6*dxy4;

   // Don't let the pt drop below some minimum
   //if (fabs(1./S(state_q_over_pt))<PT_MIN) {
   //  S(state_q_over_pt)=(1./PT_MIN)*(S(state_q_over_pt)>0.0?1.:-1.);
   //}
   // Don't let tanl exceed some maximum
   if (fabs(S(state_tanl))>TAN_MAX){
      S(state_tanl)=TAN_MAX*(S(state_tanl)>0.0?1.:-1.);
   }

   return NOERROR;
}

// Assuming that the magnetic field is constant over the step, use a helical
// model to step directly to the next point along the trajectory.
void DTrackFitterKalmanSIMD::FastStep(double &z,double ds, double dEdx,
				      DMatrix5x1 &S){
  
  // Compute convenience terms involving Bx, By, Bz
  double one_over_p=fabs(S(state_q_over_p));
  double p=1./one_over_p;
  double tx=S(state_tx),ty=S(state_ty);
  double denom=sqrt(1.+tx*tx+ty*ty);
  double px=p*tx/denom;
  double py=p*ty/denom;
  double pz=p/denom;
  double q=S(state_q_over_p)>0?1.:-1.;			 
  double k_q=qBr2p*q;
  double ds_over_p=ds*one_over_p;
  double factor=k_q*(0.25*ds_over_p);
  double Ax=factor*Bx,Ay=factor*By,Az=factor*Bz;
  double Ax2=Ax*Ax,Ay2=Ay*Ay,Az2=Az*Az;
  double AxAy=Ax*Ay,AxAz=Ax*Az,AyAz=Ay*Az;
  double one_plus_Ax2=1.+Ax2;
  double scale=ds_over_p/(one_plus_Ax2+Ay2+Az2);
  
  // Compute new position 
  double dx=scale*(px*one_plus_Ax2+py*(AxAy+Az)+pz*(AxAz-Ay));
  double dy=scale*(px*(AxAy-Az)+py*(1.+Ay2)+pz*(AyAz+Ax));
  double dz=scale*(px*(AxAz+Ay)+py*(AyAz-Ax)+pz*(1.+Az2));
  S(state_x)+=dx;
  S(state_y)+=dy;
  z+=dz;
      
  // Compute new momentum
  px+=k_q*(Bz*dy-By*dz);
  py+=k_q*(Bx*dz-Bz*dx);
  pz+=k_q*(By*dx-Bx*dy); 
  S(state_tx)=px/pz;
  S(state_ty)=py/pz;
  if (fabs(dEdx)>EPS){
    double one_over_p_sq=one_over_p*one_over_p;
    double E=sqrt(1./one_over_p_sq+mass2); 
    S(state_q_over_p)-=S(state_q_over_p)*one_over_p_sq*E*dEdx*ds;    
  }
}
// Assuming that the magnetic field is constant over the step, use a helical
// model to step directly to the next point along the trajectory.
void DTrackFitterKalmanSIMD::FastStep(DVector2 &xy,double ds, double dEdx,
				      DMatrix5x1 &S){
  
  // Compute convenience terms involving Bx, By, Bz
  double pt=fabs(1./S(state_q_over_pt)); 
  double one_over_p=cos(atan(S(state_tanl)))/pt;
  double px=pt*cos(S(state_phi));
  double py=pt*sin(S(state_phi));
  double pz=pt*S(state_tanl);
  double q=S(state_q_over_pt)>0?1.:-1.; 
  double k_q=qBr2p*q;
  double ds_over_p=ds*one_over_p;
  double factor=k_q*(0.25*ds_over_p);
  double Ax=factor*Bx,Ay=factor*By,Az=factor*Bz;
  double Ax2=Ax*Ax,Ay2=Ay*Ay,Az2=Az*Az;
  double AxAy=Ax*Ay,AxAz=Ax*Az,AyAz=Ay*Az;
  double one_plus_Ax2=1.+Ax2;
  double scale=ds_over_p/(one_plus_Ax2+Ay2+Az2);
  
  // Compute new position 
  double dx=scale*(px*one_plus_Ax2+py*(AxAy+Az)+pz*(AxAz-Ay));
  double dy=scale*(px*(AxAy-Az)+py*(1.+Ay2)+pz*(AyAz+Ax));
  double dz=scale*(px*(AxAz+Ay)+py*(AyAz-Ax)+pz*(1.+Az2));
  xy.Set(xy.X()+dx,xy.Y()+dy);
  S(state_z)+=dz;
      
  // Compute new momentum
  px+=k_q*(Bz*dy-By*dz);
  py+=k_q*(Bx*dz-Bz*dx);
  pz+=k_q*(By*dx-Bx*dy);
  pt=sqrt(px*px+py*py); 
  S(state_q_over_pt)=q/pt;
  S(state_phi)=atan2(py,px);
  S(state_tanl)=pz/pt;
  if (fabs(dEdx)>EPS){
    double one_over_p_sq=one_over_p*one_over_p;
    double E=sqrt(1./one_over_p_sq+mass2); 
    S(state_q_over_p)-=S(state_q_over_pt)*one_over_p_sq*E*dEdx*ds;    
  }
}

// Calculate the jacobian matrix for the alternate parameter set 
// {q/pT,phi,tanl(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::StepJacobian(const DVector2 &xy,
      const DVector2 &dxy,
      double ds,const DMatrix5x1 &S,
      double dEdx,DMatrix5x5 &J){
   // Initialize the Jacobian matrix
   //J.Zero();
   //for (int i=0;i<5;i++) J(i,i)=1.;
   J=I5x5;

   if (fabs(ds)<EPS) return NOERROR; // break out if ds is too small
   // B-field and gradient at current (x,y,z)
   //bfield->GetFieldAndGradient(xy.X(),xy.Y(),S(state_z),Bx,By,Bz,
   //dBxdx,dBxdy,dBxdz,dBydx,
   //dBydy,dBydz,dBzdx,dBzdy,dBzdz);

   // Charge
   double q=(S(state_q_over_pt)>0.0)?1.:-1.;

   //kinematic quantities
   double q_over_pt=S(state_q_over_pt);
   double sinphi=sin(S(state_phi));
   double cosphi=cos(S(state_phi));
   double D=S(state_D);
   double lambda=atan(S(state_tanl));
   double sinl=sin(lambda);
   double cosl=cos(lambda);
   double cosl2=cosl*cosl;
   double cosl3=cosl*cosl2;
   double one_over_cosl=1./cosl;
   double pt=fabs(1./q_over_pt);

   // Turn off dEdx if pt drops below some minimum
   if (pt<PT_MIN) {
      dEdx=0.;
   }
   double kds=qBr2p*ds;
   double kq_ds_over_pt=kds*q_over_pt;
   double By_cosphi_minus_Bx_sinphi=By*cosphi-Bx*sinphi;
   double By_sinphi_plus_Bx_cosphi=By*sinphi+Bx*cosphi;

   // Jacobian matrix elements
   J(state_phi,state_phi)+=kq_ds_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
   J(state_phi,state_q_over_pt)=kds*(By_sinphi_plus_Bx_cosphi*sinl-Bz*cosl);
   J(state_phi,state_tanl)=kq_ds_over_pt*(By_sinphi_plus_Bx_cosphi*cosl
         +Bz*sinl)*cosl2;
   J(state_phi,state_z)
      =kq_ds_over_pt*(dBxdz*cosphi*sinl+dBydz*sinphi*sinl-dBzdz*cosl);

   J(state_tanl,state_phi)=-kq_ds_over_pt*By_sinphi_plus_Bx_cosphi*one_over_cosl;
   J(state_tanl,state_q_over_pt)=kds*By_cosphi_minus_Bx_sinphi*one_over_cosl;
   J(state_tanl,state_tanl)+=kq_ds_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
   J(state_tanl,state_z)=kq_ds_over_pt*(dBydz*cosphi-dBxdz*sinphi)*one_over_cosl;  
   J(state_q_over_pt,state_phi)
      =-kq_ds_over_pt*q_over_pt*sinl*By_sinphi_plus_Bx_cosphi;  
   J(state_q_over_pt,state_q_over_pt)
      +=2.*kq_ds_over_pt*sinl*By_cosphi_minus_Bx_sinphi;
   J(state_q_over_pt,state_tanl)
      =kq_ds_over_pt*q_over_pt*cosl3*By_cosphi_minus_Bx_sinphi;
   if (CORRECT_FOR_ELOSS && fabs(dEdx)>EPS){  
      double p=pt*one_over_cosl;
      double p_sq=p*p;
      double m2_over_p2=mass2/p_sq;
      double E=sqrt(p_sq+mass2);
      double dE_over_E=dEdx*ds/E;

      J(state_q_over_pt,state_q_over_pt)-=dE_over_E*(2.+3.*m2_over_p2);
      J(state_q_over_pt,state_tanl)+=q*dE_over_E*sinl*(1.+2.*m2_over_p2)/p;
   }
   J(state_q_over_pt,state_z)
      =kq_ds_over_pt*q_over_pt*sinl*(dBydz*cosphi-dBxdz*sinphi);
   J(state_z,state_tanl)=cosl3*ds;

   // Deal with changes in D
   double B=sqrt(Bx*Bx+By*By+Bz*Bz);
   //double qrc_old=qpt/(qBr2p*fabs(Bz));
   double qpt=FactorForSenseOfRotation/q_over_pt;
   double qrc_old=qpt/(qBr2p*B);
   double qrc_plus_D=D+qrc_old;
   double dx=dxy.X();
   double dy=dxy.Y();
   double dx_sinphi_minus_dy_cosphi=dx*sinphi-dy*cosphi;
   double rc=sqrt(dxy.Mod2()
         +2.*qrc_plus_D*(dx_sinphi_minus_dy_cosphi)
         +qrc_plus_D*qrc_plus_D);
   double q_over_rc=FactorForSenseOfRotation*q/rc;

   J(state_D,state_D)=q_over_rc*(dx_sinphi_minus_dy_cosphi+qrc_plus_D);
   J(state_D,state_q_over_pt)=qpt*qrc_old*(J(state_D,state_D)-1.);
   J(state_D,state_phi)=q_over_rc*qrc_plus_D*(dx*cosphi+dy*sinphi);

   return NOERROR;
}




// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
jerror_t DTrackFitterKalmanSIMD::StepJacobian(const DVector2 &xy,
      double ds,const DMatrix5x1 &S,
      double dEdx,DMatrix5x5 &J){
   // Initialize the Jacobian matrix
   //J.Zero();
   //for (int i=0;i<5;i++) J(i,i)=1.;
   J=I5x5;

   if (fabs(ds)<EPS) return NOERROR; // break out if ds is too small

   // Matrices for intermediate steps
   DMatrix5x5 J1;
   DMatrix5x1 D1;
   DVector2 dxy1;

   // charge
   double q=(S(state_q_over_pt)>0.0)?1.:-1.;
   q*=FactorForSenseOfRotation;

   //kinematic quantities
   double qpt=1./S(state_q_over_pt) * FactorForSenseOfRotation;
   double sinphi=sin(S(state_phi));
   double cosphi=cos(S(state_phi));
   double D=S(state_D);

   CalcDerivAndJacobian(xy,dxy1,S,dEdx,J1,D1);
   // double Bz_=fabs(Bz); // needed for computing D

   // New Jacobian matrix
   J+=ds*J1;

   // change in position
   DVector2 dxy =ds*dxy1;

   // Deal with changes in D
   double B=sqrt(Bx*Bx+By*By+Bz*Bz);
   //double qrc_old=qpt/(qBr2p*Bz_);
   double qrc_old=qpt/(qBr2p*B);
   double qrc_plus_D=D+qrc_old;
   double dx=dxy.X();
   double dy=dxy.Y();
   double dx_sinphi_minus_dy_cosphi=dx*sinphi-dy*cosphi;
   double rc=sqrt(dxy.Mod2()
         +2.*qrc_plus_D*(dx_sinphi_minus_dy_cosphi)
         +qrc_plus_D*qrc_plus_D);
   double q_over_rc=q/rc;

   J(state_D,state_D)=q_over_rc*(dx_sinphi_minus_dy_cosphi+qrc_plus_D);
   J(state_D,state_q_over_pt)=qpt*qrc_old*(J(state_D,state_D)-1.);
   J(state_D,state_phi)=q_over_rc*qrc_plus_D*(dx*cosphi+dy*sinphi);

   return NOERROR;
}

// Compute contributions to the covariance matrix due to multiple scattering
// using the Lynch/Dahl empirical formulas
jerror_t DTrackFitterKalmanSIMD::GetProcessNoiseCentral(double ds,
      double chi2c_factor,
      double chi2a_factor,
      double chi2a_corr,
      const DMatrix5x1 &Sc,
      DMatrix5x5 &Q){
   Q.Zero();
   //return NOERROR;
   if (USE_MULS_COVARIANCE && chi2c_factor>0. && fabs(ds)>EPS){
      double tanl=Sc(state_tanl);
      double tanl2=tanl*tanl;
      double one_plus_tanl2=1.+tanl2;
      double one_over_pt=fabs(Sc(state_q_over_pt)); 
      double my_ds=fabs(ds);
      double my_ds_2=0.5*my_ds;

      Q(state_phi,state_phi)=one_plus_tanl2;
      Q(state_tanl,state_tanl)=one_plus_tanl2*one_plus_tanl2;
      Q(state_q_over_pt,state_q_over_pt)=one_over_pt*one_over_pt*tanl2;
      Q(state_q_over_pt,state_tanl)=Q(state_tanl,state_q_over_pt)
         =Sc(state_q_over_pt)*tanl*one_plus_tanl2;
      Q(state_D,state_D)=ONE_THIRD*ds*ds;

      // I am not sure the following is correct...
      double sign_D=(Sc(state_D)>0.0)?1.:-1.;
      Q(state_D,state_phi)=Q(state_phi,state_D)=sign_D*my_ds_2*sqrt(one_plus_tanl2);
      Q(state_D,state_tanl)=Q(state_tanl,state_D)=sign_D*my_ds_2*one_plus_tanl2;
      Q(state_D,state_q_over_pt)=Q(state_q_over_pt,state_D)=sign_D*my_ds_2*tanl*Sc(state_q_over_pt);

      double one_over_p_sq=one_over_pt*one_over_pt/one_plus_tanl2;
      double one_over_beta2=1.+mass2*one_over_p_sq;
      double chi2c_p_sq=chi2c_factor*my_ds*one_over_beta2;
      double chi2a_p_sq=chi2a_factor*(1.+chi2a_corr*one_over_beta2);
      // F=Fraction of Moliere distribution to be taken into account
      // nu=0.5*chi2c/(chi2a*(1.-F));
      double nu=MOLIERE_RATIO1*chi2c_p_sq/chi2a_p_sq;
      double one_plus_nu=1.+nu;
      // sig2_ms=chi2c*1e-6/(1.+F*F)*((one_plus_nu)/nu*log(one_plus_nu)-1.);
      double sig2_ms=chi2c_p_sq*one_over_p_sq*MOLIERE_RATIO2
         *(one_plus_nu/nu*log(one_plus_nu)-1.);

      Q*=sig2_ms;
   }

   return NOERROR;
}

// Compute contributions to the covariance matrix due to multiple scattering
// using the Lynch/Dahl empirical formulas
jerror_t DTrackFitterKalmanSIMD::GetProcessNoise(double z, double ds,
      double chi2c_factor,
      double chi2a_factor,
      double chi2a_corr,
      const DMatrix5x1 &S,
      DMatrix5x5 &Q){

   Q.Zero();
   //return NOERROR;
   if (USE_MULS_COVARIANCE && chi2c_factor>0. && fabs(ds)>EPS){
      double tx=S(state_tx),ty=S(state_ty);
      double one_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
      double my_ds=fabs(ds);
      double my_ds_2=0.5*my_ds;
      double tx2=tx*tx;
      double ty2=ty*ty;
      double one_plus_tx2=1.+tx2;
      double one_plus_ty2=1.+ty2;
      double tsquare=tx2+ty2;
      double one_plus_tsquare=1.+tsquare;

      Q(state_tx,state_tx)=one_plus_tx2*one_plus_tsquare;
      Q(state_ty,state_ty)=one_plus_ty2*one_plus_tsquare;
      Q(state_tx,state_ty)=Q(state_ty,state_tx)=tx*ty*one_plus_tsquare;

      Q(state_x,state_x)=ONE_THIRD*ds*ds;
      Q(state_y,state_y)=Q(state_x,state_x);
      Q(state_y,state_ty)=Q(state_ty,state_y)
         = my_ds_2*sqrt(one_plus_tsquare*one_plus_ty2);
      Q(state_x,state_tx)=Q(state_tx,state_x)
         = my_ds_2*sqrt(one_plus_tsquare*one_plus_tx2);

      double one_over_beta2=1.+one_over_p_sq*mass2;
      double chi2c_p_sq=chi2c_factor*my_ds*one_over_beta2;  
      double chi2a_p_sq=chi2a_factor*(1.+chi2a_corr*one_over_beta2);
      // F=MOLIERE_FRACTION =Fraction of Moliere distribution to be taken into account
      // nu=0.5*chi2c/(chi2a*(1.-F));
      double nu=MOLIERE_RATIO1*chi2c_p_sq/chi2a_p_sq;
      double one_plus_nu=1.+nu;
      double sig2_ms=chi2c_p_sq*one_over_p_sq*MOLIERE_RATIO2
         *(one_plus_nu/nu*log(one_plus_nu)-1.);

      //   printf("fac %f %f %f\n",chi2c_factor,chi2a_factor,chi2a_corr);
      //printf("omega %f\n",chi2c/chi2a);


      Q*=sig2_ms;
   }

   return NOERROR;
}

// Calculate the energy loss per unit length given properties of the material
// through which a particle of momentum p is passing
double DTrackFitterKalmanSIMD::GetdEdx(double q_over_p,double K_rho_Z_over_A,
      double rho_Z_over_A,double LnI,double Z){
   if (rho_Z_over_A<=0.) return 0.;
   //return 0.;

   double p=fabs(1./q_over_p);
   double betagamma=p/MASS;
   double betagamma2=betagamma*betagamma;
   double gamma2=1.+betagamma2;
   double beta2=betagamma2/gamma2;
   if (beta2<EPS) beta2=EPS;

   // density effect
   double delta=CalcDensityEffect(betagamma,rho_Z_over_A,LnI);

   double dEdx=0.;
   // For particles heavier than electrons:
   if (IsHadron){
      double two_Me_betagamma_sq=two_m_e*betagamma2;
      double Tmax
         =two_Me_betagamma_sq/(1.+2.*sqrt(gamma2)*m_ratio+m_ratio_sq);

      dEdx=K_rho_Z_over_A/beta2*(-log(two_Me_betagamma_sq*Tmax)
            +2.*(LnI + beta2)+delta);
   }
   else{
      // relativistic kinetic energy in units of M_e c^2
      double tau=sqrt(gamma2)-1.;
      double tau_sq=tau*tau;
      double tau_plus_1=tau+1.;
      double tau_plus_2=tau+2.;
      double tau_plus_2_sq=tau_plus_2*tau_plus_2;
      double f=0.; // function that depends on tau; see Leo (2nd ed.), p. 38.
      if (IsElectron){
         f=1.-beta2+(0.125*tau_sq-(2.*tau+1.)*log(2.))/(tau_plus_1*tau_plus_1);
      }
      else{
         f=2.*log(2.)-(beta2/12.)*(23.+14./tau_plus_2+10./tau_plus_2_sq
               +4./(tau_plus_2*tau_plus_2_sq));
      }

      // collision loss (Leo eq. 2.66)
      double dEdx_coll
         =-K_rho_Z_over_A/beta2*(log(0.5*tau_sq*tau_plus_2*m_e_sq)-LnI+f-delta);

      // radiation loss (Leo eqs. 2.74, 2.76 with Z^2 -> Z(Z+1)
      double a=Z*ALPHA;
      double a2=a*a;
      double a4=a2*a2;
      double epsilon=1.-PHOTON_ENERGY_CUTOFF;
      double epsilon2=epsilon*epsilon;
      double f_Z=a2*(1./(1.+a2)+0.20206-0.0369*a2+0.0083*a4-0.002*a2*a4);
      // The expression below is the integral of the photon energy weighted 
      // by the bremsstrahlung cross section up to a maximum photon energy 
      // expressed as a fraction of the incident electron energy. 
      double dEdx_rad=-K_rho_Z_over_A*tau_plus_1*(2.*ALPHA/M_PI)*(Z+1.)
         *((log(183.*pow(Z,-1./3.))-f_Z)
               *(1.-epsilon-(1./3.)*(epsilon2-epsilon*epsilon2))
               +1./18.*(1.-epsilon2));


      // dEdx_rad=0.;

      dEdx=dEdx_coll+dEdx_rad;
   }

   return dEdx;
}

// Calculate the variance in the energy loss in a Gaussian approximation.
// The standard deviation of the energy loss distribution is
//      var=0.1535*density*(Z/A)*x*(1-0.5*beta^2)*Tmax  [MeV]
// where Tmax is the maximum energy transfer.
// (derived from Leo (2nd ed.), eq. 2.100.  Note that I think there is a typo 
// in this formula in the text...)
double DTrackFitterKalmanSIMD::GetEnergyVariance(double ds,
						 double one_over_beta2,
						 double K_rho_Z_over_A){
   if (K_rho_Z_over_A<=0.) return 0.;
   
   double betagamma2=1./(one_over_beta2-1.);
   double gamma2=betagamma2*one_over_beta2;
   double two_Me_betagamma_sq=two_m_e*betagamma2;
   double Tmax=two_Me_betagamma_sq/(1.+2.*sqrt(gamma2)*m_ratio+m_ratio_sq);
   double var=K_rho_Z_over_A*one_over_beta2*fabs(ds)*Tmax*(1.-0.5/one_over_beta2);
   return var;
}

// Interface routine for Kalman filter
jerror_t DTrackFitterKalmanSIMD::KalmanLoop(void){
   if (z_<Z_MIN) return VALUE_OUT_OF_RANGE;

   // Vector to store the list of hits used in the fit for the forward parametrization
   vector<const DCDCTrackHit*>forward_cdc_used_in_fit;

   // State vector and initial guess for covariance matrix
   DMatrix5x1 S0;
   DMatrix5x5 C0;

   chisq_=-1.;

   // Angle with respect to beam line
   double theta_deg=(180/M_PI)*input_params.momentum().Theta();
   //double theta_deg_sq=theta_deg*theta_deg;
   double tanl0=tanl_=tan(M_PI_2-input_params.momentum().Theta());

   // Azimuthal angle
   double phi0=phi_=input_params.momentum().Phi();

   // Guess for momentum error
   double dpt_over_pt=0.1;
   /*
      if (theta_deg<15){
      dpt_over_pt=0.107-0.0178*theta_deg+0.000966*theta_deg_sq;
      }
      else {
      dpt_over_pt=0.0288+0.00579*theta_deg-2.77e-5*theta_deg_sq;
      }
      */
   /* 
      if (theta_deg<28.){
      theta_deg=28.;
      theta_deg_sq=theta_deg*theta_deg;
      }
      else if (theta_deg>125.){          
      theta_deg=125.;
      theta_deg_sq=theta_deg*theta_deg;
      }
      */
   double sig_lambda=0.02;
   double dp_over_p_sq
      =dpt_over_pt*dpt_over_pt+tanl_*tanl_*sig_lambda*sig_lambda;

   // Input charge
   double q=input_params.charge();

   // Input momentum 
   DVector3 pvec=input_params.momentum();
   double p_mag=pvec.Mag();
   double pz=pvec.z();
   double q_over_p0=q_over_p_=q/p_mag;
   double q_over_pt0=q_over_pt_=q/pvec.Perp();

   // Initial position
   double x0=x_=input_params.position().x();
   double y0=y_=input_params.position().y();
   double z0=z_=input_params.position().z();
   
   // Check integrity of input parameters
   if (!isfinite(x0) || !isfinite(y0) || !isfinite(q_over_p0)){
      if (DEBUG_LEVEL>0) _DBG_ << "Invalid input parameters!" <<endl;
      return UNRECOVERABLE_ERROR;
   }

   // Initial direction tangents
   double tx0=tx_=pvec.x()/pz;
   double ty0=ty_=pvec.y()/pz;
   double one_plus_tsquare=1.+tx_*tx_+ty_*ty_;

   // deal with hits in FDC
   double fdc_prob=0.,fdc_chisq=-1.;
   unsigned int fdc_ndf=0;
   if (my_fdchits.size()>0 
         && // Make sure that these parameters are valid for forward-going tracks
         (isfinite(tx0) && isfinite(ty0))
      ){
      if (DEBUG_LEVEL>0){
         _DBG_ << "Using forward parameterization." <<endl;
      }

      // Initial guess for the state vector
      S0(state_x)=x_;
      S0(state_y)=y_;
      S0(state_tx)=tx_;
      S0(state_ty)=ty_;
      S0(state_q_over_p)=q_over_p_;

      // Initial guess for forward representation covariance matrix   
      C0(state_x,state_x)=2.0;
      C0(state_y,state_y)=2.0;  
      C0(state_tx,state_tx)=0.001;
      C0(state_ty,state_ty)=0.001;
      if (theta_deg>12.35)
      {
         double temp=sig_lambda*one_plus_tsquare;
         C0(state_tx,state_tx)=C0(state_ty,state_ty)=temp*temp;
      }
      C0(state_q_over_p,state_q_over_p)=dp_over_p_sq*q_over_p_*q_over_p_;
      C0*=COVARIANCE_SCALE_FACTOR_FORWARD;

      if (my_cdchits.size()>0){
	mCDCInternalStepSize=0.25;
      }

      // The position from the track candidate is reported just outside the 
      // start counter for tracks containing cdc hits. Propagate to the distance
      // of closest approach to the beam line
      if (fit_type==kWireBased) ExtrapolateToVertex(S0);

      kalman_error_t error=ForwardFit(S0,C0); 
      if (error!=FIT_FAILED){
	if (fit_type==kWireBased) return NOERROR;

	if (my_cdchits.size()<6){
	  if (ndf_==0) return UNRECOVERABLE_ERROR;
	  return NOERROR;
	}
	fdc_prob=TMath::Prob(chisq_,ndf_);
	if (fdc_prob>0.001 && error==FIT_SUCCEEDED) return NOERROR;
	fdc_ndf=ndf_;
	fdc_chisq=chisq_;
      }
      if (my_cdchits.size()<6) return UNRECOVERABLE_ERROR;
   }

   // Deal with hits in the CDC 
   if (my_cdchits.size()>5){
      kalman_error_t error=FIT_NOT_DONE;
      kalman_error_t cdc_error=FIT_NOT_DONE;

      // Save the current state of the extrapolation vector if it exists
      map<DetectorSystem_t,vector<Extrapolation_t> >saved_extrapolations;
      if (!extrapolations.empty()){
	saved_extrapolations=extrapolations;
	ClearExtrapolations();
      }
      bool save_IsSmoothed=IsSmoothed;

      // Chi-squared, degrees of freedom, and probability
      double forward_prob=0.;
      double chisq_forward=-1.;
      unsigned int ndof_forward=0;

      // Parameters at "vertex"
      double phi=phi_,q_over_pt=q_over_pt_,tanl=tanl_,x=x_,y=y_,z=z_;
      vector< vector <double> > fcov_save;
      vector<pull_t>pulls_save;
      pulls_save.assign(pulls.begin(),pulls.end());     
      if (!fcov.empty()){
         fcov_save.assign(fcov.begin(),fcov.end());
      }
      if (my_fdchits.size()>0){
         if (error==INVALID_FIT) _DBG_<< "Invalid fit " << fcov.size() << " " << fdc_ndf <<endl;
      }

      // Use forward parameters for CDC-only tracks with theta<THETA_CUT degrees
      if (theta_deg<THETA_CUT){
         if (DEBUG_LEVEL>0){
            _DBG_ << "Using forward parameterization." <<endl;
         }

         // Step size
         mStepSizeS=mCentralStepSize;

         // Initialize the state vector
         S0(state_x)=x_=x0;
         S0(state_y)=y_=y0;
         S0(state_tx)=tx_=tx0;
         S0(state_ty)=ty_=ty0;
         S0(state_q_over_p)=q_over_p_=q_over_p0; 
         z_=z0;

         // Initial guess for forward representation covariance matrix  
         double temp=sig_lambda*one_plus_tsquare;
         C0(state_x,state_x)=4.0;
         C0(state_y,state_y)=4.0;
         C0(state_tx,state_tx)=C0(state_ty,state_ty)=temp*temp;
         C0(state_q_over_p,state_q_over_p)=dp_over_p_sq*q_over_p_*q_over_p_;
	 C0*=COVARIANCE_SCALE_FACTOR_FORWARD;

         //C0*=1.+1./tsquare;

         // The position from the track candidate is reported just outside the 
         // start counter for tracks containing cdc hits. Propagate to the 
         // distance of closest approach to the beam line   
         if (fit_type==kWireBased) ExtrapolateToVertex(S0);

         error=ForwardCDCFit(S0,C0);

         if (error!=FIT_FAILED && error!=EXTRAPOLATION_FAILED){
            // Find the CL of the fit
            forward_prob=TMath::Prob(chisq_,ndf_);
            if (my_fdchits.size()>0){
               if (fdc_ndf==0 || forward_prob>fdc_prob){
                  // We did not end up using the fdc hits after all...
                  fdchits_used_in_fit.clear();
               }
               else{
                  chisq_=fdc_chisq;
                  ndf_=fdc_ndf;
                  x_=x;
                  y_=y;
                  z_=z;
                  phi_=phi;
                  tanl_=tanl;
                  q_over_pt_=q_over_pt;
                  if (!fcov_save.empty()){
                     fcov.assign(fcov_save.begin(),fcov_save.end());
                  }
		  if (!saved_extrapolations.empty()){
		    extrapolations=saved_extrapolations;
		  }
		  IsSmoothed=save_IsSmoothed;
		  pulls.assign(pulls_save.begin(),pulls_save.end());
		  
                  //                         _DBG_ << endl;
                  return NOERROR;
               }
            }
            if (forward_prob>0.001 
                  && error==FIT_SUCCEEDED) return NOERROR;

            // Save the best values for the parameters and chi2 for now
            chisq_forward=chisq_;
            ndof_forward=ndf_;
            x=x_;
            y=y_;
            z=z_;
            phi=phi_;
            tanl=tanl_;
            q_over_pt=q_over_pt_;
            fcov_save.assign(fcov.begin(),fcov.end());	
	    pulls_save.assign(pulls.begin(),pulls.end());
	    save_IsSmoothed=IsSmoothed;
	    if (!extrapolations.empty()){
	      saved_extrapolations=extrapolations;
	      ClearExtrapolations();
	    }

            // Save the list of hits used in the fit
            forward_cdc_used_in_fit.assign(cdchits_used_in_fit.begin(),cdchits_used_in_fit.end());

         }
      }

      // Attempt to fit the track using the central parametrization.
      if (DEBUG_LEVEL>0){
         _DBG_ << "Using central parameterization." <<endl;
      }

      // Step size
      mStepSizeS=mCentralStepSize;

      // Initialize the state vector
      S0(state_q_over_pt)=q_over_pt_=q_over_pt0;
      S0(state_phi)=phi_=phi0;
      S0(state_tanl)=tanl_=tanl0;
      S0(state_z)=z_=z0;  
      S0(state_D)=0.;

      // Initialize the covariance matrix
      double dz=1.0;
      C0(state_z,state_z)=dz*dz;
      C0(state_q_over_pt,state_q_over_pt)
         =q_over_pt_*q_over_pt_*dpt_over_pt*dpt_over_pt;
      double dphi=0.02;
      C0(state_phi,state_phi)=dphi*dphi;
      C0(state_D,state_D)=1.0;
      double tanl2=tanl_*tanl_;
      double one_plus_tanl2=1.+tanl2;
      C0(state_tanl,state_tanl)=(one_plus_tanl2)*(one_plus_tanl2)
         *sig_lambda*sig_lambda;
      C0*=COVARIANCE_SCALE_FACTOR_CENTRAL;

      //if (theta_deg>90.) C0*=1.+5.*tanl2;
      //else C0*=1.+5.*tanl2*tanl2;
      
      mCentralStepSize=0.4;
      mCDCInternalStepSize=0.2;

      // The position from the track candidate is reported just outside the 
      // start counter for tracks containing cdc hits. Propagate to the 
      // distance of closest approach to the beam line     
      DVector2 xy(x0,y0);  
      if (fit_type==kWireBased){  
         ExtrapolateToVertex(xy,S0);
      }
    
      cdc_error=CentralFit(xy,S0,C0);
      if (cdc_error==FIT_SUCCEEDED){
         // if the result of the fit using the forward parameterization succeeded
         // but the chi2 was too high, it still may provide the best estimate for 
         // the track parameters... 
         double central_prob=TMath::Prob(chisq_,ndf_);

         if (central_prob<forward_prob){
            phi_=phi;
            q_over_pt_=q_over_pt;
            tanl_=tanl;
            x_=x;
            y_=y;
            z_=z;
            chisq_=chisq_forward;
            ndf_= ndof_forward;
            fcov.assign(fcov_save.begin(),fcov_save.end());
	    pulls.assign(pulls_save.begin(),pulls_save.end());
	    IsSmoothed=save_IsSmoothed;
	    if (!saved_extrapolations.empty()){
	      extrapolations=saved_extrapolations;
	    }
	    
            cdchits_used_in_fit.assign(forward_cdc_used_in_fit.begin(),forward_cdc_used_in_fit.end());

            // We did not end up using any fdc hits...
            fdchits_used_in_fit.clear();

         }
         return NOERROR;

      }
      // otherwise if the fit using the forward parametrization worked, return that 
      else if (error==FIT_SUCCEEDED || error==LOW_CL_FIT){
         phi_=phi;
         q_over_pt_=q_over_pt;
         tanl_=tanl;
         x_=x;
         y_=y;
         z_=z;
         chisq_=chisq_forward;
         ndf_= ndof_forward;

	 if (!saved_extrapolations.empty()){
	   extrapolations=saved_extrapolations;
	 }
	 IsSmoothed=save_IsSmoothed; 
	 fcov.assign(fcov_save.begin(),fcov_save.end());
	 pulls.assign(pulls_save.begin(),pulls_save.end());
         cdchits_used_in_fit.assign(forward_cdc_used_in_fit.begin(),forward_cdc_used_in_fit.end());	

         // We did not end up using any fdc hits...
         fdchits_used_in_fit.clear();

         return NOERROR;
      }
      else return UNRECOVERABLE_ERROR;
   }

   if (ndf_==0) return UNRECOVERABLE_ERROR;

   return NOERROR;
}

#define ITMAX 20
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b)>=0.0?fabs(a):-fabs(a))
// Routine for finding the minimum of a function bracketed between two values
// (see Numerical Recipes in C, pp. 404-405).
jerror_t DTrackFitterKalmanSIMD::BrentsAlgorithm(double ds1,double ds2,
      double dedx,DVector2 &pos,
      const double z0wire,
      const DVector2 &origin,
      const DVector2 &dir,  
      DMatrix5x1 &Sc,
      double &ds_out){
   double d=0.;
   double e=0.0; // will be distance moved on step before last 
   double ax=0.;
   double bx=-ds1;
   double cx=-ds1-ds2;

   double a=(ax<cx?ax:cx);
   double b=(ax>cx?ax:cx);
   double x=bx,w=bx,v=bx;

   //  printf("ds1 %f ds2 %f\n",ds1,ds2);
   // Initialize return step size
   ds_out=0.;

   // Save the starting position 
   // DVector3 pos0=pos;
   // DMatrix S0(Sc);

   // Step to intermediate point
   Step(pos,x,Sc,dedx);
   // Bail if the transverse momentum has dropped below some minimum
   if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
      if (DEBUG_LEVEL>2)
      {
         _DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
            << endl;
      }
      return VALUE_OUT_OF_RANGE;
   }

   DVector2 wirepos=origin+(Sc(state_z)-z0wire)*dir;
   double u_old=x;
   double u=0.;

   // initialization
   double fw=(pos-wirepos).Mod2();
   double fv=fw,fx=fw;

   // main loop
   for (unsigned int iter=1;iter<=ITMAX;iter++){
      double xm=0.5*(a+b);
      double tol1=EPS2*fabs(x)+ZEPS;
      double tol2=2.0*tol1;

      if (fabs(x-xm)<=(tol2-0.5*(b-a))){
         if (Sc(state_z)<=cdc_origin[2]){
            unsigned int iter2=0;
            double ds_temp=0.;
            while (fabs(Sc(state_z)-cdc_origin[2])>EPS2 && iter2<20){
               u=x-(cdc_origin[2]-Sc(state_z))*sin(atan(Sc(state_tanl)));
               x=u;
               ds_temp+=u_old-u;
               // Bail if the transverse momentum has dropped below some minimum
               if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
                  if (DEBUG_LEVEL>2)
                  {
                     _DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
                        << endl;
                  }
                  return VALUE_OUT_OF_RANGE;
               }

               // Function evaluation
               Step(pos,u_old-u,Sc,dedx);
               u_old=u;
               iter2++;
            }
            //printf("new z %f ds %f \n",pos.z(),x);	
            ds_out=ds_temp;
            return NOERROR;
         }
         else if (Sc(state_z)>=endplate_z){
            unsigned int iter2=0;
            double ds_temp=0.;
            while (fabs(Sc(state_z)-endplate_z)>EPS2 && iter2<20){
               u=x-(endplate_z-Sc(state_z))*sin(atan(Sc(state_tanl)));
               x=u;
               ds_temp+=u_old-u;

               // Bail if the transverse momentum has dropped below some minimum
               if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
                  if (DEBUG_LEVEL>2)
                  {
                     _DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
                        << endl;
                  }
                  return VALUE_OUT_OF_RANGE;
               }

               // Function evaluation
               Step(pos,u_old-u,Sc,dedx);
               u_old=u;
               iter2++;
            }
            //printf("new z %f ds %f \n",pos.z(),x);
            ds_out=ds_temp;
            return NOERROR;	
         }
         ds_out=cx-x;
         return NOERROR;
      }
      // trial parabolic fit
      if (fabs(e)>tol1){
         double x_minus_w=x-w;
         double x_minus_v=x-v;
         double r=x_minus_w*(fx-fv);
         double q=x_minus_v*(fx-fw);
         double p=x_minus_v*q-x_minus_w*r;
         q=2.0*(q-r);
         if (q>0.0) p=-p;
         q=fabs(q);
         double etemp=e;
         e=d;
         if (fabs(p)>=fabs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x))
            // fall back on the Golden Section technique
            d=CGOLD*(e=(x>=xm?a-x:b-x));
         else{
            // parabolic step
            d=p/q;
            u=x+d;
            if (u-a<tol2 || b-u <tol2)
               d=SIGN(tol1,xm-x);
         }						
      } else{
         d=CGOLD*(e=(x>=xm?a-x:b-x));
      }
      u=(fabs(d)>=tol1 ? x+d: x+SIGN(tol1,d));

      // Bail if the transverse momentum has dropped below some minimum
      if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
         if (DEBUG_LEVEL>2)
         {
            _DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
               << endl;
         }
         return VALUE_OUT_OF_RANGE;
      }

      // Function evaluation
      Step(pos,u_old-u,Sc,dedx);
      u_old=u;

      wirepos=origin;
      wirepos+=(Sc(state_z)-z0wire)*dir;
      double fu=(pos-wirepos).Mod2();

      //cout << "Brent: z="<<Sc(state_z) << " d="<<sqrt(fu) << endl;

      if (fu<=fx){
         if (u>=x) a=x; else b=x;
         SHFT(v,w,x,u);
         SHFT(fv,fw,fx,fu);      
      }
      else {
         if (u<x) a=u; else b=u;
         if (fu<=fw || w==x){
            v=w;
            w=u;
            fv=fw;
            fw=fu;
         }
         else if (fu<=fv || v==x || v==w){
            v=u;
            fv=fu;
         }
      }
   }  
   ds_out=cx-x;

   return NOERROR;
}

// Routine for finding the minimum of a function bracketed between two values
// (see Numerical Recipes in C, pp. 404-405).
jerror_t DTrackFitterKalmanSIMD::BrentsAlgorithm(double z,double dz,
      double dedx,
      const double z0wire,
      const DVector2 &origin,
      const DVector2 &dir,
      DMatrix5x1 &S,
      double &dz_out){
   double d=0.,u=0.;
   double e=0.0; // will be distance moved on step before last 
   double ax=0.;
   double bx=-dz;
   double cx=-2.*dz;

   double a=(ax<cx?ax:cx);
   double b=(ax>cx?ax:cx);
   double x=bx,w=bx,v=bx;

   // Initialize dz_out
   dz_out=0.;

   // Step to intermediate point
   double z_new=z+x;
   Step(z,z_new,dedx,S); 
   // Bail if the momentum has dropped below some minimum
   if (fabs(S(state_q_over_p))>Q_OVER_P_MAX){
      if (DEBUG_LEVEL>2)
      {
         _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p))
            << endl;
      }
      return VALUE_OUT_OF_RANGE;
   }

   double dz0wire=z-z0wire;
   DVector2 wirepos=origin+(dz0wire+x)*dir;
   DVector2 pos(S(state_x),S(state_y));
   double z_old=z_new;

   // initialization
   double fw=(pos-wirepos).Mod2();
   double fv=fw;
   double fx=fw;

   // main loop
   for (unsigned int iter=1;iter<=ITMAX;iter++){
      double xm=0.5*(a+b);
      double tol1=EPS2*fabs(x)+ZEPS;
      double tol2=2.0*tol1;
      if (fabs(x-xm)<=(tol2-0.5*(b-a))){
         if (z_new>=endplate_z){
            x=endplate_z-z_new;

            // Bail if the momentum has dropped below some minimum
            if (fabs(S(state_q_over_p))>Q_OVER_P_MAX){
               if (DEBUG_LEVEL>2)
               {
                  _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p))
                     << endl;
               }   
               return VALUE_OUT_OF_RANGE;
            }
            if (!isfinite(S(state_x)) || !isfinite(S(state_y))){
               _DBG_ <<endl;
               return VALUE_OUT_OF_RANGE;    
            }
            Step(z_new,endplate_z,dedx,S);
         }
         dz_out=x;
         return NOERROR;
      }
      // trial parabolic fit
      if (fabs(e)>tol1){
         double x_minus_w=x-w;
         double x_minus_v=x-v;
         double r=x_minus_w*(fx-fv);
         double q=x_minus_v*(fx-fw);
         double p=x_minus_v*q-x_minus_w*r;
         q=2.0*(q-r);
         if (q>0.0) p=-p;
         q=fabs(q);
         double etemp=e;
         e=d;
         if (fabs(p)>=fabs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x))
            // fall back on the Golden Section technique
            d=CGOLD*(e=(x>=xm?a-x:b-x));
         else{
            // parabolic step
            d=p/q;
            u=x+d;
            if (u-a<tol2 || b-u <tol2)
               d=SIGN(tol1,xm-x);
         }						
      } else{
         d=CGOLD*(e=(x>=xm?a-x:b-x));
      }
      u=(fabs(d)>=tol1 ? x+d: x+SIGN(tol1,d));

      // Function evaluation
      //S=S0;
      z_new=z+u;
      // Bail if the momentum has dropped below some minimum
      if (fabs(S(state_q_over_p))>Q_OVER_P_MAX){
         if (DEBUG_LEVEL>2)
         {
            _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p))
               << endl;
         }
         return VALUE_OUT_OF_RANGE;
      }

      Step(z_old,z_new,dedx,S);
      z_old=z_new;

      wirepos=origin;
      wirepos+=(dz0wire+u)*dir;
      pos.Set(S(state_x),S(state_y));
      double fu=(pos-wirepos).Mod2();

      // _DBG_ << "Brent: z="<< z+u << " d^2="<<fu << endl;

      if (fu<=fx){
         if (u>=x) a=x; else b=x;
         SHFT(v,w,x,u);
         SHFT(fv,fw,fx,fu);      
      }
      else {
         if (u<x) a=u; else b=u;
         if (fu<=fw || w==x){
            v=w;
            w=u;
            fv=fw;
            fw=fu;
         }
         else if (fu<=fv || v==x || v==w){
            v=u;
            fv=fu;
         }
      }
   }
   dz_out=x;
   return NOERROR;
}

// Kalman engine for Central tracks; updates the position on the trajectory
// after the last hit (closest to the target) is added
kalman_error_t DTrackFitterKalmanSIMD::KalmanCentral(double anneal_factor,
      DMatrix5x1 &Sc,DMatrix5x5 &Cc,
      DVector2 &xy,double &chisq,
      unsigned int &my_ndf){
   DMatrix1x5 H;  // Track projection matrix
   DMatrix5x1 H_T; // Transpose of track projection matrix
   DMatrix5x5 J;  // State vector Jacobian matrix
   DMatrix5x5 Q;  // Process noise covariance matrix
   DMatrix5x1 K;  // KalmanSIMD gain matrix
   DMatrix5x5 Ctest; // covariance matrix
   // double V=0.2028; //1.56*1.56/12.;  // Measurement variance
   double V=0.0507;
   double InvV; // inverse of variance
   //DMatrix5x1 dS;  // perturbation in state vector
   DMatrix5x1 S0,S0_; // state vector

   // set the used_in_fit flags to false for cdc hits 
   unsigned int num_cdc=cdc_used_in_fit.size();
   for (unsigned int i=0;i<num_cdc;i++) cdc_used_in_fit[i]=false;
   for (unsigned int i=0;i<central_traj.size();i++){
      central_traj[i].h_id=0;
   }

   // Initialize the chi2 for this part of the track
   chisq=0.;
   my_ndf=0;
   double var_cut=NUM_CDC_SIGMA_CUT*NUM_CDC_SIGMA_CUT;
   double my_anneal=anneal_factor*anneal_factor;
   double chi2cut=my_anneal*var_cut;

   // path length increment
   double ds2=0.;

   //printf(">>>>>>>>>>>>>>>>\n");   

   // beginning position
   double phi=Sc(state_phi);
   xy.Set(central_traj[break_point_step_index].xy.X()-Sc(state_D)*sin(phi),
         central_traj[break_point_step_index].xy.Y()+Sc(state_D)*cos(phi));

   // Wire origin and direction
   //  unsigned int cdc_index=my_cdchits.size()-1;
   unsigned int cdc_index=break_point_cdc_index;
   if (break_point_cdc_index<num_cdc-1){
     num_cdc=break_point_cdc_index+1;
   }

   if (cdc_index<MIN_HITS_FOR_REFIT) chi2cut=1000.0;

   // Wire origin and direction
   DVector2 origin=my_cdchits[cdc_index]->origin;
   double z0w=my_cdchits[cdc_index]->z0wire;
   DVector2 dir=my_cdchits[cdc_index]->dir;
   DVector2 wirexy=origin+(Sc(state_z)-z0w)*dir;

   // Save the starting values for C and S in the deque
   central_traj[break_point_step_index].Skk=Sc;
   central_traj[break_point_step_index].Ckk=Cc;

   // doca variables
   double doca2,old_doca2=(xy-wirexy).Mod2();

   // energy loss
   double dedx=0.;

   // Boolean for flagging when we are done with measurements
   bool more_measurements=true;

   // Initialize S0_ and perform the loop over the trajectory
   S0_=central_traj[break_point_step_index].S;

   for (unsigned int k=break_point_step_index+1;k<central_traj.size();k++){
      unsigned int k_minus_1=k-1;

      // Check that C matrix is positive definite
      if (!Cc.IsPosDef()){
         if (DEBUG_LEVEL>0) _DBG_ << "Broken covariance matrix!" <<endl;
         return BROKEN_COVARIANCE_MATRIX;
      }

      // Get the state vector, jacobian matrix, and multiple scattering matrix 
      // from reference trajectory
      S0=central_traj[k].S;
      J=central_traj[k].J;
      Q=central_traj[k].Q;

      //Q.Print();
      //J.Print();

      // State S is perturbation about a seed S0
      //dS=Sc-S0_;
      //dS.Zero();

      // Update the actual state vector and covariance matrix
      Sc=S0+J*(Sc-S0_);
      // Cc=J*(Cc*JT)+Q;   
      // Cc=Q.AddSym(Cc.SandwichMultiply(J));
      Cc=Q.AddSym(J*Cc*J.Transpose());

      // Save the current state and covariance matrix in the deque
      if (fit_type==kTimeBased){
         central_traj[k].Skk=Sc;
         central_traj[k].Ckk=Cc;
      }

      // update position based on new doca to reference trajectory
      xy.Set(central_traj[k].xy.X()-Sc(state_D)*sin(Sc(state_phi)),
            central_traj[k].xy.Y()+Sc(state_D)*cos(Sc(state_phi)));

      // Bail if the position is grossly outside of the tracking volume
      if (xy.Mod2()>R2_MAX || Sc(state_z)<Z_MIN || Sc(state_z)>Z_MAX){
         if (DEBUG_LEVEL>2)
         {
            _DBG_<< "Went outside of tracking volume at z="<<Sc(state_z)
               << " r="<<xy.Mod()<<endl;
            _DBG_ << " Break indexes:  " << break_point_cdc_index <<","
               << break_point_step_index << endl;
         }
         return POSITION_OUT_OF_RANGE;
      }
      // Bail if the transverse momentum has dropped below some minimum
      if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
         if (DEBUG_LEVEL>2)
         {
            _DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
               << " at step " << k 
               << endl;
         }
         return MOMENTUM_OUT_OF_RANGE;
      }


      // Save the current state of the reference trajectory
      S0_=S0;

      // new wire position
      wirexy=origin;
      wirexy+=(Sc(state_z)-z0w)*dir;

      // new doca
      doca2=(xy-wirexy).Mod2();

      // Check if the doca is no longer decreasing
      if (more_measurements && (doca2>old_doca2 && Sc(state_z)>cdc_origin[2])){
         if (my_cdchits[cdc_index]->status==good_hit){
            if (DEBUG_LEVEL>9) {
               _DBG_ << " Good Hit Ring " << my_cdchits[cdc_index]->hit->wire->ring << " Straw " << my_cdchits[cdc_index]->hit->wire->straw << endl;
               _DBG_ << " doca " << sqrt(doca2) << endl;
            }

            // Save values at end of current step
            DVector2 xy0=central_traj[k].xy;

            // Variables for the computation of D at the doca to the wire
            double D=Sc(state_D);
            double q=(Sc(state_q_over_pt)>0.0)?1.:-1.;

            q*=FactorForSenseOfRotation;

            double qpt=1./Sc(state_q_over_pt) * FactorForSenseOfRotation;
            double sinphi=sin(Sc(state_phi));
            double cosphi=cos(Sc(state_phi));
            //double qrc_old=qpt/fabs(qBr2p*bfield->GetBz(pos.x(),pos.y(),pos.z()));
            double qrc_old=qpt/fabs(qBr2p*central_traj[k].B);
            double qrc_plus_D=D+qrc_old;

            // wire direction variable
            double ux=dir.X();
            double uy=dir.Y();
            double cosstereo=my_cdchits[cdc_index]->cosstereo;
            // Variables relating wire direction and track direction
            //double my_ux=ux*sinl-cosl*cosphi;
            //double my_uy=uy*sinl-cosl*sinphi;
            //double denom=my_ux*my_ux+my_uy*my_uy;
            // distance variables
            DVector2 diff,dxy1;

	    // use Brent's algorithm to find the poca to the wire
	    // See Numerical Recipes in C, pp 404-405

	    // dEdx for current position along trajectory
	    double q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));            
	    if (CORRECT_FOR_ELOSS){
	      dedx=GetdEdx(q_over_p, central_traj[k].K_rho_Z_over_A,
			   central_traj[k].rho_Z_over_A,
			   central_traj[k].LnI,central_traj[k].Z);
	    }

	    if (BrentCentral(dedx,xy,z0w,origin,dir,Sc,ds2)!=NOERROR) return MOMENTUM_OUT_OF_RANGE;

	    //Step along the reference trajectory and compute the new covariance matrix
	    StepStateAndCovariance(xy0,ds2,dedx,S0,J,Cc);

	    // Compute the value of D (signed distance to the reference trajectory)
	    // at the doca to the wire
	    dxy1=xy0-central_traj[k].xy;
	    double rc=sqrt(dxy1.Mod2()
			   +2.*qrc_plus_D*(dxy1.X()*sinphi-dxy1.Y()*cosphi)
			   +qrc_plus_D*qrc_plus_D);
	    Sc(state_D)=q*rc-qrc_old;

	    // wire position
	    wirexy=origin;
	    wirexy+=(Sc(state_z)-z0w)*dir;
	    diff=xy-wirexy;

            // prediction for measurement  
            double doca=diff.Mod()+EPS;
            double prediction=doca*cosstereo;

            // Measurement
            double measurement=0.39,tdrift=0.,tcorr=0.,dDdt0=0.;
            if (fit_type==kTimeBased || USE_PASS1_TIME_MODE){	
               // Find offset of wire with respect to the center of the
               // straw at this z position
               const DCDCWire *mywire=my_cdchits[cdc_index]->hit->wire;
               int ring_index=mywire->ring-1;
               int straw_index=mywire->straw-1;
               double dz=Sc(state_z)-z0w;
               double phi_d=diff.Phi();
               double delta
                  =max_sag[ring_index][straw_index]*(1.-dz*dz/5625.)
                  *cos(phi_d + sag_phi_offset[ring_index][straw_index]);
               double dphi=phi_d-mywire->origin.Phi();
               while (dphi>M_PI) dphi-=2*M_PI;
               while (dphi<-M_PI) dphi+=2*M_PI;
               if (mywire->origin.Y()<0) dphi*=-1.;

               tdrift=my_cdchits[cdc_index]->tdrift-mT0
                  -central_traj[k_minus_1].t*TIME_UNIT_CONVERSION;
               double B=central_traj[k_minus_1].B;
               ComputeCDCDrift(dphi,delta,tdrift,B,measurement,V,tcorr);
               if (ALIGNMENT_CENTRAL){
                  double myV=0.;
                  double mytcorr=0.;
                  double d_shifted;
                  double dt=2.0;
                  ComputeCDCDrift(dphi,delta,tdrift+dt,B,d_shifted,myV,mytcorr);
                  dDdt0=(d_shifted-measurement)/dt;
               }

               //_DBG_ << tcorr << " " << dphi << " " << dm << endl;

            }

            // Projection matrix        
            sinphi=sin(Sc(state_phi));
            cosphi=cos(Sc(state_phi));
            double dx=diff.X();
            double dy=diff.Y();
            double cosstereo_over_doca=cosstereo/doca;
            H_T(state_D)=(dy*cosphi-dx*sinphi)*cosstereo_over_doca;
            H_T(state_phi)
               =-Sc(state_D)*cosstereo_over_doca*(dx*cosphi+dy*sinphi);
            H_T(state_z)=-cosstereo_over_doca*(dx*ux+dy*uy);
	    H(state_tanl)=0.;
	    H_T(state_tanl)=0.;
            H(state_D)=H_T(state_D);
            H(state_z)=H_T(state_z);
            H(state_phi)=H_T(state_phi);

            // Difference and inverse of variance
            //InvV=1./(V+H*(Cc*H_T));
            //double Vproj=Cc.SandwichMultiply(H_T);
	    double Vproj=H*Cc*H_T;
            InvV=1./(V+Vproj);
            double dm=measurement-prediction;

            if (InvV<0.){
               if (DEBUG_LEVEL>1)
                  _DBG_ << k <<" "<< central_traj.size()-1<<" Negative variance??? " << V << " " << H*(Cc*H_T) <<endl;

               break_point_cdc_index=(3*num_cdc)/4;
               return NEGATIVE_VARIANCE;
            }
            /*
               if (fabs(cosstereo)<1.){
               printf("t %f delta %f sigma %f V %f chi2 %f\n",my_cdchits[cdc_index]->hit->tdrift-mT0,dm,sqrt(V),1./InvV,dm*dm*InvV);
               }
               */

            // Check how far this hit is from the expected position
            double chi2check=dm*dm*InvV;
            if (DEBUG_LEVEL>9) _DBG_ << " Prediction " << prediction << " Measurement " << measurement << " Chi2 " << chi2check << endl;
            if (chi2check<chi2cut)
            {
               if (DEBUG_LEVEL>9) _DBG_ << " Passed Chi^2 check Ring " << my_cdchits[cdc_index]->hit->wire->ring << " Straw " << my_cdchits[cdc_index]->hit->wire->straw << endl;
               /*
                  if (chi2check>var_cut){
               // Give hits that satisfy the wide cut but are still pretty far
               // from the projected position less weight

               // ad hoc correction 
               double diff = chi2check-var_cut;    
               V*=1.+my_anneal*diff;
               InvV=1./(V+Vproj);
               }
               */
               // Compute Kalman gain matrix
               K=InvV*(Cc*H_T);

               // Update state vector covariance matrix
               //Cc=Cc-(K*(H*Cc));  
               Ctest=Cc.SubSym(K*(H*Cc)); 
               // Joseph form
               // C = (I-KH)C(I-KH)^T + KVK^T
               //Ctest=Cc.SandwichMultiply(I5x5-K*H)+V*MultiplyTranspose(K);
               // Check that Ctest is positive definite
               if (!Ctest.IsPosDef()){
                  if (DEBUG_LEVEL>0) _DBG_ << "Broken covariance matrix!" <<endl;
                  return BROKEN_COVARIANCE_MATRIX;
               }
               bool skip_ring
                  =(my_cdchits[cdc_index]->hit->wire->ring==RING_TO_SKIP);
               //Update covariance matrix  and state vector
               if (skip_ring==false && tdrift >= CDC_T_DRIFT_MIN){
                  Cc=Ctest;
                  Sc+=dm*K;
               }

               // Mark point on ref trajectory with a hit id for the straw
               central_traj[k].h_id=cdc_index+1;
               if (DEBUG_LEVEL>9) _DBG_ << " Marked Trajectory central_traj[k].h_id=cdc_index+1 (k cdc_index)" << k << " " << cdc_index << endl;

               // Save some updated information for this hit   
               double scale=(skip_ring)?1.:(1.-H*K);		   
               cdc_updates[cdc_index].tcorr=tcorr;
               cdc_updates[cdc_index].tdrift=tdrift;
               cdc_updates[cdc_index].doca=measurement;
               cdc_updates[cdc_index].variance=V;
               cdc_updates[cdc_index].dDdt0=dDdt0;
               cdc_used_in_fit[cdc_index]=true;
               if (tdrift < CDC_T_DRIFT_MIN) cdc_used_in_fit[cdc_index]=false;

               // Update chi2 for this hit
               if (skip_ring==false && tdrift >= CDC_T_DRIFT_MIN){
                  chisq+=scale*dm*dm/V;      
                  my_ndf++;
               }
               if (DEBUG_LEVEL>10) 
                  cout 
                     << "ring " << my_cdchits[cdc_index]->hit->wire->ring
                     << " t " << my_cdchits[cdc_index]->hit->tdrift 
                     << " Dm-Dpred " << dm
                     << " chi2 " << (1.-H*K)*dm*dm/V
                     << endl;

               break_point_cdc_index=cdc_index;
               break_point_step_index=k_minus_1;

               //else printf("Negative variance!!!\n");


            }

            // Move back to the right step along the reference trajectory.
	    StepStateAndCovariance(xy,-ds2,dedx,Sc,J,Cc);

            //  Save state and covariance matrix to update vector
            cdc_updates[cdc_index].S=Sc;
            cdc_updates[cdc_index].C=Cc;

	    //Sc.Print();
	    //Cc.Print();

            // update position on current trajectory based on corrected doca to 
            // reference trajectory
            xy.Set(central_traj[k].xy.X()-Sc(state_D)*sin(Sc(state_phi)),
                  central_traj[k].xy.Y()+Sc(state_D)*cos(Sc(state_phi)));

         }

         // new wire origin and direction
         if (cdc_index>0){
            cdc_index--;
            origin=my_cdchits[cdc_index]->origin;
            z0w=my_cdchits[cdc_index]->z0wire;
            dir=my_cdchits[cdc_index]->dir;
         }
         else{
            more_measurements=false;
         }

         // Update the wire position
         wirexy=origin+(Sc(state_z)-z0w)*dir;

         //s+=ds2;
         // new doca
         doca2=(xy-wirexy).Mod2();
      }

      old_doca2=doca2;

   } 

   // If there are not enough degrees of freedom, something went wrong...
   if (my_ndf<6){
      chisq=-1.;
      my_ndf=0;
      return PRUNED_TOO_MANY_HITS;
   }
   else my_ndf-=5;

   // Check if the momentum is unphysically large
   double p=cos(atan(Sc(state_tanl)))/fabs(Sc(state_q_over_pt));
   if (p>12.0){
      if (DEBUG_LEVEL>2)
      {
         _DBG_ << "Unphysical momentum: P = " << p <<endl;
      }
      return MOMENTUM_OUT_OF_RANGE;
   }

   // Check if we have a kink in the track or threw away too many cdc hits
   if (num_cdc>=MIN_HITS_FOR_REFIT){
      if (break_point_cdc_index>1){
         if (break_point_cdc_index<num_cdc/2){
            break_point_cdc_index=(3*num_cdc)/4;
         }
         return BREAK_POINT_FOUND;
      }

      unsigned int num_good=0; 
      for (unsigned int j=0;j<num_cdc;j++){
         if (cdc_used_in_fit[j]) num_good++;
      }
      if (double(num_good)/double(num_cdc)<MINIMUM_HIT_FRACTION){
         return PRUNED_TOO_MANY_HITS;
      }
   }

   return FIT_SUCCEEDED;
}

// Kalman engine for forward tracks
kalman_error_t DTrackFitterKalmanSIMD::KalmanForward(double fdc_anneal_factor,
      double cdc_anneal_factor,
      DMatrix5x1 &S, 
      DMatrix5x5 &C,
      double &chisq, 
      unsigned int &numdof){
   DMatrix2x1 Mdiff; // difference between measurement and prediction 
   DMatrix2x5 H;  // Track projection matrix
   DMatrix5x2 H_T; // Transpose of track projection matrix 
   DMatrix1x5 Hc;  // Track projection matrix for cdc hits
   DMatrix5x1 Hc_T; // Transpose of track projection matrix for cdc hits
   DMatrix5x5 J;  // State vector Jacobian matrix
   //DMatrix5x5 J_T; // transpose of this matrix
   DMatrix5x5 Q;  // Process noise covariance matrix
   DMatrix5x2 K;  // Kalman gain matrix
   DMatrix5x1 Kc;  // Kalman gain matrix for cdc hits
   DMatrix2x2 V(0.0833,0.,0.,FDC_CATHODE_VARIANCE);  // Measurement covariance matrix
   DMatrix2x1 R;  // Filtered residual
   DMatrix2x2 RC;  // Covariance of filtered residual
   DMatrix5x1 S0,S0_; //State vector 
   DMatrix5x5 Ctest; // Covariance matrix
   DMatrix2x2 InvV; // Inverse of error matrix

   double Vc=0.0507;

   // Vectors for cdc wires
   DVector2 origin,dir,wirepos;
   double z0w=0.; // origin in z for wire

   // Set used_in_fit flags to false for fdc and cdc hits
   unsigned int num_cdc=cdc_used_in_fit.size();
   unsigned int num_fdc=fdc_used_in_fit.size();
   for (unsigned int i=0;i<num_cdc;i++) cdc_used_in_fit[i]=false;
   for (unsigned int i=0;i<num_fdc;i++) fdc_used_in_fit[i]=false;
   for (unsigned int i=0;i<forward_traj.size();i++){
      if (forward_traj[i].h_id>999)
         forward_traj[i].h_id=0;
   }

   // Save the starting values for C and S in the deque
   forward_traj[break_point_step_index].Skk=S;
   forward_traj[break_point_step_index].Ckk=C;

   // Initialize chi squared
   chisq=0;

   // Initialize number of degrees of freedom
   numdof=0;

   double my_cdc_anneal=cdc_anneal_factor*cdc_anneal_factor;
   double my_fdc_anneal=fdc_anneal_factor*fdc_anneal_factor;

   double var_fdc_cut=NUM_FDC_SIGMA_CUT*NUM_FDC_SIGMA_CUT;
   double fdc_chi2cut=my_fdc_anneal*var_fdc_cut;

   double var_cdc_cut=NUM_CDC_SIGMA_CUT*NUM_CDC_SIGMA_CUT;
   double cdc_chi2cut=my_cdc_anneal*var_cdc_cut;

   unsigned int num_fdc_hits=break_point_fdc_index+1;
   unsigned int max_num_fdc_used_in_fit=num_fdc_hits;
   unsigned int num_cdc_hits=my_cdchits.size(); 
   unsigned int cdc_index=0;
   if (num_cdc_hits>0) cdc_index=num_cdc_hits-1;
   bool more_cdc_measurements=(num_cdc_hits>0);
   double old_doca2=1e6;

   if (num_fdc_hits+num_cdc_hits<MIN_HITS_FOR_REFIT){
      cdc_chi2cut=1000.0;
      fdc_chi2cut=1000.0;
   }

   if (more_cdc_measurements){
      origin=my_cdchits[cdc_index]->origin;  
      dir=my_cdchits[cdc_index]->dir;   
      z0w=my_cdchits[cdc_index]->z0wire;
      wirepos=origin+(forward_traj[break_point_step_index].z-z0w)*dir;
   }

   S0_=(forward_traj[break_point_step_index].S);

   if (DEBUG_LEVEL > 25){
      jout << "Entering DTrackFitterKalmanSIMD::KalmanForward ================================" << endl;
      jout << " We have the following starting parameters for our fit. S = State vector, C = Covariance matrix" << endl;
      S.Print();
      C.Print();
      jout << setprecision(6); 
      jout << " There are " << num_cdc <<  " CDC Updates and " << num_fdc << " FDC Updates on this track" << endl;
      jout << " There are " << num_cdc_hits << " CDC Hits and " << num_fdc_hits << " FDC Hits on this track" << endl;
      jout << " With NUM_FDC_SIGMA_CUT = " << NUM_FDC_SIGMA_CUT << " and NUM_CDC_SIGMA_CUT = " << NUM_CDC_SIGMA_CUT << endl;
      jout << "      fdc_anneal_factor = " << fdc_anneal_factor << "     cdc_anneal_factor = " << cdc_anneal_factor << endl;
      jout << " yields     fdc_chi2cut = " << fdc_chi2cut       << "           cdc_chi2cut = " << cdc_chi2cut       << endl; 
      jout << " Starting from break_point_step_index " << break_point_step_index << endl;
      jout << " S0_ which is the state vector of the reference trajectory at the break point step:" << endl;
      S0_.Print();
      jout << " ===== Beginning pass over reference trajectory ======== " << endl;
   }

   for (unsigned int k=break_point_step_index+1;k<forward_traj.size();k++){
      unsigned int k_minus_1=k-1;

      // Check that C matrix is positive definite
      if (!C.IsPosDef()){
         if (DEBUG_LEVEL>0) _DBG_ << "Broken covariance matrix!" <<endl;
         return BROKEN_COVARIANCE_MATRIX;
      }

      // Get the state vector, jacobian matrix, and multiple scattering matrix 
      // from reference trajectory
      S0=(forward_traj[k].S);
      J=(forward_traj[k].J);
      Q=(forward_traj[k].Q);

      // State S is perturbation about a seed S0
      //dS=S-S0_;

      // Update the actual state vector and covariance matrix
      S=S0+J*(S-S0_);

      // Bail if the momentum has dropped below some minimum
      if (fabs(S(state_q_over_p))>=Q_OVER_P_MAX){
         if (DEBUG_LEVEL>2)
         {
            _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p)) << endl;
         }
         break_point_fdc_index=(3*num_fdc)/4;
         return MOMENTUM_OUT_OF_RANGE;
      }


      //C=J*(C*J_T)+Q;   
      //C=Q.AddSym(C.SandwichMultiply(J));
      C=Q.AddSym(J*C*J.Transpose());

      // Save the current state and covariance matrix in the deque
      forward_traj[k].Skk=S;
      forward_traj[k].Ckk=C;

      // Save the current state of the reference trajectory
      S0_=S0;

      // Z position along the trajectory 
      double z=forward_traj[k].z;

      if (DEBUG_LEVEL > 25){
         jout << " At reference trajectory index " << k << " at z=" << z << endl;
         jout << " The State vector from the reference trajectory" << endl;
         S0.Print();
         jout << " The Jacobian matrix " << endl;
         J.Print();
         jout << " The Q matrix "<< endl;
         Q.Print();
         jout << " The updated State vector S=S0+J*(S-S0_)" << endl;
         S.Print();
         jout << " The updated Covariance matrix C=J*(C*J_T)+Q;" << endl;
         C.Print();
      }

      // Add the hit
      if (num_fdc_hits>0){
         if (forward_traj[k].h_id>0 && forward_traj[k].h_id<1000){
            unsigned int id=forward_traj[k].h_id-1;

            // Make the small alignment rotations
            // Use small-angle form.

            // Position and direction from state vector
            double x=S(state_x) + my_fdchits[id]->phiZ*S(state_y);
            double y=S(state_y) - my_fdchits[id]->phiZ*S(state_x);
	    double tx = (S(state_tx) + my_fdchits[id]->phiZ*S(state_ty) - my_fdchits[id]->phiY) ;
            double ty = (S(state_ty) - my_fdchits[id]->phiZ*S(state_tx) + my_fdchits[id]->phiX) ;

            double cosa=my_fdchits[id]->cosa;
            double sina=my_fdchits[id]->sina;
            double u=my_fdchits[id]->uwire;
            double v=my_fdchits[id]->vstrip;

            // Projected position along the wire without doca-dependent corrections
            double vpred_uncorrected=x*sina+y*cosa;

            // Projected postion in the plane of the wires transverse to the wires
            double upred=x*cosa-y*sina;

            // Direction tangent in the u-z plane
            double tu=tx*cosa-ty*sina;
            double alpha=atan(tu);
            double cosalpha=cos(alpha);
            double cosalpha2=cosalpha*cosalpha;
            double sinalpha=sin(alpha);

            // (signed) distance of closest approach to wire
            double du=upred-u;
            double doca=du*cosalpha;

            // Correction for lorentz effect
            double nz=my_fdchits[id]->nz;
            double nr=my_fdchits[id]->nr;
            double nz_sinalpha_plus_nr_cosalpha=nz*sinalpha+nr*cosalpha;

            // Variance in coordinate along wire
            V(1,1)=my_fdchits[id]->vvar;

            // Difference between measurement and projection
            double tv=tx*sina+ty*cosa;
            Mdiff(1)=v-(vpred_uncorrected+doca*(nz_sinalpha_plus_nr_cosalpha
                     -tv*sinalpha
                     ));
            Mdiff(0)=-doca;

            if (fit_type==kTimeBased && USE_FDC_DRIFT_TIMES){
               double drift_time=my_fdchits[id]->t-mT0
                  -forward_traj[k].t*TIME_UNIT_CONVERSION;
               //double drift=DRIFT_SPEED*drift_time*(du>0?1.:-1.); 
               double drift=(du>0.0?1.:-1.)*fdc_drift_distance(drift_time,forward_traj[k].B);

               Mdiff(0)=drift-doca;

               // Variance in drift distance
               V(0,0)=fdc_drift_variance(drift_time);

            }

            // To transform from (x,y) to (u,v), need to do a rotation:
            //   u = x*cosa-y*sina
            //   v = y*cosa+x*sina
            double temp2=nz_sinalpha_plus_nr_cosalpha-tv*sinalpha;
            H_T(state_x,1)=sina+cosa*cosalpha*temp2;	
            H_T(state_y,1)=cosa-sina*cosalpha*temp2;	

            double cos2_minus_sin2=cosalpha2-sinalpha*sinalpha;
            double fac=nz*cos2_minus_sin2-2.*nr*cosalpha*sinalpha;
            double doca_cosalpha=doca*cosalpha;
            double temp=doca_cosalpha*fac;	
            H_T(state_tx,1)=cosa*temp
               -doca_cosalpha*(tu*sina+tv*cosa*cos2_minus_sin2)
               ;
            H_T(state_ty,1)=-sina*temp
               -doca_cosalpha*(tu*cosa-tv*sina*cos2_minus_sin2)
               ;

            H_T(state_x,0)=cosa*cosalpha;
            H_T(state_y,0)=-sina*cosalpha;
            double one_plus_tu2=1.+tu*tu;
            double factor=du*tu/sqrt(one_plus_tu2)/one_plus_tu2;
            H_T(state_ty,0)=sina*factor;
            H_T(state_tx,0)=-cosa*factor;

            // Matrix transpose H_T -> H
            H=Transpose(H_T);

            // Check to see if we have multiple hits in the same plane
            if (!ALIGNMENT_FORWARD && forward_traj[k].num_hits>1){ 
               // If we do have multiple hits, then all of the hits within some
               // validation region are included with weights determined by how
               // close the hits are to the track projection of the state to the
               // "hit space".
               vector<DMatrix5x2> Klist;
               vector<DMatrix2x1> Mlist;
               vector<DMatrix2x5> Hlist;
               vector<DMatrix5x2> HTlist;
               vector<DMatrix2x2> Vlist;
               vector<double>probs;
               vector<unsigned int>used_ids;

               // Deal with the first hit:
               //double Vtemp=V+H*C*H_T;
               DMatrix2x2 Vtemp=V+H*C*H_T;
               InvV=Vtemp.Invert();

               //probability
               double chi2_hit=Vtemp.Chi2(Mdiff);
               double prob_hit=exp(-0.5*chi2_hit)
                  /(M_TWO_PI*sqrt(Vtemp.Determinant()));

               if (DEBUG_LEVEL > 25) jout << " == There are multiple (" << forward_traj[k].num_hits << ") FDC hits" << endl;

               // Cut out outliers
               if (chi2_hit<fdc_chi2cut && my_fdchits[id]->status==good_hit){
                  probs.push_back(prob_hit);
                  Vlist.push_back(V);
                  Hlist.push_back(H);
                  HTlist.push_back(H_T);
                  Mlist.push_back(Mdiff);
                  Klist.push_back(C*H_T*InvV); // Kalman gain

                  used_ids.push_back(id);
                  fdc_used_in_fit[id]=true;
               }

               // loop over the remaining hits
               for (unsigned int m=1;m<forward_traj[k].num_hits;m++){
                  unsigned int my_id=id-m;
                  if (my_fdchits[my_id]->status==good_hit){
                     u=my_fdchits[my_id]->uwire;
                     v=my_fdchits[my_id]->vstrip;
		     
                     // Doca to this wire
                     du=upred-u;
                     doca=du*cosalpha;

                     // variance for coordinate along the wire
                     V(1,1)=my_fdchits[my_id]->vvar;

                     // Difference between measurement and projection
                     Mdiff(1)=v-(vpred_uncorrected+doca*(nz_sinalpha_plus_nr_cosalpha
                              -tv*sinalpha
                              ));
                     Mdiff(0)=-doca;
                     if (fit_type==kTimeBased && USE_FDC_DRIFT_TIMES){
                        double drift_time=my_fdchits[id]->t-mT0
                           -forward_traj[k].t*TIME_UNIT_CONVERSION;
                        //double drift=DRIFT_SPEED*drift_time*(du>0?1.:-1.); 
                        double drift=(du>0.0?1.:-1.)*fdc_drift_distance(drift_time,forward_traj[k].B);

                        Mdiff(0)=drift-doca;

                        // Variance in drift distance
                        V(0,0)=fdc_drift_variance(drift_time);

                     }

                     // Update the terms in H/H_T that depend on the particular hit    
                     doca_cosalpha=doca*cosalpha;
                     temp=doca_cosalpha*fac;	
                     H_T(state_tx,1)=cosa*temp	 
                        -doca_cosalpha*(tu*sina+tv*cosa*cos2_minus_sin2)
                        ;
                     H_T(state_ty,1)=-sina*temp
                        -doca_cosalpha*(tu*cosa-tv*sina*cos2_minus_sin2)
                        ;
                     factor=du*tu/sqrt(one_plus_tu2)/one_plus_tu2;
                     H_T(state_ty,0)=sina*factor;	     
                     H_T(state_tx,0)=-cosa*factor;   

                     // Matrix transpose H_T -> H
                     H(1,state_tx)=H_T(state_tx,1);
                     H(1,state_ty)=H_T(state_ty,1);
                     H(0,state_ty)=H_T(state_ty,0);
                     H(0,state_tx)=H_T(state_tx,0);

                     // Calculate the kalman gain for this hit 
                     ///Vtemp=V+H*C*H_T;
                     Vtemp=V+H*C*H_T;
                     InvV=Vtemp.Invert();

                     // probability
                     double chi2_hit=Vtemp.Chi2(Mdiff);
                     double prob_hit=exp(-0.5*chi2_hit)
                        /(M_TWO_PI*sqrt(Vtemp.Determinant()));

                     // Cut out outliers
                     if(chi2_hit<fdc_chi2cut){	      
                        probs.push_back(prob_hit);	
                        Mlist.push_back(Mdiff);
                        Vlist.push_back(V);
                        Hlist.push_back(H);   
                        HTlist.push_back(H_T);
                        Klist.push_back(C*H_T*InvV);	  

                        used_ids.push_back(my_id);
                        fdc_used_in_fit[my_id]=true;

                     }
                  }
               }
               double prob_tot=1e-100;
               for (unsigned int m=0;m<probs.size();m++){
                  prob_tot+=probs[m];
               }

               // Adjust the state vector and the covariance using the hit 
               //information
               bool skip_plane=(my_fdchits[id]->hit->wire->layer==PLANE_TO_SKIP);
               if (skip_plane==false){
                  DMatrix5x5 sum=I5x5;
                  DMatrix5x5 sum2;
                  for (unsigned int m=0;m<Klist.size();m++){
                     double my_prob=probs[m]/prob_tot;
                     S+=my_prob*(Klist[m]*Mlist[m]);
                     sum-=my_prob*(Klist[m]*Hlist[m]);
                     sum2+=(my_prob*my_prob)*(Klist[m]*Vlist[m]*Transpose(Klist[m]));

                     // Update chi2
                     DMatrix2x2 HK=Hlist[m]*Klist[m];
                     R=Mlist[m]-HK*Mlist[m];
                     RC=Vlist[m]-HK*Vlist[m];
                     chisq+=my_prob*RC.Chi2(R);

                     unsigned int my_id=used_ids[m];  
                     fdc_updates[my_id].V=RC;

                     if (DEBUG_LEVEL > 25) {
                        jout << " Adjusting state vector for FDC hit " << m << " with prob " << my_prob << " S:" << endl;
                        S.Print();
                     }
                  }
		  // C=C.SandwichMultiply(sum)+sum2;  
		  C=sum2.AddSym(sum*C*sum.Transpose());

                  if (DEBUG_LEVEL > 25) { jout << " C: " << endl; C.Print();}
               }

               for (unsigned int m=0;m<Hlist.size();m++){
                  unsigned int my_id=used_ids[m];  
                  fdc_updates[my_id].S=S;
                  fdc_updates[my_id].C=C; 
                  fdc_updates[my_id].tdrift
                     =my_fdchits[my_id]->t-forward_traj[k].t*TIME_UNIT_CONVERSION-mT0;
                  fdc_updates[my_id].tcorr=fdc_updates[my_id].tdrift; // temporary!
                  fdc_updates[my_id].doca=doca;

                  if (skip_plane){
                     fdc_updates[my_id].V=Vlist[m];
                  }
               }

               // update number of degrees of freedom
               if (skip_plane==false){
                  numdof+=2;
               }
            }
            else{
               if (DEBUG_LEVEL > 25) jout << " == There is a single FDC hit on this plane" << endl;

               // Variance for this hit
               DMatrix2x2 Vtemp=V+H*C*H_T;
               InvV=Vtemp.Invert();

               // Check if this hit is an outlier
               double chi2_hit=Vtemp.Chi2(Mdiff);
               if (chi2_hit<fdc_chi2cut){
                  // Compute Kalman gain matrix
                  K=C*H_T*InvV;

                  bool skip_plane=(my_fdchits[id]->hit->wire->layer==PLANE_TO_SKIP);
                  if (skip_plane==false){
                     // Update the state vector 
                     S+=K*Mdiff;

                     // Update state vector covariance matrix
                     //C=C-K*(H*C);    
                     C=C.SubSym(K*(H*C));

                     if (DEBUG_LEVEL > 25) {
                        jout << "S Update: " << endl; S.Print();
                        jout << "C Uodate: " << endl; C.Print();
                     }
                  }

                  // Store the "improved" values for the state vector and covariance
                  fdc_updates[id].S=S;
                  fdc_updates[id].C=C;
                  fdc_updates[id].tdrift
                     =my_fdchits[id]->t-forward_traj[k].t*TIME_UNIT_CONVERSION-mT0;
                  fdc_updates[id].tcorr=fdc_updates[id].tdrift; // temporary!
                  fdc_updates[id].doca=doca;
                  fdc_used_in_fit[id]=true;

                  if (skip_plane==false){  
                     // Filtered residual and covariance of filtered residual
                     R=Mdiff-H*K*Mdiff;   
                     RC=V-H*(C*H_T);

                     fdc_updates[id].V=RC;

                     // Update chi2 for this segment
                     chisq+=RC.Chi2(R);

                     // update number of degrees of freedom
                     numdof+=2;


                     if (DEBUG_LEVEL>20)
                     {
                        printf("hit %d p %5.2f t %f dm %5.2f sig %f chi2 %5.2f z %5.2f\n",
                              id,1./S(state_q_over_p),fdc_updates[id].tdrift,Mdiff(1),
                              sqrt(V(1,1)),RC.Chi2(R),
                              forward_traj[k].z);

                     }
                  }
                  else{
                     fdc_updates[id].V=V;
                  }
		  
                  break_point_fdc_index=id;
                  break_point_step_index=k;		 
               }
            }
            if (num_fdc_hits>=forward_traj[k].num_hits)
               num_fdc_hits-=forward_traj[k].num_hits;
         }
      }
      else if (more_cdc_measurements /* && z<endplate_z*/){   
         // new wire position
         wirepos=origin;
         wirepos+=(z-z0w)*dir;

         // doca variables
         double dx=S(state_x)-wirepos.X();
         double dy=S(state_y)-wirepos.Y();
         double doca2=dx*dx+dy*dy;

         // Check if the doca is no longer decreasing
         if (doca2>old_doca2 /* && z<endplate_z */){
            if(my_cdchits[cdc_index]->status==good_hit){
               double newz=z;

               // energy loss 
               double dedx=0.;

               // track direction variables
               double tx=S(state_tx);
               double ty=S(state_ty);	
               double tanl=1./sqrt(tx*tx+ty*ty);
               double sinl=sin(atan(tanl));

               // Wire direction variables
               double ux=dir.X();
               double uy=dir.Y();
               // Variables relating wire direction and track direction
               double my_ux=tx-ux;
               double my_uy=ty-uy;
               double denom=my_ux*my_ux+my_uy*my_uy;
               double dz=0.;


               // variables for dealing with propagation of S and C if we 
               // need to use Brent's algorithm to find the doca to the wire
               int num_steps=0;
               double dz3=0.;
               double my_dz=0.;

               // if the path length increment is small relative to the radius 
               // of curvature, use a linear approximation to find dz	
               bool do_brent=false;
               double step1=mStepSizeZ;
               double step2=mStepSizeZ;
               if (k>=2){
                  step1=-forward_traj[k].z+forward_traj[k_minus_1].z;
                  step2=-forward_traj[k_minus_1].z+forward_traj[k-2].z;
               }
               //printf("step1 %f step 2 %f \n",step1,step2);
               double two_step=step1+step2;
               if (fabs(qBr2p*S(state_q_over_p)
                        *bfield->GetBz(S(state_x),S(state_y),z)
                        *two_step/sinl)<0.05
                     && denom>EPS)
               {
                  double dzw=z-z0w;
                  dz=-((S(state_x)-origin.X()-ux*dzw)*my_ux
                        +(S(state_y)-origin.Y()-uy*dzw)*my_uy)/denom;

                  if (fabs(dz)>two_step || dz<0){
                     do_brent=true;
                  }
                  else{
                     newz=z+dz;
                     // Check for exiting the straw
                     if (newz>endplate_z){
                        dz=endplate_z-z;
                     }
                  }
               }
               else do_brent=true;
               if (do_brent){
                  if (CORRECT_FOR_ELOSS){
                     dedx=GetdEdx(S(state_q_over_p), 
                           forward_traj[k].K_rho_Z_over_A,
                           forward_traj[k].rho_Z_over_A,
                           forward_traj[k].LnI,forward_traj[k].Z);
                  }

                  // We have bracketed the minimum doca:  use Brent's agorithm
                  if (BrentForward(z,dedx,z0w,origin,dir,S,dz)!=NOERROR){
                     break_point_fdc_index=(3*num_fdc)/4;
                     return MOMENTUM_OUT_OF_RANGE;
                  }

                  // Step the state and covariance through the field
                  if (fabs(dz)>mStepSizeZ){
                     my_dz=(dz>0?1.0:-1.)*mStepSizeZ;
                     num_steps=int(fabs(dz/my_dz));
                     dz3=dz-num_steps*my_dz;

                     double my_z=z;
                     for (int m=0;m<num_steps;m++){
                        newz=my_z+my_dz;

                        // Step current state by my_dz
                        //Step(z,newz,dedx,S);

                        // propagate error matrix to z-position of hit
                        StepJacobian(my_z,newz,S0,dedx,J);
                        C=J*C*J.Transpose();
                        //C=C.SandwichMultiply(J);

                        // Step reference trajectory by my_dz
                        Step(my_z,newz,dedx,S0); 

                        my_z=newz;
                     }

                     newz=my_z+dz3;

                     // Step current state by dz3
                     //Step(my_z,newz,dedx,S);	  

                     // propagate error matrix to z-position of hit
                     StepJacobian(my_z,newz,S0,dedx,J);
                     C=J*C*J.Transpose();
                     //C=C.SandwichMultiply(J);

                     // Step reference trajectory by dz3
                     Step(my_z,newz,dedx,S0); 	    
                  }
                  else{
                     newz = z + dz;

                     // Step current state by dz
                     //Step(z,newz,dedx,S);

                     // propagate error matrix to z-position of hit
                     StepJacobian(z,newz,S0,dedx,J);
                     C=J*C*J.Transpose();
                     //C=C.SandwichMultiply(J);

                     // Step reference trajectory by dz
                     Step(z,newz,dedx,S0); 
                  }
               }

               // Wire position at current z
               wirepos=origin;
               wirepos+=(newz-z0w)*dir;

               double xw=wirepos.X();
               double yw=wirepos.Y();

               // predicted doca taking into account the orientation of the wire
               if (do_brent==false){
                  // In this case we did not have to swim to find the doca and 
                  // the transformation from the state vector space to the 
                  // measurement space is straight-forward.
                  dy=S(state_y)+S(state_ty)*dz-yw;
                  dx=S(state_x)+S(state_tx)*dz-xw;   
               }
               else{ 
                  // In this case we swam the state vector to the position of 
                  // the doca
                  dy=S(state_y)-yw;
                  dx=S(state_x)-xw; 
               }
               double cosstereo=my_cdchits[cdc_index]->cosstereo;
               double d=sqrt(dx*dx+dy*dy)*cosstereo+EPS;

               // Track projection
               double cosstereo2_over_d=cosstereo*cosstereo/d;
               Hc_T(state_x)=dx*cosstereo2_over_d; 
               Hc(state_x)=Hc_T(state_x);
               Hc_T(state_y)=dy*cosstereo2_over_d;	  
               Hc(state_y)=Hc_T(state_y);
               if (do_brent==false){
                  Hc_T(state_ty)=Hc_T(state_y)*dz;
                  Hc(state_ty)=Hc_T(state_ty);	  
                  Hc_T(state_tx)=Hc_T(state_x)*dz;
                  Hc(state_tx)=Hc_T(state_tx);
               }
               else{
                  Hc_T(state_ty)=0.;
                  Hc(state_ty)=0.;
                  Hc_T(state_tx)=0.;
                  Hc(state_tx)=0.;
               }

               //H.Print();

               // The next measurement
               double dm=0.39,tdrift=0.,tcorr=0.,dDdt0=0.;
               if (fit_type==kTimeBased)
               {
                  // Find offset of wire with respect to the center of the
                  // straw at this z position
                  const DCDCWire *mywire=my_cdchits[cdc_index]->hit->wire;
                  int ring_index=mywire->ring-1;
                  int straw_index=mywire->straw-1;
                  double dz=newz-z0w;
                  double phi_d=atan2(dy,dx);
                  double delta
                     =max_sag[ring_index][straw_index]*(1.-dz*dz/5625.)
                     *cos(phi_d + sag_phi_offset[ring_index][straw_index]);
                  double dphi=phi_d-mywire->origin.Phi();
                  while (dphi>M_PI) dphi-=2*M_PI;
                  while (dphi<-M_PI) dphi+=2*M_PI;
                  if (mywire->origin.Y()<0) dphi*=-1.;

                  tdrift=my_cdchits[cdc_index]->tdrift-mT0
                     -forward_traj[k_minus_1].t*TIME_UNIT_CONVERSION;
                  double B=forward_traj[k_minus_1].B;
                  ComputeCDCDrift(dphi,delta,tdrift,B,dm,Vc,tcorr);
                  if (ALIGNMENT_FORWARD){
                     double myV=0.;
                     double mytcorr=0.;
                     double d_shifted;
                     double dt=5.0;
                     // Dont compute this for negative drift times
                     if (tdrift < 0.) d_shifted = dm;
                     else ComputeCDCDrift(dphi,delta,tdrift+dt,B,d_shifted,myV,mytcorr);
                     dDdt0=(d_shifted-dm)/dt;
                  }

		  if (max_num_fdc_used_in_fit>4)
		    {
		    Vc*=CDC_VAR_SCALE_FACTOR;  //de-weight CDC hits 
		  }
                  //_DBG_ << "t " << tdrift << " d " << d << " delta " << delta << " dphi " << atan2(dy,dx)-mywire->origin.Phi() << endl;

                  //_DBG_ << tcorr << " " << dphi << " " << dm << endl;
               }

               // Residual
               double res=dm-d;

               // inverse variance including prediction
               //double InvV1=1./(Vc+H*(C*H_T));
               //double Vproj=C.SandwichMultiply(Hc_T);
	       double Vproj=Hc*C*Hc_T;
               double InvV1=1./(Vc+Vproj);
               if (InvV1<0.){
                  if (DEBUG_LEVEL>0)
                     _DBG_ << "Negative variance???" << endl;
                  return NEGATIVE_VARIANCE;
               }

               // Check if this hit is an outlier
               double chi2_hit=res*res*InvV1;
               if (chi2_hit<cdc_chi2cut){
                  // Compute KalmanSIMD gain matrix
                  Kc=InvV1*(C*Hc_T);


                  // Update state vector covariance matrix
                  //C=C-K*(H*C);    
                  Ctest=C.SubSym(Kc*(Hc*C));
                  //Ctest=C.SandwichMultiply(I5x5-K*H)+Vc*MultiplyTranspose(K);	 
                  // Check that Ctest is positive definite
                  if (!Ctest.IsPosDef()){
                     if (DEBUG_LEVEL>0) _DBG_ << "Broken covariance matrix!" <<endl;
                     return BROKEN_COVARIANCE_MATRIX;
                  }
                  bool skip_ring
                     =(my_cdchits[cdc_index]->hit->wire->ring==RING_TO_SKIP);
                  // update covariance matrix and state vector
                  if (my_cdchits[cdc_index]->hit->wire->ring!=RING_TO_SKIP && tdrift >= CDC_T_DRIFT_MIN){
                     C=Ctest;
                     S+=res*Kc;

                     if (DEBUG_LEVEL > 25) {
                        jout << " == Adding CDC Hit in Ring " << my_cdchits[cdc_index]->hit->wire->ring << endl;
                        jout << " New S: " << endl; S.Print();
                        jout << " New C: " << endl; C.Print();
                     }
                  }

                  // Flag the place along the reference trajectory with hit id
                  forward_traj[k].h_id=1000+cdc_index;

                  // Store updated info related to this hit
                  double scale=(skip_ring)?1.:(1.-Hc*Kc); 
                  cdc_updates[cdc_index].tdrift=tdrift;
                  cdc_updates[cdc_index].tcorr=tcorr;
                  cdc_updates[cdc_index].variance=Vc;
                  cdc_updates[cdc_index].doca=dm;
                  cdc_updates[cdc_index].dDdt0=dDdt0;
                  cdc_used_in_fit[cdc_index]=true;
                  if(tdrift < CDC_T_DRIFT_MIN){
                     //_DBG_ << tdrift << endl;
                     cdc_used_in_fit[cdc_index]=false;
                  }

                  // Update chi2 and number of degrees of freedom for this hit
                  if (skip_ring==false && tdrift >= CDC_T_DRIFT_MIN){
                     chisq+=scale*res*res/Vc;
                     numdof++;
                  }

                  if (DEBUG_LEVEL>10)
                     jout << "Ring " <<  my_cdchits[cdc_index]->hit->wire->ring
                        << " Straw " <<  my_cdchits[cdc_index]->hit->wire->straw
                        << " Pred " << d << " Meas " << dm
                        << " Sigma meas " << sqrt(Vc)
                        << " t " << tcorr
                        << " Chi2 " << (1.-Hc*Kc)*res*res/Vc << endl;

                  break_point_cdc_index=cdc_index;
                  break_point_step_index=k_minus_1;

               }

               // If we had to use Brent's algorithm to find the true doca, 
               // we need to swim the state vector and covariance matrix back to 
               // the appropriate position along the reference trajectory.
               if (do_brent){
                  if (num_steps==0){
                     // Step C back to the z-position on the reference trajectory
                     StepJacobian(newz,z,S0,dedx,J);
                     C=J*C*J.Transpose();
                     //C=C.SandwichMultiply(J);

                     // Step S to current position on the reference trajectory
                     Step(newz,z,dedx,S);

                     // Step S0 to current position on the reference trajectory
                     Step(newz,z,dedx,S0);
                  }
                  else{
                     double my_z=newz;
                     for (int m=0;m<num_steps;m++){
                        z=my_z-my_dz;

                        // Step C along z
                        StepJacobian(my_z,z,S0,dedx,J);
                        C=J*C*J.Transpose();
                        //C=C.SandwichMultiply(J);

                        // Step S along z
                        Step(my_z,z,dedx,S); 

                        // Step S0 along z
                        Step(my_z,z,dedx,S0);

                        my_z=z;
                     }
                     z=my_z-dz3;

                     // Step C back to the z-position on the reference trajectory
                     StepJacobian(my_z,z,S0,dedx,J);
                     C=J*C*J.Transpose();
                     //C=C.SandwichMultiply(J);

                     // Step S to current position on the reference trajectory
                     Step(my_z,z,dedx,S);

                     // Step S to current position on the reference trajectory
                     Step(my_z,z,dedx,S0);
                  }
               }
               cdc_updates[cdc_index].S=S;
               cdc_updates[cdc_index].C=C;	  
            }

            // new wire origin and direction
            if (cdc_index>0){
               cdc_index--;
               origin=my_cdchits[cdc_index]->origin;
               z0w=my_cdchits[cdc_index]->z0wire;
               dir=my_cdchits[cdc_index]->dir;
            }
            else more_cdc_measurements=false;

            // Update the wire position
            wirepos=origin+(z-z0w)*dir;

            // new doca
            dx=S(state_x)-wirepos.X();
            dy=S(state_y)-wirepos.Y();
            doca2=dx*dx+dy*dy;
         }
         old_doca2=doca2;
      }
   }
   // Save final z position
   z_=forward_traj[forward_traj.size()-1].z;

   // The following code segment addes a fake point at a well-defined z position
   // that would correspond to a thin foil target.  It should not be turned on
   // for an extended target.
   if (ADD_VERTEX_POINT){
      double dz_to_target=TARGET_Z-z_;
      double my_dz=mStepSizeZ*(dz_to_target>0?1.:-1.);
      int num_steps=int(fabs(dz_to_target/my_dz));

      for (int k=0;k<num_steps;k++){
         double newz=z_+my_dz;
         // Step C along z
         StepJacobian(z_,newz,S,0.,J);
         C=J*C*J.Transpose();
         //C=C.SandwichMultiply(J);

         // Step S along z
         Step(z_,newz,0.,S);

         z_=newz;
      }

      // Step C along z
      StepJacobian(z_,TARGET_Z,S,0.,J);
      C=J*C*J.Transpose();
      //C=C.SandwichMultiply(J);

      // Step S along z
      Step(z_,TARGET_Z,0.,S);

      z_=TARGET_Z;

      // predicted doca taking into account the orientation of the wire
      double dy=S(state_y);
      double dx=S(state_x);      
      double d=sqrt(dx*dx+dy*dy);

      // Track projection
      double one_over_d=1./d;
      Hc_T(state_x)=dx*one_over_d; 
      Hc(state_x)=Hc_T(state_x);
      Hc_T(state_y)=dy*one_over_d;	  
      Hc(state_y)=Hc_T(state_y);

      // Variance of target point
      // Variance is for average beam spot size assuming triangular distribution
      // out to 2.2 mm from the beam line.
      //   sigma_r = 2.2 mm/ sqrt(18)
      Vc=0.002689;

      // inverse variance including prediction
      double InvV1=1./(Vc+Hc*(C*Hc_T));
      //double InvV1=1./(Vc+C.SandwichMultiply(H_T));
      if (InvV1<0.){
         if (DEBUG_LEVEL>0)
            _DBG_ << "Negative variance???" << endl;
         return NEGATIVE_VARIANCE;
      }
      // Compute KalmanSIMD gain matrix
      Kc=InvV1*(C*Hc_T);

      // Update the state vector with the target point
      // "Measurement" is average of expected beam spot size
      double res=0.1466666667-d;
      S+=res*Kc;  
      // Update state vector covariance matrix
      //C=C-K*(H*C);    
      C=C.SubSym(Kc*(Hc*C));

      // Update chi2 for this segment
      chisq+=(1.-Hc*Kc)*res*res/Vc;
      numdof++;
   }

   // Check that there were enough hits to make this a valid fit
   if (numdof<6){
      chisq=-1.;
      numdof=0;

      if (num_cdc==0){
         unsigned int new_index=(3*num_fdc)/4;
         break_point_fdc_index=(new_index>=MIN_HITS_FOR_REFIT)?new_index:(MIN_HITS_FOR_REFIT-1);
      }
      else{
         unsigned int new_index=num_fdc/2;
         if (new_index+num_cdc>=MIN_HITS_FOR_REFIT){
            break_point_fdc_index=new_index;
         }
         else{
            break_point_fdc_index=MIN_HITS_FOR_REFIT-num_cdc;
         }
      }
      return PRUNED_TOO_MANY_HITS;
   }

   //  chisq*=anneal_factor;
   numdof-=5;

   // Final positions in x and y for this leg
   x_=S(state_x);
   y_=S(state_y);

   if (DEBUG_LEVEL>1){
      cout << "Position after forward filter: " << x_ << ", " << y_ << ", " << z_ <<endl;
      cout << "Momentum " << 1./S(state_q_over_p) <<endl;
   }

   if (!S.IsFinite()) return FIT_FAILED;
  
   // Check if we have a kink in the track or threw away too many hits
   if (num_cdc>0 && break_point_fdc_index>0 && break_point_cdc_index>2){ 
      if (break_point_fdc_index+num_cdc<MIN_HITS_FOR_REFIT){
         //_DBG_ << endl;
         unsigned int new_index=num_fdc/2;
         if (new_index+num_cdc>=MIN_HITS_FOR_REFIT){
            break_point_fdc_index=new_index;
         }
         else{
            break_point_fdc_index=MIN_HITS_FOR_REFIT-num_cdc;
         }
      }
      return BREAK_POINT_FOUND;
   }
   if (num_cdc==0 && break_point_fdc_index>2){
      //_DBG_ << endl;
      if (break_point_fdc_index<num_fdc/2){
         break_point_fdc_index=(3*num_fdc)/4;
      }
      if (break_point_fdc_index<MIN_HITS_FOR_REFIT-1){
         break_point_fdc_index=MIN_HITS_FOR_REFIT-1;
      }
      return BREAK_POINT_FOUND;
   }
   if (num_cdc>5 && break_point_cdc_index>2){
      //_DBG_ << endl;  
      unsigned int new_index=num_fdc/2;
      if (new_index+num_cdc>=MIN_HITS_FOR_REFIT){
         break_point_fdc_index=new_index;
      }
      else{
         break_point_fdc_index=MIN_HITS_FOR_REFIT-num_cdc;
      }
      return BREAK_POINT_FOUND;
   }
   unsigned int num_good=0; 
   unsigned int num_hits=num_cdc+max_num_fdc_used_in_fit;
   for (unsigned int j=0;j<num_cdc;j++){
      if (cdc_used_in_fit[j]) num_good++;
   }
   for (unsigned int j=0;j<num_fdc;j++){
      if (fdc_used_in_fit[j]) num_good++;
   }
   if (double(num_good)/double(num_hits)<MINIMUM_HIT_FRACTION){
      //_DBG_ <<endl;
      if (num_cdc==0){
         unsigned int new_index=(3*num_fdc)/4;
         break_point_fdc_index=(new_index>=MIN_HITS_FOR_REFIT)?new_index:(MIN_HITS_FOR_REFIT-1);
      }
      else{
         unsigned int new_index=num_fdc/2;
         if (new_index+num_cdc>=MIN_HITS_FOR_REFIT){
            break_point_fdc_index=new_index;
         }
         else{
            break_point_fdc_index=MIN_HITS_FOR_REFIT-num_cdc;
         }
      }
      return PRUNED_TOO_MANY_HITS;
   }

   return FIT_SUCCEEDED;


}



// Kalman engine for forward tracks -- this routine adds CDC hits
kalman_error_t DTrackFitterKalmanSIMD::KalmanForwardCDC(double anneal,DMatrix5x1 &S, 
      DMatrix5x5 &C,double &chisq,
      unsigned int &numdof){
   DMatrix1x5 H;  // Track projection matrix
   DMatrix5x1 H_T; // Transpose of track projection matrix
   DMatrix5x5 J;  // State vector Jacobian matrix
   //DMatrix5x5 J_T; // transpose of this matrix
   DMatrix5x5 Q;  // Process noise covariance matrix
   DMatrix5x1 K;  // KalmanSIMD gain matrix
   DMatrix5x1 S0,S0_,Stest; //State vector
   DMatrix5x5 Ctest; // covariance matrix
   //DMatrix5x1 dS;  // perturbation in state vector
   double V=0.0507;

   // set used_in_fit flags to false for cdc hits
   unsigned int num_cdc=cdc_used_in_fit.size();
   for (unsigned int i=0;i<num_cdc;i++) cdc_used_in_fit[i]=false;
   for (unsigned int i=0;i<forward_traj.size();i++){
      forward_traj[i].h_id=0;
   }

   // initialize chi2 info
   chisq=0.;
   numdof=0;
   double var_cut=NUM_CDC_SIGMA_CUT*NUM_CDC_SIGMA_CUT;
   double my_anneal=anneal*anneal;
   double chi2cut=my_anneal*var_cut;

   // Save the starting values for C and S in the deque
   forward_traj[break_point_step_index].Skk=S;
   forward_traj[break_point_step_index].Ckk=C;

   // z-position
   double z=forward_traj[break_point_step_index].z;

   // Step size variables
   double step1=mStepSizeZ;
   double step2=mStepSizeZ;

   // wire information  
   unsigned int cdc_index=break_point_cdc_index;
   if (cdc_index<num_cdc-1){
     num_cdc=cdc_index+1;
   }

   if (cdc_index<MIN_HITS_FOR_REFIT) chi2cut=100.0;

   DVector2 origin=my_cdchits[cdc_index]->origin;
   double z0w=my_cdchits[cdc_index]->z0wire;
   DVector2 dir=my_cdchits[cdc_index]->dir;
   DVector2 wirepos=origin+(z-z0w)*dir;
   bool more_measurements=true;

   // doca variables
   double dx=S(state_x)-wirepos.X();
   double dy=S(state_y)-wirepos.Y();
   double doca2=0.,old_doca2=dx*dx+dy*dy;

   // loop over entries in the trajectory
   S0_=(forward_traj[break_point_step_index].S);
   for (unsigned int k=break_point_step_index+1;k<forward_traj.size()/*-1*/;k++){
      unsigned int k_minus_1=k-1;

      // Check that C matrix is positive definite
      if (!C.IsPosDef()){
         if (DEBUG_LEVEL>0) _DBG_ << "Broken covariance matrix!" <<endl;
         return BROKEN_COVARIANCE_MATRIX;
      }

      z=forward_traj[k].z;

      // Get the state vector, jacobian matrix, and multiple scattering matrix 
      // from reference trajectory
      S0=(forward_traj[k].S);
      J=(forward_traj[k].J);
      Q=(forward_traj[k].Q);

      // State S is perturbation about a seed S0
      //dS=S-S0_;
      /*
         dS.Print();
         J.Print();
         */

      // Update the actual state vector and covariance matrix
      S=S0+J*(S-S0_);

      // Bail if the position is grossly outside of the tracking volume
      if (S(state_x)*S(state_x)+S(state_y)*S(state_y)>R2_MAX){
         if (DEBUG_LEVEL>2)
         {
            _DBG_<< "Went outside of tracking volume at x=" << S(state_x)
               << " y=" << S(state_y) <<" z="<<z<<endl;
         }
         return POSITION_OUT_OF_RANGE;
      }
      // Bail if the momentum has dropped below some minimum
      if (fabs(S(state_q_over_p))>=Q_OVER_P_MAX){
         if (DEBUG_LEVEL>2)
         {
            _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p)) << endl;
         }
         return MOMENTUM_OUT_OF_RANGE;
      }



      //C=J*(C*J_T)+Q;   
      C=Q.AddSym(J*C*J.Transpose());
      //C=Q.AddSym(C.SandwichMultiply(J));

      // Save the current state of the reference trajectory
      S0_=S0;

      // new wire position
      wirepos=origin;
      wirepos+=(z-z0w)*dir;

      // new doca
      dx=S(state_x)-wirepos.X();
      dy=S(state_y)-wirepos.Y();
      doca2=dx*dx+dy*dy;

      // Save the current state and covariance matrix in the deque
      if (fit_type==kTimeBased){
         forward_traj[k].Skk=S;
         forward_traj[k].Ckk=C;
      }

      // Check if the doca is no longer decreasing
      if (more_measurements && doca2>old_doca2 && z<endplate_z){	
         if (my_cdchits[cdc_index]->status==good_hit){
            double dz=0.,newz=z;

            // energy loss 
            double dedx=0.;

            // Last 2 step sizes
            if (k>=2){
               double z1=forward_traj[k_minus_1].z;
               step1=-forward_traj[k].z+z1;
               step2=-z1+forward_traj[k-2].z;
            }

            // Track direction variables
            double tx=S(state_tx);
            double ty=S(state_ty);	
            double tanl=1./sqrt(tx*tx+ty*ty);
            double sinl=sin(atan(tanl));

            // Wire direction variables
            double ux=dir.X();
            double uy=dir.Y();
            // Variables relating wire direction and track direction
            double my_ux=tx-ux;
            double my_uy=ty-uy;
            double denom=my_ux*my_ux+my_uy*my_uy;

            // variables for dealing with propagation of S and C if we 
            // need to use Brent's algorithm to find the doca to the wire
            int num_steps=0;
            double dz3=0.;
            double my_dz=0.;

            // if the path length increment is small relative to the radius 
            // of curvature, use a linear approximation to find dz	
            bool do_brent=false;
            //printf("step1 %f step 2 %f \n",step1,step2);
            double two_step=step1+step2;
            if (fabs(qBr2p*S(state_q_over_p)
                     //*bfield->GetBz(S(state_x),S(state_y),z)
                     *forward_traj[k].B
                     *two_step/sinl)<0.05
                  && denom>EPS){
               double dzw=z-z0w;
               dz=-((S(state_x)-origin.X()-ux*dzw)*my_ux
                     +(S(state_y)-origin.Y()-uy*dzw)*my_uy)/denom;

               if (fabs(dz)>two_step || dz<0.0){
                  do_brent=true;
               }
               else{
                  newz=z+dz;
                  // Check for exiting the straw
                  if (newz>endplate_z){
                     dz=endplate_z-z;
                  }
               }
            }
            else{
               do_brent=true;
            }
            if (do_brent){
               if (CORRECT_FOR_ELOSS){
                  dedx=GetdEdx(S(state_q_over_p), 
                        forward_traj[k].K_rho_Z_over_A,
                        forward_traj[k].rho_Z_over_A,
                        forward_traj[k].LnI,forward_traj[k].Z);
               }

               // We have bracketed the minimum doca:  use Brent's agorithm
               if (BrentForward(z,dedx,z0w,origin,dir,S,dz)
                     !=NOERROR){
                  return MOMENTUM_OUT_OF_RANGE;
               }

               // Step the state and covariance through the field	  
               if (fabs(dz)>mStepSizeZ){
                  my_dz=(dz>0.0?1.0:-1.)*mStepSizeZ;
                  num_steps=int(fabs(dz/my_dz));
                  dz3=dz-num_steps*my_dz;

                  double my_z=z;
                  for (int m=0;m<num_steps;m++){
                     newz=my_z+my_dz;

                     // Step current state by my_dz
                     // Step(z,newz,dedx,S);

                     // propagate error matrix to z-position of hit
                     StepJacobian(my_z,newz,S0,dedx,J);
                     C=J*C*J.Transpose();
                     //C=C.SandwichMultiply(J);

                     // Step reference trajectory by my_dz
                     Step(my_z,newz,dedx,S0); 

                     my_z=newz;
                  }

                  newz=my_z+dz3;

                  // Step current state by dz3
                  //Step(my_z,newz,dedx,S);

                  // propagate error matrix to z-position of hit
                  StepJacobian(my_z,newz,S0,dedx,J);
                  C=J*C*J.Transpose();
                  //C=C.SandwichMultiply(J);	  

                  // Step reference trajectory by dz3
                  Step(my_z,newz,dedx,S0); 
               }
               else{
                  newz = z + dz;
                  // Step current state by dz
                  //Step(z,newz,dedx,S);

                  // propagate error matrix to z-position of hit
                  StepJacobian(z,newz,S0,dedx,J);
                  C=J*C*J.Transpose();
                  //C=C.SandwichMultiply(J);  

                  // Step reference trajectory by dz
                  Step(z,newz,dedx,S0); 
               }
            }

            // Wire position at current z
            wirepos=origin;
            wirepos+=(newz-z0w)*dir;

            double xw=wirepos.X();
            double yw=wirepos.Y();

            // predicted doca taking into account the orientation of the wire
            if (do_brent==false){
               // In this case we did not have to swim to find the doca and 
               // the transformation from the state vector space to the 
               // measurement space is straight-forward.
               dy=S(state_y)+S(state_ty)*dz-yw;
               dx=S(state_x)+S(state_tx)*dz-xw;      
            }
            else{
               // In this case we swam the state vector to the position of 
               // the doca
               dy=S(state_y)-yw;
               dx=S(state_x)-xw;
            }
            double cosstereo=my_cdchits[cdc_index]->cosstereo;
            double d=sqrt(dx*dx+dy*dy)*cosstereo+EPS;

            //printf("z %f d %f z-1 %f\n",newz,d,forward_traj[k_minus_1].z);

            // Track projection
            double cosstereo2_over_d=cosstereo*cosstereo/d;
            H_T(state_x)=dx*cosstereo2_over_d;
            H(state_x)=H_T(state_x);
            H_T(state_y)=dy*cosstereo2_over_d;
            H(state_y)=H_T(state_y);
            if (do_brent==false){
               H_T(state_ty)=H_T(state_y)*dz;
               H(state_ty)=H_T(state_ty);	  
               H_T(state_tx)=H_T(state_x)*dz;
               H(state_tx)=H_T(state_tx);
            }
            else{
               H_T(state_ty)=0.;
               H(state_ty)=0.;
               H_T(state_tx)=0.;
               H(state_tx)=0.;
            }

            //H.Print();

            // The next measurement
            double dm=0.39,tdrift=0.,tcorr=0.,dDdt0=0.;
            if (fit_type==kTimeBased || USE_PASS1_TIME_MODE){
               // Find offset of wire with respect to the center of the
               // straw at this z position
               const DCDCWire *mywire=my_cdchits[cdc_index]->hit->wire;
               int ring_index=mywire->ring-1;
               int straw_index=mywire->straw-1;
               double dz=newz-z0w;
               double phi_d=atan2(dy,dx);
               double delta
                  =max_sag[ring_index][straw_index]*(1.-dz*dz/5625.)
                  *cos(phi_d + sag_phi_offset[ring_index][straw_index]);
               double dphi=phi_d-mywire->origin.Phi();
               while (dphi>M_PI) dphi-=2*M_PI;
               while (dphi<-M_PI) dphi+=2*M_PI;
               if (mywire->origin.Y()<0) dphi*=-1.;

               tdrift=my_cdchits[cdc_index]->tdrift-mT0
                  -forward_traj[k_minus_1].t*TIME_UNIT_CONVERSION;
               double B=forward_traj[k_minus_1].B;
               ComputeCDCDrift(dphi,delta,tdrift,B,dm,V,tcorr);
               if (ALIGNMENT_FORWARD){
                  double myV=0.;
                  double mytcorr=0.;
                  double d_shifted;
                  double dt=2.0;
                  if (tdrift < 0.) d_shifted = dm;
                  else ComputeCDCDrift(dphi,delta,tdrift+dt,B,d_shifted,myV,mytcorr);
                  dDdt0=(d_shifted-dm)/dt;
               }
               //_DBG_ << tcorr << " " << dphi << " " << dm << endl;

            }
            // residual
            double res=dm-d;

            // inverse of variance including prediction
            //InvV=1./(V+H*(C*H_T));
            //double Vproj=C.SandwichMultiply(H_T);
	    double Vproj=H*C*H_T;
            double InvV=1./(V+Vproj);
            if (InvV<0.){
               if (DEBUG_LEVEL>0)
                  _DBG_ << "Negative variance???" << endl;
               break_point_cdc_index=(3*num_cdc)/4;
               return NEGATIVE_VARIANCE;
            }

            // Check how far this hit is from the expected position
            double chi2check=res*res*InvV;
            //if (sqrt(chi2check)>NUM_CDC_SIGMA_CUT) InvV*=0.8;
            if (chi2check<chi2cut){	  
               /*
                  if (chi2check>var_cut){
               // Give hits that satisfy the wide cut but are still pretty far
               // from the projected position less weight
               //_DBG_ << my_anneal << endl;

               // ad hoc correction 
               double diff = chi2check-var_cut;    
               V*=1.+my_anneal*diff*diff;
               InvV=1./(V+Vproj);
               }
               */

               // Compute KalmanSIMD gain matrix
               K=InvV*(C*H_T);

               // Update state vector covariance matrix
               Ctest=C.SubSym(K*(H*C));
               // Joseph form
               // C = (I-KH)C(I-KH)^T + KVK^T
               //Ctest=C.SandwichMultiply(I5x5-K*H)+V*MultiplyTranspose(K);
               // Check that Ctest is positive definite
               if (!Ctest.IsPosDef()){
                  if (DEBUG_LEVEL>0) _DBG_ << "Broken covariance matrix!" <<endl;
                  return BROKEN_COVARIANCE_MATRIX;
               }

               bool skip_ring
                  =(my_cdchits[cdc_index]->hit->wire->ring==RING_TO_SKIP);
               // update covariance matrix and state vector
               if (skip_ring==false && tdrift >= CDC_T_DRIFT_MIN){
                  C=Ctest;		       
                  S+=res*K;
               }
               // Mark point on ref trajectory with a hit id for the straw
               forward_traj[k].h_id=cdc_index+1000;

               // Store some updated values related to the hit
               double scale=(skip_ring)?1.:(1.-H*K);
               cdc_updates[cdc_index].tcorr=tcorr;
               cdc_updates[cdc_index].tdrift=tdrift;
               cdc_updates[cdc_index].doca=dm;
               cdc_updates[cdc_index].variance=V;
               cdc_updates[cdc_index].dDdt0=dDdt0;
               cdc_used_in_fit[cdc_index]=true;
               if(tdrift < CDC_T_DRIFT_MIN) cdc_used_in_fit[cdc_index]=false;

               // Update chi2 for this segment
               if (skip_ring==false && tdrift >= CDC_T_DRIFT_MIN){
                  chisq+=scale*res*res/V;
                  numdof++;	
               }
               break_point_cdc_index=cdc_index;
               break_point_step_index=k_minus_1;

               if (DEBUG_LEVEL>9)
                  printf("Ring %d straw %d pred %f meas %f chi2 %f useBrent %i \n",
                        my_cdchits[cdc_index]->hit->wire->ring,
                        my_cdchits[cdc_index]->hit->wire->straw,
                        d,dm,(1.-H*K)*res*res/V,do_brent);

            }

            // If we had to use Brent's algorithm to find the true doca, 
            // we need to swim the state vector and covariance matrix back 
            // to the appropriate position along the reference trajectory.
            if (do_brent){
               if (num_steps==0){
                  // Step C back to the z-position on the reference trajectory
                  StepJacobian(newz,z,S0,dedx,J);
                  C=J*C*J.Transpose();
                  //C=C.SandwichMultiply(J);

                  // Step S to current position on the reference trajectory=
                  Step(newz,z,dedx,S);

                  // Step S0 to current position on the reference trajectory=
                  Step(newz,z,dedx,S0);
               }
               else{
                  double my_z=newz;
                  for (int m=0;m<num_steps;m++){
                     z=my_z-my_dz;

                     // Step C along z
                     StepJacobian(my_z,z,S0,dedx,J);
                     C=J*C*J.Transpose();
                     //C=C.SandwichMultiply(J);

                     // Step S along z
                     Step(my_z,z,dedx,S);

                     // Step S0 along z
                     Step(my_z,z,dedx,S0);

                     my_z=z;
                  }
                  z=my_z-dz3;

                  // Step C back to the z-position on the reference trajectory
                  StepJacobian(my_z,z,S0,dedx,J);
                  C=J*C*J.Transpose();
                  //C=C.SandwichMultiply(J);

                  // Step S to current position on the reference trajectory
                  Step(my_z,z,dedx,S);

                  // Step S0 to current position on the reference trajectory
                  Step(my_z,z,dedx,S0);

               }
            }
            cdc_updates[cdc_index].S=S;
            cdc_updates[cdc_index].C=C;	
         }
         else {
            if (cdc_index>0) cdc_index--;
            else cdc_index=0;

         }

         // new wire origin and direction
         if (cdc_index>0){
            cdc_index--;
            origin=my_cdchits[cdc_index]->origin;
            z0w=my_cdchits[cdc_index]->z0wire;
            dir=my_cdchits[cdc_index]->dir;
         }
         else{
            more_measurements=false;
         }

         // Update the wire position
         wirepos=origin;
         wirepos+=(z-z0w)*dir;

         // new doca
         dx=S(state_x)-wirepos.X();
         dy=S(state_y)-wirepos.Y();
         doca2=dx*dx+dy*dy;
      }
      old_doca2=doca2;

   }

   // Check that there were enough hits to make this a valid fit
   if (numdof<6){
      chisq=-1.;
      numdof=0;
      return PRUNED_TOO_MANY_HITS;
   }
   numdof-=5;

   // Final position for this leg
   x_=S(state_x);
   y_=S(state_y);
   z_=forward_traj[forward_traj.size()-1].z;

   if (!S.IsFinite()) return FIT_FAILED;

   // Check if the momentum is unphysically large
   if (1./fabs(S(state_q_over_p))>12.0){
      if (DEBUG_LEVEL>2)
      {
         _DBG_ << "Unphysical momentum: P = " << 1./fabs(S(state_q_over_p))
            <<endl;
      }
      return MOMENTUM_OUT_OF_RANGE;
   }

   // Check if we have a kink in the track or threw away too many cdc hits
   if (num_cdc>=MIN_HITS_FOR_REFIT){
      if (break_point_cdc_index>1){
         if (break_point_cdc_index<num_cdc/2){
            break_point_cdc_index=(3*num_cdc)/4;
         }
         return BREAK_POINT_FOUND;
      }

      unsigned int num_good=0; 
      for (unsigned int j=0;j<num_cdc;j++){
         if (cdc_used_in_fit[j]) num_good++;
      }
      if (double(num_good)/double(num_cdc)<MINIMUM_HIT_FRACTION){
         return PRUNED_TOO_MANY_HITS;
      }
   }

   return FIT_SUCCEEDED;
}

// Extrapolate to the point along z of closest approach to the beam line using 
// the forward track state vector parameterization.  Converts to the central
// track representation at the end.
jerror_t DTrackFitterKalmanSIMD::ExtrapolateToVertex(DMatrix5x1 &S,
      DMatrix5x5 &C){
   DMatrix5x5 J;  // Jacobian matrix
   DMatrix5x5 Q;  // multiple scattering matrix
   DMatrix5x1 S1(S);  // copy of S

   // position variables
   double z=z_,newz=z_;

   DVector2 beam_pos=beam_center+(z-beam_z0)*beam_dir;
   double r2_old=(DVector2(S(state_x),S(state_y))-beam_pos).Mod2();
   double dz_old=0.;
   double dEdx=0.;
   double sign=1.;

   // material properties
   double rho_Z_over_A=0.,LnI=0.,K_rho_Z_over_A=0.,Z=0.;
   double chi2c_factor=0.,chi2a_factor=0.,chi2a_corr=0.;

   //  if (fit_type==kTimeBased)printf("------Extrapolating\n");

   //  printf("-----------\n");
   // Current position
   DVector3 pos(S(state_x),S(state_y),z);

   // get material properties from the Root Geometry
   if (geom->FindMatKalman(pos,K_rho_Z_over_A,rho_Z_over_A,LnI,Z,
            chi2c_factor,chi2a_factor,chi2a_corr,
            last_material_map)
         !=NOERROR){
      if (DEBUG_LEVEL>1)
         _DBG_ << "Material error in ExtrapolateToVertex at (x,y,z)=("
            << pos.X() <<"," << pos.y()<<","<< pos.z()<<")"<< endl;
      return UNRECOVERABLE_ERROR;
   }

   // Adjust the step size
   double ds_dz=sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));
   double dz=-mStepSizeS/ds_dz;
   if (fabs(dEdx)>EPS){      
      dz=(-1.)*DE_PER_STEP/fabs(dEdx)/ds_dz;
   }
   if(fabs(dz)>mStepSizeZ) dz=-mStepSizeZ;
   if(fabs(dz)<MIN_STEP_SIZE)dz=-MIN_STEP_SIZE;

   // Get dEdx for the upcoming step
   if (CORRECT_FOR_ELOSS){
      dEdx=GetdEdx(S(state_q_over_p),K_rho_Z_over_A,rho_Z_over_A,LnI,Z); 
   }


   double ztest=endplate_z;
   if (my_fdchits.size()>0){ 
      ztest =my_fdchits[0]->z-1.;
   }  
   if (z<ztest)
   {
      // Check direction of propagation	
      DMatrix5x1 S2(S); // second copy of S

      // Step test states through the field and compute squared radii
      Step(z,z-dz,dEdx,S1);	
      // Bail if the momentum has dropped below some minimum
      if (fabs(S1(state_q_over_p))>Q_OVER_P_MAX){
         if (DEBUG_LEVEL>2)
         {
            _DBG_ << "Bailing: P = " << 1./fabs(S1(state_q_over_p))
               << endl;
         }
         return UNRECOVERABLE_ERROR;
      } 
      beam_pos=beam_center+(z-dz-beam_z0)*beam_dir;
      double r2minus=(DVector2(S1(state_x),S1(state_y))-beam_pos).Mod2();

      Step(z,z+dz,dEdx,S2);	
      // Bail if the momentum has dropped below some minimum
      if (fabs(S2(state_q_over_p))>Q_OVER_P_MAX){
         if (DEBUG_LEVEL>2)
         {
            _DBG_ << "Bailing: P = " << 1./fabs(S2(state_q_over_p))
               << endl;
         }
         return UNRECOVERABLE_ERROR;
      }
      beam_pos=beam_center+(z+dz-beam_z0)*beam_dir;
      double r2plus=(DVector2(S2(state_x),S2(state_y))-beam_pos).Mod2();
      // Check to see if we have already bracketed the minimum
      if (r2plus>r2_old && r2minus>r2_old){
         newz=z+dz;  
         double dz2=0.;
         if (BrentsAlgorithm(newz,dz,dEdx,0.,beam_center,beam_dir,S2,dz2)!=NOERROR){
            if (DEBUG_LEVEL>2)
            {
               _DBG_ << "Bailing: P = " << 1./fabs(S2(state_q_over_p))
                  << endl;
            }
            return UNRECOVERABLE_ERROR;
         }
         z_=newz+dz2;

         // Compute the Jacobian matrix
         StepJacobian(z,z_,S,dEdx,J);  

         // Propagate the covariance matrix
         C=J*C*J.Transpose();
         //C=C.SandwichMultiply(J);

         // Step to the position of the doca
         Step(z,z_,dEdx,S);

         // update internal variables
         x_=S(state_x);
         y_=S(state_y);

         return NOERROR;
      }

      // Find direction to propagate toward minimum doca
      if (r2minus<r2_old && r2_old<r2plus){ 
         newz=z-dz;

         // Compute the Jacobian matrix
         StepJacobian(z,newz,S,dEdx,J);  

         // Propagate the covariance matrix
         C=J*C*J.Transpose();
         //C=C.SandwichMultiply(J);

         S2=S;
         S=S1;
         S1=S2;
         dz*=-1.;
         sign=-1.;
         dz_old=dz;
         r2_old=r2minus;
         z=z+dz;
      }
      if (r2minus>r2_old && r2_old>r2plus){
         newz=z+dz;

         // Compute the Jacobian matrix
         StepJacobian(z,newz,S,dEdx,J);  

         // Propagate the covariance matrix
         C=J*C*J.Transpose();
         //C=C.SandwichMultiply(J);

         S1=S;
         S=S2;
         dz_old=dz;
         r2_old=r2plus;
         z=z+dz;
      }
   }

   double r2=r2_old;
   while (z>Z_MIN && r2<R2_MAX && z<ztest && r2>EPS){   
      // Bail if the momentum has dropped below some minimum
      if (fabs(S(state_q_over_p))>Q_OVER_P_MAX){
         if (DEBUG_LEVEL>2)
         {
            _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p))
               << endl;
         }
         return UNRECOVERABLE_ERROR;
      }

      // Relationship between arc length and z
      double dz_ds=1./sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));

      // get material properties from the Root Geometry
      pos.SetXYZ(S(state_x),S(state_y),z);
      double s_to_boundary=1.e6;
      if (ENABLE_BOUNDARY_CHECK && fit_type==kTimeBased){
         DVector3 mom(S(state_tx),S(state_ty),1.);
         if (geom->FindMatKalman(pos,mom,K_rho_Z_over_A,rho_Z_over_A,LnI,Z,
                  chi2c_factor,chi2a_factor,chi2a_corr,
                  last_material_map,&s_to_boundary)
               !=NOERROR){
            _DBG_ << "Material error in ExtrapolateToVertex! " << endl;
            return UNRECOVERABLE_ERROR;      
         }  
      }
      else{
         if (geom->FindMatKalman(pos,K_rho_Z_over_A,rho_Z_over_A,LnI,Z,
                  chi2c_factor,chi2a_factor,chi2a_corr,
                  last_material_map)
               !=NOERROR){
            _DBG_ << "Material error in ExtrapolateToVertex! " << endl;
            break;
         }
      }

      // Get dEdx for the upcoming step
      if (CORRECT_FOR_ELOSS){
         dEdx=GetdEdx(S(state_q_over_p),K_rho_Z_over_A,rho_Z_over_A,LnI,Z); 
      }

      // Adjust the step size 
      //dz=-sign*mStepSizeS*dz_ds;
      double ds=mStepSizeS;
      if (fabs(dEdx)>EPS){     
         ds=DE_PER_STEP/fabs(dEdx);
      }
      /*
         if(fabs(dz)>mStepSizeZ) dz=-sign*mStepSizeZ;
         if (fabs(dz)<z_to_boundary) dz=-sign*z_to_boundary;
         if(fabs(dz)<MIN_STEP_SIZE)dz=-sign*MIN_STEP_SIZE;
         */
      if (ds>mStepSizeS) ds=mStepSizeS;
      if (ds>s_to_boundary) ds=s_to_boundary;
      if (ds<MIN_STEP_SIZE) ds=MIN_STEP_SIZE;
      dz=-sign*ds*dz_ds;

      // Get the contribution to the covariance matrix due to multiple 
      // scattering
      GetProcessNoise(z,ds,chi2c_factor,chi2a_factor,chi2a_corr,S,Q);

      if (CORRECT_FOR_ELOSS){
         double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
         double one_over_beta2=1.+mass2*q_over_p_sq;
         double varE=GetEnergyVariance(ds,one_over_beta2,K_rho_Z_over_A);
         Q(state_q_over_p,state_q_over_p)=varE*q_over_p_sq*q_over_p_sq*one_over_beta2;
      }


      newz=z+dz;
      // Compute the Jacobian matrix
      StepJacobian(z,newz,S,dEdx,J);  

      // Propagate the covariance matrix
      C=Q.AddSym(J*C*J.Transpose());
      //C=Q.AddSym(C.SandwichMultiply(J));

      // Step through field
      Step(z,newz,dEdx,S);

      // Check if we passed the minimum doca to the beam line
      beam_pos=beam_center+(newz-beam_z0)*beam_dir;
      r2=(DVector2(S(state_x),S(state_y))-beam_pos).Mod2();
      //r2=S(state_x)*S(state_x)+S(state_y)*S(state_y);
      if (r2>r2_old){
         double two_step=dz+dz_old;

         // Find the increment/decrement in z to get to the minimum doca to the
         // beam line   
         S1=S;
         if (BrentsAlgorithm(newz,0.5*two_step,dEdx,0.,beam_center,beam_dir,S,dz)!=NOERROR){
            //_DBG_<<endl;
            return UNRECOVERABLE_ERROR;
         }

         // Compute the Jacobian matrix
         z_=newz+dz;
         StepJacobian(newz,z_,S1,dEdx,J);

         // Propagate the covariance matrix
         //C=J*C*J.Transpose()+(dz/(newz-z))*Q;
         //C=((dz/newz-z)*Q).AddSym(C.SandwichMultiply(J));
         //C=C.SandwichMultiply(J);
	 C=J*C*J.Transpose();

         // update internal variables
         x_=S(state_x);
         y_=S(state_y);

         return NOERROR;
      }
      r2_old=r2;
      dz_old=dz;
      S1=S;
      z=newz;
   }
   // update internal variables
   x_=S(state_x);
   y_=S(state_y);
   z_=newz;

   return NOERROR;
}


// Extrapolate to the point along z of closest approach to the beam line using 
// the forward track state vector parameterization.
jerror_t DTrackFitterKalmanSIMD::ExtrapolateToVertex(DMatrix5x1 &S){
   DMatrix5x5 J;  // Jacobian matrix
   DMatrix5x1 S1(S);  // copy of S

   // position variables
   double z=z_,newz=z_;
   DVector2 beam_pos=beam_center+(z-beam_z0)*beam_dir;
   double r2_old=(DVector2(S(state_x),S(state_y))-beam_pos).Mod2();
   double dz_old=0.;
   double dEdx=0.;

   // Direction and origin for beam line
   DVector2 dir;
   DVector2 origin;

   // material properties
   double rho_Z_over_A=0.,LnI=0.,K_rho_Z_over_A=0.,Z=0.;
   double chi2c_factor=0.,chi2a_factor=0.,chi2a_corr=0.;
   DVector3 pos;  // current position along trajectory

   double r2=r2_old;
   while (z>Z_MIN && r2<R2_MAX && z<Z_MAX && r2>EPS){
      // Bail if the momentum has dropped below some minimum
      if (fabs(S(state_q_over_p))>Q_OVER_P_MAX){
         if (DEBUG_LEVEL>2)
         {
            _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p))
               << endl;
         }
         return UNRECOVERABLE_ERROR;
      }

      // Relationship between arc length and z
      double dz_ds=1./sqrt(1.+S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty));

      // get material properties from the Root Geometry
      pos.SetXYZ(S(state_x),S(state_y),z);
      if (geom->FindMatKalman(pos,K_rho_Z_over_A,rho_Z_over_A,LnI,Z,
               chi2c_factor,chi2a_factor,chi2a_corr,
               last_material_map)
            !=NOERROR){
         _DBG_ << "Material error in ExtrapolateToVertex! " << endl;
         break;
      }

      // Get dEdx for the upcoming step
      if (CORRECT_FOR_ELOSS){
         dEdx=GetdEdx(S(state_q_over_p),K_rho_Z_over_A,rho_Z_over_A,LnI,Z); 
      }

      // Adjust the step size 
      double ds=mStepSizeS;
      if (fabs(dEdx)>EPS){     
         ds=DE_PER_STEP/fabs(dEdx);
      }
      if (ds>mStepSizeS) ds=mStepSizeS;
      if (ds<MIN_STEP_SIZE) ds=MIN_STEP_SIZE;
      double dz=-ds*dz_ds;

      newz=z+dz;


      // Step through field
      Step(z,newz,dEdx,S);

      // Check if we passed the minimum doca to the beam line
      beam_pos=beam_center+(newz-beam_z0)*beam_dir;
      r2=(DVector2(S(state_x),S(state_y))-beam_pos).Mod2();

      if (r2>r2_old && newz<endplate_z){
         double two_step=dz+dz_old;

         // Find the increment/decrement in z to get to the minimum doca to the
         // beam line   
         if (BrentsAlgorithm(newz,0.5*two_step,dEdx,0.,beam_center,beam_dir,S,dz)!=NOERROR){
            return UNRECOVERABLE_ERROR;
         }
         // update internal variables
         x_=S(state_x);
         y_=S(state_y); 
         z_=newz+dz;

         return NOERROR;
      }
      r2_old=r2;
      dz_old=dz;
      S1=S;
      z=newz;
   }
   // update internal variables
   x_=S(state_x);
   y_=S(state_y);
   z_=newz;


   return NOERROR;
}




// Propagate track to point of distance of closest approach to origin
jerror_t DTrackFitterKalmanSIMD::ExtrapolateToVertex(DVector2 &xy,
      DMatrix5x1 &Sc,DMatrix5x5 &Cc){
   DMatrix5x5 Jc=I5x5;  //Jacobian matrix
   DMatrix5x5 Q; // multiple scattering matrix

   // Position and step variables
   DVector2 beam_pos=beam_center+(Sc(state_z)-beam_z0)*beam_dir;
   double r2=(xy-beam_pos).Mod2();
   double ds=-mStepSizeS; // step along path in cm
   double r2_old=r2;

   // Energy loss
   double dedx=0.;

   // Check direction of propagation
   DMatrix5x1 S0;
   S0=Sc; 
   DVector2 xy0=xy;
   DVector2 xy1=xy;
   Step(xy0,ds,S0,dedx);
   // Bail if the transverse momentum has dropped below some minimum
   if (fabs(S0(state_q_over_pt))>Q_OVER_PT_MAX){
      if (DEBUG_LEVEL>2)
      {
         _DBG_ << "Bailing: PT = " << 1./fabs(S0(state_q_over_pt))
            << endl;
      }
      return UNRECOVERABLE_ERROR;
   }
   beam_pos=beam_center+(S0(state_z)-beam_z0)*beam_dir;
   r2=(xy0-beam_pos).Mod2();
   if (r2>r2_old) ds*=-1.;
   double ds_old=ds;

   //  if (fit_type==kTimeBased)printf("------Extrapolating\n");

   if (central_traj.size()==0){
      if (DEBUG_LEVEL>1) _DBG_ << "Central trajectory size==0!" << endl;
      return UNRECOVERABLE_ERROR;
   }

   // D is now on the actual trajectory itself
   Sc(state_D)=0.;

   // Track propagation loop
   while (Sc(state_z)>Z_MIN && Sc(state_z)<Z_MAX  
         && r2<R2_MAX){  
      // Bail if the transverse momentum has dropped below some minimum
      if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
         if (DEBUG_LEVEL>2)
         {
            _DBG_ << "Bailing: PT = " << 1./fabs(Sc(state_q_over_pt))
               << endl;
         }
         return UNRECOVERABLE_ERROR;
      }

      // get material properties from the Root Geometry
      double rho_Z_over_A=0.,LnI=0.,K_rho_Z_over_A=0.,Z=0.;
      double chi2c_factor=0.,chi2a_factor=0.,chi2a_corr=0.;
      DVector3 pos3d(xy.X(),xy.Y(),Sc(state_z));
      if (geom->FindMatKalman(pos3d,K_rho_Z_over_A,rho_Z_over_A,LnI,Z,
               chi2c_factor,chi2a_factor,chi2a_corr,
               last_material_map)
            !=NOERROR){
         _DBG_ << "Material error in ExtrapolateToVertex! " << endl;
         break;
      }

      // Get dEdx for the upcoming step
      double q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
      if (CORRECT_FOR_ELOSS){
         dedx=-GetdEdx(q_over_p,K_rho_Z_over_A,rho_Z_over_A,LnI,Z); 
      }
      // Adjust the step size
      double sign=(ds>0.0)?1.:-1.;
      if (fabs(dedx)>EPS){
         ds=sign*DE_PER_STEP/fabs(dedx);
      }
      if(fabs(ds)>mStepSizeS) ds=sign*mStepSizeS;
      if(fabs(ds)<MIN_STEP_SIZE)ds=sign*MIN_STEP_SIZE;

      // Multiple scattering
      GetProcessNoiseCentral(ds,chi2c_factor,chi2a_factor,chi2a_corr,Sc,Q);

      if (CORRECT_FOR_ELOSS){
         double q_over_p_sq=q_over_p*q_over_p;
         double one_over_beta2=1.+mass2*q_over_p*q_over_p;
         double varE=GetEnergyVariance(ds,one_over_beta2,K_rho_Z_over_A);
	 Q(state_q_over_pt,state_q_over_pt)
	   +=varE*Sc(state_q_over_pt)*Sc(state_q_over_pt)*one_over_beta2
	   *q_over_p_sq;
      }

      // Propagate the state and covariance through the field
      S0=Sc;
      DVector2 old_xy=xy;  
      StepStateAndCovariance(xy,ds,dedx,Sc,Jc,Cc);

      // Add contribution due to multiple scattering
      Cc=(sign*Q).AddSym(Cc);

      beam_pos=beam_center+(Sc(state_z)-beam_z0)*beam_dir;
      r2=(xy-beam_pos).Mod2();
      //printf("r %f r_old %f \n",sqrt(r2),sqrt(r2_old));
      if (r2>r2_old) {
         // We've passed the true minimum; backtrack to find the "vertex" 
         // position
         double my_ds=0.;
	 if (BrentsAlgorithm(ds,ds_old,dedx,xy,0.,beam_center,beam_dir,Sc,my_ds)!=NOERROR){
               //_DBG_ <<endl;
	   return UNRECOVERABLE_ERROR;
	 }
	 //printf ("Brent min r %f\n",xy.Mod());
         
         // Find the field and gradient at (old_x,old_y,old_z)
         bfield->GetFieldAndGradient(old_xy.X(),old_xy.Y(),S0(state_z),Bx,By,Bz,
               dBxdx,dBxdy,dBxdz,dBydx,
               dBydy,dBydz,dBzdx,dBzdy,dBzdz);

         // Compute the Jacobian matrix
         my_ds-=ds_old;
         StepJacobian(old_xy,xy-old_xy,my_ds,S0,dedx,Jc);

         // Propagate the covariance matrix
         //Cc=Jc*Cc*Jc.Transpose()+(my_ds/ds_old)*Q;
         //Cc=((my_ds/ds_old)*Q).AddSym(Cc.SandwichMultiply(Jc));
	 Cc=((my_ds/ds_old)*Q).AddSym(Jc*Cc*Jc.Transpose());

         break;
      }
      r2_old=r2;
      ds_old=ds;
   }   

   return NOERROR;
}

// Propagate track to point of distance of closest approach to origin
jerror_t DTrackFitterKalmanSIMD::ExtrapolateToVertex(DVector2 &xy,
      DMatrix5x1 &Sc){

   // Initialize the beam position = center of target, and the direction
   DVector2 origin;  
   DVector2 dir;

   // Position and step variables
   DVector2 beam_pos=beam_center+(Sc(state_z)-beam_z0)*beam_dir;
   double r2=(xy-beam_pos).Mod2();
   double ds=-mStepSizeS; // step along path in cm
   double r2_old=r2;

   // Energy loss
   double dedx=0.;

   // Check direction of propagation
   DMatrix5x1 S0;
   S0=Sc; 
   DVector2 xy0=xy;
   DVector2 xy1=xy;
   Step(xy0,ds,S0,dedx);
   beam_pos=beam_center+(S0(state_z)-beam_z0)*beam_dir;
   r2=(xy0-beam_pos).Mod2();
   if (r2>r2_old) ds*=-1.;
   double ds_old=ds;

   // Track propagation loop
   while (Sc(state_z)>Z_MIN && Sc(state_z)<Z_MAX  
         && r2<R2_MAX){  
      // get material properties from the Root Geometry
      double rho_Z_over_A=0.,LnI=0.,K_rho_Z_over_A=0.,Z=0;
      double chi2c_factor=0.,chi2a_factor=0.,chi2a_corr=0.;
      DVector3 pos3d(xy.X(),xy.Y(),Sc(state_z));
      if (geom->FindMatKalman(pos3d,K_rho_Z_over_A,rho_Z_over_A,LnI,Z,
               chi2c_factor,chi2a_factor,chi2a_corr,
               last_material_map)
            !=NOERROR){
         _DBG_ << "Material error in ExtrapolateToVertex! " << endl;
         break;
      }

      // Get dEdx for the upcoming step
      double q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
      if (CORRECT_FOR_ELOSS){
         dedx=GetdEdx(q_over_p,K_rho_Z_over_A,rho_Z_over_A,LnI,Z); 
      }
      // Adjust the step size
      double sign=(ds>0.0)?1.:-1.;
      if (fabs(dedx)>EPS){
         ds=sign*DE_PER_STEP/fabs(dedx);
      }
      if(fabs(ds)>mStepSizeS) ds=sign*mStepSizeS;
      if(fabs(ds)<MIN_STEP_SIZE)ds=sign*MIN_STEP_SIZE;

      // Propagate the state through the field
      Step(xy,ds,Sc,dedx);

      beam_pos=beam_center+(Sc(state_z)-beam_z0)*beam_dir;
      r2=(xy-beam_pos).Mod2();
      //printf("r %f r_old %f \n",r,r_old);
      if (r2>r2_old) {
         // We've passed the true minimum; backtrack to find the "vertex" 
         // position
         double my_ds=0.;
	 BrentsAlgorithm(ds,ds_old,dedx,xy,0.,beam_center,beam_dir,Sc,my_ds);
	 //printf ("Brent min r %f\n",pos.Perp());
         break;
      }
      r2_old=r2;
      ds_old=ds;
   }   

   return NOERROR;
}




// Transform the 5x5 tracking error matrix into a 7x7 error matrix in cartesian
// coordinates
shared_ptr<TMatrixFSym> DTrackFitterKalmanSIMD::Get7x7ErrorMatrixForward(DMatrixDSym C){
   auto C7x7 = dResourcePool_TMatrixFSym->Get_SharedResource();
   C7x7->ResizeTo(7, 7);
   DMatrix J(7,5);

   double p=1./fabs(q_over_p_);
   double tanl=1./sqrt(tx_*tx_+ty_*ty_);
   double tanl2=tanl*tanl;
   double lambda=atan(tanl);
   double sinl=sin(lambda);
   double sinl3=sinl*sinl*sinl;

   J(state_X,state_x)=J(state_Y,state_y)=1.;
   J(state_Pz,state_ty)=-p*ty_*sinl3;
   J(state_Pz,state_tx)=-p*tx_*sinl3;
   J(state_Px,state_ty)=J(state_Py,state_tx)=-p*tx_*ty_*sinl3;
   J(state_Px,state_tx)=p*(1.+ty_*ty_)*sinl3;
   J(state_Py,state_ty)=p*(1.+tx_*tx_)*sinl3;
   J(state_Pz,state_q_over_p)=-p*sinl/q_over_p_;
   J(state_Px,state_q_over_p)=tx_*J(state_Pz,state_q_over_p);
   J(state_Py,state_q_over_p)=ty_*J(state_Pz,state_q_over_p); 
   J(state_Z,state_x)=tx_*tanl2;
   J(state_Z,state_y)=ty_*tanl2;
   double diff=tx_*tx_-ty_*ty_;
   double frac=tanl2*tanl2;
   J(state_Z,state_tx)=-(x_*diff+2.*tx_*ty_*y_)*frac;
   J(state_Z,state_ty)=-(2.*tx_*ty_*x_-y_*diff)*frac;

   // C'= JCJ^T
   *C7x7=C.Similarity(J);

   return C7x7;
}



// Transform the 5x5 tracking error matrix into a 7x7 error matrix in cartesian
// coordinates
shared_ptr<TMatrixFSym> DTrackFitterKalmanSIMD::Get7x7ErrorMatrix(DMatrixDSym C){
   auto C7x7 = dResourcePool_TMatrixFSym->Get_SharedResource();
   C7x7->ResizeTo(7, 7);
   DMatrix J(7,5);
   //double cosl=cos(atan(tanl_));
   double pt=1./fabs(q_over_pt_);
   //double p=pt/cosl;
   // double p_sq=p*p;
   //  double E=sqrt(mass2+p_sq);
   double pt_sq=1./(q_over_pt_*q_over_pt_);
   double cosphi=cos(phi_);
   double sinphi=sin(phi_);
   double q=(q_over_pt_>0.0)?1.:-1.;

   J(state_Px,state_q_over_pt)=-q*pt_sq*cosphi;
   J(state_Px,state_phi)=-pt*sinphi;

   J(state_Py,state_q_over_pt)=-q*pt_sq*sinphi;
   J(state_Py,state_phi)=pt*cosphi;

   J(state_Pz,state_q_over_pt)=-q*pt_sq*tanl_;
   J(state_Pz,state_tanl)=pt;

   J(state_X,state_phi)=-D_*cosphi;
   J(state_X,state_D)=-sinphi;

   J(state_Y,state_phi)=-D_*sinphi;
   J(state_Y,state_D)=cosphi;

   J(state_Z,state_z)=1.;

   // C'= JCJ^T
   *C7x7=C.Similarity(J);

   return C7x7;
}

// Track recovery for Central tracks
//-----------------------------------
// This code attempts to recover tracks that are "broken".  Sometimes the fit fails because too many hits were pruned 
// such that the number of surviving hits is less than the minimum number of degrees of freedom for a five-parameter fit.
// This condition is flagged as PRUNED_TOO_MANY_HITS.  It may also be the case that even if we used enough hits for the fit to
// be valid (i.e., make sense in terms of the number of degrees of freedom for the number of parameters (5) we are trying 
// to determine), we throw away a number of hits near the target because the projected trajectory deviates too far away from 
// the actual hits.  This may be an indication of a kink in the trajectory.  This condition is flagged as BREAK_POINT_FOUND.
// Sometimes the fit may succeed and no break point is found, but the pruning threw out a large number of hits.  This 
// condition is flagged as PRUNED_TOO_MANY_HITS.  The recovery code proceeds by truncating the reference trajectory to 
// to a point just after the hit where a break point was found and all hits are ignored beyond this point subject to a 
// minimum number of hits set by MIN_HITS_FOR_REFIT.  The recovery code always attempts to use the hits closest to the 
// target.  The code is allowed to iterate; with each iteration the trajectory and list of useable hits is further truncated.
kalman_error_t DTrackFitterKalmanSIMD::RecoverBrokenTracks(double anneal_factor, 
      DMatrix5x1 &S, 
      DMatrix5x5 &C,
      const DMatrix5x5 &C0,
      DVector2 &xy,
      double &chisq, 
      unsigned int &numdof){
   if (DEBUG_LEVEL>1) _DBG_  << "Attempting to recover broken track ... " <<endl;

   // Initialize degrees of freedom and chi^2
   double refit_chisq=-1.;
   unsigned int refit_ndf=0;
   // State vector and covariance matrix
   DMatrix5x5 C1;
   DMatrix5x1 S1;
   // Position vector
   DVector2 refit_xy;

   // save the status of the hits used in the fit
   unsigned int num_hits=cdc_used_in_fit.size();
   vector<bool>old_cdc_used_status(num_hits);  
   for (unsigned int j=0;j<num_hits;j++){
      old_cdc_used_status[j]=cdc_used_in_fit[j];
   }

   // Truncate the reference trajectory to just beyond the break point (or the minimum number of hits)
   unsigned int min_cdc_index_for_refit=MIN_HITS_FOR_REFIT-1;
   //if (break_point_cdc_index<num_hits/2)
   //  break_point_cdc_index=num_hits/2;
   if (break_point_cdc_index<min_cdc_index_for_refit){
      break_point_cdc_index=min_cdc_index_for_refit;
   }
   // Next determine where we need to truncate the trajectory  
   DVector2 origin=my_cdchits[break_point_cdc_index]->origin;
   DVector2 dir=my_cdchits[break_point_cdc_index]->dir;
   double z0=my_cdchits[break_point_cdc_index]->z0wire;
   unsigned int k=0;
   for (k=central_traj.size()-1;k>1;k--){
      double r2=central_traj[k].xy.Mod2();
      double r2next=central_traj[k-1].xy.Mod2();
      double r2_cdc=(origin+(central_traj[k].S(state_z)-z0)*dir).Mod2();
      if (r2next>r2 && r2>r2_cdc) break;
   }
   break_point_step_index=k;

   if (break_point_step_index==central_traj.size()-1){
      if (DEBUG_LEVEL>0) _DBG_ << "Invalid reference trajectory in track recovery..." << endl;
      return FIT_FAILED;
   }

   kalman_error_t refit_error=FIT_NOT_DONE;
   unsigned int old_cdc_index=break_point_cdc_index;
   unsigned int old_step_index=break_point_step_index;
   unsigned int refit_iter=0;
   unsigned int num_cdchits=my_cdchits.size();
   while (break_point_cdc_index>4 && break_point_step_index>0 
         && break_point_step_index<central_traj.size()){
      refit_iter++;

      // Flag the cdc hits within the radius of the break point cdc index
      // as good, the rest as bad.
      for (unsigned int j=0;j<=break_point_cdc_index;j++){
         if (my_cdchits[j]->status!=late_hit)my_cdchits[j]->status=good_hit;
      }
      for (unsigned int j=break_point_cdc_index+1;j<num_cdchits;j++){
         my_cdchits[j]->status=bad_hit;
      }

      // Now refit with the truncated trajectory and list of hits
      //C1=C0;
      //C1=4.0*C0;
      C1=1.0*C0;
      S1=central_traj[break_point_step_index].S;
      refit_chisq=-1.;
      refit_ndf=0; 
      refit_error=KalmanCentral(anneal_factor,S1,C1,refit_xy,
            refit_chisq,refit_ndf);

      if (refit_error==FIT_SUCCEEDED
            || (refit_error==BREAK_POINT_FOUND 
               && break_point_cdc_index==1
               )
            //|| refit_error==PRUNED_TOO_MANY_HITS
         ){
         C=C1;
         S=S1;
         xy=refit_xy;
         chisq=refit_chisq;
         numdof=refit_ndf;

         if (DEBUG_LEVEL>0) _DBG_ << "Fit recovery succeeded..." << endl;
         return FIT_SUCCEEDED;
      }

      break_point_cdc_index=old_cdc_index-refit_iter;
      break_point_step_index=old_step_index-refit_iter;
   }

   // If the refit did not work, restore the old list hits used in the fit 
   // before the track recovery was attempted.
   for (unsigned int k=0;k<num_hits;k++){
      cdc_used_in_fit[k]=old_cdc_used_status[k];
   }

   if (DEBUG_LEVEL>0) _DBG_ << "Fit recovery failed..." << endl;
   return FIT_FAILED;
}

// Track recovery for forward tracks
//-----------------------------------
// This code attempts to recover tracks that are "broken".  Sometimes the fit fails because too many hits were pruned 
// such that the number of surviving hits is less than the minimum number of degrees of freedom for a five-parameter fit.
// This condition is flagged as PRUNED_TOO_MANY_HITS.  It may also be the case that even if we used enough hits for the fit to
// be valid (i.e., make sense in terms of the number of degrees of freedom for the number of parameters (5) we are trying 
// to determine), we throw away a number of hits near the target because the projected trajectory deviates too far away from 
// the actual hits.  This may be an indication of a kink in the trajectory.  This condition is flagged as BREAK_POINT_FOUND.
// Sometimes the fit may succeed and no break point is found, but the pruning threw out a large number of hits.  This 
// condition is flagged as PRUNED_TOO_MANY_HITS.  The recovery code proceeds by truncating the reference trajectory to 
// to a point just after the hit where a break point was found and all hits are ignored beyond this point subject to a 
// minimum number of hits.  The recovery code always attempts to use the hits closest to the target.  The code is allowed to 
// iterate; with each iteration the trajectory and list of useable hits is further truncated.
kalman_error_t DTrackFitterKalmanSIMD::RecoverBrokenTracks(double anneal_factor, 
      DMatrix5x1 &S, 
      DMatrix5x5 &C,
      const DMatrix5x5 &C0,
      double &chisq, 
      unsigned int &numdof){
   if (DEBUG_LEVEL>1) _DBG_  << "Attempting to recover broken track ... " <<endl;

   unsigned int num_cdchits=my_cdchits.size();
   unsigned int num_fdchits=my_fdchits.size();

   // Initialize degrees of freedom and chi^2
   double refit_chisq=-1.;
   unsigned int refit_ndf=0;
   // State vector and covariance matrix
   DMatrix5x5 C1;
   DMatrix5x1 S1;

   // save the status of the hits used in the fit
   vector<bool>old_cdc_used_status(num_cdchits);
   vector<bool>old_fdc_used_status(num_fdchits);
   for (unsigned int j=0;j<num_fdchits;j++){
      old_fdc_used_status[j]=fdc_used_in_fit[j];
   }
   for (unsigned int j=0;j<num_cdchits;j++){
      old_cdc_used_status[j]=cdc_used_in_fit[j];
   }

   unsigned int min_cdc_index=MIN_HITS_FOR_REFIT-1;
   if (my_fdchits.size()>0){
      if (break_point_cdc_index<5){
         break_point_cdc_index=0;
         min_cdc_index=5;
      }
   }
   /*else{
     unsigned int half_num_cdchits=num_cdchits/2;
     if (break_point_cdc_index<half_num_cdchits
     && half_num_cdchits>min_cdc_index)
     break_point_cdc_index=half_num_cdchits;
     }
     */
   if (break_point_cdc_index<min_cdc_index){
      break_point_cdc_index=min_cdc_index;
   }

   // Find the index at which to truncate the reference trajectory
   DVector2 origin=my_cdchits[break_point_cdc_index]->origin;
   DVector2 dir=my_cdchits[break_point_cdc_index]->dir;
   double z0=my_cdchits[break_point_cdc_index]->z0wire;
   unsigned int k=forward_traj.size()-1;
   for (;k>1;k--){
      DMatrix5x1 S1=forward_traj[k].S;
      double x1=S1(state_x);
      double y1=S1(state_y);
      double r2=x1*x1+y1*y1;
      DMatrix5x1 S2=forward_traj[k-1].S;
      double x2=S2(state_x);
      double y2=S2(state_y);
      double r2next=x2*x2+y2*y2;
      double r2cdc=(origin+(forward_traj[k].z-z0)*dir).Mod2();

      if (r2next>r2 && r2>r2cdc) break;
   }
   break_point_step_index=k;

   if (break_point_step_index==forward_traj.size()-1){
      if (DEBUG_LEVEL>0) _DBG_ << "Invalid reference trajectory in track recovery..." << endl;
      return FIT_FAILED;
   }

   // Attemp to refit the track using the abreviated list of hits and the truncated
   // reference trajectory.  Iterates if the fit fails.
   kalman_error_t refit_error=FIT_NOT_DONE;
   unsigned int old_cdc_index=break_point_cdc_index;
   unsigned int old_step_index=break_point_step_index;
   unsigned int refit_iter=0;
   while (break_point_cdc_index>4 && break_point_step_index>0 
         && break_point_step_index<forward_traj.size()){
      refit_iter++;

      // Flag the cdc hits within the radius of the break point cdc index
      // as good, the rest as bad.
      for (unsigned int j=0;j<=break_point_cdc_index;j++){
         if (my_cdchits[j]->status!=late_hit) my_cdchits[j]->status=good_hit;
      }
      for (unsigned int j=break_point_cdc_index+1;j<num_cdchits;j++){
         my_cdchits[j]->status=bad_hit;
      }

      // Re-initialize the state vector, covariance, chisq and number of degrees of freedom    
      //C1=4.0*C0;
      C1=1.0*C0;
      S1=forward_traj[break_point_step_index].S;
      refit_chisq=-1.;
      refit_ndf=0;
      // Now refit with the truncated trajectory and list of hits
      refit_error=KalmanForwardCDC(anneal_factor,S1,C1,
            refit_chisq,refit_ndf);   
      if (refit_error==FIT_SUCCEEDED 
            || (refit_error==BREAK_POINT_FOUND
               && break_point_cdc_index==1
               )
            //|| refit_error==PRUNED_TOO_MANY_HITS
         ){
         C=C1;
         S=S1;
         chisq=refit_chisq;
         numdof=refit_ndf;
         return FIT_SUCCEEDED;
      }
      break_point_cdc_index=old_cdc_index-refit_iter;
      break_point_step_index=old_step_index-refit_iter;
   }
   // If the refit did not work, restore the old list hits used in the fit 
   // before the track recovery was attempted.
   for (unsigned int k=0;k<num_cdchits;k++){
      cdc_used_in_fit[k]=old_cdc_used_status[k];
   }
   for (unsigned int k=0;k<num_fdchits;k++){
      fdc_used_in_fit[k]=old_fdc_used_status[k];
   }   

   return FIT_FAILED;
}


// Track recovery for forward-going tracks with hits in the FDC 
kalman_error_t
DTrackFitterKalmanSIMD::RecoverBrokenForwardTracks(double fdc_anneal,
      double cdc_anneal,
      DMatrix5x1 &S, 
      DMatrix5x5 &C,
      const DMatrix5x5 &C0,
      double &chisq, 
      unsigned int &numdof){
   if (DEBUG_LEVEL>1)
      _DBG_  << "Attempting to recover broken track ... " <<endl;
   unsigned int num_cdchits=my_cdchits.size();
   unsigned int num_fdchits=fdc_updates.size();

   // Initialize degrees of freedom and chi^2
   double refit_chisq=-1.;
   unsigned int refit_ndf=0;
   // State vector and covariance matrix
   DMatrix5x5 C1;
   DMatrix5x1 S1;

   // save the status of the hits used in the fit
   vector<int>old_cdc_used_status(num_cdchits);
   vector<int>old_fdc_used_status(num_fdchits);
   for (unsigned int j=0;j<num_fdchits;j++){
      old_fdc_used_status[j]=fdc_used_in_fit[j];
   }
   for (unsigned int j=0;j<num_cdchits;j++){     
      old_cdc_used_status[j]=cdc_used_in_fit[j];
   }

   // Truncate the trajectory
   double zhit=my_fdchits[break_point_fdc_index]->z;   
   unsigned int k=0;
   for (;k<forward_traj.size();k++){
      double z=forward_traj[k].z;
      if (z<zhit) break;
   }
   for (unsigned int j=0;j<=break_point_fdc_index;j++){
      my_fdchits[j]->status=good_hit;
   }
   for (unsigned int j=break_point_fdc_index+1;j<num_fdchits;j++){
      my_fdchits[j]->status=bad_hit;
   }

   if (k==forward_traj.size()) return FIT_NOT_DONE;

   break_point_step_index=k;

   // Attemp to refit the track using the abreviated list of hits and the truncated
   // reference trajectory. 
   kalman_error_t refit_error=FIT_NOT_DONE;
   int refit_iter=0;
   unsigned int break_id=break_point_fdc_index;
   while (break_id+num_cdchits>=MIN_HITS_FOR_REFIT && break_id>0 
         && break_point_step_index<forward_traj.size()
         && break_point_step_index>1
         && refit_iter<10){
      refit_iter++;

      // Mark the hits as bad if they are not included
      if (break_id >= 0){
         for (unsigned int j=0;j<num_cdchits;j++){
            if (my_cdchits[j]->status!=late_hit)my_cdchits[j]->status=good_hit;
         }
         for (unsigned int j=0;j<=break_id;j++){
            my_fdchits[j]->status=good_hit;
         }
         for (unsigned int j=break_id+1;j<num_fdchits;j++){
            my_fdchits[j]->status=bad_hit;
         }
      }
      else{
         // BreakID should always be 0 or positive, so this should never happen, but could be investigated in the future.
         for (unsigned int j=0;j<num_cdchits+break_id;j++){
            if (my_cdchits[j]->status!=late_hit) my_cdchits[j]->status=good_hit;
         }
         for (unsigned int j=num_cdchits+break_id;j<num_cdchits;j++){
            my_cdchits[j]->status=bad_hit;
         }
         for (unsigned int j=0;j<num_fdchits;j++){
            my_fdchits[j]->status=bad_hit;
         }
      }

      // Re-initialize the state vector, covariance, chisq and number of degrees of freedom    
      //C1=4.0*C0;
      C1=1.0*C0;
      S1=forward_traj[break_point_step_index].S;
      refit_chisq=-1.;
      refit_ndf=0;

      // Now refit with the truncated trajectory and list of hits
      refit_error=KalmanForward(fdc_anneal,cdc_anneal,S1,C1,refit_chisq,refit_ndf);   
      if (refit_error==FIT_SUCCEEDED
            //|| (refit_error==PRUNED_TOO_MANY_HITS)
         ){
         C=C1;
         S=S1;
         chisq=refit_chisq;
         numdof=refit_ndf;

         if (DEBUG_LEVEL>1)  _DBG_ << "Refit succeeded" << endl;
         return FIT_SUCCEEDED;
      }
      // Truncate the trajectory
      if (break_id>0) break_id--;
      else break;
      zhit=my_fdchits[break_id]->z;  
      k=0;
      for (;k<forward_traj.size();k++){
         double z=forward_traj[k].z;
         if (z<zhit) break;
      }
      break_point_step_index=k;

   }

   // If the refit did not work, restore the old list hits used in the fit 
   // before the track recovery was attempted.
   for (unsigned int k=0;k<num_cdchits;k++){
      cdc_used_in_fit[k]=old_cdc_used_status[k];
   }
   for (unsigned int k=0;k<num_fdchits;k++){
      fdc_used_in_fit[k]=old_fdc_used_status[k];
   }   

   return FIT_FAILED;
}



// Routine to fit hits in the FDC and the CDC using the forward parametrization
kalman_error_t DTrackFitterKalmanSIMD::ForwardFit(const DMatrix5x1 &S0,const DMatrix5x5 &C0){   
   unsigned int num_cdchits=my_cdchits.size();
   unsigned int num_fdchits=my_fdchits.size();
   unsigned int max_fdc_index=num_fdchits-1;

   // Vectors to keep track of updated state vectors and covariance matrices (after
   // adding the hit information)
   vector<bool>last_fdc_used_in_fit(num_fdchits);
   vector<bool>last_cdc_used_in_fit(num_cdchits);
   vector<pull_t>forward_pulls;
   vector<pull_t>last_forward_pulls;

   // Charge
   // double q=input_params.charge();

   // Covariance matrix and state vector
   DMatrix5x5 C;
   DMatrix5x1 S=S0;

   // Create matrices to store results from previous iteration
   DMatrix5x1 Slast(S);
   DMatrix5x5 Clast(C0); 
   // last z position
   double last_z=z_;

   double fdc_anneal=FORWARD_ANNEAL_SCALE+1.;  // variable for scaling cut for hit pruning
   double my_fdc_anneal_const=FORWARD_ANNEAL_POW_CONST;
   //  if (fit_type==kTimeBased && fabs(1./S(state_q_over_p))<1.0
   // && my_anneal_const>=2.0) my_anneal_const*=0.5;
   double cdc_anneal=(fit_type==kTimeBased?ANNEAL_SCALE+1.:2.);  // variable for scaling cut for hit pruning
   double my_cdc_anneal_const=ANNEAL_POW_CONST;


   // Chi-squared and degrees of freedom
   double chisq=-1.,chisq_forward=-1.;
   unsigned int my_ndf=0;
   unsigned int last_ndf=1;  

   // Iterate over reference trajectories
   for (int iter=0;
         iter<(fit_type==kTimeBased?MAX_TB_PASSES:MAX_WB_PASSES);
         iter++) {      
      // These variables provide the approximate location along the trajectory
      // where there is an indication of a kink in the track      
      break_point_fdc_index=max_fdc_index;
      break_point_cdc_index=0;
      break_point_step_index=0;

      // Reset material map index
      last_material_map=0;

      // Abort if momentum is too low
      if (fabs(S(state_q_over_p))>Q_OVER_P_MAX) break;

      // Initialize path length variable and flight time
      len=0;
      ftime=0.;
      var_ftime=0.;

      // Scale cut for pruning hits according to the iteration number
      if (fit_type==kTimeBased)
      {
         fdc_anneal=FORWARD_ANNEAL_SCALE/pow(my_fdc_anneal_const,iter)+1.;
         cdc_anneal=ANNEAL_SCALE/pow(my_cdc_anneal_const,iter)+1.;
      }

      // Swim through the field out to the most upstream FDC hit
      jerror_t ref_track_error=SetReferenceTrajectory(S);
      if (ref_track_error!=NOERROR){
	if (iter==0) return FIT_FAILED;  
	break;
      }

      // Reset the status of the cdc hits 
      for (unsigned int j=0;j<num_cdchits;j++){
	if (my_cdchits[j]->status!=late_hit)my_cdchits[j]->status=good_hit;
      }	
      
      // perform the kalman filter 
      C=C0;
      kalman_error_t error=KalmanForward(fdc_anneal,cdc_anneal,S,C,chisq,my_ndf);
      if (DEBUG_LEVEL>1){
	_DBG_ << "Iter: " << iter+1 << " Chi2=" << chisq << " Ndf=" << my_ndf << " Error code: " << error << endl; 
      }
      
      // Try to recover tracks that failed the first attempt at fitting
      if (error!=FIT_SUCCEEDED && RECOVER_BROKEN_TRACKS
	  && num_fdchits>2  // some minimum to make this worthwhile...
	  && break_point_fdc_index<num_fdchits
	  && break_point_fdc_index+num_cdchits>=MIN_HITS_FOR_REFIT
	  && forward_traj.size()>2*MIN_HITS_FOR_REFIT // avoid small track stubs
	  ){
	DMatrix5x5 Ctemp=C;
	DMatrix5x1 Stemp=S;
	unsigned int temp_ndf=my_ndf;
	double temp_chi2=chisq;
	double x=x_,y=y_,z=z_;
	
	kalman_error_t refit_error=RecoverBrokenForwardTracks(fdc_anneal,
							      cdc_anneal,
							      S,C,C0,chisq,
							      my_ndf);
	if (refit_error!=FIT_SUCCEEDED){
	  if (error==PRUNED_TOO_MANY_HITS || error==BREAK_POINT_FOUND){
	    C=Ctemp;
	    S=Stemp;
	    my_ndf=temp_ndf;
	    chisq=temp_chi2;
	    x_=x,y_=y,z_=z;
	    
	    if (num_cdchits<6) error=FIT_SUCCEEDED;
	  }
	  else error=FIT_FAILED;
	}
	else error=FIT_SUCCEEDED;
      }
      if ((error==POSITION_OUT_OF_RANGE || error==MOMENTUM_OUT_OF_RANGE)
	  && iter==0){ 
	return FIT_FAILED;
      }
      if (error==FIT_FAILED || error==INVALID_FIT  || my_ndf==0){
	if (iter==0) return FIT_FAILED; // first iteration failed
	break;
      }
      
      if (iter>0){
	double new_reduced_chisq=chisq/my_ndf;
	double old_reduced_chisq=chisq_forward/last_ndf;
	if (new_reduced_chisq>old_reduced_chisq 
	    || fabs(new_reduced_chisq-old_reduced_chisq)<CHISQ_DELTA) break;
      }
      
      chisq_forward=chisq; 
      last_ndf=my_ndf;
      Slast=S;
      Clast=C;	 
      last_z=z_;
      
      IsSmoothed=false;
      if(fit_type==kTimeBased){
	forward_pulls.clear();
	if (SmoothForward(forward_pulls) == NOERROR){
	  IsSmoothed = true;
	}
	last_forward_pulls.assign(forward_pulls.begin(),forward_pulls.end());
      }
      
      last_fdc_used_in_fit=fdc_used_in_fit;
      last_cdc_used_in_fit=cdc_used_in_fit;
   } //iteration
   
   // total chisq and ndf
   chisq_=chisq_forward;
   ndf_=last_ndf;

   // output lists of hits used in the fit and fill pull vector
   cdchits_used_in_fit.clear();
   for (unsigned int m=0;m<last_cdc_used_in_fit.size();m++){
      if (last_cdc_used_in_fit[m]){
         cdchits_used_in_fit.push_back(my_cdchits[m]->hit);
      }
   }
   fdchits_used_in_fit.clear();
   for (unsigned int m=0;m<last_fdc_used_in_fit.size();m++){
      if (last_fdc_used_in_fit[m]){
         fdchits_used_in_fit.push_back(my_fdchits[m]->hit);
      }
   }
   // fill pull vector
   pulls.assign(last_forward_pulls.begin(),last_forward_pulls.end());

   // fill vector of extrapolations
   ClearExtrapolations();
   ExtrapolateForwardToOtherDetectors(); 
   if (extrapolations.at(SYS_BCAL).size()==1){
     // There needs to be some steps inside the the volume of the BCAL for 
     // the extrapolation to be useful.  If this is not the case, clear 
     // the extrolation vector.
     extrapolations[SYS_BCAL].clear();
   }
   
   // Extrapolate to the point of closest approach to the beam line
   z_=last_z;
   if (sqrt(Slast(state_x)*Slast(state_x)+Slast(state_y)*Slast(state_y))
         >EPS2){  
      DMatrix5x5 Ctemp=Clast;
      DMatrix5x1 Stemp=Slast; 
      double ztemp=z_;
      if (ExtrapolateToVertex(Stemp,Ctemp)==NOERROR){
         Clast=Ctemp;
         Slast=Stemp;
      }
      else{
         //_DBG_ << endl;
         z_=ztemp;
      }
   }

   // Final momentum, positions and tangents
   x_=Slast(state_x), y_=Slast(state_y);
   tx_=Slast(state_tx),ty_=Slast(state_ty);
   q_over_p_=Slast(state_q_over_p);

   // Convert from forward rep. to central rep.
   double tsquare=tx_*tx_+ty_*ty_;
   tanl_=1./sqrt(tsquare);
   double cosl=cos(atan(tanl_));
   q_over_pt_=q_over_p_/cosl;
   phi_=atan2(ty_,tx_);
   if (FORWARD_PARMS_COV==false){
     DVector2 beam_pos=beam_center+(z_-beam_z0)*beam_dir;
     double dx=x_-beam_pos.X();
     double dy=y_-beam_pos.Y();
     D_=sqrt(dx*dx+dy*dy)+EPS;
     x_ = dx; y_ = dy;
     double cosphi=cos(phi_);
     double sinphi=sin(phi_);
     if ((dx>0.0 && sinphi>0.0) || (dy<0.0 && cosphi>0.0) 
	 || (dy>0.0 && cosphi<0.0) || (dx<0.0 && sinphi<0.0)) D_*=-1.;
     TransformCovariance(Clast);      
   }
   // Covariance matrix  
   vector<double>dummy;
   for (unsigned int i=0;i<5;i++){
      dummy.clear();
      for(unsigned int j=0;j<5;j++){
         dummy.push_back(Clast(i,j));
      }
      fcov.push_back(dummy);
   }

   return FIT_SUCCEEDED;
}

// Routine to fit hits in the CDC using the forward parametrization
kalman_error_t DTrackFitterKalmanSIMD::ForwardCDCFit(const DMatrix5x1 &S0,const DMatrix5x5 &C0){   
   unsigned int num_cdchits=my_cdchits.size();
   unsigned int max_cdc_index=num_cdchits-1;
   unsigned int min_cdc_index_for_refit=MIN_HITS_FOR_REFIT-1;

   // Charge
   //  double q=input_params.charge();

   // Covariance matrix and state vector
   DMatrix5x5 C;
   DMatrix5x1 S=S0;

   // Create matrices to store results from previous iteration
   DMatrix5x1 Slast(S);
   DMatrix5x5 Clast(C0); 

   // Vectors to keep track of updated state vectors and covariance matrices (after
   // adding the hit information)
   vector<pull_t>cdc_pulls;
   vector<pull_t>last_cdc_pulls;
   vector<bool>last_cdc_used_in_fit;

   double anneal_scale=ANNEAL_SCALE; // variable for scaling cut for hit pruning
   double my_anneal_const=ANNEAL_POW_CONST;
   /* 
      double tsquare=S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty);
      double tanl=1./sqrt(tsquare);
      if (tanl>2.5){
      anneal_scale=FORWARD_ANNEAL_SCALE;
      my_anneal_const=FORWARD_ANNEAL_POW_CONST;
      }
      */
   double anneal_factor=anneal_scale+1.;
   /*
      if (fit_type==kTimeBased && fabs(1./S(state_q_over_p))<1.0 
      && my_anneal_const>=2.0){ 
      my_anneal_const*=0.5;
      }
      */
   // Chi-squared and degrees of freedom
   double chisq=-1.,chisq_forward=-1.;
   unsigned int my_ndf=0;
   unsigned int last_ndf=1;
   // last z position
   double zlast=0.;
   //  unsigned int last_break_point_index=0,last_break_point_step_index=0;

   // Iterate over reference trajectories
   for (int iter2=0;
         iter2<(fit_type==kTimeBased?MAX_TB_PASSES:MAX_WB_PASSES);
         iter2++){   
      if (DEBUG_LEVEL>1){
         _DBG_ <<"-------- iteration " << iter2+1 << " -----------" <<endl;
      } 

      // These variables provide the approximate location along the trajectory
      // where there is an indication of a kink in the track      
      break_point_cdc_index=max_cdc_index;
      break_point_step_index=0;

      // Reset material map index
      last_material_map=0;

      // Abort if momentum is too low
      if (fabs(S(state_q_over_p))>Q_OVER_P_MAX){
         //printf("Too low momentum? %f\n",1/S(state_q_over_p));
         break;
      }

      // Scale cut for pruning hits according to the iteration number
      if (fit_type==kTimeBased)
      {   
         anneal_factor=anneal_scale/pow(my_anneal_const,iter2)+1.;
      }

      // Initialize path length variable and flight time
      len=0;
      ftime=0.;
      var_ftime=0.;

      // Swim to create the reference trajectory
      jerror_t ref_track_error=SetCDCForwardReferenceTrajectory(S);
      if (ref_track_error!=NOERROR){
	if (iter2==0) return FIT_FAILED;  
	break;
      }

      // Reset the status of the cdc hits 
      for (unsigned int j=0;j<num_cdchits;j++){
	if (my_cdchits[j]->status!=late_hit)my_cdchits[j]->status=good_hit;
      }
      
      // perform the filter 
      C=C0;
      kalman_error_t error=KalmanForwardCDC(anneal_factor,S,C,chisq,my_ndf); 
      
      // Try to recover tracks that failed the first attempt at fitting
      if (error!=FIT_SUCCEEDED && RECOVER_BROKEN_TRACKS 
	  && num_cdchits>=MIN_HITS_FOR_REFIT){
	DMatrix5x5 Ctemp=C;
	DMatrix5x1 Stemp=S;
	unsigned int temp_ndf=my_ndf;
	double temp_chi2=chisq;
	double x=x_,y=y_,z=z_;
	
	if (error==MOMENTUM_OUT_OF_RANGE){
	  // _DBG_ << "Momentum out of range" <<endl;	
	  unsigned int new_index=(3*max_cdc_index)/4;
	  break_point_cdc_index=(new_index>min_cdc_index_for_refit)?new_index:min_cdc_index_for_refit;
	}
	
	if (error==BROKEN_COVARIANCE_MATRIX){
	  break_point_cdc_index=min_cdc_index_for_refit;
	  //_DBG_ << "Bad Cov" <<endl;
	}
	if (error==POSITION_OUT_OF_RANGE){
	  //_DBG_ << "Bad position" << endl;
	  unsigned int new_index=(max_cdc_index)/2;
	  break_point_cdc_index=(new_index>min_cdc_index_for_refit)?new_index:min_cdc_index_for_refit;
	}
	if (error==PRUNED_TOO_MANY_HITS){
	  //_DBG_ << "Prune" << endl;
	  unsigned int new_index=(3*max_cdc_index)/4;
	  break_point_cdc_index=(new_index>min_cdc_index_for_refit)?new_index:min_cdc_index_for_refit;
	  // anneal_factor*=10.;
	}
	
	kalman_error_t refit_error=RecoverBrokenTracks(anneal_factor,S,C,C0,chisq,my_ndf);    
	
	if (refit_error!=FIT_SUCCEEDED){
	  if (error==PRUNED_TOO_MANY_HITS || error==BREAK_POINT_FOUND){
	    C=Ctemp;
	    S=Stemp;
	    my_ndf=temp_ndf;
	    chisq=temp_chi2;
	    x_=x,y_=y,z_=z;
	    
	    error=FIT_SUCCEEDED;
	  }
	  else error=FIT_FAILED;
	}
	else error=FIT_SUCCEEDED;
      }
      if ((error==POSITION_OUT_OF_RANGE || error==MOMENTUM_OUT_OF_RANGE)
	  && iter2==0){ 
	return FIT_FAILED;
      }
      if (error==FIT_FAILED || error==INVALID_FIT || my_ndf==0){  
	if (iter2==0) return error;
	break;
      }

      if (DEBUG_LEVEL>1)  _DBG_ << "--> Chisq " << chisq << " NDF " 
				<< my_ndf 
				<< " Prob: " << TMath::Prob(chisq,my_ndf)
				<< endl;
      
      if (iter2>0){
	double new_reduced_chisq=chisq/my_ndf;
	double old_reduced_chisq=chisq_forward/last_ndf;
	if (new_reduced_chisq>old_reduced_chisq 
	    || fabs(new_reduced_chisq-old_reduced_chisq)<CHISQ_DELTA) break;
      }
      
      // Run the smoother
      IsSmoothed=false;
      if(fit_type==kTimeBased){
	cdc_pulls.clear();
	if (SmoothForwardCDC(cdc_pulls) == NOERROR){
	  IsSmoothed = true;
	}
	last_cdc_pulls.assign(cdc_pulls.begin(),cdc_pulls.end());
      }
      
      chisq_forward=chisq;
      Slast=S;
      Clast=C;
      last_ndf=my_ndf;
      zlast=z_;
      
      last_cdc_used_in_fit=cdc_used_in_fit;
   } //iteration
 
   // total chisq and ndf
   chisq_=chisq_forward;
   ndf_=last_ndf;

   // output lists of hits used in the fit and fill the pull vector
   cdchits_used_in_fit.clear();
   for (unsigned int m=0;m<last_cdc_used_in_fit.size();m++){
      if (last_cdc_used_in_fit[m]){
         cdchits_used_in_fit.push_back(my_cdchits[m]->hit);
      }
   }  
   // output pulls vector
   pulls.assign(last_cdc_pulls.begin(),last_cdc_pulls.end());

   // Fill extrapolation vector
   ClearExtrapolations();
   ExtrapolateForwardToOtherDetectors();  
   if (extrapolations.at(SYS_BCAL).size()==1){
     // There needs to be some steps inside the the volume of the BCAL for 
     // the extrapolation to be useful.  If this is not the case, clear 
     // the extrolation vector.
     extrapolations[SYS_BCAL].clear();
   }

   // Extrapolate to the point of closest approach to the beam line
   z_=zlast;
   if (sqrt(Slast(state_x)*Slast(state_x)+Slast(state_y)*Slast(state_y))
         >EPS2) 
      if (ExtrapolateToVertex(Slast,Clast)!=NOERROR) return EXTRAPOLATION_FAILED;

   // Final momentum, positions and tangents
   x_=Slast(state_x), y_=Slast(state_y);
   tx_=Slast(state_tx),ty_=Slast(state_ty);
   q_over_p_=Slast(state_q_over_p);

   // Convert from forward rep. to central rep.
   double tsquare=tx_*tx_+ty_*ty_;
   tanl_=1./sqrt(tsquare);
   double cosl=cos(atan(tanl_));
   q_over_pt_=q_over_p_/cosl;
   phi_=atan2(ty_,tx_);
   if (FORWARD_PARMS_COV==false){
     DVector2 beam_pos=beam_center+(z_-beam_z0)*beam_dir;
     double dx=x_-beam_pos.X();
     double dy=y_-beam_pos.Y();
     D_=sqrt(dx*dx+dy*dy)+EPS;
     x_ = dx; y_ = dy;
     double cosphi=cos(phi_);
     double sinphi=sin(phi_);
     if ((dx>0.0 && sinphi>0.0) || (dy<0.0 && cosphi>0.0) 
	 || (dy>0.0 && cosphi<0.0) || (dx<0.0 && sinphi<0.0)) D_*=-1.;
     TransformCovariance(Clast);
   }
   // Covariance matrix  
   vector<double>dummy;
   // ... forward parameterization
   fcov.clear();
   for (unsigned int i=0;i<5;i++){
      dummy.clear();
      for(unsigned int j=0;j<5;j++){
         dummy.push_back(Clast(i,j));
      }
      fcov.push_back(dummy);
   }  

   return FIT_SUCCEEDED;
}

// Routine to fit hits in the CDC using the central parametrization
kalman_error_t DTrackFitterKalmanSIMD::CentralFit(const DVector2 &startpos,
      const DMatrix5x1 &S0,
      const DMatrix5x5 &C0){   
   // Initial position in x and y
   DVector2 pos(startpos);

   // Charge
   //  double q=input_params.charge();

   // Covariance matrix and state vector
   DMatrix5x5 Cc;
   DMatrix5x1 Sc=S0;

   // Variables to store values from previous iterations
   DMatrix5x1 Sclast(Sc);
   DMatrix5x5 Cclast(C0);
   DVector2 last_pos=pos;

   unsigned int num_cdchits=my_cdchits.size();
   unsigned int max_cdc_index=num_cdchits-1;
   unsigned int min_cdc_index_for_refit=MIN_HITS_FOR_REFIT-1;

   // Vectors to keep track of updated state vectors and covariance matrices (after
   // adding the hit information)
   vector<pull_t>last_cdc_pulls;
   vector<pull_t>cdc_pulls;
   vector<bool>last_cdc_used_in_fit(num_cdchits);

   double anneal_factor=ANNEAL_SCALE+1.; // variable for scaling cut for hit pruning
   double my_anneal_const=ANNEAL_POW_CONST;
   //if (fit_type==kTimeBased && fabs(1./Sc(state_q_over_p))<1.0) my_anneal_const*=0.5;

   //Initialization of chisq, ndf, and error status
   double chisq_iter=-1.,chisq=-1.;
   unsigned int my_ndf=0;
   ndf_=0.;
   unsigned int last_ndf=1;
   kalman_error_t error=FIT_NOT_DONE;

   // Iterate over reference trajectories
   int iter2=0;
   for (;iter2<(fit_type==kTimeBased?MAX_TB_PASSES:MAX_WB_PASSES);
         iter2++){     
      if (DEBUG_LEVEL>1){
         _DBG_ <<"-------- iteration " << iter2+1 << " -----------" <<endl;
      }  

      // These variables provide the approximate location along the trajectory
      // where there is an indication of a kink in the track
      break_point_cdc_index=max_cdc_index;
      break_point_step_index=0;

      // Reset material map index
      last_material_map=0;

      // Break out of loop if p is too small
      double q_over_p=Sc(state_q_over_pt)*cos(atan(Sc(state_tanl)));
      if (fabs(q_over_p)>Q_OVER_P_MAX) break;

      // Initialize path length variable and flight time
      len=0.;
      ftime=0.;
      var_ftime=0.;

      // Scale cut for pruning hits according to the iteration number
      if (fit_type==kTimeBased)
      {
         anneal_factor=ANNEAL_SCALE/pow(my_anneal_const,iter2)+1.;
      }

      // Initialize trajectory deque
      jerror_t ref_track_error=SetCDCReferenceTrajectory(last_pos,Sc);
      if (ref_track_error!=NOERROR){
	if (iter2==0) return FIT_FAILED;  
	break;
      }

      // Reset the status of the cdc hits 
      for (unsigned int j=0;j<num_cdchits;j++){
	if (my_cdchits[j]->status!=late_hit)my_cdchits[j]->status=good_hit;
      }
      
      // perform the fit
      Cc=C0;
      error=KalmanCentral(anneal_factor,Sc,Cc,pos,chisq,my_ndf);
      // Try to recover tracks that failed the first attempt at fitting
      if (error!=FIT_SUCCEEDED && RECOVER_BROKEN_TRACKS
	  && num_cdchits>=MIN_HITS_FOR_REFIT){
	DVector2 temp_pos=pos;
	DMatrix5x1 Stemp=Sc;
	DMatrix5x5 Ctemp=Cc;
	unsigned int temp_ndf=my_ndf;
	double temp_chi2=chisq;
	
	if (error==MOMENTUM_OUT_OF_RANGE){
	  //_DBG_ << "Momentum out of range" <<endl; 
	  unsigned int new_index=(3*max_cdc_index)/4;
	  break_point_cdc_index=(new_index>min_cdc_index_for_refit)?new_index:min_cdc_index_for_refit;
	}
	
	if (error==BROKEN_COVARIANCE_MATRIX){ 
	  break_point_cdc_index=min_cdc_index_for_refit;  
	  //_DBG_ << "Bad Cov" <<endl;
	}
	if (error==POSITION_OUT_OF_RANGE){
	  //_DBG_ << "Bad position" << endl;
	  unsigned int new_index=(max_cdc_index)/2;
	  break_point_cdc_index=(new_index>min_cdc_index_for_refit)?new_index:min_cdc_index_for_refit;
	}
	if (error==PRUNED_TOO_MANY_HITS){	 
	  unsigned int new_index=(3*max_cdc_index)/4;
	  break_point_cdc_index=(new_index>min_cdc_index_for_refit)?new_index:min_cdc_index_for_refit;
	  //anneal_factor*=10.;
	  //_DBG_ << "Prune" << endl;	  
	}

	kalman_error_t refit_error=RecoverBrokenTracks(anneal_factor,Sc,Cc,C0,pos,chisq,my_ndf);  
	if (refit_error!=FIT_SUCCEEDED){
	  if (error==PRUNED_TOO_MANY_HITS || error==BREAK_POINT_FOUND){
	    Cc=Ctemp;
	    Sc=Stemp;
	    my_ndf=temp_ndf;
	    chisq=temp_chi2;
	    pos=temp_pos;
	    if (DEBUG_LEVEL > 1) _DBG_ << " Refit did not succeed, but restoring old values" << endl;
	    
	    error=FIT_SUCCEEDED;
	  }
	  else error=FIT_FAILED;
	}
	else error=FIT_SUCCEEDED;
      }
      if ((error==POSITION_OUT_OF_RANGE || error==MOMENTUM_OUT_OF_RANGE)
	  && iter2==0){ 
	return FIT_FAILED;
      }
      if (error==FIT_FAILED || error==INVALID_FIT || my_ndf==0){
	if (iter2==0) return error;
	break;
      } 

      if (DEBUG_LEVEL>1) _DBG_ << "--> Chisq " << chisq << " Ndof " << my_ndf 
			       << " Prob: " << TMath::Prob(chisq,my_ndf)
			       << endl;
      
      if (iter2>0){
	double new_reduced_chisq=chisq/my_ndf;
	double old_reduced_chisq=chisq_iter/last_ndf;
	if (new_reduced_chisq>old_reduced_chisq 
	    || fabs(new_reduced_chisq-old_reduced_chisq)<CHISQ_DELTA) break;
      }
      
      // Save the current state vector and covariance matrix
      Cclast=Cc;
      Sclast=Sc;
      last_pos=pos;
      chisq_iter=chisq;
      last_ndf=my_ndf;

      // Run smoother and fill pulls vector
      IsSmoothed=false;
      if(fit_type==kTimeBased){
	cdc_pulls.clear();
	if (SmoothCentral(cdc_pulls) == NOERROR){
	  IsSmoothed = true;
	}
	last_cdc_pulls.assign(cdc_pulls.begin(),cdc_pulls.end()); 
      }
      
      last_cdc_used_in_fit=cdc_used_in_fit;
   }
	 
   // Fill extrapolations vector
   ClearExtrapolations();
   ExtrapolateCentralToOtherDetectors();
   if (extrapolations.at(SYS_BCAL).size()==1){
     // There needs to be some steps inside the the volume of the BCAL for 
     // the extrapolation to be useful.  If this is not the case, clear 
     // the extrolation vector.
     extrapolations[SYS_BCAL].clear();
   }
   if (last_pos.Mod()>0.001){ // in cm
      if (ExtrapolateToVertex(last_pos,Sclast,Cclast)!=NOERROR) return EXTRAPOLATION_FAILED; 
   }

   // output lists of hits used in the fit 
   cdchits_used_in_fit.clear();
   for (unsigned int m=0;m<last_cdc_used_in_fit.size();m++){
      if (last_cdc_used_in_fit[m]){
         cdchits_used_in_fit.push_back(my_cdchits[m]->hit);
      }
   }
   // output the pull information
   pulls.assign(last_cdc_pulls.begin(),last_cdc_pulls.end());

   // Track Parameters at "vertex"
   phi_=Sclast(state_phi);
   q_over_pt_=Sclast(state_q_over_pt);
   tanl_=Sclast(state_tanl);
   x_=last_pos.X();
   y_=last_pos.Y();
   z_=Sclast(state_z);  
   // Find the (signed) distance of closest approach to the beam line
   DVector2 beam_pos=beam_center+(z_-beam_z0)*beam_dir;
   double dx=x_-beam_pos.X();
   double dy=y_-beam_pos.Y();
   D_=sqrt(dx*dx+dy*dy)+EPS;
   x_ = dx; y_ = dy;
   double cosphi=cos(phi_);
   double sinphi=sin(phi_);
   if ((dx>0.0 && sinphi>0.0) || (dy <0.0 && cosphi>0.0) 
         || (dy>0.0 && cosphi<0.0) || (dx<0.0 && sinphi<0.0)) D_*=-1.; 
   // Rotate covariance matrix to coordinate system centered on x=0,y=0 in the 
   // lab 
   DMatrix5x5 Jc=I5x5;
   Jc(state_D,state_D)=(dy*cosphi-dx*sinphi)/D_;
   //Cclast=Cclast.SandwichMultiply(Jc);
   Cclast=Jc*Cclast*Jc.Transpose();

   if (!isfinite(x_) || !isfinite(y_) || !isfinite(z_) || !isfinite(phi_) 
         || !isfinite(q_over_pt_) || !isfinite(tanl_)){
      if (DEBUG_LEVEL>0){
         _DBG_ << "At least one parameter is NaN or +-inf!!" <<endl;
         _DBG_ << "x " << x_ << " y " << y_ << " z " << z_ << " phi " << phi_
            << " q/pt " << q_over_pt_ << " tanl " << tanl_ << endl;
      }
      return INVALID_FIT;	       
   }

   // Covariance matrix at vertex
   fcov.clear();
   vector<double>dummy;
   for (unsigned int i=0;i<5;i++){
      dummy.clear();
      for(unsigned int j=0;j<5;j++){
         dummy.push_back(Cclast(i,j));
      }
      cov.push_back(dummy);
   }

   // total chisq and ndf
   chisq_=chisq_iter;
   ndf_=last_ndf;
   //printf("NDof %d\n",ndf);

   return FIT_SUCCEEDED;
}

// Smoothing algorithm for the forward trajectory.  Updates the state vector
// at each step (going in the reverse direction to the filter) based on the 
// information from all the steps and outputs the state vector at the
// outermost step.
jerror_t DTrackFitterKalmanSIMD::SmoothForward(vector<pull_t>&forward_pulls){ 
   if (forward_traj.size()<2) return RESOURCE_UNAVAILABLE;

   unsigned int max=forward_traj.size()-1;
   DMatrix5x1 S=(forward_traj[max].Skk);
   DMatrix5x5 C=(forward_traj[max].Ckk);
   DMatrix5x5 JT=forward_traj[max].J.Transpose();
   DMatrix5x1 Ss=S;
   DMatrix5x5 Cs=C;
   DMatrix5x5 A,dC;

   if (DEBUG_LEVEL>19){
      jout << "---- Smoothed residuals ----" <<endl;
      jout << setprecision(4); 
   }

   for (unsigned int m=max-1;m>0;m--){
      if (forward_traj[m].h_id>0){
         if (forward_traj[m].h_id<1000){
            unsigned int id=forward_traj[m].h_id-1;
            if (DEBUG_LEVEL>1) _DBG_ << " Smoothing FDC ID " << id << endl;
            if (fdc_used_in_fit[id] && my_fdchits[id]->status==good_hit){
               if (DEBUG_LEVEL>1) _DBG_ << " Used in fit " << endl;
               A=fdc_updates[id].C*JT*C.InvertSym();
               Ss=fdc_updates[id].S+A*(Ss-S);

               if (!Ss.IsFinite()){ 
                  if (DEBUG_LEVEL>1) _DBG_ << "Invalid values for smoothed parameters..." << endl;
                  return VALUE_OUT_OF_RANGE;
               }
               dC=A*(Cs-C)*A.Transpose();
               Cs=fdc_updates[id].C+dC;
               if (!Cs.IsPosDef()){
                  if (DEBUG_LEVEL>1)
                     _DBG_ << "Covariance Matrix not PosDef..." << endl;
                  return VALUE_OUT_OF_RANGE;
               }

               // Position and direction from state vector with small angle
	       // alignment correction
               double x=Ss(state_x) + my_fdchits[id]->phiZ*Ss(state_y);
               double y=Ss(state_y) - my_fdchits[id]->phiZ*Ss(state_x);
               double tx=Ss(state_tx)+ my_fdchits[id]->phiZ*Ss(state_ty) 
		 - my_fdchits[id]->phiY;
               double ty=Ss(state_ty) - my_fdchits[id]->phiZ*Ss(state_tx) 
		 + my_fdchits[id]->phiX;

               double cosa=my_fdchits[id]->cosa;
               double sina=my_fdchits[id]->sina;
               double u=my_fdchits[id]->uwire;
               double v=my_fdchits[id]->vstrip;

               // Projected position along the wire without doca-dependent corrections
               double vpred_uncorrected=x*sina+y*cosa;

               // Projected position in the plane of the wires transverse to the wires
               double upred=x*cosa-y*sina;

               // Direction tangent in the u-z plane
               double tu=tx*cosa-ty*sina;
               double one_plus_tu2=1.+tu*tu;
               double alpha=atan(tu);
               double cosalpha=cos(alpha);
               //double cosalpha2=cosalpha*cosalpha;
               double sinalpha=sin(alpha);

               // (signed) distance of closest approach to wire
               double du=upred-u;
               double doca=du*cosalpha;

               // Correction for lorentz effect
               double nz=my_fdchits[id]->nz;
               double nr=my_fdchits[id]->nr;
               double nz_sinalpha_plus_nr_cosalpha=nz*sinalpha+nr*cosalpha;

               // Difference between measurement and projection
               double tv=tx*sina+ty*cosa;
               double resi=v-(vpred_uncorrected+doca*(nz_sinalpha_plus_nr_cosalpha
                        -tv*sinalpha));	
               double drift_time=my_fdchits[id]->t-mT0
                  -forward_traj[m].t*TIME_UNIT_CONVERSION;
               double drift = 0.0;
               if (USE_FDC_DRIFT_TIMES){
                  drift=(du>0.0?1.:-1.)*fdc_drift_distance(drift_time,forward_traj[m].B);
               }

               double resi_a=drift-doca;

               // Variance from filter step
               // This V is really "R" in Fruhwirths notation, in the case that the track is used in the fit.
               DMatrix2x2 V=fdc_updates[id].V;
               // Compute projection matrix and find the variance for the residual
               DMatrix5x2 H_T;
               double temp2=nz_sinalpha_plus_nr_cosalpha-tv*sinalpha;
               H_T(state_x,1)=sina+cosa*cosalpha*temp2;
               H_T(state_y,1)=cosa-sina*cosalpha*temp2;

               double cos2_minus_sin2=cosalpha*cosalpha-sinalpha*sinalpha;
               double fac=nz*cos2_minus_sin2-2.*nr*cosalpha*sinalpha;
               double doca_cosalpha=doca*cosalpha;
               double temp=doca_cosalpha*fac;
               H_T(state_tx,1)=cosa*temp
                  -doca_cosalpha*(tu*sina+tv*cosa*cos2_minus_sin2)
                  ;
               H_T(state_ty,1)=-sina*temp
                  -doca_cosalpha*(tu*cosa-tv*sina*cos2_minus_sin2)
                  ;

               H_T(state_x,0)=cosa*cosalpha;
               H_T(state_y,0)=-sina*cosalpha;
               double factor=du*tu/sqrt(one_plus_tu2)/one_plus_tu2;
               H_T(state_ty,0)=sina*factor;
               H_T(state_tx,0)=-cosa*factor;

               // Matrix transpose H_T -> H
               DMatrix2x5 H=Transpose(H_T);

               if (my_fdchits[id]->hit->wire->layer==PLANE_TO_SKIP){
                  //V+=Cs.SandwichMultiply(H_T);
                  V=V+H*Cs*H_T;
               }
               else{
                  //V-=dC.SandwichMultiply(H_T);
                  V=V-H*dC*H_T;
               }


               vector<double> alignmentDerivatives;
               if (ALIGNMENT_FORWARD){
                  alignmentDerivatives.resize(FDCTrackD::size);
                  // Let's get the alignment derivatives

                  // Things are assumed to be linear near the wire, derivatives can be determined analytically.
                  // First for the wires

                  //dDOCAW/ddeltax
                  alignmentDerivatives[FDCTrackD::dDOCAW_dDeltaX] = -(1/sqrt(1 + pow(tx*cosa - ty*sina,2)));

                  //dDOCAW/ddeltaz
                  alignmentDerivatives[FDCTrackD::dDOCAW_dDeltaZ] = (tx*cosa - ty*sina)/sqrt(1 + pow(tx*cosa - ty*sina,2));

                  //dDOCAW/ddeltaPhiX
                  double cos2a = cos(2.*my_fdchits[id]->hit->wire->angle);
                  double sin2a = sin(2.*my_fdchits[id]->hit->wire->angle);
                  double cos3a = cos(3.*my_fdchits[id]->hit->wire->angle);
                  double sin3a = sin(3.*my_fdchits[id]->hit->wire->angle);
                  //double tx2 = tx*tx;
                  //double ty2 = ty*ty;
                  alignmentDerivatives[FDCTrackD::dDOCAW_dDeltaPhiX] = (sina*(-(tx*cosa) + ty*sina)*(u - x*cosa + y*sina))/
                     pow(1 + pow(tx*cosa - ty*sina,2),1.5);
                  alignmentDerivatives[FDCTrackD::dDOCAW_dDeltaPhiY] = -((cosa*(tx*cosa - ty*sina)*(u - x*cosa + y*sina))/
                        pow(1 + pow(tx*cosa - ty*sina,2),1.5));
                  alignmentDerivatives[FDCTrackD::dDOCAW_dDeltaPhiZ] = (tx*ty*u*cos2a + (x + pow(ty,2)*x - tx*ty*y)*sina + 
                        cosa*(-(tx*ty*x) + y + pow(tx,2)*y + (tx - ty)*(tx + ty)*u*sina))/
                     pow(1 + pow(tx*cosa - ty*sina,2),1.5);

                  // dDOCAW/dt0
                  double t0shift=4.;//ns
                  double drift_shift = 0.0;
                  if(USE_FDC_DRIFT_TIMES){
                     if (drift_time < 0.) drift_shift = drift;
                     else drift_shift = (du>0.0?1.:-1.)*fdc_drift_distance(drift_time+t0shift,forward_traj[m].B);
                  }
                  alignmentDerivatives[FDCTrackD::dW_dt0]= (drift_shift-drift)/t0shift;

                  // dDOCAW/dx
                  alignmentDerivatives[FDCTrackD::dDOCAW_dx] = cosa/sqrt(1 + pow(tx*cosa - ty*sina,2));

                  // dDOCAW/dy
                  alignmentDerivatives[FDCTrackD::dDOCAW_dy] = -(sina/sqrt(1 + pow(tx*cosa - ty*sina,2)));

                  // dDOCAW/dtx
                  alignmentDerivatives[FDCTrackD::dDOCAW_dtx] = (cosa*(tx*cosa - ty*sina)*(u - x*cosa + y*sina))/pow(1 + pow(tx*cosa - ty*sina,2),1.5);

                  // dDOCAW/dty
                  alignmentDerivatives[FDCTrackD::dDOCAW_dty] = (sina*(-(tx*cosa) + ty*sina)*(u - x*cosa + y*sina))/
                     pow(1 + pow(tx*cosa - ty*sina,2),1.5);

                  // Then for the cathodes. The magnetic field correction now correlates the alignment constants for the wires and cathodes.

                  //dDOCAC/ddeltax
                  alignmentDerivatives[FDCTrackD::dDOCAC_dDeltaX] =
                     (-nr + (-nz + ty*cosa + tx*sina)*(tx*cosa - ty*sina))/(1 + pow(tx*cosa - ty*sina,2));

                  //dDOCAC/ddeltaz
                  alignmentDerivatives[FDCTrackD::dDOCAC_dDeltaZ] =
                     nz + (-nz + (nr*tx + ty)*cosa + (tx - nr*ty)*sina)/(1 + pow(tx*cosa - ty*sina,2));

                  //dDOCAC/ddeltaPhiX
                  alignmentDerivatives[FDCTrackD::dDOCAC_dDeltaPhiX] =
                     (-2*y*cosa*sina*(tx*cosa - ty*sina) - 2*x*pow(sina,2)*(tx*cosa - ty*sina) - 
                      (u - x*cosa + y*sina)*(-(nz*sina) + sina*(ty*cosa + tx*sina) - 
                         cosa*(tx*cosa - ty*sina)))/(1 + pow(tx*cosa - ty*sina,2)) + 
                     (2*sina*(tx*cosa - ty*sina)*(-((u - x*cosa + y*sina)*
                                                    (nr + nz*(tx*cosa - ty*sina) - (ty*cosa + tx*sina)*(tx*cosa - ty*sina))) + 
                                                  y*cosa*(1 + pow(tx*cosa - ty*sina,2)) + x*sina*(1 + pow(tx*cosa - ty*sina,2))))/
                     pow(1 + pow(tx*cosa - ty*sina,2),2); 


                  //dDOCAC/ddeltaPhiY
                  alignmentDerivatives[FDCTrackD::dDOCAC_dDeltaPhiY] = (-2*y*pow(cosa,2)*(tx*cosa - ty*sina) - 2*x*cosa*sina*(tx*cosa - ty*sina) - 
                        (u - x*cosa + y*sina)*(-(nz*cosa) + cosa*(ty*cosa + tx*sina) + 
                           sina*(tx*cosa - ty*sina)))/(1 + pow(tx*cosa - ty*sina,2)) + 
                     (2*cosa*(tx*cosa - ty*sina)*(-((u - x*cosa + y*sina)*
                                                    (nr + nz*(tx*cosa - ty*sina) - (ty*cosa + tx*sina)*(tx*cosa - ty*sina))) + 
                                                  y*cosa*(1 + pow(tx*cosa - ty*sina,2)) + x*sina*(1 + pow(tx*cosa - ty*sina,2))))/
                     pow(1 + pow(tx*cosa - ty*sina,2),2);

                  //dDOCAC/ddeltaPhiZ
                  alignmentDerivatives[FDCTrackD::dDOCAC_dDeltaPhiZ] = (-2*(ty*cosa + tx*sina)*(tx*cosa - ty*sina)*
                        (-((u - x*cosa + y*sina)*(nr + nz*(tx*cosa - ty*sina) - 
                                                  (ty*cosa + tx*sina)*(tx*cosa - ty*sina))) + 
                         y*cosa*(1 + pow(tx*cosa - ty*sina,2)) + x*sina*(1 + pow(tx*cosa - ty*sina,2))))/
                     pow(1 + pow(tx*cosa - ty*sina,2),2) + 
                     (2*y*cosa*(ty*cosa + tx*sina)*(tx*cosa - ty*sina) + 
                      2*x*sina*(ty*cosa + tx*sina)*(tx*cosa - ty*sina) - 
                      (-(y*cosa) - x*sina)*(nr + nz*(tx*cosa - ty*sina) - 
                         (ty*cosa + tx*sina)*(tx*cosa - ty*sina)) - 
                      x*cosa*(1 + pow(tx*cosa - ty*sina,2)) + y*sina*(1 + pow(tx*cosa - ty*sina,2)) - 
                      (u - x*cosa + y*sina)*(nz*(ty*cosa + tx*sina) - pow(ty*cosa + tx*sina,2) - 
                         (tx*cosa - ty*sina)*(-(tx*cosa) + ty*sina)))/(1 + pow(tx*cosa - ty*sina,2));

                  //dDOCAC/dx
                  alignmentDerivatives[FDCTrackD::dDOCAC_dx] = (cosa*(nr - tx*ty + nz*tx*cosa) + sina + ty*(ty - nz*cosa)*sina)/
                     (1 + pow(tx*cosa - ty*sina,2));

                  //dDOCAC/dy
                  alignmentDerivatives[FDCTrackD::dDOCAC_dy] = ((1 + pow(tx,2))*cosa - (nr + tx*ty + nz*tx*cosa)*sina + nz*ty*pow(sina,2))/
                     (1 + pow(tx*cosa - ty*sina,2));

                  //dDOCAC/dtx
                  alignmentDerivatives[FDCTrackD::dDOCAC_dtx] = ((u - x*cosa + y*sina)*(4*nr*tx - 2*ty*(pow(tx,2) + pow(ty,2)) + nz*(-4 + 3*pow(tx,2) + pow(ty,2))*cosa + 
                           2*(2*nr*tx + ty*(2 - pow(tx,2) + pow(ty,2)))*cos2a + nz*(tx - ty)*(tx + ty)*cos3a - 2*nz*tx*ty*sina + 
                           4*(tx - nr*ty + tx*pow(ty,2))*sin2a - 2*nz*tx*ty*sin3a))/
                     pow(2 + pow(tx,2) + pow(ty,2) + (tx - ty)*(tx + ty)*cos2a - 2*tx*ty*sin2a,2);

                  //dDOCAC/dty
                  alignmentDerivatives[FDCTrackD::dDOCAC_dty] = -(((u - x*cosa + y*sina)*(-2*(pow(tx,3) + 2*nr*ty + tx*pow(ty,2)) - 2*nz*tx*ty*cosa - 
                              2*(2*tx + pow(tx,3) - 2*nr*ty - tx*pow(ty,2))*cos2a + 2*nz*tx*ty*cos3a + 
                              nz*(-4 + pow(tx,2) + 3*pow(ty,2))*sina + 4*(ty + tx*(nr + tx*ty))*sin2a + nz*(tx - ty)*(tx + ty)*sin3a))
                        /pow(2 + pow(tx,2) + pow(ty,2) + (tx - ty)*(tx + ty)*cos2a - 2*tx*ty*sin2a,2));

               }

               if (DEBUG_LEVEL>19){
                  jout << "Layer " << my_fdchits[id]->hit->wire->layer
                     <<":   t " << drift_time << " x "<< x << " y " << y 
                     << " coordinate along wire " << v << " resi_c " <<resi
                     << " coordinate transverse to wire " << drift 
                     <<" resi_a " << resi_a
                     <<endl;
               }

	       double scale=1./sqrt(1.+tx*tx+ty*ty);
	       double cosThetaRel=my_fdchits[id]->hit->wire->udir.Dot(DVector3(scale*tx,scale*ty,scale));
               DTrackFitter::pull_t thisPull = pull_t(resi_a,sqrt(V(0,0)),
						      forward_traj[m].s,
						      fdc_updates[id].tdrift,
						      fdc_updates[id].doca,
						      NULL,my_fdchits[id]->hit,
						      0.,
						      forward_traj[m].z,
						      cosThetaRel,0.,
						      resi,sqrt(V(1,1)));
               thisPull.AddTrackDerivatives(alignmentDerivatives);
               forward_pulls.push_back(thisPull);
            }
            else{
               A=forward_traj[m].Ckk*JT*C.InvertSym();
               Ss=forward_traj[m].Skk+A*(Ss-S);
               Cs=forward_traj[m].Ckk+A*(Cs-C)*A.Transpose();
            }

         }
         else{
            unsigned int id=forward_traj[m].h_id-1000;
            if (DEBUG_LEVEL>1) _DBG_ << " Smoothing CDC ID " << id << endl;
            if (cdc_used_in_fit[id]&&my_cdchits[id]->status==good_hit){
               if (DEBUG_LEVEL>1) _DBG_ << " Used in fit " << endl;
               A=cdc_updates[id].C*JT*C.InvertSym();
               Ss=cdc_updates[id].S+A*(Ss-S);
               Cs=cdc_updates[id].C+A*(Cs-C)*A.Transpose();
               if (!Cs.IsPosDef()){
                  if (DEBUG_LEVEL>1)
                     _DBG_ << "Covariance Matrix not PosDef..." << endl;
                  return VALUE_OUT_OF_RANGE;
               } 
               if (!Ss.IsFinite()){
                  if (DEBUG_LEVEL>5) _DBG_ << "Invalid values for smoothed parameters..." << endl;
                  return VALUE_OUT_OF_RANGE;
               }

               // Fill in pulls information for cdc hits
               if(FillPullsVectorEntry(Ss,Cs,forward_traj[m],my_cdchits[id],
				       cdc_updates[id],forward_pulls) != NOERROR) return VALUE_OUT_OF_RANGE;
            }
            else{
               A=forward_traj[m].Ckk*JT*C.InvertSym();
               Ss=forward_traj[m].Skk+A*(Ss-S);
               Cs=forward_traj[m].Ckk+A*(Cs-C)*A.Transpose();
            }
         }
      }
      else{
         A=forward_traj[m].Ckk*JT*C.InvertSym();
         Ss=forward_traj[m].Skk+A*(Ss-S);
         Cs=forward_traj[m].Ckk+A*(Cs-C)*A.Transpose();
      }

      S=forward_traj[m].Skk;
      C=forward_traj[m].Ckk;
      JT=forward_traj[m].J.Transpose();
   }

   return NOERROR;
}

// at each step (going in the reverse direction to the filter) based on the 
// information from all the steps.
jerror_t DTrackFitterKalmanSIMD::SmoothCentral(vector<pull_t>&cdc_pulls){ 
   if (central_traj.size()<2) return RESOURCE_UNAVAILABLE;

   unsigned int max = central_traj.size()-1;
   DMatrix5x1 S=(central_traj[max].Skk);
   DMatrix5x5 C=(central_traj[max].Ckk);
   DMatrix5x5 JT=central_traj[max].J.Transpose();
   DMatrix5x1 Ss=S;
   DMatrix5x5 Cs=C;
   DMatrix5x5 A,dC;

   if (DEBUG_LEVEL>1) {
      _DBG_ << " S C JT at start of smoothing " << endl;
      S.Print(); C.Print(); JT.Print();
   }

   for (unsigned int m=max-1;m>0;m--){    
     if (central_traj[m].h_id>0){
         unsigned int id=central_traj[m].h_id-1;
         if (DEBUG_LEVEL>1) _DBG_ << " Encountered Hit ID " << id << " At trajectory position " << m << "/" << max << endl;
         if (cdc_used_in_fit[id] && my_cdchits[id]->status == good_hit){
            if (DEBUG_LEVEL>1) _DBG_ << " SmoothCentral CDC Hit ID " << id << " used in fit " << endl;

            A=cdc_updates[id].C*JT*C.InvertSym();
	    dC=Cs-C;
            Ss=cdc_updates[id].S+A*(Ss-S);
	    Cs=cdc_updates[id].C+A*dC*A.Transpose();
	    
            if (!Ss.IsFinite()){
               if (DEBUG_LEVEL>5) 
                  _DBG_ << "Invalid values for smoothed parameters..." << endl;
               return VALUE_OUT_OF_RANGE;
            }
            if (!Cs.IsPosDef()){
               if (DEBUG_LEVEL>5){
		 _DBG_ << "Covariance Matrix not PosDef... Ckk dC A" << endl;
		 cdc_updates[id].C.Print(); dC.Print(); A.Print();
               }
               return VALUE_OUT_OF_RANGE;
            }

            // Get estimate for energy loss 
            double q_over_p=Ss(state_q_over_pt)*cos(atan(Ss(state_tanl)));
            double dEdx=GetdEdx(q_over_p,central_traj[m].K_rho_Z_over_A,
                  central_traj[m].rho_Z_over_A,
                  central_traj[m].LnI,central_traj[m].Z);

            // Use Brent's algorithm to find doca to the wire
            DVector2 xy(central_traj[m].xy.X()-Ss(state_D)*sin(Ss(state_phi)),
                  central_traj[m].xy.Y()+Ss(state_D)*cos(Ss(state_phi)));
            DVector2 old_xy=xy;
            DMatrix5x1 myS=Ss;
            double myds;  
            DVector2 origin=my_cdchits[id]->origin;
            DVector2 dir=my_cdchits[id]->dir;
            double z0wire=my_cdchits[id]->z0wire;
            //BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dEdx,xy,z0wire,origin,dir,myS,myds);
            if(BrentCentral(dEdx,xy,z0wire,origin,dir,myS,myds)!=NOERROR) return VALUE_OUT_OF_RANGE;
            if(DEBUG_HISTS) alignDerivHists[0]->Fill(myds);
            DVector2 wirepos=origin+(myS(state_z)-z0wire)*dir;
            double cosstereo=my_cdchits[id]->cosstereo;
            DVector2 diff=xy-wirepos;
            // here we add a small number to avoid division by zero errors
            double d=cosstereo*diff.Mod()+EPS; 

            // If we are doing the alignment, we need to numerically calculate the derivatives
            // wrt the wire origin, direction, and the track parameters.
            vector<double> alignmentDerivatives;
            if (ALIGNMENT_CENTRAL){
               double dscut_min=0., dscut_max=1.;
               DVector3 wireDir = my_cdchits[id]->hit->wire->udir;
               double cosstereo_shifted;
               DMatrix5x1 alignS=Ss; // We will mess with this one
               double alignds;
               alignmentDerivatives.resize(12);
               alignmentDerivatives[CDCTrackD::dDdt0]=cdc_updates[id].dDdt0;
               // Wire position shift
               double wposShift=0.025; 
               double wdirShift=0.00005;  

               // Shift each track parameter value 
               double shiftFraction=0.01;
               double shift_q_over_pt=shiftFraction*Ss(state_q_over_pt);
               double shift_phi=0.0001;
               double shift_tanl=shiftFraction*Ss(state_tanl);
               double shift_D=0.01;
               double shift_z=0.01;

               // Some data containers we don't need multiples of
               double z0_shifted;
               DVector2 shift, origin_shifted, dir_shifted, wirepos_shifted, diff_shifted, xy_shifted;

               // The DOCA for the shifted states == f(x+h)
               double d_dOriginX=0., d_dOriginY=0., d_dOriginZ=0.;
               double d_dDirX=0., d_dDirY=0., d_dDirZ=0.;
               double d_dS0=0., d_dS1=0., d_dS2=0., d_dS3=0., d_dS4=0.;
               // Let's do the wire shifts first

               //dOriginX
               alignS=Ss;
               alignds=0.;
               shift.Set(wposShift, 0.);
               origin_shifted=origin+shift;
               dir_shifted=dir;
               z0_shifted=z0wire;
               xy_shifted.Set(central_traj[m].xy.X()-alignS(state_D)*sin(alignS(state_phi)),
                     central_traj[m].xy.Y()+alignS(state_D)*cos(alignS(state_phi)));
               if (BrentCentral(dEdx,xy_shifted,z0_shifted,origin_shifted,
                        dir_shifted,alignS,alignds)!=NOERROR) return VALUE_OUT_OF_RANGE;
               if (alignds < dscut_min || alignds > dscut_max) return VALUE_OUT_OF_RANGE;
               //if (BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dEdx,xy_shifted,z0_shifted,origin_shifted,
               //         dir_shifted,alignS,alignds)!=NOERROR) return VALUE_OUT_OF_RANGE;
               wirepos_shifted=origin_shifted+(alignS(state_z)-z0_shifted)*dir_shifted;
               diff_shifted=xy_shifted-wirepos_shifted;
               d_dOriginX=cosstereo*diff_shifted.Mod()+EPS;
               alignmentDerivatives[CDCTrackD::dDOCAdOriginX] = (d_dOriginX - d)/wposShift;
               if(DEBUG_HISTS){
                  alignDerivHists[12]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdOriginX]);
                  alignDerivHists[1]->Fill(alignds);
                  brentCheckHists[1]->Fill(alignds,d_dOriginX);
               }

               //dOriginY
               alignS=Ss;
               alignds=0.;
               shift.Set(0.,wposShift);
               origin_shifted=origin+shift;
               dir_shifted=dir;
               z0_shifted=z0wire;
               xy_shifted.Set(central_traj[m].xy.X()-alignS(state_D)*sin(alignS(state_phi)),
                     central_traj[m].xy.Y()+alignS(state_D)*cos(alignS(state_phi)));
               if (BrentCentral(dEdx,xy_shifted,z0_shifted,origin_shifted,
                        dir_shifted,alignS,alignds)!=NOERROR) return VALUE_OUT_OF_RANGE;
               if (alignds < dscut_min || alignds > dscut_max) return VALUE_OUT_OF_RANGE;
               //if(BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dEdx,xy_shifted,z0_shifted,origin_shifted,
               //         dir_shifted,alignS,alignds) != NOERROR) return VALUE_OUT_OF_RANGE;
               wirepos_shifted=origin_shifted+(alignS(state_z)-z0_shifted)*dir_shifted;
               diff_shifted=xy_shifted-wirepos_shifted;
               d_dOriginY=cosstereo*diff_shifted.Mod()+EPS;
               alignmentDerivatives[CDCTrackD::dDOCAdOriginY] = (d_dOriginY - d)/wposShift;
               if(DEBUG_HISTS){
                  alignDerivHists[13]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdOriginY]);
                  alignDerivHists[2]->Fill(alignds);
                  brentCheckHists[1]->Fill(alignds,d_dOriginY);
               }

               //dOriginZ
               alignS=Ss;
               alignds=0.;
               origin_shifted=origin;
               dir_shifted=dir;
               z0_shifted=z0wire+wposShift;
               xy_shifted.Set(central_traj[m].xy.X()-alignS(state_D)*sin(alignS(state_phi)),
                     central_traj[m].xy.Y()+alignS(state_D)*cos(alignS(state_phi)));
               if (BrentCentral(dEdx,xy_shifted,z0_shifted,origin_shifted,
                        dir_shifted,alignS,alignds)!=NOERROR) return VALUE_OUT_OF_RANGE;
               if (alignds < dscut_min || alignds > dscut_max) return VALUE_OUT_OF_RANGE;
               //if(BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dEdx,xy_shifted,z0_shifted,origin_shifted,
               //         dir_shifted,alignS,alignds) != NOERROR) return VALUE_OUT_OF_RANGE;
               wirepos_shifted=origin_shifted+(alignS(state_z)-z0_shifted)*dir_shifted;
               diff_shifted=xy_shifted-wirepos_shifted;
               d_dOriginZ=cosstereo*diff_shifted.Mod()+EPS;
               alignmentDerivatives[CDCTrackD::dDOCAdOriginZ] = (d_dOriginZ - d)/wposShift;
               if(DEBUG_HISTS){
                  alignDerivHists[14]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdOriginZ]);
                  alignDerivHists[3]->Fill(alignds);
                  brentCheckHists[1]->Fill(alignds,d_dOriginZ);
               }

               //dDirX
               alignS=Ss;
               alignds=0.;
               shift.Set(wdirShift,0.);
               origin_shifted=origin;
               z0_shifted=z0wire;
               xy_shifted.Set(central_traj[m].xy.X()-alignS(state_D)*sin(alignS(state_phi)),
                     central_traj[m].xy.Y()+alignS(state_D)*cos(alignS(state_phi)));
               dir_shifted=dir+shift;
               cosstereo_shifted = cos((wireDir+DVector3(wdirShift,0.,0.)).Angle(DVector3(0.,0.,1.)));
               if (BrentCentral(dEdx,xy_shifted,z0_shifted,origin_shifted,
                        dir_shifted,alignS,alignds)!=NOERROR) return VALUE_OUT_OF_RANGE;
               if (alignds < dscut_min || alignds > dscut_max) return VALUE_OUT_OF_RANGE;
               //if(BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dEdx,xy_shifted,z0_shifted,origin_shifted,
               //         dir_shifted,alignS,alignds) != NOERROR) return VALUE_OUT_OF_RANGE;
               wirepos_shifted=origin_shifted+(alignS(state_z)-z0_shifted)*dir_shifted;
               diff_shifted=xy_shifted-wirepos_shifted;
               d_dDirX=cosstereo_shifted*diff_shifted.Mod()+EPS;
               alignmentDerivatives[CDCTrackD::dDOCAdDirX] = (d_dDirX - d)/wdirShift;
               if(DEBUG_HISTS){
                  alignDerivHists[15]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdDirX]);
                  alignDerivHists[4]->Fill(alignds);
               }

               //dDirY
               alignS=Ss;
               alignds=0.;
               shift.Set(0.,wdirShift);
               origin_shifted=origin;
               z0_shifted=z0wire;
               xy_shifted.Set(central_traj[m].xy.X()-alignS(state_D)*sin(alignS(state_phi)),
                     central_traj[m].xy.Y()+alignS(state_D)*cos(alignS(state_phi)));
               dir_shifted=dir+shift;
               cosstereo_shifted = cos((wireDir+DVector3(0.,wdirShift,0.)).Angle(DVector3(0.,0.,1.)));
               if (BrentCentral(dEdx,xy_shifted,z0_shifted,origin_shifted,
                        dir_shifted,alignS,alignds)!=NOERROR) return VALUE_OUT_OF_RANGE;
               if (alignds < dscut_min || alignds > dscut_max) return VALUE_OUT_OF_RANGE;
               //if(BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dEdx,xy_shifted,z0_shifted,origin_shifted,
               //         dir_shifted,alignS,alignds) != NOERROR) return VALUE_OUT_OF_RANGE;
               wirepos_shifted=origin_shifted+(alignS(state_z)-z0_shifted)*dir_shifted;
               diff_shifted=xy_shifted-wirepos_shifted;
               d_dDirY=cosstereo_shifted*diff_shifted.Mod()+EPS;
               alignmentDerivatives[CDCTrackD::dDOCAdDirY] = (d_dDirY - d)/wdirShift;
               if(DEBUG_HISTS){
                  alignDerivHists[16]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdDirY]);
                  alignDerivHists[5]->Fill(alignds);
               }

               //dDirZ
               alignS=Ss;
               alignds=0.;
               origin_shifted=origin;
               dir_shifted.Set(wireDir.X()/(wireDir.Z()+wdirShift), wireDir.Y()/(wireDir.Z()+wdirShift));
               z0_shifted=z0wire;
               xy_shifted.Set(central_traj[m].xy.X()-alignS(state_D)*sin(alignS(state_phi)),
                     central_traj[m].xy.Y()+alignS(state_D)*cos(alignS(state_phi)));
               cosstereo_shifted = cos((wireDir+DVector3(0.,0.,wdirShift)).Angle(DVector3(0.,0.,1.)));
               if (BrentCentral(dEdx,xy_shifted,z0_shifted,origin_shifted,
                        dir_shifted,alignS,alignds)!=NOERROR) return VALUE_OUT_OF_RANGE;
               if (alignds < dscut_min || alignds > dscut_max) return VALUE_OUT_OF_RANGE;
               //if(BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dEdx,xy_shifted,z0_shifted,origin_shifted,
               //         dir_shifted,alignS,alignds)!=NOERROR) return VALUE_OUT_OF_RANGE;
               wirepos_shifted=origin_shifted+(alignS(state_z)-z0_shifted)*dir_shifted;
               diff_shifted=xy_shifted-wirepos_shifted;
               d_dDirZ=cosstereo_shifted*diff_shifted.Mod()+EPS;
               alignmentDerivatives[CDCTrackD::dDOCAdDirZ] = (d_dDirZ - d)/wdirShift;
               if(DEBUG_HISTS){
                  alignDerivHists[17]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdDirZ]);
                  alignDerivHists[6]->Fill(alignds);
               }

               // And now the derivatives wrt the track parameters
               //DMatrix5x1 trackShift(shift_q_over_pt, shift_phi, shift_tanl, shift_D, shift_z);

               DMatrix5x1 trackShiftS0(shift_q_over_pt, 0., 0., 0., 0.);
               DMatrix5x1 trackShiftS1(0., shift_phi, 0., 0., 0.);
               DMatrix5x1 trackShiftS2(0., 0., shift_tanl, 0., 0.);
               DMatrix5x1 trackShiftS3(0., 0., 0., shift_D, 0.);
               DMatrix5x1 trackShiftS4(0., 0., 0., 0., shift_z);

               // dS0
               alignS=Ss+trackShiftS0;
               alignds=0.;
               xy_shifted.Set(central_traj[m].xy.X()-alignS(state_D)*sin(alignS(state_phi)),
                     central_traj[m].xy.Y()+alignS(state_D)*cos(alignS(state_phi)));
               if(BrentCentral(dEdx,xy_shifted,z0wire,origin,dir,alignS,alignds) != NOERROR) return VALUE_OUT_OF_RANGE;
               if (alignds < dscut_min || alignds > dscut_max) return VALUE_OUT_OF_RANGE;
               //if(BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dEdx,xy_shifted,z0wire,origin,dir,alignS,alignds) != NOERROR) return VALUE_OUT_OF_RANGE;
               wirepos_shifted=origin+(alignS(state_z)-z0wire)*dir;
               diff_shifted=xy_shifted-wirepos_shifted;
               d_dS0=cosstereo*diff_shifted.Mod()+EPS;
               alignmentDerivatives[CDCTrackD::dDOCAdS0] = (d_dS0 - d)/shift_q_over_pt;
               if(DEBUG_HISTS){
                  alignDerivHists[18]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdS0]);
                  alignDerivHists[7]->Fill(alignds);
               }

               // dS1
               alignS=Ss+trackShiftS1;
               alignds=0.;
               xy_shifted.Set(central_traj[m].xy.X()-alignS(state_D)*sin(alignS(state_phi)),
                     central_traj[m].xy.Y()+alignS(state_D)*cos(alignS(state_phi)));
               if(BrentCentral(dEdx,xy_shifted,z0wire,origin,dir,alignS,alignds) != NOERROR) return VALUE_OUT_OF_RANGE;
               if (alignds < dscut_min || alignds > dscut_max) return VALUE_OUT_OF_RANGE;
               //if(BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dEdx,xy_shifted,z0wire,origin,dir,alignS,alignds) != NOERROR) return VALUE_OUT_OF_RANGE;
               wirepos_shifted=origin+(alignS(state_z)-z0wire)*dir;
               diff_shifted=xy_shifted-wirepos_shifted;
               d_dS1=cosstereo*diff_shifted.Mod()+EPS;
               alignmentDerivatives[CDCTrackD::dDOCAdS1] = (d_dS1 - d)/shift_phi;
               if(DEBUG_HISTS){
                  alignDerivHists[19]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdS1]);
                  alignDerivHists[8]->Fill(alignds);
               }

               // dS2
               alignS=Ss+trackShiftS2;
               alignds=0.;
               xy_shifted.Set(central_traj[m].xy.X()-alignS(state_D)*sin(alignS(state_phi)),
                     central_traj[m].xy.Y()+alignS(state_D)*cos(alignS(state_phi)));
               if(BrentCentral(dEdx,xy_shifted,z0wire,origin,dir,alignS,alignds) != NOERROR) return VALUE_OUT_OF_RANGE;
               if (alignds < dscut_min || alignds > dscut_max) return VALUE_OUT_OF_RANGE;
               //if(BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dEdx,xy_shifted,z0wire,origin,dir,alignS,alignds) != NOERROR) return VALUE_OUT_OF_RANGE;
               wirepos_shifted=origin+(alignS(state_z)-z0wire)*dir;
               diff_shifted=xy_shifted-wirepos_shifted;
               d_dS2=cosstereo*diff_shifted.Mod()+EPS;
               alignmentDerivatives[CDCTrackD::dDOCAdS2] = (d_dS2 - d)/shift_tanl;
               if(DEBUG_HISTS){
                  alignDerivHists[20]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdS2]);
                  alignDerivHists[9]->Fill(alignds);
               }

               // dS3
               alignS=Ss+trackShiftS3;
               alignds=0.;
               xy_shifted.Set(central_traj[m].xy.X()-alignS(state_D)*sin(alignS(state_phi)),
                     central_traj[m].xy.Y()+alignS(state_D)*cos(alignS(state_phi)));
               if(BrentCentral(dEdx,xy_shifted,z0wire,origin,dir,alignS,alignds) != NOERROR) return VALUE_OUT_OF_RANGE;
               if (alignds < dscut_min || alignds > dscut_max) return VALUE_OUT_OF_RANGE;
               //if(BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dEdx,xy_shifted,z0wire,origin,dir,alignS,alignds) != NOERROR) return VALUE_OUT_OF_RANGE;
               wirepos_shifted=origin+(alignS(state_z)-z0wire)*dir;
               diff_shifted=xy_shifted-wirepos_shifted;
               d_dS3=cosstereo*diff_shifted.Mod()+EPS;
               alignmentDerivatives[CDCTrackD::dDOCAdS3] = (d_dS3 - d)/shift_D;
               if(DEBUG_HISTS){
                  alignDerivHists[21]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdS3]);
                  alignDerivHists[10]->Fill(alignds);
               }

               // dS4
               alignS=Ss+trackShiftS4;
               alignds=0.;
               xy_shifted.Set(central_traj[m].xy.X()-alignS(state_D)*sin(alignS(state_phi)),
                     central_traj[m].xy.Y()+alignS(state_D)*cos(alignS(state_phi)));
               if(BrentCentral(dEdx,xy_shifted,z0wire,origin,dir,alignS,alignds) != NOERROR) return VALUE_OUT_OF_RANGE;
               if (alignds < dscut_min || alignds > dscut_max) return VALUE_OUT_OF_RANGE;
               //if(BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dEdx,xy_shifted,z0wire,origin,dir,alignS,alignds) != NOERROR) return VALUE_OUT_OF_RANGE;
               wirepos_shifted=origin+(alignS(state_z)-z0wire)*dir;
               diff_shifted=xy_shifted-wirepos_shifted;
               d_dS4=cosstereo*diff_shifted.Mod()+EPS;
               alignmentDerivatives[CDCTrackD::dDOCAdS4] = (d_dS4 - d)/shift_z;
               if(DEBUG_HISTS){
                  alignDerivHists[22]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdS4]);
                  alignDerivHists[11]->Fill(alignds);
               }
            }

            // Compute the Jacobian matrix
            // Find the field and gradient at (old_x,old_y,old_z)
            bfield->GetFieldAndGradient(old_xy.X(),old_xy.Y(),Ss(state_z),
                  Bx,By,Bz,
                  dBxdx,dBxdy,dBxdz,dBydx,
                  dBydy,dBydz,dBzdx,dBzdy,dBzdz);
            DMatrix5x5 Jc;
            StepJacobian(old_xy,xy-old_xy,myds,Ss,dEdx,Jc);

            // Projection matrix        
            DMatrix5x1 H_T;
            double sinphi=sin(myS(state_phi));
            double cosphi=cos(myS(state_phi));
            double dx=diff.X();
            double dy=diff.Y();
            double cosstereo2_over_doca=cosstereo*cosstereo/d;
            H_T(state_D)=(dy*cosphi-dx*sinphi)*cosstereo2_over_doca;
            H_T(state_phi)
               =-myS(state_D)*cosstereo2_over_doca*(dx*cosphi+dy*sinphi);
            H_T(state_z)=-cosstereo2_over_doca*(dx*dir.X()+dy*dir.Y());
	    DMatrix1x5 H;
	    H(state_D)=H_T(state_D);
	    H(state_phi)=H_T(state_phi);
	    H(state_z)=H_T(state_z);	    

            double Vhit=cdc_updates[id].variance;
	    Cs=Jc*Cs*Jc.Transpose();
            //double Vtrack = Cs.SandwichMultiply(Jc*H_T);
	    double Vtrack=H*Cs*H_T;
            double VRes;

            bool skip_ring=(my_cdchits[id]->hit->wire->ring==RING_TO_SKIP);
            if (skip_ring) VRes = Vhit + Vtrack;
            else VRes = Vhit - Vtrack;

            if (DEBUG_LEVEL>1 && (!isfinite(VRes) || VRes < 0.0) ) _DBG_ << " SmoothCentral Problem: VRes is " << VRes << " = " << Vhit << " - " << Vtrack << endl;

	    double lambda=atan(Ss(state_tanl));
	    double sinl=sin(lambda);
	    double cosl=cos(lambda);
	    double cosThetaRel=my_cdchits[id]->hit->wire->udir.Dot(DVector3(cosphi*cosl,
								 sinphi*cosl,
								 sinl));
            pull_t thisPull(cdc_updates[id].doca-d,sqrt(VRes),
			    central_traj[m].s,cdc_updates[id].tdrift,
			    d,my_cdchits[id]->hit,NULL,
			    diff.Phi(),myS(state_z),cosThetaRel,
			    cdc_updates[id].tcorr);

            thisPull.AddTrackDerivatives(alignmentDerivatives);
            cdc_pulls.push_back(thisPull);
	 }
         else{
	   A=central_traj[m].Ckk*JT*C.InvertSym();
	   Ss=central_traj[m].Skk+A*(Ss-S);
	   Cs=central_traj[m].Ckk+A*(Cs-C)*A.Transpose();      
         }
      }
      else{
	A=central_traj[m].Ckk*JT*C.InvertSym();
	Ss=central_traj[m].Skk+A*(Ss-S);
	Cs=central_traj[m].Ckk+A*(Cs-C)*A.Transpose();      
      }
      S=central_traj[m].Skk;
      C=central_traj[m].Ckk;
      JT=central_traj[m].J.Transpose();
   }

   // ... last entries?
   // Don't really need since we should have already encountered all of the hits

   return NOERROR; 

}

// Smoothing algorithm for the forward_traj_cdc trajectory.  
// Updates the state vector
// at each step (going in the reverse direction to the filter) based on the 
// information from all the steps and outputs the state vector at the
// outermost step.
jerror_t DTrackFitterKalmanSIMD::SmoothForwardCDC(vector<pull_t>&cdc_pulls){  
   if (forward_traj.size()<2) return RESOURCE_UNAVAILABLE;

   unsigned int max=forward_traj.size()-1;
   DMatrix5x1 S=(forward_traj[max].Skk);
   DMatrix5x5 C=(forward_traj[max].Ckk);
   DMatrix5x5 JT=forward_traj[max].J.Transpose();
   DMatrix5x1 Ss=S;
   DMatrix5x5 Cs=C;
   DMatrix5x5 A;
   for (unsigned int m=max-1;m>0;m--){
      if (forward_traj[m].h_id>999){ 
         unsigned int cdc_index=forward_traj[m].h_id-1000; 	
         if(cdc_used_in_fit[cdc_index] && my_cdchits[cdc_index]->status == good_hit){
            if (DEBUG_LEVEL > 5)  {
               _DBG_ << " Smoothing CDC index " << cdc_index << " ring " << my_cdchits[cdc_index]->hit->wire->ring
                  << " straw " << my_cdchits[cdc_index]->hit->wire->straw << endl;
            }

            A=cdc_updates[cdc_index].C*JT*C.InvertSym();
            Ss=cdc_updates[cdc_index].S+A*(Ss-S);
            if (!Ss.IsFinite()){
               if (DEBUG_LEVEL>5) 
                  _DBG_ << "Invalid values for smoothed parameters..." << endl;
               return VALUE_OUT_OF_RANGE;
            }

            Cs=cdc_updates[cdc_index].C+A*(Cs-C)*A.Transpose();

            if (!Cs.IsPosDef()){
               if (DEBUG_LEVEL>5){
                  _DBG_ << "Covariance Matrix not Pos Def..." << endl;
                  _DBG_ << " cdc_updates[cdc_index].C A C_ Cs " << endl;
                  cdc_updates[cdc_index].C.Print();
                  A.Print();
                  C.Print();
                  Cs.Print();
               }
               return VALUE_OUT_OF_RANGE;
            }
            if(FillPullsVectorEntry(Ss,Cs,forward_traj[m],my_cdchits[cdc_index],
				    cdc_updates[cdc_index],cdc_pulls) != NOERROR) return VALUE_OUT_OF_RANGE;

         }
         else{
            A=forward_traj[m].Ckk*JT*C.InvertSym();
            Ss=forward_traj[m].Skk+A*(Ss-S);
            Cs=forward_traj[m].Ckk+A*(Cs-C)*A.Transpose();
         }
      }
      else{
         A=forward_traj[m].Ckk*JT*C.InvertSym();
         Ss=forward_traj[m].Skk+A*(Ss-S);
         Cs=forward_traj[m].Ckk+A*(Cs-C)*A.Transpose();
      }

      S=forward_traj[m].Skk;
      C=forward_traj[m].Ckk;
      JT=forward_traj[m].J.Transpose();
   }

   return NOERROR;
}

// Fill the pulls vector with the best residual information using the smoothed
// filter results.  Uses Brent's algorithm to find the distance of closest 
// approach to the wire hit.
jerror_t DTrackFitterKalmanSIMD::FillPullsVectorEntry(const DMatrix5x1 &Ss,
      const DMatrix5x5 &Cs,
						      const DKalmanForwardTrajectory_t &traj,const DKalmanSIMDCDCHit_t *hit,const DKalmanUpdate_t &update,
						      vector<pull_t>&my_pulls){

   // Get estimate for energy loss
   double dEdx=GetdEdx(Ss(state_q_over_p),traj.K_rho_Z_over_A,traj.rho_Z_over_A,
         traj.LnI,traj.Z);

   // Use Brent's algorithm to find the doca to the wire
   DMatrix5x1 myS=Ss;
   DMatrix5x1 myS_temp=Ss;
   DMatrix5x5 myC=Cs;
   double mydz;
   double z=traj.z;
   DVector2 origin=hit->origin;
   DVector2 dir=hit->dir;
   double z0wire=hit->z0wire;
   if(BrentForward(z,dEdx,z0wire,origin,dir,myS,mydz) != NOERROR) return VALUE_OUT_OF_RANGE;
   if(DEBUG_HISTS)alignDerivHists[23]->Fill(mydz);
   double new_z=z+mydz;
   DVector2 wirepos=origin+(new_z-z0wire)*dir;
   double cosstereo=hit->cosstereo;
   DVector2 xy(myS(state_x),myS(state_y));

   DVector2 diff=xy-wirepos;
   double d=cosstereo*diff.Mod();

   // If we are doing the alignment, we need to numerically calculate the derivatives
   // wrt the wire origin, direction, and the track parameters.
   vector<double> alignmentDerivatives;
   if (ALIGNMENT_FORWARD){
      double dzcut_min=0., dzcut_max=1.;
      DMatrix5x1 alignS=Ss; // We will mess with this one
      DVector3 wireDir = hit->hit->wire->udir;
      double cosstereo_shifted;
      double aligndz;
      alignmentDerivatives.resize(12);

      // Set t0 derivative
      alignmentDerivatives[CDCTrackD::dDdt0]=update.dDdt0;

      // Wire position shift
      double wposShift=0.025; 
      double wdirShift=0.00005; 

      // Shift each track parameter
      double shiftFraction=0.01;
      double shift_x=0.01;
      double shift_y=0.01;
      double shift_tx=shiftFraction*Ss(state_tx);
      double shift_ty=shiftFraction*Ss(state_ty);;
      double shift_q_over_p=shiftFraction*Ss(state_q_over_p);

      // Some data containers we don't need multiples of
      double z0_shifted, new_z_shifted;
      DVector2 shift, origin_shifted, dir_shifted, wirepos_shifted, diff_shifted, xy_shifted;

      // The DOCA for the shifted states == f(x+h)
      double d_dOriginX=0., d_dOriginY=0., d_dOriginZ=0.;
      double d_dDirX=0., d_dDirY=0., d_dDirZ=0.;
      double d_dS0=0., d_dS1=0., d_dS2=0., d_dS3=0., d_dS4=0.;
      // Let's do the wire shifts first

      //dOriginX
      alignS=Ss;
      aligndz=0.;
      shift.Set(wposShift, 0.);
      origin_shifted=origin+shift;
      dir_shifted=dir;
      z0_shifted=z0wire;
      if(BrentForward(z,dEdx,z0_shifted,origin_shifted,dir_shifted,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      if(aligndz < dzcut_min || aligndz > dzcut_max) return VALUE_OUT_OF_RANGE; 
      //if(BrentsAlgorithm(z,-mStepSizeZ,dEdx,z0_shifted,origin_shifted,dir_shifted,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      new_z_shifted=z+aligndz;
      wirepos_shifted=origin_shifted+(new_z_shifted-z0_shifted)*dir_shifted;
      xy_shifted.Set(alignS(state_x),alignS(state_y));
      diff_shifted=xy_shifted-wirepos_shifted;
      d_dOriginX=cosstereo*diff_shifted.Mod()+EPS;
      alignmentDerivatives[CDCTrackD::dDOCAdOriginX] = (d_dOriginX - d)/wposShift;
      if(DEBUG_HISTS){
         alignDerivHists[24]->Fill(aligndz);
         alignDerivHists[35]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdOriginX]);
         brentCheckHists[0]->Fill(aligndz,d_dOriginX);
      }

      //dOriginY
      alignS=Ss;
      aligndz=0.;
      shift.Set(0.,wposShift);
      origin_shifted=origin+shift;
      dir_shifted=dir;
      z0_shifted=z0wire;
      if(BrentForward(z,dEdx,z0_shifted,origin_shifted,dir_shifted,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      if(aligndz < dzcut_min || aligndz > dzcut_max) return VALUE_OUT_OF_RANGE;
      //if(BrentsAlgorithm(z,-mStepSizeZ,dEdx,z0_shifted,origin_shifted,dir_shifted,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      new_z_shifted=z+aligndz;
      wirepos_shifted=origin_shifted+(new_z_shifted-z0_shifted)*dir_shifted;
      xy_shifted.Set(alignS(state_x),alignS(state_y));
      diff_shifted=xy_shifted-wirepos_shifted;
      d_dOriginY=cosstereo*diff_shifted.Mod()+EPS;
      alignmentDerivatives[CDCTrackD::dDOCAdOriginY] = (d_dOriginY - d)/wposShift;
      if(DEBUG_HISTS){
         alignDerivHists[25]->Fill(aligndz);
         alignDerivHists[36]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdOriginY]);
         brentCheckHists[0]->Fill(aligndz,d_dOriginY);
      }

      //dOriginZ
      alignS=Ss;
      aligndz=0.;
      origin_shifted=origin;
      dir_shifted=dir;
      z0_shifted=z0wire+wposShift;
      if(BrentForward(z,dEdx,z0_shifted,origin_shifted,dir_shifted,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      if(aligndz < dzcut_min || aligndz > dzcut_max) return VALUE_OUT_OF_RANGE;
      //if(BrentsAlgorithm(z,-mStepSizeZ,dEdx,z0_shifted,origin_shifted,dir_shifted,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      new_z_shifted=z+aligndz;
      wirepos_shifted=origin_shifted+(new_z_shifted-z0_shifted)*dir_shifted;
      xy_shifted.Set(alignS(state_x),alignS(state_y));
      diff_shifted=xy_shifted-wirepos_shifted;
      d_dOriginZ=cosstereo*diff_shifted.Mod()+EPS;
      alignmentDerivatives[CDCTrackD::dDOCAdOriginZ] = (d_dOriginZ - d)/wposShift;
      if(DEBUG_HISTS){
         alignDerivHists[26]->Fill(aligndz);
         alignDerivHists[37]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdOriginZ]);
         brentCheckHists[0]->Fill(aligndz,d_dOriginZ);
      }

      //dDirX
      alignS=Ss;
      aligndz=0.;
      shift.Set(wdirShift,0.);
      origin_shifted=origin;
      z0_shifted=z0wire;
      dir_shifted=dir+shift;
      cosstereo_shifted = cos((wireDir+DVector3(wdirShift,0.,0.)).Angle(DVector3(0.,0.,1.)));
      if(BrentForward(z,dEdx,z0_shifted,origin_shifted,dir_shifted,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      if(aligndz < dzcut_min || aligndz > dzcut_max) return VALUE_OUT_OF_RANGE;
      //if(BrentsAlgorithm(z,-mStepSizeZ,dEdx,z0_shifted,origin_shifted,dir_shifted,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      new_z_shifted=z+aligndz;
      wirepos_shifted=origin_shifted+(new_z_shifted-z0_shifted)*dir_shifted;
      xy_shifted.Set(alignS(state_x),alignS(state_y));
      diff_shifted=xy_shifted-wirepos_shifted;
      d_dDirX=cosstereo_shifted*diff_shifted.Mod()+EPS;
      alignmentDerivatives[CDCTrackD::dDOCAdDirX] = (d_dDirX - d)/wdirShift;
      if(DEBUG_HISTS){
         alignDerivHists[27]->Fill(aligndz);
         alignDerivHists[38]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdDirX]);
      }

      //dDirY
      alignS=Ss;
      aligndz=0.;
      shift.Set(0.,wdirShift);
      origin_shifted=origin;
      z0_shifted=z0wire;
      dir_shifted=dir+shift;
      cosstereo_shifted = cos((wireDir+DVector3(0.,wdirShift,0.)).Angle(DVector3(0.,0.,1.)));
      if(BrentForward(z,dEdx,z0_shifted,origin_shifted,dir_shifted,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      if(aligndz < dzcut_min || aligndz > dzcut_max) return VALUE_OUT_OF_RANGE;
      //if(BrentsAlgorithm(z,-mStepSizeZ,dEdx,z0_shifted,origin_shifted,dir_shifted,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      new_z_shifted=z+aligndz;
      wirepos_shifted=origin_shifted+(new_z_shifted-z0_shifted)*dir_shifted;
      xy_shifted.Set(alignS(state_x),alignS(state_y));
      diff_shifted=xy_shifted-wirepos_shifted;
      d_dDirY=cosstereo_shifted*diff_shifted.Mod()+EPS;
      alignmentDerivatives[CDCTrackD::dDOCAdDirY] = (d_dDirY - d)/wdirShift;
      if(DEBUG_HISTS){
         alignDerivHists[28]->Fill(aligndz);
         alignDerivHists[39]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdDirY]);
      }

      //dDirZ - This is divided out in this code
      alignS=Ss;
      aligndz=0.;
      origin_shifted=origin;
      dir_shifted.Set(wireDir.X()/(wireDir.Z()+wdirShift), wireDir.Y()/(wireDir.Z()+wdirShift));
      z0_shifted=z0wire;
      cosstereo_shifted = cos((wireDir+DVector3(0.,0.,wdirShift)).Angle(DVector3(0.,0.,1.)));
      if(BrentForward(z,dEdx,z0_shifted,origin_shifted,dir_shifted,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      if(aligndz < dzcut_min || aligndz > dzcut_max) return VALUE_OUT_OF_RANGE;
      //if(BrentsAlgorithm(z,-mStepSizeZ,dEdx,z0_shifted,origin_shifted,dir_shifted,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      new_z_shifted=z+aligndz;
      wirepos_shifted=origin_shifted+(new_z_shifted-z0_shifted)*dir_shifted;
      xy_shifted.Set(alignS(state_x),alignS(state_y));
      diff_shifted=xy_shifted-wirepos_shifted;
      d_dDirZ=cosstereo_shifted*diff_shifted.Mod()+EPS;
      alignmentDerivatives[CDCTrackD::dDOCAdDirZ] = (d_dDirZ - d)/wdirShift;
      if(DEBUG_HISTS){
         alignDerivHists[29]->Fill(aligndz);
         alignDerivHists[40]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdDirZ]);
      }

      // And now the derivatives wrt the track parameters

      DMatrix5x1 trackShiftS0(shift_x, 0., 0., 0., 0.);
      DMatrix5x1 trackShiftS1(0., shift_y, 0., 0., 0.);
      DMatrix5x1 trackShiftS2(0., 0., shift_tx, 0., 0.);
      DMatrix5x1 trackShiftS3(0., 0., 0., shift_ty, 0.);
      DMatrix5x1 trackShiftS4(0., 0., 0., 0., shift_q_over_p);

      // dS0
      alignS=Ss+trackShiftS0;
      aligndz=0.;
      if(BrentForward(z,dEdx,z0wire,origin,dir,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      if(aligndz < dzcut_min || aligndz > dzcut_max) return VALUE_OUT_OF_RANGE;
      //if(BrentsAlgorithm(z,-mStepSizeZ,dEdx,z0wire,origin,dir,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      new_z_shifted=z+aligndz;
      wirepos_shifted=origin+(new_z_shifted-z0wire)*dir;
      xy_shifted.Set(alignS(state_x),alignS(state_y));
      diff_shifted=xy_shifted-wirepos_shifted;
      d_dS0=cosstereo*diff_shifted.Mod()+EPS;
      alignmentDerivatives[CDCTrackD::dDOCAdS0] = (d_dS0 - d)/shift_x;
      if(DEBUG_HISTS){
         alignDerivHists[30]->Fill(aligndz);
         alignDerivHists[41]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdS0]);
      }

      // dS1
      alignS=Ss+trackShiftS1;
      aligndz=0.;
      if(BrentForward(z,dEdx,z0wire,origin,dir,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      if(aligndz < dzcut_min || aligndz > dzcut_max) return VALUE_OUT_OF_RANGE;
      //if(BrentsAlgorithm(z,-mStepSizeZ,dEdx,z0wire,origin,dir,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      new_z_shifted=z+aligndz;
      wirepos_shifted=origin+(new_z_shifted-z0wire)*dir;
      xy_shifted.Set(alignS(state_x),alignS(state_y));
      diff_shifted=xy_shifted-wirepos_shifted;
      d_dS1=cosstereo*diff_shifted.Mod()+EPS;
      alignmentDerivatives[CDCTrackD::dDOCAdS1] = (d_dS1 - d)/shift_y;
      if(DEBUG_HISTS){
         alignDerivHists[31]->Fill(aligndz);
         alignDerivHists[42]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdS1]);
      }

      // dS2
      alignS=Ss+trackShiftS2;
      aligndz=0.;
      if(BrentForward(z,dEdx,z0wire,origin,dir,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      if(aligndz < dzcut_min || aligndz > dzcut_max) return VALUE_OUT_OF_RANGE;
      //if(BrentsAlgorithm(z,-mStepSizeZ,dEdx,z0wire,origin,dir,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      new_z_shifted=z+aligndz;
      wirepos_shifted=origin+(new_z_shifted-z0wire)*dir;
      xy_shifted.Set(alignS(state_x),alignS(state_y));
      diff_shifted=xy_shifted-wirepos_shifted;
      d_dS2=cosstereo*diff_shifted.Mod()+EPS;
      alignmentDerivatives[CDCTrackD::dDOCAdS2] = (d_dS2 - d)/shift_tx;
      if(DEBUG_HISTS){
         alignDerivHists[32]->Fill(aligndz);
         alignDerivHists[43]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdS2]);
      }

      // dS3
      alignS=Ss+trackShiftS3;
      aligndz=0.;
      if(BrentForward(z,dEdx,z0wire,origin,dir,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      if(aligndz < dzcut_min || aligndz > dzcut_max) return VALUE_OUT_OF_RANGE;
      //if(BrentsAlgorithm(z,-mStepSizeZ,dEdx,z0wire,origin,dir,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      new_z_shifted=z+aligndz;
      wirepos_shifted=origin+(new_z_shifted-z0wire)*dir;
      xy_shifted.Set(alignS(state_x),alignS(state_y));
      diff_shifted=xy_shifted-wirepos_shifted;
      d_dS3=cosstereo*diff_shifted.Mod()+EPS;
      alignmentDerivatives[CDCTrackD::dDOCAdS3] = (d_dS3 - d)/shift_ty;
      if(DEBUG_HISTS){
         alignDerivHists[33]->Fill(aligndz);
         alignDerivHists[44]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdS3]);
      }

      // dS4
      alignS=Ss+trackShiftS4;
      aligndz=0.;
      if(BrentForward(z,dEdx,z0wire,origin,dir,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      if(aligndz < dzcut_min || aligndz > dzcut_max) return VALUE_OUT_OF_RANGE;
      //if(BrentsAlgorithm(z,-mStepSizeZ,dEdx,z0wire,origin,dir,alignS,aligndz) != NOERROR) return VALUE_OUT_OF_RANGE;
      new_z_shifted=z+aligndz;
      wirepos_shifted=origin+(new_z_shifted-z0wire)*dir;
      xy_shifted.Set(alignS(state_x),alignS(state_y));
      diff_shifted=xy_shifted-wirepos_shifted;
      d_dS4=cosstereo*diff_shifted.Mod()+EPS;
      alignmentDerivatives[CDCTrackD::dDOCAdS4] = (d_dS4 - d)/shift_q_over_p;
      if(DEBUG_HISTS){
         alignDerivHists[34]->Fill(aligndz);
         alignDerivHists[45]->Fill(alignmentDerivatives[CDCTrackD::dDOCAdS4]);
      }
   }

   // Find the field and gradient at (old_x,old_y,old_z) and compute the 
   // Jacobian matrix for transforming from S to myS
   bfield->GetFieldAndGradient(Ss(state_x),Ss(state_y),z,
         Bx,By,Bz,dBxdx,dBxdy,dBxdz,dBydx,
         dBydy,dBydz,dBzdx,dBzdy,dBzdz);
   DMatrix5x5 Jc;
   StepJacobian(z,new_z,Ss,dEdx,Jc);

   // Find the projection matrix
   DMatrix5x1 H_T;
   double cosstereo2_over_d=cosstereo*cosstereo/d;
   H_T(state_x)=diff.X()*cosstereo2_over_d;
   H_T(state_y)=diff.Y()*cosstereo2_over_d;
   DMatrix1x5 H;
   H(state_x)=H_T(state_x);
   H(state_y)=H_T(state_y);

   // Find the variance for this hit

   bool skip_ring=(hit->hit->wire->ring==RING_TO_SKIP);

   double V=update.variance;
   myC=Jc*myC*Jc.Transpose();
   if (skip_ring) V+=H*myC*H_T;
   else V-=H*myC*H_T;

   if (DEBUG_LEVEL>1 && (!isfinite(V) || V < 0.0) ) _DBG_ << " Problem: V is " << V << endl;

   double tx=Ss(state_tx);
   double ty=Ss(state_ty);
   double scale=1./sqrt(1.+tx*tx+ty*ty);
   double cosThetaRel=hit->hit->wire->udir.Dot(DVector3(scale*tx,scale*ty,scale));

   pull_t thisPull(update.doca-d,sqrt(V),traj.s,update.tdrift,d,hit->hit,
		   NULL,diff.Phi(),new_z,cosThetaRel,update.tcorr);
   thisPull.AddTrackDerivatives(alignmentDerivatives);
   my_pulls.push_back(thisPull);
   return NOERROR;
}

// Transform the 5x5 covariance matrix from the forward parametrization to the 
// central parametrization.
void DTrackFitterKalmanSIMD::TransformCovariance(DMatrix5x5 &C){
   DMatrix5x5 J;
   double tsquare=tx_*tx_+ty_*ty_;
   double cosl=cos(atan(tanl_));
   double tanl2=tanl_*tanl_;
   double tanl3=tanl2*tanl_;
   double factor=1./sqrt(1.+tsquare);
   J(state_z,state_x)=tx_/tsquare;
   J(state_z,state_y)=ty_/tsquare;
   double diff=tx_*tx_-ty_*ty_;
   double frac=1./(tsquare*tsquare);
   J(state_z,state_tx)=-(x_*diff+2.*tx_*ty_*y_)*frac;
   J(state_z,state_ty)=-(2.*tx_*ty_*x_-y_*diff)*frac;
   J(state_tanl,state_tx)=-tx_*tanl3;
   J(state_tanl,state_ty)=-ty_*tanl3;
   J(state_q_over_pt,state_q_over_p)=1./cosl;
   J(state_q_over_pt,state_tx)=-q_over_p_*tx_*tanl3*factor;
   J(state_q_over_pt,state_ty)=-q_over_p_*ty_*tanl3*factor;
   J(state_phi,state_tx)=-ty_*tanl2;
   J(state_phi,state_ty)=tx_*tanl2;
   J(state_D,state_x)=x_/D_;
   J(state_D,state_y)=y_/D_;

   C=J*C*J.Transpose();      

}

jerror_t DTrackFitterKalmanSIMD::BrentForward(double z, double dedx, const double z0w,
      const DVector2 &origin, const DVector2 &dir, DMatrix5x1 &S, double &dz){

   DVector2 wirepos=origin;
   wirepos+=(z-z0w)*dir;
   double dx=S(state_x)-wirepos.X();
   double dy=S(state_y)-wirepos.Y();
   double doca2 = dx*dx+dy*dy;

   if (BrentsAlgorithm(z,-mStepSizeZ,dedx,z0w,origin,dir,S,dz)!=NOERROR){
      return VALUE_OUT_OF_RANGE;
   }

   double newz = z+dz;
   unsigned int maxSteps=5;
   unsigned int stepCounter=0;
   if (fabs(dz)<EPS3){
      // doca
      double old_doca2=doca2;

      double ztemp=newz;
      newz=ztemp-mStepSizeZ;
      Step(ztemp,newz,dedx,S);
      // new wire position
      wirepos=origin;
      wirepos+=(newz-z0w)*dir;

      dx=S(state_x)-wirepos.X();
      dy=S(state_y)-wirepos.Y();
      doca2=dx*dx+dy*dy;
      ztemp=newz;

      while(doca2<old_doca2 && stepCounter<maxSteps){
         newz=ztemp-mStepSizeZ;
         old_doca2=doca2;

         // Step to the new z position
         Step(ztemp,newz,dedx,S);
         stepCounter++;

         // find the new distance to the wire
         wirepos=origin;
         wirepos+=(newz-z0w)*dir;

         dx=S(state_x)-wirepos.X();
         dy=S(state_y)-wirepos.Y();
         doca2=dx*dx+dy*dy;

         ztemp=newz;
      }

      // Find the true doca
      double dz2=0.;
      if (BrentsAlgorithm(newz,-mStepSizeZ,dedx,z0w,origin,dir,S,dz2)!=NOERROR){
         return VALUE_OUT_OF_RANGE;
      }
      newz=ztemp+dz2;

      // Change in z relative to where we started for this wire
      dz=newz-z;
   }
   else if (fabs(dz)>2.*mStepSizeZ-EPS3){    
      // whoops, looks like we didn't actually bracket the minimum 
      // after all.  Swim to make sure we pass the minimum doca.

      double ztemp=newz;
      // new wire position
      wirepos=origin;
      wirepos+=(ztemp-z0w)*dir;

      // doca
      double old_doca2=doca2;

      dx=S(state_x)-wirepos.X();
      dy=S(state_y)-wirepos.Y();
      doca2=dx*dx+dy*dy;

      while(doca2<old_doca2 && stepCounter<10*maxSteps){
         newz=ztemp+mStepSizeZ;
         old_doca2=doca2;

         // Step to the new z position
         Step(ztemp,newz,dedx,S);
         stepCounter++;

         // find the new distance to the wire
         wirepos=origin;
         wirepos+=(newz-z0w)*dir;

         dx=S(state_x)-wirepos.X();
         dy=S(state_y)-wirepos.Y();
         doca2=dx*dx+dy*dy;

         ztemp=newz;
      }

      // Find the true doca
      double dz2=0.;
      if (BrentsAlgorithm(newz,mStepSizeZ,dedx,z0w,origin,dir,S,dz2)!=NOERROR){
         return VALUE_OUT_OF_RANGE;
      }
      newz=ztemp+dz2;

      // Change in z relative to where we started for this wire
      dz=newz-z;
   }
   return NOERROR;
}

jerror_t DTrackFitterKalmanSIMD::BrentCentral(double dedx, DVector2 &xy, const double z0w, const DVector2 &origin, const DVector2 &dir, DMatrix5x1 &Sc, double &ds){

   DVector2 wirexy=origin;
   wirexy+=(Sc(state_z)-z0w)*dir;

   // new doca
   double doca2=(xy-wirexy).Mod2();
   double old_doca2=doca2;

   if (BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dedx,xy,z0w,
            origin,dir,Sc,ds)!=NOERROR){
      return VALUE_OUT_OF_RANGE;
   }

   unsigned int maxSteps=3;
   unsigned int stepCounter=0;

   if (fabs(ds)<EPS3){
      double my_ds=ds;
      old_doca2=doca2;
      Step(xy,-mStepSizeS,Sc,dedx);
      my_ds-=mStepSizeS;
      wirexy=origin;
      wirexy+=(Sc(state_z)-z0w)*dir;
      doca2=(xy-wirexy).Mod2();
      while(doca2<old_doca2 && stepCounter<maxSteps){
         old_doca2=doca2;
         // Bail if the transverse momentum has dropped below some minimum
         if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
            return VALUE_OUT_OF_RANGE;
         }

         // Step through the field
         Step(xy,-mStepSizeS,Sc,dedx);
         stepCounter++;

         wirexy=origin;
         wirexy+=(Sc(state_z)-z0w)*dir;
         doca2=(xy-wirexy).Mod2();

         my_ds-=mStepSizeS;
      }
      // Swim to the "true" doca
      double ds2=0.;
      if (BrentsAlgorithm(-mStepSizeS,-mStepSizeS,dedx,xy,z0w,
               origin,dir,Sc,ds2)!=NOERROR){
         return VALUE_OUT_OF_RANGE;
      }
      ds=my_ds+ds2;
   }
   else if (fabs(ds)>2*mStepSizeS-EPS3){
      double my_ds=ds;

      // new wire position
      wirexy=origin;
      wirexy+=(Sc(state_z)-z0w)*dir;

      // doca
      old_doca2=doca2;
      doca2=(xy-wirexy).Mod2();

      while(doca2<old_doca2 && stepCounter<maxSteps){
         old_doca2=doca2;

         // Bail if the transverse momentum has dropped below some minimum
         if (fabs(Sc(state_q_over_pt))>Q_OVER_PT_MAX){
            return VALUE_OUT_OF_RANGE;
         }

         // Step through the field
         Step(xy,mStepSizeS,Sc,dedx);
         stepCounter++;

         // Find the new distance to the wire
         wirexy=origin;
         wirexy+=(Sc(state_z)-z0w)*dir;
         doca2=(xy-wirexy).Mod2();

         my_ds+=mStepSizeS;
      }
      // Swim to the "true" doca
      double ds2=0.;
      if (BrentsAlgorithm(mStepSizeS,mStepSizeS,dedx,xy,z0w,
               origin,dir,Sc,ds2)!=NOERROR){
         return VALUE_OUT_OF_RANGE;
      }
      ds=my_ds+ds2;
   }
   return NOERROR;
}

// Routine to find intersections with surfaces useful at a later stage for track
// matching
jerror_t DTrackFitterKalmanSIMD::ExtrapolateForwardToOtherDetectors(){
  if (forward_traj.size()<2) return RESOURCE_UNAVAILABLE;

  // First deal with start counter.  Only do this if the track has a chance
  // to intersect with the start counter volume.
  unsigned int inner_index=forward_traj.size()-1; 
  unsigned int index_beyond_start_counter=inner_index;
  DMatrix5x1 S=forward_traj[inner_index].S;
  bool intersected_start_counter=false;
  if (sc_norm.empty()==false
      && S(state_x)*S(state_x)+S(state_y)*S(state_y)<SC_BARREL_R2
      && forward_traj[inner_index].z<SC_END_NOSE_Z){
    double d_old=1000.,d=1000.,z=0.;
    unsigned int index=0;
    for (unsigned int m=0;m<12;m++){
      unsigned int k=inner_index;
      for (;k>1;k--){ 
	S=forward_traj[k].S;
	z=forward_traj[k].z;

	double dphi=atan2(S(state_y),S(state_x))-SC_PHI_SECTOR1;
	if (dphi<0) dphi+=2.*M_PI;
	index=int(floor(dphi/(2.*M_PI/30.)));
	if (index>29) index=0;
	d=sc_norm[index][m].Dot(DVector3(S(state_x),S(state_y),z)
				-sc_pos[index][m]);
	if (d*d_old<0){ // break if we cross the current plane
	  if (m==0) index_beyond_start_counter=k;
	  break;
	}
	d_old=d;
      }
      // if the z position would be beyond the current segment along z of 
      // the start counter, move to the next plane
      if (z>sc_pos[index][m+1].z()&&m<11){
	continue;
      }   
      // allow for a little slop at the end of the nose
      else if (z<sc_pos[index][sc_pos[0].size()-1].z()+0.1){
	// Hone in on intersection with the appropriate segment of the start 
	// counter
	int count=0;
	DMatrix5x1 bestS=S;
	double dmin=d;
	double bestz=z;
	double t=forward_traj[k].t;
	double s=forward_traj[k].s;  
	// Magnetic field 
	bfield->GetField(S(state_x),S(state_y),z,Bx,By,Bz); 

	while (fabs(d)>0.05 && count<20){
	  // track direction
	  DVector3 phat(S(state_tx),S(state_ty),1);
	  phat.SetMag(1.);

	  // Step to the start counter plane
	  double ds=d/sc_norm[index][m].Dot(phat);
	  FastStep(z,-ds,0.,S);
	
	  // Flight time
	  double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
	  double one_over_beta2=1.+mass2*q_over_p_sq;
	  if (one_over_beta2>BIG) one_over_beta2=BIG;
	  t-=ds*sqrt(one_over_beta2); // in units where c=1
	  s-=ds;

	  // Find the index for the nearest start counter paddle
         double dphi=atan2(S(state_y),S(state_x))-SC_PHI_SECTOR1;
         if (dphi<0) dphi+=2.*M_PI;
         index=int(floor(dphi/(2.*M_PI/30.)));

	 // Find the new distance to the start counter (which could now be to
         // a plane in the one adjacent to the one before the step...)
	 d=sc_norm[index][m].Dot(DVector3(S(state_x),S(state_y),z)
                                 -sc_pos[index][m]);
	 if (fabs(d)<fabs(dmin)){
           bestS=S;
           dmin=d;
	   bestz=z;
	 }
	 count++;
	}

	// New position and momentum
	double tsquare=bestS(state_tx)*bestS(state_tx)
	  +bestS(state_ty)*bestS(state_ty);
	double tanl=1./sqrt(tsquare);
	double cosl=cos(atan(tanl));
	double pt=cosl/fabs(bestS(state_q_over_p));
	double phi=atan2(bestS(state_ty),bestS(state_tx));
	DVector3 position(bestS(state_x),bestS(state_y),bestz);
	DVector3 momentum(pt*cos(phi),pt*sin(phi),pt*tanl);
	extrapolations[SYS_START].push_back(Extrapolation_t(position,momentum,
						    t*TIME_UNIT_CONVERSION,s));

	//printf("forward track:\n");
	//position.Print();
	intersected_start_counter=true;
	break;
      }
    }   
  }
  // Accumulate multiple-scattering terms for use in matching routines
  double s_theta_ms_sum=0.;
  double theta2ms_sum=0.;
  if (intersected_start_counter){
    for (unsigned int k=inner_index;k>index_beyond_start_counter;k--){
      double ds_theta_ms_sq=3.*fabs(forward_traj[k].Q(state_x,state_x));
      s_theta_ms_sum+=sqrt(ds_theta_ms_sq);
      double ds=forward_traj[k].s-forward_traj[k-1].s;
      theta2ms_sum+=ds_theta_ms_sq/(ds*ds);
    }
  }

  // Deal with points within fiducial volume of chambers
  unsigned int fdc_plane=0;
  mT0Detector=SYS_NULL;
  mT0MinimumDriftTime=1e6;
  for (int k=intersected_start_counter?index_beyond_start_counter:inner_index;k>=0;k--){
    double z=forward_traj[k].z;
    double t=forward_traj[k].t*TIME_UNIT_CONVERSION;
    double s=forward_traj[k].s;
    DMatrix5x1 S=forward_traj[k].S;
    double tsquare=S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty);
    double tanl=1./sqrt(tsquare);
    double cosl=cos(atan(tanl));
    double pt=cosl/fabs(S(state_q_over_p));
    double phi=atan2(S(state_ty),S(state_tx)); 

    // Find estimate for t0 using earliest drift time
    if (forward_traj[k].h_id>999){
      unsigned int index=forward_traj[k].h_id-1000;
      double dt=my_cdchits[index]->tdrift-t;
      if (dt<mT0MinimumDriftTime){
	mT0MinimumDriftTime=dt;
	mT0Detector=SYS_CDC;
      }
    }
    else if (forward_traj[k].h_id>0){
      unsigned int index=forward_traj[k].h_id-1;
      double dt=my_fdchits[index]->t-t;  
      if (dt<mT0MinimumDriftTime){
	mT0MinimumDriftTime=dt;
	mT0Detector=SYS_FDC;
      }
    }

    //multiple scattering terms
    if (k>0){
      double ds_theta_ms_sq=3.*fabs(forward_traj[k].Q(state_x,state_x));
      s_theta_ms_sum+=sqrt(ds_theta_ms_sq);  
      double ds=forward_traj[k].s-forward_traj[k-1].s;
      theta2ms_sum+=ds_theta_ms_sq/(ds*ds);
    }
    // Extrapolations in CDC region
    if (z<endplate_z){
      DVector3 position(S(state_x),S(state_y),z);
      DVector3 momentum(pt*cos(phi),pt*sin(phi),pt*tanl);
      extrapolations[SYS_CDC].push_back(Extrapolation_t(position,momentum,
						t*TIME_UNIT_CONVERSION,s,
						s_theta_ms_sum,theta2ms_sum));

    }
    else{ // extrapolations in FDC region
      if (fdc_plane==24){
	break;	
      }
    
      // output step near wire plane  
      if (z>fdc_z_wires[fdc_plane]-0.25){	
	double dz=z-fdc_z_wires[fdc_plane];
	//printf("extrp dz %f\n",dz);
	if (fabs(dz)>EPS2){
	  Step(z,fdc_z_wires[fdc_plane],0.,S);
	  tsquare=S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty);
	  tanl=1./sqrt(tsquare);
	  cosl=cos(atan(tanl));
	  pt=cosl/fabs(S(state_q_over_p));
	  phi=atan2(S(state_ty),S(state_tx)); 
	}
	DVector3 position(S(state_x),S(state_y),fdc_z_wires[fdc_plane]);
	DVector3 momentum(pt*cos(phi),pt*sin(phi),pt*tanl);
	extrapolations[SYS_FDC].push_back(Extrapolation_t(position,momentum,
						t*TIME_UNIT_CONVERSION,s,
						s_theta_ms_sum,theta2ms_sum));

	fdc_plane++;
      }
      
    }

  }
  

  //--------------------------------
  // Next swim to outer detectors...
  //--------------------------------
  DMatrix5x5 Q;  // multiple scattering matrix

  // Direction and origin of beam line
  DVector2 dir;
  DVector2 origin;

  // Energy loss
  double dEdx=0.;
  
  // material properties
  double rho_Z_over_A=0.,LnI=0.,K_rho_Z_over_A=0.,Z=0.;
  double chi2c_factor=0.,chi2a_factor=0.,chi2a_corr=0.;

  // Position variables
  double z=forward_traj[0].z;
  double newz=z,dz=0.;
  S=forward_traj[0].S;

  // Current time and path length
  double t=forward_traj[0].t;
  double s=forward_traj[0].s;
  
  // Loop to propagate track to outer detectors
  const double z_outer_max=650.;
  const double x_max=130.;
  const double y_max=130.;
  bool hit_tof=false; 
  bool hit_dirc=false;
  while (z>Z_MIN && z<z_outer_max && fabs(S(state_x))<x_max 
	 && fabs(S(state_y))<y_max){   
    // Bail if the momentum has dropped below some minimum
    if (fabs(S(state_q_over_p))>Q_OVER_P_MAX){
      if (DEBUG_LEVEL>2)
	{
	  _DBG_ << "Bailing: P = " << 1./fabs(S(state_q_over_p))
		<< endl;
	}
      return VALUE_OUT_OF_RANGE;
    }

    // Check if we have passed into the BCAL
    double r2=S(state_x)*S(state_x)+S(state_y)*S(state_y);
    if (r2>89.*89. && z<400.) return VALUE_OUT_OF_RANGE;
    if (r2>64.9*64.9 && r2<89.*89.){
      if (extrapolations.at(SYS_BCAL).size()>299){
	return VALUE_OUT_OF_RANGE;
      }
      if (z<406.){
	double tsquare=S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty);
	double tanl=1./sqrt(tsquare);
	double cosl=cos(atan(tanl));
	double pt=cosl/fabs(S(state_q_over_p));
	double phi=atan2(S(state_ty),S(state_tx));
	DVector3 position(S(state_x),S(state_y),z);
	DVector3 momentum(pt*cos(phi),pt*sin(phi),pt*tanl);
	extrapolations[SYS_BCAL].push_back(Extrapolation_t(position,momentum,
							   t*TIME_UNIT_CONVERSION,s));
      }
      else if (extrapolations.at(SYS_BCAL).size()<5){
	// There needs to be some steps inside the the volume of the BCAL for 
	// the extrapolation to be useful.  If this is not the case, clear 
	// the extrapolation vector.
	extrapolations[SYS_BCAL].clear();
      }
    }    
   
    // Relationship between arc length and z
    double dz_ds=1./sqrt(1.+S(state_tx)*S(state_tx)
			   +S(state_ty)*S(state_ty));
    
    // get material properties from the Root Geometry
    DVector3 pos(S(state_x),S(state_y),z); 
    DVector3 dir(S(state_tx),S(state_ty),1.);
    double s_to_boundary=0.;
    if (geom->FindMatKalman(pos,dir,K_rho_Z_over_A,rho_Z_over_A,LnI,Z,
			    chi2c_factor,chi2a_factor,chi2a_corr,
			    last_material_map,&s_to_boundary)
	  !=NOERROR){
      if (DEBUG_LEVEL>0)
	{
	_DBG_ << "Material error in ExtrapolateForwardToOuterDetectors!"<< endl;
	_DBG_ << " Position (x,y,z)=("<<pos.x()<<","<<pos.y()<<","<<pos.z()<<")"
	      <<endl;
      }
      return VALUE_OUT_OF_RANGE;
    }
    
    // Get dEdx for the upcoming step
    if (CORRECT_FOR_ELOSS){
      dEdx=GetdEdx(S(state_q_over_p),K_rho_Z_over_A,rho_Z_over_A,LnI,Z); 
    }
    
    // Adjust the step size 
    double ds=mStepSizeS;
    if (fabs(dEdx)>EPS){     
      ds=DE_PER_STEP/fabs(dEdx);
    }
    if (ds>mStepSizeS) ds=mStepSizeS;
    if (s_to_boundary<ds) ds=s_to_boundary;
    if (ds<MIN_STEP_SIZE) ds=MIN_STEP_SIZE;
    if (ds<0.5 && z<406. && r2>65.*65.) ds=0.5;
    dz=ds*dz_ds;
    newz=z+dz;
    if (hit_tof==false && newz>dTOFz){
      newz=dTOFz+EPS;
      ds=(newz-z)/dz_ds;
    }
    if (hit_dirc==false && newz>dDIRCz){
      newz=dDIRCz+EPS;
      ds=(newz-z)/dz_ds;
    }
    if (newz>dFCALz){
      newz=dFCALz+EPS;
      ds=(newz-z)/dz_ds;
    }
    bool got_fdc_hit=false;
    if (fdc_plane<24 && newz>fdc_z_wires[fdc_plane]){
      newz=fdc_z_wires[fdc_plane];
      ds=(newz-z)/dz_ds;
      got_fdc_hit=true;
    }
    s+=ds;

    // Flight time
    double q_over_p_sq=S(state_q_over_p)*S(state_q_over_p);
    double one_over_beta2=1.+mass2*q_over_p_sq;
    if (one_over_beta2>BIG) one_over_beta2=BIG;
    t+=ds*sqrt(one_over_beta2); // in units where c=1
    
    // Get the contribution to the covariance matrix due to multiple 
    // scattering
    GetProcessNoise(z,ds,chi2c_factor,chi2a_factor,chi2a_corr,S,Q);
    double ds_theta_ms_sq=3.*fabs(Q(state_x,state_x));
    s_theta_ms_sum+=sqrt(ds_theta_ms_sq);
    theta2ms_sum+=ds_theta_ms_sq/(ds*ds);

    // Step through field
    Step(z,newz,dEdx,S); 

    if (got_fdc_hit){
      double tsquare=S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty);
      double tanl=1./sqrt(tsquare);
      double cosl=cos(atan(tanl));
      double pt=cosl/fabs(S(state_q_over_p));
      double phi=atan2(S(state_ty),S(state_tx));
      DVector3 position(S(state_x),S(state_y),z);
      DVector3 momentum(pt*cos(phi),pt*sin(phi),pt*tanl);
      extrapolations[SYS_FDC].push_back(Extrapolation_t(position,momentum,
						t*TIME_UNIT_CONVERSION,s,
						s_theta_ms_sum,theta2ms_sum));

      fdc_plane++;
    }        
    if (hit_dirc==false && newz>dDIRCz){
      hit_dirc=true;

      double tsquare=S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty);
      double tanl=1./sqrt(tsquare);
      double cosl=cos(atan(tanl));
      double pt=cosl/fabs(S(state_q_over_p));
      double phi=atan2(S(state_ty),S(state_tx));
      DVector3 position(S(state_x),S(state_y),z);
      DVector3 momentum(pt*cos(phi),pt*sin(phi),pt*tanl);
      extrapolations[SYS_DIRC].push_back(Extrapolation_t(position,momentum,
							t*TIME_UNIT_CONVERSION,s));
    }  
    if (hit_tof==false && newz>dTOFz){
      hit_tof=true;

      double tsquare=S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty);
      double tanl=1./sqrt(tsquare);
      double cosl=cos(atan(tanl));
      double pt=cosl/fabs(S(state_q_over_p));
      double phi=atan2(S(state_ty),S(state_tx));
      DVector3 position(S(state_x),S(state_y),z);
      DVector3 momentum(pt*cos(phi),pt*sin(phi),pt*tanl);
      extrapolations[SYS_TOF].push_back(Extrapolation_t(position,momentum,
							t*TIME_UNIT_CONVERSION,s));
    }  
    if (newz>dFCALz){
      double tsquare=S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty);
      double tanl=1./sqrt(tsquare);
      double cosl=cos(atan(tanl));
      double pt=cosl/fabs(S(state_q_over_p));
      double phi=atan2(S(state_ty),S(state_tx));
      DVector3 position(S(state_x),S(state_y),z);
      DVector3 momentum(pt*cos(phi),pt*sin(phi),pt*tanl);
      extrapolations[SYS_FCAL].push_back(Extrapolation_t(position,momentum,
					    t*TIME_UNIT_CONVERSION,s)); 

      // add another extrapolation point at downstream end of FCAL
      double zend=newz+45.;
      Step(newz,zend,dEdx,S); 
      ds=(zend-newz)/dz_ds;
      t+=ds*sqrt(one_over_beta2); // in units where c=1 
      s+=ds;
      tsquare=S(state_tx)*S(state_tx)+S(state_ty)*S(state_ty);
      tanl=1./sqrt(tsquare); 
      cosl=cos(atan(tanl));
      pt=cosl/fabs(S(state_q_over_p));
      phi=atan2(S(state_ty),S(state_tx));
      position.SetXYZ(S(state_x),S(state_y),zend);
      momentum.SetXYZ(pt*cos(phi),pt*sin(phi),pt*tanl);
      extrapolations[SYS_FCAL].push_back(Extrapolation_t(position,momentum,
					    t*TIME_UNIT_CONVERSION,s)); 

      return NOERROR;
    }
    z=newz;
  }
  return NOERROR;
}

// Routine to find intersections with surfaces useful at a later stage for track
// matching
jerror_t DTrackFitterKalmanSIMD::ExtrapolateCentralToOtherDetectors(){
  if (central_traj.size()<2) return RESOURCE_UNAVAILABLE;

  // First deal with start counter.  Only do this if the track has a chance
  // to intersect with the start counter volume.
  unsigned int inner_index=central_traj.size()-1;  
  unsigned int index_beyond_start_counter=inner_index;
  DVector2 xy=central_traj[inner_index].xy;
  DMatrix5x1 S=central_traj[inner_index].S;
  if (sc_norm.empty()==false
      &&xy.Mod2()<SC_BARREL_R2&& S(state_z)<SC_END_NOSE_Z){ 
    double d_old=1000.,d=1000.,z=0.;
    unsigned int index=0;
    for (unsigned int m=0;m<12;m++){
      unsigned int k=inner_index;
      for (;k>1;k--){ 
	S=central_traj[k].S;
	z=S(state_z);
	xy=central_traj[k].xy;

	double dphi=xy.Phi()-SC_PHI_SECTOR1;
	if (dphi<0) dphi+=2.*M_PI;
	index=int(floor(dphi/(2.*M_PI/30.)));
	if (index>29) index=0;
	//cout << "dphi " << dphi << " " << index << endl;
	
	d=sc_norm[index][m].Dot(DVector3(xy.X(),xy.Y(),z)
				       -sc_pos[index][m]);
	
	if (d*d_old<0){ // break if we cross the current plane
	  if (m==0) index_beyond_start_counter=k;
	  break;
	}
	d_old=d;
      }
      // if the z position would be beyond the current segment along z of 
      // the start counter, move to the next plane
      if (z>sc_pos[index][m+1].z()&&m<11){
	continue;
      } 
      // allow for a little slop at the end of the nose
      else if (z<sc_pos[index][sc_pos[0].size()-1].z()+0.1){ 
	// Propagate the state and covariance through the field
	// using a straight-line approximation for each step to zero in on the 
	// start counter paddle
	int count=0;
	DMatrix5x1 bestS=S;
	double dmin=d;
	DVector2 bestXY=central_traj[k].xy;
	double t=central_traj[k].t;
	double s=central_traj[k].s;
	// Magnetic field 
	bfield->GetField(xy.X(),xy.Y(),S(state_z),Bx,By,Bz); 
	
	while (fabs(d)>0.05 && count<20){
	  // track direction
	  DVector3 phat(cos(S(state_phi)),sin(S(state_phi)),S(state_tanl));
	  phat.SetMag(1.);

	  // path length increment
	  double ds=d/sc_norm[index][m].Dot(phat);
	  s-=ds;

	  // Flight time   
	  double q_over_p=S(state_q_over_pt)*cos(atan(S(state_tanl)));
	  double q_over_p_sq=q_over_p*q_over_p;
	  double one_over_beta2=1.+mass2*q_over_p_sq;
	  if (one_over_beta2>BIG) one_over_beta2=BIG;
	  t-=ds*sqrt(one_over_beta2); // in units where c=1
	  
	  // Step along the trajectory using d to estimate path length 
	  FastStep(xy,-ds,0.,S);
	  // Find the index for the nearest start counter paddle
	  double dphi=xy.Phi()-SC_PHI_SECTOR1;
	  if (dphi<0) dphi+=2.*M_PI;
	  index=int(floor(dphi/(2.*M_PI/30.)));
	  if (index>29) index=0;  

	  // Find the new distance to the start counter (which could now be to
	  // a plane in the one adjacent to the one before the step...)
	  d=sc_norm[index][m].Dot(DVector3(xy.X(),xy.Y(),S(state_z))
				  -sc_pos[index][m]);
	  if (fabs(d)<fabs(dmin)){
	    bestS=S;
	    dmin=d;
	    bestXY=xy;
	  }
	  count++;
	}

	if (bestS(state_z)>sc_pos[0][0].z()-0.1){
	  double tanl=bestS(state_tanl);
	  double pt=1/fabs(bestS(state_q_over_pt));
	  double phi=bestS(state_phi);
	  DVector3 position(bestXY.X(),bestXY.Y(),bestS(state_z));
	  DVector3 momentum(pt*cos(phi),pt*sin(phi),pt*tanl); 
	  extrapolations[SYS_START].push_back(Extrapolation_t(position,momentum,
	  					   t*TIME_UNIT_CONVERSION,s));
	  //printf("Central track:\n");
	  //position.Print();
	}  
	break;
      }
    }
  }

   // Accumulate multiple-scattering terms for use in matching routines
  double s_theta_ms_sum=0.,theta2ms_sum=0.;
  for (unsigned int k=inner_index;k>index_beyond_start_counter;k--){
    double ds_theta_ms_sq=3.*fabs(central_traj[k].Q(state_D,state_D));
    s_theta_ms_sum+=sqrt(ds_theta_ms_sq);  
    double ds=central_traj[k].s-central_traj[k-1].s;
    theta2ms_sum+=ds_theta_ms_sq/(ds*ds);
  }
  
  // Deal with points within fiducial volume of chambers
  mT0Detector=SYS_NULL;
  mT0MinimumDriftTime=1e6;
  for (int k=index_beyond_start_counter;k>=0;k--){ 
    S=central_traj[k].S;
    xy=central_traj[k].xy;
    double t=central_traj[k].t*TIME_UNIT_CONVERSION; // convert to ns
    double s=central_traj[k].s;
    double tanl=S(state_tanl);
    double pt=1/fabs(S(state_q_over_pt));
    double phi=S(state_phi); 

    // Find estimate for t0 using earliest drift time
    if (central_traj[k].h_id>0){
      unsigned int index=central_traj[k].h_id-1;
      double dt=my_cdchits[index]->tdrift-t;  
      if (dt<mT0MinimumDriftTime){
	mT0MinimumDriftTime=dt;
	mT0Detector=SYS_CDC;
      }
    }

    //multiple scattering terms
    if (k>0){
      double ds_theta_ms_sq=3.*fabs(central_traj[k].Q(state_D,state_D));
      s_theta_ms_sum+=sqrt(ds_theta_ms_sq);  
      double ds=central_traj[k].s-central_traj[k-1].s;
      theta2ms_sum+=ds_theta_ms_sq/(ds*ds);
    }
    if (S(state_z)<endplate_z){
      DVector3 position(xy.X(),xy.Y(),S(state_z));
      DVector3 momentum(pt*cos(phi),pt*sin(phi),pt*tanl); 
      extrapolations[SYS_CDC].push_back(Extrapolation_t(position,momentum,t,s,
					       s_theta_ms_sum,theta2ms_sum));
      
    }
  }

  //------------------------------
  // Next swim to outer detectors
  //------------------------------
  S=central_traj[0].S;
 
  // Position and step variables 
  xy=central_traj[0].xy;
  double r2=xy.Mod2();
  double ds=mStepSizeS; // step along path in cm
  
  // Energy loss
  double dedx=0.;
  
  // Current time and path length
  double t=central_traj[0].t;
  double s=central_traj[0].s;

  // Matrix for multiple scattering covariance terms
  DMatrix5x5 Q;

  // Track propagation loop
  //if (false)
  while (S(state_z)>0. && S(state_z)<Z_MAX  
            && r2<89.*89.){  
    // Bail if the transverse momentum has dropped below some minimum
    if (fabs(S(state_q_over_pt))>Q_OVER_PT_MAX){
      if (DEBUG_LEVEL>2)
	{
	  _DBG_ << "Bailing: PT = " << 1./fabs(S(state_q_over_pt))
                    << endl;
	}
      return VALUE_OUT_OF_RANGE;
    }
    
    // get material properties from the Root Geometry
    double rho_Z_over_A=0.,LnI=0.,K_rho_Z_over_A=0.,Z=0.;
    double chi2c_factor=0.,chi2a_factor=0.,chi2a_corr=0.;
    DVector3 pos3d(xy.X(),xy.Y(),S(state_z));
    double s_to_boundary=0.;
    DVector3 dir(cos(S(state_phi)),sin(S(state_phi)),S(state_tanl));
    if (geom->FindMatKalman(pos3d,dir,K_rho_Z_over_A,rho_Z_over_A,LnI,Z,
			    chi2c_factor,chi2a_factor,chi2a_corr,
			    last_material_map,&s_to_boundary)
	!=NOERROR){
      _DBG_ << "Material error in ExtrapolateToVertex! " << endl;
      _DBG_ << " Position (x,y,z)=("<<pos3d.x()<<","<<pos3d.y()<<","
	    << pos3d.z()<<")"
	    <<endl;
      break;
    }
    
    // Get dEdx for the upcoming step
    double q_over_p=S(state_q_over_pt)*cos(atan(S(state_tanl)));
    if (CORRECT_FOR_ELOSS){
      dedx=GetdEdx(q_over_p,K_rho_Z_over_A,rho_Z_over_A,LnI,Z); 
    }
    // Adjust the step size
    if (fabs(dedx)>EPS){
      ds=DE_PER_STEP/fabs(dedx);
    }

    if (ds>mStepSizeS) ds=mStepSizeS;
    if (s_to_boundary<ds) ds=s_to_boundary;
    if (ds<MIN_STEP_SIZE)ds=MIN_STEP_SIZE; 
    if (ds<0.5 && S(state_z)<400. && pos3d.Perp()>65.) ds=0.5;
    s+=ds;
    
    // Flight time
    double q_over_p_sq=q_over_p*q_over_p;
    double one_over_beta2=1.+mass2*q_over_p_sq;
    if (one_over_beta2>BIG) one_over_beta2=BIG;
    t+=ds*sqrt(one_over_beta2); // in units where c=1
    
    // Multiple scattering
    GetProcessNoiseCentral(ds,chi2c_factor,chi2a_factor,chi2a_corr,S,Q);
   
    double ds_theta_ms_sq=3.*fabs(Q(state_D,state_D));
    s_theta_ms_sum+=sqrt(fabs(Q(state_D,state_D)));
    theta2ms_sum+=ds_theta_ms_sq/(ds*ds);

    // Propagate the state through the field
    Step(xy,ds,S,dedx);
     
    r2=xy.Mod2(); 
    // Check if we have passed into the BCAL
    if (r2>64.9*64.9){  
      if (extrapolations.at(SYS_BCAL).size()>299){
	return VALUE_OUT_OF_RANGE;
      }
      if (S(state_z)<406.&&S(state_z)>17.0){
	double tanl=S(state_tanl);
	double pt=1/fabs(S(state_q_over_pt));
	double phi=S(state_phi);
	DVector3 position(xy.X(),xy.Y(),S(state_z));   
	DVector3 momentum(pt*cos(phi),pt*sin(phi),pt*tanl);
	extrapolations[SYS_BCAL].push_back(Extrapolation_t(position,momentum,
							   t*TIME_UNIT_CONVERSION,s));
      }
      else if (extrapolations.at(SYS_BCAL).size()<5){
	// There needs to be some steps inside the the volume of the BCAL for 
	// the extrapolation to be useful.  If this is not the case, clear 
	// the extrolation vector.
	extrapolations[SYS_BCAL].clear();
      }
    }
  }   

  return NOERROR;
}

/*---------------------------------------------------------------------------*/


