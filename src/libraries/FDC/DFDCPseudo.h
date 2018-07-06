//***********************************************************************
// DFDCPseudo.h : definition for a set of FDCHits that have gone 
// through first-order reconstruction. 
// Author: Craig Bookwalter (craigb at jlab.org)
// Date:	March 2006
//***********************************************************************

#ifndef DFDCPSEUDO_H
#define DFDCPSEUDO_H

#include <JANA/JObject.h>
using namespace jana;

#include "DFDCWire.h"
#include "DFDCCathodeCluster.h"
#include <sstream>
#include <vector>
#include <DVector2.h>
#include <DMatrix.h>
#include <DMatrixSIMD.h>

typedef struct {
   double pos;
   double q;
   double q_from_pulse_height;
   int numstrips;
   double t; // mean time of strips in peak
   double t_rms; // rms of strips in peak
   unsigned int cluster; // index for cluster from which this centroid was generated
   DMatrix3x1 X,N,NRaw,index;
}centroid_t;

enum FDCTrackD {
   dDOCAW_dDeltaX=0,
   dDOCAW_dDeltaZ,
   dDOCAW_dDeltaPhiX,
   dDOCAW_dDeltaPhiY,
   dDOCAW_dDeltaPhiZ,
   dW_dt0,
   dDOCAW_dx,
   dDOCAW_dy,
   dDOCAW_dtx,
   dDOCAW_dty,
   dDOCAC_dDeltaX,
   dDOCAC_dDeltaZ,
   dDOCAC_dDeltaPhiX,
   dDOCAC_dDeltaPhiY,
   dDOCAC_dDeltaPhiZ,
   dDOCAC_dx,
   dDOCAC_dy,
   dDOCAC_dtx,
   dDOCAC_dty,
   size
};

enum FDCPseudoD {
   dWcddeltaU=0,
   dWcddeltaV,
   dWcddeltaPhiU,
   dWcddeltaPhiV,
   dSddeltaU,
   dSddeltaV,
   dSddeltaPhiU,
   dSddeltaPhiV,
   dWdX,
   dSdX,
   dWddeltaPhiWz,
   dWddeltaPhiWy
};

///
/// class DFDCPseudo: definition for a reconstructed point in the FDC
/// 
class DFDCPseudo : public JObject {
   public :
      JOBJECT_PUBLIC(DFDCPseudo);			/// DANA identifier

      /// 
      /// DFDCPseudo::DFDCPseudo():
      /// Default constructor-- provide the X, Y, global layer #, and resolution
      ///
      DFDCPseudo(){}


      double u,v; ///< centroid positions in the two cathode views
      double t_u,t_v; ///< time of the two cathode clusters
      double phi_u,phi_v; ///< rotation angles for cathode planes
      centroid_t cluster_u, cluster_v; ///< Cathode cluster centroids, Useful for gain/strip pitch calibration.
      double w,dw; ///< local coordinate of pseudopoint in the direction perpendicular to the wires and its uncertainty
      double w_c; /// < wire position computed from cathode data, assuming the avalanche occurs at the wire
      double s,ds; ///< local coordinate of pseudopoint in the direction along the wire and its uncertainty
      const DFDCWire* wire; ///< DFDCWire for this wire 
      double time; ///< time corresponding to this pseudopoint.
      int status; ///< status word for pseudopoint
      double covxx,covxy,covyy; ///< Covariance terms for (x,y) 
      double dE; ///< energy deposition, from integral 
      double dE_amp; /// < energy deposition, from pulse height
      double q; ///< anode charge deduced from cathode strips
      int itrack;
      DVector2 xy; ///< rough x,y coordinates in lab coordinate system

      void toStrings(vector<pair<string,string> > &items)const{ 
         AddString(items,"u","%3.2f",u);
         AddString(items,"v","%3.2f",v);
         AddString(items,"t_u","%3.2f",t_u);
         AddString(items,"t_v","%3.2f",t_v);
         AddString(items,"phi_u","%3.2f",phi_u);
         AddString(items,"phi_v","%3.2f",phi_v);
         AddString(items, "w", "%3.4f", w);
         AddString(items, "w_c", "%3.4f", w_c);
         AddString(items, "s", "%3.4f", s);
         AddString(items, "layer", "%d", wire->layer);
         AddString(items, "wire", "%d", wire->wire);
         AddString(items, "time", "%3.1f", time);
         AddString(items, "status", "%d", status);
         AddString(items, "x", "%.4f", xy.X());
         AddString(items, "y", "%.4f", xy.Y());
         AddString(items, "dE", "%3.1f", dE);
      }

      // For alignment purposes the residuals wrt the alignment parameters are needed.
      // These routines calculate the derivatives of the "s" and "w" variables wrt the alignment parameters
      // Currently implemented alignment parameters are dU, dV, dPhiU, dPhiV, dX (Acts in direction w), dY (acts in direction s)

      vector<double> GetFDCPseudoAlignmentDerivatives(){
         // Create the storage vector for each of the derivatives
         size_t nDerivatives = 12;
         vector<double> derivatives(nDerivatives);

         // Useful numbers...
         double sinPhiU = sin(phi_u);
         double cosPhiU = cos(phi_u);
         double sinPhiV = sin(phi_v);
         double cosPhiV = cos(phi_v);
         double sinPhiUmPhiV = sin(phi_u-phi_v);
         double sinPhiUmPhiV2 = sinPhiUmPhiV*sinPhiUmPhiV;
         double cosPhiUmPhiV = cos(phi_u-phi_v);

         // Calculate the derivatives
         derivatives[dWcddeltaU] = sinPhiV/sinPhiUmPhiV;
         derivatives[dWcddeltaV] = -sinPhiU/sinPhiUmPhiV;
         derivatives[dWcddeltaPhiU] = (v-u*cosPhiUmPhiV)*sinPhiV/sinPhiUmPhiV2;
         derivatives[dWcddeltaPhiV] = (u-v*cosPhiUmPhiV)*sinPhiU/sinPhiUmPhiV2;

         derivatives[dSddeltaU] = -cosPhiV/sinPhiUmPhiV;
         derivatives[dSddeltaV] = cosPhiU/sinPhiUmPhiV;
         derivatives[dSddeltaPhiU] = -(v-u*cosPhiUmPhiV)*cosPhiV/sinPhiUmPhiV2;
         derivatives[dSddeltaPhiV] = -(u-v*cosPhiUmPhiV)*cosPhiU/sinPhiUmPhiV2;

         derivatives[dWdX]=1.0;
         derivatives[dSdX]=1.0;

         derivatives[dWddeltaPhiWz]=-s; // For small angles. Rotations of the wire plane about the z axis
         derivatives[dWddeltaPhiWy]=-w; // Rotations about the wire axis

         return derivatives;

      }

      vector<double> GetFDCStripGainDerivatives() {
         // Numerically calculate the derivatives of the cathode position wrt the strip gains
         // Numerical calculation necessary since it comes from the solution of a nonlinear set of equations
         vector<double> derivatives; // u and v

         // Loop over the cluster and find the three hit section used
         double positionNominal, positionShifted;
         FindCentroid(cluster_u.N, cluster_u.X, positionNominal);

         double delta = 0.05;

         DMatrix3x1 deltaVectU0(delta*cluster_u.NRaw(0), 0.0, 0.0);
         DMatrix3x1 deltaVectU1(0.0,delta*cluster_u.NRaw(1), 0.0);
         DMatrix3x1 deltaVectU2(0.0,0.0,delta*cluster_u.NRaw(2));

         DMatrix3x1 deltaVectV0(delta*cluster_v.NRaw(0), 0.0, 0.0);
         DMatrix3x1 deltaVectV1(0.0,delta*cluster_v.NRaw(1), 0.0);
         DMatrix3x1 deltaVectV2(0.0,0.0,delta*cluster_v.NRaw(2));

         FindCentroid(cluster_u.N+deltaVectU0, cluster_u.X, positionShifted);
         derivatives.push_back((positionShifted-positionNominal)/delta);

         FindCentroid(cluster_u.N+deltaVectU1, cluster_u.X, positionShifted);
         derivatives.push_back((positionShifted-positionNominal)/delta);

         FindCentroid(cluster_u.N+deltaVectU2, cluster_u.X, positionShifted);
         derivatives.push_back((positionShifted-positionNominal)/delta);

         FindCentroid(cluster_v.N, cluster_v.X, positionNominal);

         FindCentroid(cluster_v.N+deltaVectV0, cluster_v.X, positionShifted);
         derivatives.push_back((positionShifted-positionNominal)/delta);

         FindCentroid(cluster_v.N+deltaVectV1, cluster_v.X, positionShifted);
         derivatives.push_back((positionShifted-positionNominal)/delta);

         FindCentroid(cluster_v.N+deltaVectV2, cluster_v.X, positionShifted);
         derivatives.push_back((positionShifted-positionNominal)/delta);

         //jout << " New Hit" << endl;
         //for (size_t i=0 ; i < derivatives.size(); i++){
         //   jout << "i " << i << " der " << derivatives[i] << endl;
         //}

         return derivatives;
      }

      const vector<double> GetFDCStripPitchDerivatives(){
         // Numerically calculate the derivatives of the cathode position wrt the strip gains
         // Numerical calculation necessary since it comes from the solution of a nonlinear set of equations
         vector<double> derivatives; // u and v

         double positionNominal, positionShifted;
         FindCentroid(cluster_u.N, cluster_u.X, positionNominal);

         double delta = 0.005;
         DMatrix3x1 deltaVect(0.0, 0.0, 0.0);
         // delta_U_SP_1
         for (int i=0 ; i < 3; i++){
            if(cluster_u.index(i) < 48) deltaVect(i) = (cluster_u.index(i)-47.)*delta;
            else deltaVect(i) = 0.0;
         }  

         FindCentroid(cluster_u.N, cluster_u.X + deltaVect, positionShifted);
         derivatives.push_back((positionShifted-positionNominal)/delta);

         // delta_U_G_1
         for (int i=0 ; i < 3; i++){
            if(cluster_u.index(i) < 48) deltaVect(i) = -delta;
            else deltaVect(i) = 0.0;
         }

         FindCentroid(cluster_u.N, cluster_u.X + deltaVect, positionShifted);
         derivatives.push_back((positionShifted-positionNominal)/delta);

         // delta_U_SP_2
         for (int i=0 ; i < 3; i++){
            if(cluster_u.index(i) < 48) deltaVect(i) = -47.5*delta;
            else if (cluster_u.index(i) < 144) deltaVect(i) = (cluster_u.index(i)-95.5)*delta;
            else deltaVect(i) = 47.5*delta;
         }

         FindCentroid(cluster_u.N, cluster_u.X + deltaVect, positionShifted);
         derivatives.push_back((positionShifted-positionNominal)/delta);

         // delta_U_G_2
         for (int i=0 ; i < 3; i++){
            if(cluster_u.index(i) >= 144) deltaVect(i) = delta;
            else deltaVect(i) = 0.0;
         }

         FindCentroid(cluster_u.N, cluster_u.X + deltaVect, positionShifted);
         derivatives.push_back((positionShifted-positionNominal)/delta);

         // delta_U_SP_3
         for (int i=0 ; i < 3; i++){
            if(cluster_u.index(i) >= 144) deltaVect(i) = (cluster_u.index(i)-144.)*delta;
            else deltaVect(i) = 0.0;
         }

         FindCentroid(cluster_u.N, cluster_u.X + deltaVect, positionShifted);
         derivatives.push_back((positionShifted-positionNominal)/delta);

         //==========
         // Get the nominal V shift
         FindCentroid(cluster_v.N, cluster_v.X, positionNominal);

         // delta_V_SP_1
         for (int i=0 ; i < 3; i++){
            if(cluster_v.index(i) < 48) deltaVect(i) = (cluster_v.index(i)-47.)*delta;
            else deltaVect(i) = 0.0;
         }

         FindCentroid(cluster_v.N, cluster_v.X + deltaVect, positionShifted);
         derivatives.push_back((positionShifted-positionNominal)/delta);

         // delta_V_G_1
         for (int i=0 ; i < 3; i++){
            if(cluster_v.index(i) < 48) deltaVect(i) = -delta;
            else deltaVect(i) = 0.0;
         }

         FindCentroid(cluster_v.N, cluster_v.X + deltaVect, positionShifted);
         derivatives.push_back((positionShifted-positionNominal)/delta);

         // delta_V_SP_2
         for (int i=0 ; i < 3; i++){
            if(cluster_v.index(i) < 48) deltaVect(i) = -47.5*delta;
            else if (cluster_v.index(i) < 144) deltaVect(i) = (cluster_v.index(i)-95.5)*delta;
            else deltaVect(i) = 47.5*delta;
         }

         FindCentroid(cluster_v.N, cluster_v.X + deltaVect, positionShifted);
         derivatives.push_back((positionShifted-positionNominal)/delta);

         // delta_V_G_2
         for (int i=0 ; i < 3; i++){
            if(cluster_v.index(i) >= 144) deltaVect(i) = delta;
            else deltaVect(i) = 0.0;
         }

         FindCentroid(cluster_v.N, cluster_v.X + deltaVect, positionShifted);
         derivatives.push_back((positionShifted-positionNominal)/delta);

         // delta_V_SP_3
         for (int i=0 ; i < 3; i++){
            if(cluster_v.index(i) >= 144) deltaVect(i) = (cluster_v.index(i)-144.)*delta;
            else deltaVect(i) = 0.0;
         }

         FindCentroid(cluster_v.N, cluster_v.X + deltaVect, positionShifted);
         derivatives.push_back((positionShifted-positionNominal)/delta);

         return derivatives;

      }

      jerror_t FindCentroid(DMatrix3x1 N, DMatrix3x1 X, double &pos){

         unsigned int X0=0;
         unsigned int QA=1;
         unsigned int K2=2;
         int ITER_MAX=100;
         float TOLX=1e-4;
         float TOLF=1e-4;
         float A_OVER_H=0.4;
         float ONE_OVER_H=2.0;

         // Define some matrices for use in the Newton-Raphson iteration
         DMatrix3x3 J;  //Jacobean matrix
         DMatrix3x1 F,par,newpar,f;

         double sum=0.;

         for (int j=0; j<3; j++){
            sum+=N(j);
         }

         // Starting values for the minimization
         par(X0)=X(1); // X0
         par(QA)=2.*sum; // QA
         par(K2)=1.; // K2
         newpar=par;

         // Newton-Raphson procedure
         double errf=0.,errx=0;
         for (int iter=1;iter<=ITER_MAX;iter++){
            errf=0.;
            errx=0.;

            // Compute Jacobian matrix: J_ij = dF_i/dx_j.
            for (int i=0;i<3;i++){
               double dx=(par(X0)-X(i))*ONE_OVER_H;
               double argp=par(K2)*(dx+A_OVER_H);
               double argm=par(K2)*(dx-A_OVER_H);
               double tanh_p=tanh(argp);
               double tanh_m=tanh(argm);
               double tanh2_p=tanh_p*tanh_p;
               double tanh2_m=tanh_m*tanh_m;
               double q_over_4=0.25*par(QA);

               f(i)=tanh_p-tanh_m;
               J(i,QA)=-0.25*f(i);
               J(i,K2)=-q_over_4*(argp/par(K2)*(1.-tanh2_p)
                     -argm/par(K2)*(1.-tanh2_m));
               J(i,X0)=-q_over_4*par(K2)*(tanh2_m-tanh2_p); 
               F(i)=N(i)-q_over_4*f(i);
               double new_errf=fabs(F(i));
               if (new_errf>errf) errf=new_errf;
            }
            // Check for convergence
            if (errf<TOLF){
               pos = par(X0);
               return NOERROR;
            }

            // Find the new set of parameters
            FindNewParmVec(N,X,F,J,par,newpar);

            //Check for convergence
            for (int i=0;i<3;i++){
               double new_err=fabs(par(i)-newpar(i));
               if (new_err>errx) errx=new_err;
            }
            if (errx<TOLX){
               pos = par(X0);
               return NOERROR;
            }
            par=newpar;
         } // iterations

         return INFINITE_RECURSION; // error placeholder
      }

      jerror_t FindNewParmVec(const DMatrix3x1 &N,
            const DMatrix3x1 &X,
            const DMatrix3x1 &F,
            const DMatrix3x3 &J,
            const DMatrix3x1 &par,
            DMatrix3x1 &newpar){

         unsigned int X0=0;
         unsigned int QA=1;
         unsigned int K2=2;
         float A_OVER_H=0.4;
         float ONE_OVER_H=2.0;
         float ALPHA=1e-4; // rate parameter for Newton step backtracking algorithm

         // Invert the J matrix
         DMatrix3x3 InvJ=J.Invert();

         // Find the full Newton step
         DMatrix3x1 fullstep=(-1.)*(InvJ*F);

         // find the rate of decrease for the Newton-Raphson step
         double slope=(-1.0)*F.Mag2(); //dot product

         // This should be a negative number...
         if (slope>=0){
            return VALUE_OUT_OF_RANGE;
         }

         double lambda=1.0;  // Start out with full Newton step
         double lambda_temp,lambda2=lambda;
         DMatrix3x1 newF;
         double f2=0.,newf;

         // Compute starting values for f=1/2 F.F 
         double f=-0.5*slope;

         for (;;){
            newpar=par+lambda*fullstep;

            // Compute the value of the vector F and f=1/2 F.F with the current step
            for (int i=0;i<3;i++){
               double dx=(newpar(X0)-X(i))*ONE_OVER_H;
               double argp=newpar(K2)*(dx+A_OVER_H);
               double argm=newpar(K2)*(dx-A_OVER_H);
               newF(i)=N(i)-0.25*newpar(QA)*(tanh(argp)-tanh(argm));
            }
            newf=0.5*newF.Mag2(); // dot product

            if (lambda<0.1) {  // make sure the step is not too small
               newpar=par;
               return NOERROR;
            } // Check if we have sufficient function decrease
            else if (newf<=f+ALPHA*lambda*slope){
               return NOERROR;
            }
            else{
               // g(lambda)=f(par+lambda*fullstep)
               if (lambda==1.0){//first attempt: quadratic approximation for g(lambda)
                  lambda_temp=-0.5*slope/(newf-f-slope);
               }
               else{ // cubic approximation for g(lambda)
                  double temp1=newf-f-lambda*slope;
                  double temp2=f2-f-lambda2*slope;
                  double one_over_lambda2_sq=1./(lambda2*lambda2);
                  double one_over_lambda_sq=1./(lambda*lambda);
                  double one_over_lambda_diff=1./(lambda-lambda2);
                  double a=(temp1*one_over_lambda_sq-temp2*one_over_lambda2_sq)*one_over_lambda_diff;
                  double b=(-lambda2*temp1*one_over_lambda_sq+lambda*temp2*one_over_lambda2_sq)
                     *one_over_lambda_diff;
                  if (a==0.0) lambda_temp=-0.5*slope/b;
                  else{
                     double disc=b*b-3.0*a*slope;
                     if (disc<0.0) lambda_temp=0.5*lambda;
                     else if (b<=0.0) lambda_temp=(-b+sqrt(disc))/(3.*a);
                     else lambda_temp=-slope/(b+sqrt(disc));
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

};

#endif //DFDCPSEUDO_H
