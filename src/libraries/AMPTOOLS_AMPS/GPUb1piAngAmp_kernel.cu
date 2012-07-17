/*
 *  GPUb1piAngAmp_kernel.cu
 *
 */


#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "GPUUtils/lorentzBoost.cuh"
#include "GPUUtils/threeVector.cuh"
#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"

#include "AMPTOOLS_AMPS/breakupMomentum.cuh"
#include "AMPTOOLS_AMPS/barrierFactor.cuh"


#define ADD4(a,b) { a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3] }

#define MASS(v)   (G_SQRT(v[0]*v[0]-v[1]*v[1]-v[2]*v[2]-v[3]*v[3]))

#define Nterm(J)  (G_SQRT((2*J+1)/(4*M_PI)))


// Macro to ease definition of loops
#define LOOP(INDEX,START,END,INC) for (int INDEX=START;INDEX<=END;INDEX+=INC)


static __device__ void //note: 4-vector input presumed
rotateZ( GDouble* v, GDouble phi ){
  GDouble sinphi = G_SIN(phi);
  GDouble cosphi = G_COS(phi);
  GDouble tx;
  tx   = v[1] * cosphi - v[2] * sinphi;
  v[2] = v[2] * cosphi + v[1] * sinphi;
  v[1] = tx;
}

static __device__ void //note: 4-vector input presumed
rotateY ( GDouble* v, GDouble theta) {
  double sinphi = G_SIN(theta);
  double cosphi = G_COS(theta);
  double tz;
  tz = v[3] * cosphi - v[1] * sinphi;
  v[1] = v[1] * cosphi + v[3] * sinphi;
  v[3] = tz;
}

static __device__ GDouble  //note: 3-vector input presumed
theta( GDouble* pv ){
  GDouble r= G_SQRT(pv[0]*pv[0] + pv[1]*pv[1]);
  return ( ( pv[2] == 0 ) && ( r == 0 ) ? 0 : G_ATAN2( r , pv[2] ) );
}


static __device__ void
MoveToRF(GDouble *parent, GDouble *daughter)
{
  GDouble *par3vec=parent+1;
  rotateZ( daughter , -phi(par3vec) );
  rotateY( daughter , -theta(par3vec) );

  GDouble beta[]={0,0, -G_SQRT(dot(par3vec,par3vec))/parent[0]};
  boost( daughter , beta );

}



static __device__ WCUComplex
BreitWigner_loc(GDouble m0, GDouble Gamma0, int L,
                GDouble *P1, GDouble *P2)
{
  
  GDouble  Ptot[4] = ADD4(P1, P2);  
  GDouble m  = MASS(Ptot);
  GDouble mass1 = MASS(P1);
  GDouble mass2 = MASS(P2);
  
  
  // assert positive breakup momenta     
  GDouble q0 = fabs( breakupMomentum(m0, mass1, mass2) );
  GDouble q  = fabs( breakupMomentum(m,  mass1, mass2) );
  
  GDouble F0 = L==0 ? 1.0 : barrierFactor(q0, L);
  GDouble F  = L==0 ? 1.0 : barrierFactor(q,  L);
  
  GDouble width_coef=Gamma0*(m0/m);
  //GDouble qq0=q/q0;
  //GDouble width_qdep = (L==0 ? qq0 : (L==1 ? qq0*qq0*qq0 : pow(qq0,2*L+1)))*((F*F)/(F0*F0));
  GDouble width_qdep = q/q0  * (F*F)/(F0*F0);
  //GDouble num_qdep = (L==0 ? q : (L==1 ? q*q*q : pow(q,2*L+1)))*(F*F);
  GDouble num_qdep = q*(F*F);
  
  GDouble width = width_coef * width_qdep;
  
  //complex<GDouble> bwtop(m0 * width, 0.0 );
  WCUComplex bwtop = { G_SQRT(m0*width_coef) * num_qdep, 0 };
  
  WCUComplex bwbottom =  { m0*m0 - m*m  , -1.0 * ( m0 * width ) };
  
  return  ( bwtop / bwbottom );
  
}



// 2012-07-12 JR  Remove polFrac in parameter list.
__global__ void
GPUb1piAngAmp_kernel( GPU_AMP_PROTO , int polBeam, 
  int J_X, int Par_X, int L_X, int I_X, int epsilon_R, int Iz_b1, int Iz_pi,
  GDouble u_rho_1, GDouble u_rho_3, GDouble u_omega_1, GDouble u_omega_3,
  GDouble u_b1_0, GDouble u_b1_2, 
  GDouble G0_omega, GDouble G0_b1, bool orthocheck) 
{
  // Calculate event for this thread.
  int iEvent = GPU_THIS_EVENT;

  WCUComplex CZero = { 0, 0 };
  WCUComplex i =     { 0, 1 };
  WCUComplex COne =  { 1, 0 };

  int pol=(polBeam==1 ? +1 : -1); // y and x-pol. respectively
  
  if (J_X==0 && Par_X*pol*epsilon_R==-1) {
    pcDevAmp[iEvent] = CZero;
    return;
  }



  int m_X;
  GDouble u_rho, u_omega, u_b1;
  GDouble InvSqrt2 = 1.0/G_SQRT(2.0);
  GDouble m0_rho   = 0.775;
  GDouble G0_rho   = 0.149;
  GDouble m0_omega = 0.783;
  GDouble m0_b1    = 1.223;
  bool useCutoff   = true;
  bool isZero;



  //  Copy four-vectors for this thread from global memory.
  //  2012-05-19 JR  rhos_pip0,omega0,rho0 added for use
  //                 in BreitWigner_loc() below.
  GDouble  beam     [4] = GPU_P4(0);
  GDouble  recoil   [4] = GPU_P4(1);
  GDouble  Xs_pi    [4] = GPU_P4(2);
  GDouble  b1s_pi   [4] = GPU_P4(3);
  GDouble  omegas_pi[4] = GPU_P4(4);
  GDouble  rhos_pim [4] = GPU_P4(5);
  GDouble  rhos_pip [4] = GPU_P4(6);
  GDouble  rhos_pip0[4] = GPU_P4(6);

  //  Make four-vector sums
  GDouble  rho   [4] = ADD4(rhos_pip, rhos_pim );
  GDouble  rho0  [4] = ADD4(rhos_pip, rhos_pim );
  GDouble  omega [4] = ADD4(rho,     omegas_pi);
  GDouble  omega0[4] = ADD4(rho,     omegas_pi);
  GDouble  b1    [4] = ADD4(omega,   b1s_pi);

  //  Store mass of b1; for other vectors we can calculate mass on the fly.
  GDouble b1mass = MASS(b1);

  //  Is this term zero?
  if (useCutoff) {
      isZero  = MASS(rho)+0.135            > m0_omega+3*G0_omega;
      isZero |= fabs(MASS(omega)-m0_omega) > 3*G0_omega;
      isZero |= fabs(b1mass-m0_b1)         > 3*G0_b1;
      isZero |= b1mass                     < (m0_omega - 3*G0_omega);
      // Zero amplitude
      if (isZero) {
         pcDevAmp[iEvent] = CZero;
         return;
      }
  }

  // Continue to Calculate amplitude
  GDouble  X[4] = ADD4(b1,      Xs_pi);

  GDouble q = breakupMomentum( MASS(X), b1mass, MASS(Xs_pi) );

  GDouble alpha = phi( &(recoil[1]) );

  //  NOTE: Values of beam and recoil are changed below.
  boostToRest (beam,   X);
  boostToRest (recoil, X);

  //  Define new coordinate system with 
  //    - beam parallel to z direction
  //    - recoil in the x,z plain (i.e., y is normal to recoil and beam)
  //    - y is normal to beam and recoil.
  GDouble zGJ[3] = { beam[1], beam[2], beam[3] };
  makeUnit( zGJ );
  
  GDouble yGJ[3] = { -recoil[1], -recoil[2], -recoil[3] };
  cross( yGJ, zGJ );
  makeUnit( yGJ );
  
  GDouble xGJ[3] = { yGJ[0], yGJ[1], yGJ[2] };
  cross( xGJ, zGJ );

  //particles to rest frames of their parents
  boostToRest (b1,      X);
  boostToRest (omega,   X);
  boostToRest (rho,     X);
  boostToRest (rhos_pip, X);
 
  // Note that in this form of the cascade of boosts, we are not
  // saving the 4-vecs in their intermediate RF, but going sequentially
  // straight to their immediate parent's RF.
  MoveToRF(b1,omega);
  MoveToRF(b1,rho);      MoveToRF(omega,rho);
  MoveToRF(b1,rhos_pip); MoveToRF(omega,rhos_pip); MoveToRF(rho,rhos_pip);
  
  GDouble *b1_3vec=b1+1;
  GDouble ang_b1[]={dot(b1_3vec, xGJ),
                    dot(b1_3vec, yGJ),
                    dot(b1_3vec, zGJ)};
  GDouble b1_XRF_cosTheta = cosTheta(ang_b1);
  GDouble b1_XRF_phi      = phi(ang_b1);
   
  GDouble rho_omegaRF_cosTheta = cosTheta(rho+1);
  GDouble rho_omegaRF_phi      = phi(rho+1);
  GDouble rhos_pip_rhoRF_cosTheta = cosTheta(rhos_pip+1);
  GDouble rhos_pip_rhoRF_phi      = phi(rhos_pip+1);
  GDouble omega_b1RF_cosTheta     = cosTheta(omega+1);
  GDouble omega_b1RF_phi          = phi(omega+1);


  // SUMMATION GUIDE:
  // notation meant to resemble TeX symbols in derivation
  // exception: pol = \epsilon_\gamma
  // l -> lambda, indicating helicity
  // u_[particle](q.n.) -> amplitude strength coefficient 

  int l_R_lim     = J_X + 1;
  
  //shortcut:  CB(L_X, J_b1, 0, l_b1 ; J_X, l_b1) vanishes when
  //  = CB(1, 1, 0, 0 ; 1, 0),  so omit l_b1=0 when J_X=L_X=1
  int l_b1_inc    = L_X==1 && J_X==1 ? 2 : 1;
  
  // restrict omega decay to just p wave
  int L_omega_lim = 1; // set to 3 to allow F wave
  int L_Rsign_lim;
  
  GDouble cosAlpha=G_COS(alpha), sinAlpha=G_SIN(alpha);
  WCUComplex expFact = {cosAlpha, sinAlpha};
  WCUComplex expFact_conj = {cosAlpha, -sinAlpha};
  
  WCUComplex ThelSum = { 0 , 0 };

  //  Setup dependent loop limits
  LOOP(l_gamma, -1, 1, 2) {

    
    LOOP(l_R, 0, l_R_lim, 1) {
      if(l_R==0 && epsilon_R==-1) continue;

      //summing positive and negative helicity terms of R's reflectivity state
      L_Rsign_lim = l_R > 0 ? -1 : +1;


      LOOP(l_Rsign, L_Rsign_lim, 1, 2) {

        m_X = l_gamma - l_Rsign * l_R;
        if (m_X==0) {
          //testing for cancelation in |J 0>+pol*P*epsilon_R*(-1)^J|J 0>
          if(Par_X*pol*epsilon_R == (J_X % 2 ==0 ? -1:+1)) continue;
        } else {
          //enforcing that the selected projection <= vector magnitude 
          if( abs(m_X)>J_X) continue; 
        }
        
        
        WCUComplex l_b1DepTerm = {0,0};
        LOOP(l_b1, -1,1,l_b1_inc) {
          

          WCUComplex L_b1DepTerm = {0,0};

          LOOP(L_b1,0,2,2) {
            
          
            WCUComplex l_omegaDepTerm = {0,0};

            LOOP(l_omega,-1,1,1) {
              
              WCUComplex L_omegaDepTerm = {0,0};
              LOOP(L_omega, 1, L_omega_lim, 2) { 
                
                WCUComplex J_rhoDepTerm = {0,0};
                LOOP(J_rho, 1, L_omega_lim, 2) {
                  
                  //enforces triang. ineq. betw. J_omega=1, J_rho and L_omega
                  // in effect, L_omega and J_rho take identical values
                  if( abs(J_rho-L_omega) > 1) continue; 
                  
                  
                  WCUComplex l_rhoDepTerm = {0,0};
                  LOOP(l_rho,-1,1,1) {
                    //shortcut CB(1,1,0,0;1,0)=0
                    if(L_omega==1 && J_rho==1 && l_rho==0) continue;
                    
                    l_rhoDepTerm += 
                      Conjugate(wignerD(1, l_omega, l_rho, 
                                        rho_omegaRF_cosTheta, rho_omegaRF_phi))
                      * clebsch(L_omega, 0, J_rho, l_rho, 1, l_rho)
                      * Y(J_rho, l_rho, rhos_pip_rhoRF_cosTheta, rhos_pip_rhoRF_phi);
                  }
                  
                  u_rho = J_rho==1 ? u_rho_1 : (J_rho==3 ? u_rho_3 : 0);
                  J_rhoDepTerm += u_rho * l_rhoDepTerm * 
                    BreitWigner_loc(m0_rho,G0_rho, J_rho,rhos_pip0,rhos_pim);
                }
                
                J_rhoDepTerm *= BreitWigner_loc(m0_omega, G0_omega, L_omega, omegas_pi,rho0);
                
                u_omega = L_omega==1 ? u_omega_1 : (L_omega==3 ? u_omega_3 : 0);
                L_omegaDepTerm += u_omega * J_rhoDepTerm * Nterm(L_omega);
              }
              
              l_omegaDepTerm += L_omegaDepTerm * 
                clebsch(L_b1, 0, 1, l_omega, 1, l_omega) *
                Conjugate(wignerD(1, l_b1, l_omega, 
                                  omega_b1RF_cosTheta, omega_b1RF_phi));
            }
            
            l_omegaDepTerm *= BreitWigner_loc(m0_b1, G0_b1, L_b1, b1s_pi, omega0);
            
            u_b1 = L_b1==0 ? u_b1_0 : (L_b1==2 ? u_b1_2 : 0); 
            L_b1DepTerm += u_b1 * l_omegaDepTerm * Nterm(L_b1);
          }
          
          l_b1DepTerm += L_b1DepTerm *
            Conjugate(wignerD(J_X, m_X, l_b1, b1_XRF_cosTheta, b1_XRF_phi)) *
            clebsch(L_X, 0, 1, l_b1, J_X, l_b1);
        }
        
        ThelSum += l_b1DepTerm  
          //to account for |eps_g> ~ (|1,-1>exp(-ia)-pol|1,+1>exp(ia)) 
          * (l_gamma==1 ? (-pol)*expFact : expFact_conj)
          //Assemble reflectivity eigenvector with epsilon_X=pol*epslion_R
          * (GDouble) (m_X<0 ? Par_X*pol*epsilon_R*((J_X-m_X) % 2 == 0 ? +1:-1) : 1) 
          * (GDouble) (m_X == 0 ? 1.0 : InvSqrt2 )
          // to apply th(l_R) reflectivity state prefactor: 
          // m=0: 1/2  m>0: 1/sqrt(2)  m<0: 0 (last just skipped in this sum)  
          * (GDouble) (l_R > 0 ? InvSqrt2 : 1.0 )
          //apply coefficients to the reflectivity basis terms:
          * (GDouble) (l_Rsign==1 ? 1 : epsilon_R)
          ; //v(*epsilon_R) *     
        
      }
    }
  }
  

  ThelSum *= Nterm(L_X) * 
    // barrier factor
    (GDouble)(L_X==0 ? 1.0 : (L_X==1 ? q : G_POW(q,L_X))) *
    // to apply polarization fraction weights: OBSOLETE! moved to polCoef class
    // (GDouble)G_SQRT((1.0-pol*polFrac)*0.5) * //(1+g) for x-pol, (1-g) for y-pol   
    (pol==1 ? i : COne)*InvSqrt2 * //to account for |eps_g> ~ sqrt(-eps/2)
    clebsch(1, Iz_b1, 1, Iz_pi, I_X, Iz_b1 + Iz_pi);
  pcDevAmp[iEvent] = ThelSum;
  
}


void
GPUb1piAngAmp_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                   int polBeam, //GDouble polFrac,
                   int J_X, int Par_X, int L_X, int I_X, int epsilon_R, 
                   int Iz_b1, int Iz_pi,
                   GDouble u_rho_1, GDouble u_rho_3, 
                   GDouble u_omega_1, GDouble u_omega_3,
                   GDouble u_b1_0, GDouble u_b1_2, 
                   GDouble G0_omega, GDouble G0_b1, bool orthocheck)
{  
  GPUb1piAngAmp_kernel<<< dimGrid, dimBlock >>>
    ( GPU_AMP_ARGS, polBeam,
      J_X, Par_X, L_X, I_X, epsilon_R, Iz_b1, Iz_pi,
      u_rho_1, u_rho_3, u_omega_1, u_omega_3, u_b1_0, u_b1_2, 
      G0_omega, G0_b1, orthocheck ); 

}

