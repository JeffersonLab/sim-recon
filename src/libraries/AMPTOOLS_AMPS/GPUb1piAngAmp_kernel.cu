/*
 *  GPUb1piAngAmp_kernel.cu
 *
 */

#include <stdio.h>
#include "cuda.h"

//  Original headers were scattered around file system
#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "GPUUtils/lorentzBoost.cuh"
#include "GPUUtils/threeVector.cuh"
#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"

#include "AMPTOOLS_AMPS/breakupMomentum.cuh"
#include "AMPTOOLS_AMPS/barrierFactor.cuh"


//  Test headers
#if 0
#include "GPUCustomTypes.h"
#include "CUDA-Complex.cuh"

#include "lorentzBoost.cuh"
#include "threeVector.cuh"
#include "wignerD.cuh"
#include "clebsch.cuh"

#include "breakupMomentum.cuh"
#include "barrierFactor.cuh"
#endif



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
  return  G_ATAN2( r , pv[2] );
}


static __device__ void
MoveToRF(GDouble *parent, GDouble *daughter)
{
  GDouble *par3vec=parent+1;
  rotateZ( daughter , -phi(par3vec) );
  rotateY( daughter , -theta(par3vec) );

  GDouble beta[]={0,0, -G_SQRT(dot(par3vec,par3vec))/parent[0]};
  //** (x)  Might this be bootToRest???
  // beta is defined to boost to parent's rest frame
  // I just adapted GPUUtil boost fcn with vector beta input
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
  
  //printf("BW: (%5.3f, %5.3f, %d) m=%6.4f m1=%6.4f m2=%6.4f q=%6.4f q0=%6.4f\n",
  //  m0,Gamma0,L,m,mass1,mass2,q,q0);
  
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


//  JR 2012-07-29
//  Set all Amplitudes to 0 on the Device.  This is needed now because we only
//  calculate amplitudes for those momenta sets with non-zero amplitudes.  If
//  this function were not performed, amplitudes which are supposed to be zero will
//  be undefined.
__global__ void Setzero_kernel(WCUComplex *pcDevAmp, int iNEvents) {
  int iEvent = GPU_THIS_EVENT;
  if (iEvent>=iNEvents) return;
  pcDevAmp[iEvent].m_dRe = 0.0;
  pcDevAmp[iEvent].m_dIm = 0.0;
}


//  JR 2012-07-29
//  Perform beginning of b1pi calculation, just enough to determine those
//  amplitude which will be set to zero.  Amplitudes are set to (1,0) if 
//  they are not zero.  These amplitudes will need set to their correct
//  values on the call to GPUb1piAngAmp_kernel().
__global__ void Pretest_kernel( GPU_AMP_PROTO , int polBeam, GDouble polFrac,
  int J_X, int Par_X, int L_X, int I_X, int epsilon_R, int Iz_b1, int Iz_pi,
  GDouble u_rho_1, GDouble u_rho_3, GDouble u_omega_1, GDouble u_omega_3,
  GDouble u_b1_0, GDouble u_b1_2, 
  GDouble G0_omega, GDouble G0_b1, bool orthocheck) 
{
  // Calculate event for this thread.
  int iEvent = GPU_THIS_EVENT;
  WCUComplex CZero = { 0, 0 };
  WCUComplex COne  = { 1, 0 };

  int pol=(polBeam==1 ? +1 : -1); // y and x-pol. respectively
  
  //** (x)  This statement can be evaluated at top of function?
  if (J_X==0 && Par_X*pol*epsilon_R==-1) {
    pcDevAmp[iEvent] = CZero;
    return;
  }

  GDouble m0_omega = 0.783;
  GDouble m0_b1    = 1.223;
  bool isZero;


  //  Copy four-vectors for this thread from global memory.
  GDouble  b1s_pi   [4] = GPU_P4(3);
  GDouble  omegas_pi[4] = GPU_P4(4);
  GDouble  rhos_pim [4] = GPU_P4(5);
  GDouble  rhos_pip [4] = GPU_P4(6);

  //  Make four-vector sums
  GDouble  rho   [4] = ADD4(rhos_pip, rhos_pim );
  GDouble  omega [4] = ADD4(rho,     omegas_pi);
  GDouble  b1    [4] = ADD4(omega,   b1s_pi);

  //  Store mass of b1; for other vectors we can calculate mass on the fly.
  GDouble b1mass = MASS(b1);

  //  Is this term zero?
  isZero  = MASS(rho)+0.135            > m0_omega+3*G0_omega;
  isZero |= fabs(MASS(omega)-m0_omega) > 3*G0_omega;
  isZero |= fabs(b1mass-m0_b1)         > 3*G0_b1;
  isZero |= b1mass                     < (m0_omega - 3*G0_omega);
  if    (isZero) pcDevAmp[iEvent] = CZero;
  else           pcDevAmp[iEvent] = COne;
}






//  JR 2012-07-29
//  Calculate amplitudes only for those momenta sets with known non-zero
//  amplitudes.
__global__ void
GPUb1piAngAmp_kernel( 
  int cnt,
  // GPU_AMP_PROTO , 
  GDouble* pfDevData, WCUComplex* pcDevAmp, int* piDevPerm, int iNParticles, int iNEvents,
  int polBeam, GDouble polFrac,
  int J_X, int Par_X, int L_X, int I_X, int epsilon_R, int Iz_b1, int Iz_pi,
  GDouble u_rho_1, GDouble u_rho_3, GDouble u_omega_1, GDouble u_omega_3,
  GDouble u_b1_0, GDouble u_b1_2, 
  GDouble G0_omega, GDouble G0_b1, bool orthocheck) 
{

  // Calculate event for this thread.
  // int iEvent = GPU_THIS_EVENT;

  //  JR 2012-07-29
  //  NOTE:  This vesrsion of this function is called with different settings
  //         for threadIdx, blockIdx and blockDim than for the original version.
  //         The next line relects that change.
  int iEvent = threadIdx.x + blockIdx.x * blockDim.x;

  //  Skip this event index if it overruns number of events. 
  if (iEvent>=iNEvents) return;

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
  
  //**  (x)  Be care of cross order, need to check this
  // 2012-05-19 JR - Invert yGJ to make cross come out right.
  // GDouble yGJ[3] = { recoil[1], recoil[2], recoil[3] };
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
  // Make sure to verify that the intermediares were not in fact needed
  // and that we didn't break anything with this simplification.
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

/*
   List_l_R:        0 1 
   List_J_rho:        1 
   List_l_rho:   -1   1 
   List_L_omega:      1 
   List_l_omega: -1 0 1 
   List_L_b1:       0   2 
   List_l_b1:    -1 0 1 
*/

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
    // LOOP(l_R, (1-epsilon_R)/2, l_R_lim, 1)  // if this still causes some GPU core
      // misalignment, try setting lower bound back to zero and tacking on
      //  * !(l_R==0 && epsilon_R==-1) 
      // to the long list of factors multiplying Thelsum below -IS


      //summing positive and negative helicity terms of R's reflectivity state
      L_Rsign_lim = l_R > 0 ? -1 : +1;
                        // Switch order of loop, because LOOP can only handle increasing increments
      // LOOP(l_Rsign, 1, L_Rsign_lim, -2) 
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
            // 2012-05-19 JR  Fix l_omega loop
            // LOOP(l_omega,-1,0,1) 
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
          //-- (_) understand why assignment here produces:
          // KERNEL LAUNCH ERROR [b1piAngAmp]: the launch timed out and was terminated
          // assigning/incrementing integers causes no problems
          
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
//    (GDouble)(L_X==0 ? 1.0 : (L_X==1 ? q : G_POW(q,L_X))) *
    (GDouble)(L_X==0 ? 1.0 : (L_X==1 ? q : ::pow(q,L_X))) *
    // to apply polarization fraction weights: 
    (GDouble)G_SQRT((1.0-pol*polFrac)*0.5) * //(1+g) for x-pol, (1-g) for y-pol   
    (pol==1 ? i : COne)*InvSqrt2 * //to account for |eps_g> ~ sqrt(-eps/2)
    clebsch(1, Iz_b1, 1, Iz_pi, I_X, Iz_b1 + Iz_pi);
  pcDevAmp[iEvent] = ThelSum;

  
}





#ifdef DEBUG
//   This is for debugging 
//     It reads the amplitdues and momemta vectors from the CUDA device and prints them.
void 
printCudaArrays(GDouble* pfDevData, WCUComplex* pcDevAmp, int* piDevPerm, int iNParticles, int iNEvents, int cnt) {

  //  Read amplitudes from GPU to CPU
	GDouble *amp = (GDouble *) malloc (iNEvents * 2 * sizeof(GDouble));
	cudaMemcpy (amp, pcDevAmp, iNEvents * 2 * sizeof(GDouble), cudaMemcpyDeviceToHost);

  //  Copy 4momenta from GPU to CPU - make part() big enough to hold the entire set of momenta
  GDouble *part = (GDouble *) malloc (iNEvents * 4 * iNParticles * sizeof(GDouble));
  cudaMemcpy (part, pfDevData, iNEvents * 4 * iNParticles * sizeof(GDouble), cudaMemcpyDeviceToHost);

	//  Print arrays
	int ievent, ipart, idim;
	int ndim = 4;
	for (ievent=0; ievent<iNEvents; ievent++) {
		printf ("test: CUDA: %2d %6d ", cnt, ievent);
		//  Print amplitude
		printf ("  %12.4e %12.4e", amp[2*ievent], amp[2*ievent+1]);
		for (ipart=0;ipart<iNParticles;ipart++) {
			printf (" ");
			for (idim=0;idim<4;idim++) {
				printf ( " %8.4f", part[ievent + idim*iNEvents + ipart*ndim*iNEvents ] );
			}
		}
		printf("\n");
	}

	//  Free allocations from arrays allocated withing this function
	if (amp)  free(amp);
	if (part) free(part);

}
#endif


void
GPUb1piAngAmp_exec(dim3 dimGrid, dim3 dimBlock, 
									 // GPU_AMP_PROTO,
                   GDouble* pfDevData, WCUComplex* pcDevAmp, int* piDevPerm, int iNParticles, int iNEvents,
                   int polBeam, GDouble polFrac,
                   int J_X, int Par_X, int L_X, int I_X, int epsilon_R, 
                   int Iz_b1, int Iz_pi,
                   GDouble u_rho_1, GDouble u_rho_3, 
                   GDouble u_omega_1, GDouble u_omega_3,
                   GDouble u_b1_0, GDouble u_b1_2, 
                   GDouble G0_omega, GDouble G0_b1, bool orthocheck)
{  
  int ievent, ievent1, idim, ipart, i, j, k;
  int nonZero = 0;
  int static cnt = 0;
	cnt++;
// printf("test: Call to GPUb1piAngAmp_exec: cnt %d\n", cnt);

  // Identify amplitudes which are zero
  Pretest_kernel<<< dimGrid, dimBlock >>>
    ( 
      // GPU_AMP_ARGS, 
      pfDevData, pcDevAmp, piDevPerm, iNParticles, iNEvents,
      polBeam, polFrac, 
      J_X, Par_X, L_X, I_X, epsilon_R, Iz_b1, Iz_pi,
      u_rho_1, u_rho_3, u_omega_1, u_omega_3, u_b1_0, u_b1_2, 
      G0_omega, G0_b1, orthocheck ); 
// printf("test: after call to Pretest_kernel()\n");


  //  Copy pcDevAmp from device to host  */
  GDouble *hostAmp = (GDouble *) malloc(2*iNEvents*sizeof(GDouble));
  cudaMemcpy (hostAmp, pcDevAmp, 2*iNEvents*sizeof(GDouble), cudaMemcpyDeviceToHost);

  //  Initialize all on-device amplitudes to zero
  Setzero_kernel<<< dimGrid, dimBlock >>>(pcDevAmp,iNEvents);
// printf("test: after call to Setzero_kernel()\n");
  

  //  Count number of nonZero amplitudes
  for (i=0;i<iNEvents;i++) {
    if (hostAmp[2*i]==1.0) nonZero++;
  }

  //  Allocate array to hold indices of nonZero amplitudes
  int *nonZeroIndices = (int *) malloc(nonZero * sizeof(int));
  j = 0;
  for (i=0;i<iNEvents;i++) {
    if (hostAmp[2*i]==1.0) nonZeroIndices[j++] = i;
  }

  //  Copy 4momenta from GPU to CPU - make part() big enough to hold the entire set of momenta
  GDouble *part = (GDouble *) malloc (iNEvents * 4 * iNParticles * sizeof(GDouble));
  cudaMemcpy (part, pfDevData, iNEvents * 4 * iNParticles * sizeof(GDouble), cudaMemcpyDeviceToHost);
// printf("test: after copy pfDevData to Device\n");

  //  Copy nonZero momenta in place to the start of the array part
  //   Make sure order of copying moves continuously from lower to higher indice.
  for (ipart=0;ipart<iNParticles;ipart++) {
    for (idim=0;idim<4;idim++) {
      for (ievent1=0;ievent1<nonZero;ievent1++) {
        ievent = nonZeroIndices[ievent1];
        //  Index of nonZero event in original particle array
        i = ievent  + idim * iNEvents + ipart * 4 * iNEvents;
        //  Index of nonZero event in new      particle array
        j = ievent1 + idim * nonZero  + ipart * 4 * nonZero;
        part[j] = part[i];
      }
    }
  }



  //  Copy new particles on CPU back to GPU, only need those momenta sets which were non-zero, not the size of the entire set.
	GDouble *part_dev;
  cudaMalloc(&part_dev,       nonZero * 4 * iNParticles * sizeof(GDouble) );
  cudaMemcpy( part_dev, part, nonZero * 4 * iNParticles * sizeof(GDouble), cudaMemcpyHostToDevice );
// printf("test: after copy Part to Device\n");
  

  //  Reset dimGrid and dimBlock for the value of nonZero
  int Nthreads = 32;
  dim3 dimBlock1(Nthreads);
  dim3 dimGrid1((nonZero-1)/Nthreads+1);
  

  //   Evaluate non-zero amplitudes
  // iNEvents = nonZero;
  GPUb1piAngAmp_kernel<<< dimGrid1, dimBlock1 >>>
    ( 
      cnt, 
      // GPU_AMP_ARGS, 
      // pfDevData, pcDevAmp, piDevPerm, iNParticles, nonZero,
      part_dev, pcDevAmp, piDevPerm, iNParticles, nonZero,
      polBeam, polFrac, 
      J_X, Par_X, L_X, I_X, epsilon_R, Iz_b1, Iz_pi,
      u_rho_1, u_rho_3, u_omega_1, u_omega_3, u_b1_0, u_b1_2, 
      G0_omega, G0_b1, orthocheck ); 
// printf("test: after call to GUPb1piAngAmp_kernel()\n");

  //  Read amplitudes from GPU to CPU
	GDouble *amp = (GDouble *) malloc (iNEvents * 2 * sizeof(GDouble));
	cudaMemcpy (amp, pcDevAmp, iNEvents * 2 * sizeof(GDouble), cudaMemcpyDeviceToHost);

// printf("test: after copy Amp to Host\n");


  //  Re-arrange location of amplitudes on GPU to match original distribution of vectors
	//  Progress through the index array backward.
	k = iNEvents;
	for (i=nonZero-1;i>=0;i--) {
		//  Zero those elements between this element and last.
		for (j=nonZeroIndices[i]+1;j<k;j++) {
			amp[2*j  ] = 0.0;
			amp[2*j+1] = 0.0;
		}
		k = nonZeroIndices[i];
		amp[2*k  ] = amp[2*i  ];
		amp[2*k+1] = amp[2*i+1];
	}
	//  Zero remaining elements
	for (j=0;j<nonZeroIndices[0];j++) {
			amp[2*j  ] = 0.0;
			amp[2*j+1] = 0.0;
	}

	//  Write values back to GPU so calling program will find them where they
	//  expect them.
	cudaMemcpy (pcDevAmp, amp, iNEvents * 2 * sizeof(GDouble), cudaMemcpyHostToDevice);
// printf("test: after copy Amp to Device\n");
	
	//  Free allocations
  if (part_dev) cudaFree(part_dev);
	if (amp)  free(amp);
	if (part) free(part);

// printf("test: after Free allocations\n");
//  Print Particle and Amplitude CUDA arrays
#ifdef DEBUG
printCudaArrays(pfDevData, pcDevAmp, piDevPerm, iNParticles, iNEvents, cnt);
#endif

}
