#include <ctime>
#include <stdlib.h>
#include <stdio.h>

#include <cassert>
#include <iostream>
#include <string>
#include <sstream>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/AmpParameter.h"
#include "b1piAngAmp.h"
#include "AMPTOOLS_AMPS/barrierFactor.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"
#include "AMPTOOLS_AMPS/breakupMomentum.h"

b1piAngAmp::b1piAngAmp( const vector< string >& args ):
  UserAmplitude< b1piAngAmp >( args ),
  m_fastCalc( false )
{

  const unsigned int base_arg_num=8;

  bool tweakBW_omega = args.size() == base_arg_num+1;
  bool tweakBW_omega_b1 = args.size() == base_arg_num+2;
  //accept either base number of arguments, extended set (14)
  // for orthogonality check diagnostics
  // or base number plus omega width or base number plus omega and b1 widths
  assert(args.size() == base_arg_num || args.size() == 14 ||
         tweakBW_omega || tweakBW_omega_b1);
  
  mpolBeam = atoi( args[0].c_str() );  // beam polarization component (X=0, Y=1)
  mpolFrac = AmpParameter( args[1] );  // fraction of polarization 0=0% 1=100%.
  mJ_X      = atoi( args[1].c_str() ); // total J of produced resonance
  // parity of produced resonance
  mPar_X    = ( atoi( args[2].c_str() ) == 1 ? 1 : -1 ); // for convenience let Par_X=0 --> -1 
  int L_X = atoi( args[3].c_str() ); // L between bachelor (pi) and isobar (b1)
  mL_X      = ( L_X >= 0 ? L_X : ( mJ_X == 0 ? 1 : 0) );
  mI_X      = atoi( args[4].c_str() );  // isospin of the resonance
  mepsilon_R= ( atoi( args[5].c_str() ) == 1 ? 1 : -1 );

  int Iz_b1    = atoi( args[6].c_str() );
  int Iz_pi    = atoi( args[7].c_str() );
  
  bool m_ORTHOCHECK = (args.size() == base_arg_num || tweakBW_omega || tweakBW_omega_b1);
  
  // Note, the following have no effect since L_\omega & J_\rho
  // have been restricted to value 1
  m_u_rho_1  = m_ORTHOCHECK ? sqrt(.9) : atoi( args[8].c_str());
  m_u_rho_3  = m_ORTHOCHECK ? sqrt(.1) : atoi( args[9].c_str());
  m_u_omega_1= m_ORTHOCHECK ? sqrt(.9) : atoi( args[10].c_str());
  m_u_omega_3= m_ORTHOCHECK ? sqrt(.1) : atoi( args[11].c_str());
  
  mG0_omega = 0.0085;
  mG0_b1 = 0.143;
  
  if(tweakBW_omega && args[8][0]=='F') m_fastCalc=true;
  else{
    if(tweakBW_omega || tweakBW_omega_b1) mG0_omega = atof( args[8].c_str());
    if(tweakBW_omega_b1) mG0_b1 = atof( args[9].c_str());
  }
  
  const GDouble b1DSratio2 = 0.277*0.277; //from PDG: D/S amp ratio=0.277+/-0.027
  m_u_b1_0   = m_ORTHOCHECK ? sqrt(1/(1 + b1DSratio2))
  : atoi( args[12].c_str());
  m_u_b1_2   = m_ORTHOCHECK ? sqrt(b1DSratio2/(1 + b1DSratio2))
  : atoi( args[13].c_str());

               
  assert( ( mpolBeam == 0 ) || ( mpolBeam == 1 ) );
  if(!(( mpolFrac >= 0 ) && ( mpolFrac <= 1 ))){
    cout << "ERROR: polFrac set to " << mpolFrac << endl << "Should be 0.0-1.0" << endl;
    assert(false);
  }

  assert( mJ_X >= 0 && mJ_X <=2 );
//  assert( abs( (double)mPar_X ) <= 1 );
  assert( abs( (double)mI_X )   <= 1 );
  //assert( mL_X <= mJ_X );
  assert( (mL_X+1 >= mJ_X && abs(mL_X-1) <= mJ_X) || mL_X==-1 );
  assert( abs(mepsilon_R)<=1 );
  assert( abs(Iz_b1) <= 1 );
  assert( abs(Iz_pi) <= 1 );

  registerParameter( mpolFrac );

  m_disableBW_omega = mG0_omega <= 0;
  m_disableBW_b1 = mG0_b1 <= 0;

  // create nominal Iz list: 0,0,-1,+1,0,-1,+1
  // (proton isospin irrelevant at the moment)
  mIz.assign(7, 0);
  mIz[2]=Iz_pi;
  mIz[3]=Iz_b1;
  mIz[5]=-1;
  mIz[6]=+1;
}

void PrintHEPvector(TLorentzVector &v){
  printf("(%6.3f, %6.3f, %6.3f; %6.3f m=%6.3f)\n",v.X(),v.Y(),v.Z(),v.T(),v.M());
}
void PrintArrVector(GDouble *v){
  printf("arr(%6.3f, %6.3f, %6.3f; %6.3f m=%6.3f)\n",v[1],v[2],v[3],v[0],
	 v[0]*v[0]-v[1]*v[1]-v[2]*v[2]-v[3]*v[3]);
}



inline GDouble b1piAngAmp::N(int J) const
{ return sqrt((2*J+1)/(4*M_PI)); }

TLorentzVector& MoveToRF(TLorentzVector &parent,
			  TLorentzVector &daughter)
{
  daughter.RotateZ(-parent.Phi());
  daughter.RotateY(-parent.Theta());
  daughter.Boost(0,0,-parent.Rho()/parent.E());
  return daughter;
}


inline complex <GDouble> b1piAngAmp::
BreitWigner(GDouble m0, GDouble Gamma0, int L,
	    TLorentzVector &P1, TLorentzVector &P2) const
{
  
  TLorentzVector Ptot=P1+P2;
  GDouble m  = Ptot.M();
  GDouble mass1 = P1.M();
  GDouble mass2 = P2.M();
  
  
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
  complex<GDouble> bwtop(sqrt(m0*width_coef) * num_qdep, 0.0 );
  
  complex<GDouble> bwbottom( ( m0*m0 - m*m ) ,
			     -1.0 * ( m0 * width ) );
  
  return( bwtop / bwbottom );
  
}


inline GDouble b1piAngAmp::CB(int j1, int j2, int m1, int m2, int J, int M) const
{
  if( j1*j2 == 0 ) return 1.0;
  
  if(J == 1){
    if(j1==1 && j2==1){
      if(M == -1){
	if(m1==0 && m2==-1)
	  return 1/sqrt(2.0);
	else if(m1==-1 && m2==0)
	  return -1/sqrt(2.0);
	else
	  return 0.0;
      }else if(M == 1){
	if(m1==0 && m2==1)
	  return -1.0/sqrt(2.0);
	else if(m1==1 && m2==0)
	  return 1.0/sqrt(2.0);
	else
	  return 0.0;
      }else if(M == 0){
	return m1/sqrt(2.0);
      }
    }else if( j1==2 && j2==1 ){
      if(M == -1 || M == 1){
	if(m1==0) return 1.0/sqrt(10.0);
      } else if(M==0 && m1==0)
	return -sqrt(2.0/5.0);
    }else if( j1==2 && j2==3 ){
      if( M == -1 || M == 1 ) {
	if (m1==0) return sqrt(6.0/35.0);
      }else if (M==0 && m1==0)
	return 3.0/sqrt(35.0);
    }else if( j1==3 && j2==3 ){
      if( M == -1 || M == 1 ) {
	if (m1==0) return -m2*0.5*sqrt(6.0/7.0);
      }else if (M==0 && m1==0)
	return 0.0;
    }
  }
   
  //printf("Resorting to clebschGordan(%3d,%3d,%3d,%3d,%3d,%3d)\n",
  //j1, j2, m1, m2, J, M);
  return clebschGordan(j1, j2, m1, m2, J, M);
}



GDouble b1piAngAmp::u_rho(int J_rho) const 
{
  return J_rho==1 ? m_u_rho_1 : (J_rho==3 ? m_u_rho_3 : 0);
}
GDouble b1piAngAmp::u_omega(int L_omega) const 
{
  return L_omega==1 ? m_u_omega_1 : (L_omega==3 ? m_u_omega_3 : 0);
}
GDouble b1piAngAmp::u_b1(int L_b1) const 
{
  return L_b1==0 ? m_u_b1_0 : (L_b1==2 ? m_u_b1_2 : 0);
}
/*GDouble b1piAngAmp::v(int epsilon_R) const 
{
assert( abs(epsilon_R)==1 );
return epsilon_R==-1 ? m_v_m : (epsilon_R==+1 ? m_v_p);
}*/



complex< GDouble >
b1piAngAmp::calcAmplitude( GDouble** pKin ) const
{
  int m_X,IMLnum=0;
  bool useCutoff=true;
  complex <GDouble> i(0, 1), COne(1, 0),CZero(0,0);  

  const vector< int >& perm = getCurrentPermutation();  
  int Iz_b1 = mIz[perm[2]];
  int Iz_pi = mIz[perm[3]];

  if(abs(Iz_b1+Iz_pi) > mI_X) return CZero;

  TLorentzVector beam  (pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0]); 
  TLorentzVector recoil(pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0]);
  GDouble InvSqrt2=1/sqrt(2.0);
  GDouble m0_rho=0.775,G0_rho=0.149;
  GDouble m0_omega=0.783, m0_b1=1.223;


  //Exprected particle list: 
  // pi- b1(pi+ omega(pi0 "rho"(pi- pi+)))
  //  2      3         4         5   6
                        

  TLorentzVector rhos_pip(pKin[6][1], pKin[6][2], pKin[6][3], pKin[6][0]);
  TLorentzVector rhos_pim(pKin[5][1], pKin[5][2], pKin[5][3], pKin[5][0]);
  TLorentzVector rho = rhos_pip + rhos_pim;

  if( useCutoff && rho.M()+0.135 > m0_omega+3*mG0_omega){
    //cout << "s";
    return CZero;
  }

  TLorentzVector omegas_pi(pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0]);
  TLorentzVector omega = rho + omegas_pi;

  if(useCutoff && fabs(omega.M()-m0_omega) > 3*mG0_omega){
    //cout << "s";
    return CZero; 
  }

  TLorentzVector b1s_pi(pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);
  TLorentzVector b1 = omega + b1s_pi;

  if( useCutoff && (fabs(b1.M()-m0_b1) > 3*mG0_b1 ||
		    b1.M() < (m0_omega - 3*mG0_omega)) ){
    //cout << "s";
    return CZero;
  }

  //printf("DEBUG: proceeding with b1pi amp calc. mG0_omega=%f, disbale BWomega=%d\n",mG0_omega,m_disableBW_omega);

  TLorentzVector Xs_pi(pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0]);
  TLorentzVector X = b1 + Xs_pi;

  GDouble q = breakupMomentum( X.M(), b1.M(), Xs_pi.M() );


  // orientation of production plane in lab
  GDouble alpha = recoil.Vect().Phi();
  
  //Resonance RF, Godfried-Jackson frame
  TLorentzRotation XRFboost( -X.BoostVector() );

  TLorentzVector beam_XRF   = XRFboost * beam;
  TLorentzVector recoil_XRF = XRFboost * recoil;
  
  //Define coordinate system
  TVector3 zGJ = beam_XRF.Vect().Unit();
  TVector3 yGJ = zGJ.Cross(recoil_XRF.Vect()).Unit();
  TVector3 xGJ = yGJ.Cross(zGJ);

  
  TLorentzVector b1_XRF      = XRFboost * b1;
  TLorentzVector omega_XRF   = XRFboost * omega;
  TLorentzVector rho_XRF     = XRFboost * rho;
  TLorentzVector rhos_pip_XRF= XRFboost * rhos_pip;

  TLorentzVector omega_b1RF(MoveToRF(b1_XRF, omega_XRF));
  TLorentzVector rho_omegaRF(MoveToRF(omega_b1RF, 
					MoveToRF(b1_XRF, rho_XRF)));
  TLorentzVector rhos_pip_rhoRF(MoveToRF(rho_omegaRF,
					    MoveToRF(omega_b1RF, 
						     MoveToRF(b1_XRF,rhos_pip_XRF))));

  TVector3 ang_b1( (b1_XRF.Vect()).Dot(xGJ),
		     (b1_XRF.Vect()).Dot(yGJ),
		     (b1_XRF.Vect()).Dot(zGJ) );
  
  //printf("%f %f    %f %f \n",
  // omega.M(), (rho+omegas_pi).M(),
  //	 rho.M(), (rhos_pip+rhos_pim).M());


  // SUMMATION GUIDE:
  // notation meant to resemble TeX symbols in derivation
  // exception: pol = \epsilon_\gamma
  // l -> lambda, indicating helicity
  // u_[particle](q.n.) -> amplitude strength coefficient 

  int pol=(mpolBeam==1 ? +1 : -1); // y and x-pol. respectively
  const int* epsilon_R=&mepsilon_R;
  GDouble rho_omegaRF_cosTheta=rho_omegaRF.CosTheta();
  GDouble rho_omegaRF_phi     =rho_omegaRF.Phi();
  GDouble rhos_pip_rhoRF_cosTheta=rhos_pip_rhoRF.CosTheta();
  GDouble rhos_pip_rhoRF_phi     =rhos_pip_rhoRF.Phi();



  // Prepare sets of quantum numbers over which to sum ----------

  // Full L_omega=1, J_rho=1 and L_omega=3, J_rho=3
  /*int aList_J_rho[]={1,3}; 
  vector<int> List_J_rho(aList_J_rho, aList_J_rho+2);
  int aList_L_omega[]={1,3};  
  vector<int> List_L_omega(aList_L_omega, aList_L_omega+2);
  */

  //Restricted  L_omega=1, J_rho=1 combination only
  int aList_J_rho[]={1}; 
  vector<int> List_J_rho(aList_J_rho, aList_J_rho+1);
  int aList_L_omega[]={1};  
  vector<int> List_L_omega(aList_L_omega, aList_L_omega+1);

  int aList_l_rho[]={-1,1}; //'0' term vanishes for L_omega=J_rho=1
  // because CB(1,1,0,0;1,0)=0 ...change if introducing omega F-wave
  vector<int> List_l_rho(aList_l_rho, aList_l_rho+2);

  int aList_l_omega[]={-1,0,1};   
  vector<int> List_l_omega(aList_l_omega, aList_l_omega+3);

  int aList_L_b1[]={0,2};  
  vector<int> List_L_b1(aList_L_b1,aList_L_b1+2);

  //shortcut:  CB(L_X, J_b1, 0, l_b1 ; J_X, l_b1) vanishes when
  //  = CB(1, 1, 0, 0 ; 1, 0),  so omit l_b1=0 when J_X=L_X=1
  int aList_l_b1[3]={-1};
  if(mL_X==1 && mJ_X==1)  aList_l_b1[1]=1;
  else{
    aList_l_b1[1]=0;
    aList_l_b1[2]=1; }
  vector<int> List_l_b1(aList_l_b1,aList_l_b1+
			(mL_X==1 && mJ_X==1 ? 2 : 3));
  // --------------------

  vector<int> List_l_R;
  for(int n=0; n <= mJ_X+1 ; n++) {
    List_l_R.push_back(n);
    //printf("List_l_R.push_back(%d)\n",n);
  }
  
  //int aList_epsilon_R[]={-1,1}; 
  //vector<int> List_epsilon_R(aList_epsilon_R,aList_epsilon_R+2);

  // End of q.n. set prep ----------------------------------------


  complex <GDouble> ThelSum(0,0);

  if(mJ_X==0) if(mPar_X*pol*(*epsilon_R) ==  -1 ) return CZero;

  complex <GDouble> expFact(cos(alpha), sin(alpha));
  complex <GDouble> expFact_conj(conj(expFact));
  //summing positive and negative helicity terms
  for(int l_gamma=-1; l_gamma <= +1 ; l_gamma+=2){
    
    //for(vector<int>::iterator epsilon_R=List_epsilon_R.begin();
	//epsilon_R != List_epsilon_R.end() ; epsilon_R++){
      
      for(vector<int>::iterator l_R=List_l_R.begin(); 
	  l_R != List_l_R.end() ; l_R++){
	if(*l_R==0 && *epsilon_R==-1) continue;
	
	//summing positive and negative helicity terms of R's reflectivity state
	for(int l_Rsign = +1; l_Rsign >= (*l_R>0 ? -1:+1) ; l_Rsign-=2){
	  m_X=l_gamma - l_Rsign * (*l_R);
	  if(m_X==0){
	    //testing for cancelation in |J 0>+pol*P*epsilon_R*(-1)^J|J 0>
	    if(mPar_X*pol*(*epsilon_R) == (mJ_X % 2 ==0 ? -1:+1)) continue;
	  }else
	    //enforcing that the selected projection <= vector magnitude 
	    if( abs(m_X)>mJ_X) continue; 
	  
	  
	  complex <GDouble> l_b1DepTerm(0,0);
	  for(vector<int>::iterator l_b1=List_l_b1.begin(); 
	      l_b1 != List_l_b1.end() ; l_b1++){
	    
	    complex <GDouble> L_b1DepTerm(0,0);
	    for(vector<int>::iterator L_b1=List_L_b1.begin(); 
		L_b1 != List_L_b1.end() ; L_b1++){
	      
	      complex <GDouble> l_omegaDepTerm(0,0);
	      for(vector<int>::iterator l_omega=List_l_omega.begin(); 
		  l_omega != List_l_omega.end() ; l_omega++){

		complex <GDouble> L_omegaDepTerm(0,0);
		//only odd L_omega allowed to given off J_rho to honor P_omega=-1
		for(vector<int>::iterator L_omega=List_L_omega.begin(); 
		    L_omega != List_L_omega.end() ; L_omega++){

		  complex <GDouble> J_rhoDepTerm(0,0);
		  for(vector<int>::iterator J_rho=List_J_rho.begin(); 
		      J_rho != List_J_rho.end() ; J_rho++){

		    //enforces triang. ineq. betw. J_omega=1, J_rho and L_omega
		    if( abs(*J_rho-*L_omega) > 1) continue; 
		    
		    complex <GDouble> l_rhoDepTerm(0,0);
		    for(vector<int>::iterator l_rho = List_l_rho.begin(); 
			l_rho != List_l_rho.end() ; l_rho++){
		      //shortcut CB(1,1,0,0;1,0)=0
		      if(*L_omega==1 && *J_rho==1 && *l_rho==0) continue;
		      l_rhoDepTerm+= conj(wignerD(1, *l_omega, *l_rho,
						  rho_omegaRF_cosTheta, 
						  rho_omegaRF_phi))*
			CB(*L_omega, *J_rho, 0, *l_rho, 1, *l_rho) *
			Y(*J_rho, *l_rho, rhos_pip_rhoRF_cosTheta, rhos_pip_rhoRF_phi);
		      
		      IMLnum++;
		    }

		    J_rhoDepTerm += u_rho(*J_rho) * l_rhoDepTerm *
		      BreitWigner(m0_rho,G0_rho, *J_rho,rhos_pip,rhos_pim);
		  }
		  
		  if(!m_disableBW_omega) J_rhoDepTerm*=
		    BreitWigner(m0_omega,mG0_omega, *L_omega, omegas_pi,rho);
		  
		  L_omegaDepTerm += u_omega(*L_omega)*J_rhoDepTerm*N(*L_omega);
		}
		
		l_omegaDepTerm += 
		  L_omegaDepTerm *
		  conj(wignerD(1, *l_b1, *l_omega, omega_b1RF.CosTheta(), 
			       omega_b1RF.Phi())) *
		  CB(*L_b1, 1, 0, *l_omega, 1, *l_omega);
	      }
	      
	      if(!m_disableBW_b1) l_omegaDepTerm*=
		BreitWigner(m0_b1, mG0_b1, *L_b1, b1s_pi, omega);
	      
	      L_b1DepTerm += u_b1(*L_b1)*l_omegaDepTerm * N(*L_b1);
	    }
	    
	    l_b1DepTerm += 
	      L_b1DepTerm * CB(mL_X, 1, 0, *l_b1, mJ_X, *l_b1)*
	      conj(wignerD(mJ_X, m_X, *l_b1, ang_b1.CosTheta(), ang_b1.Phi()));
	    
	    
	  }
	  
	  ThelSum += 
	    //Assemble reflectivity eigenvector with epsilon_X=pol*epslion_R
	    l_b1DepTerm*(GDouble)
	    (m_X<0 ? mPar_X*pol*(*epsilon_R)*((mJ_X-m_X) % 2 == 0 ? +1:-1) : 1) * 
	    (GDouble)(m_X == 0 ? 1.0 : InvSqrt2 ) *
	    //to account for |eps_g> ~ (|1,-1>exp(-ia)-pol|1,+1>exp(ia)) 
	    (l_gamma==1 ? (GDouble)(-pol)*expFact : expFact_conj)*
	    // to apply th(l_R) reflectivity state prefactor: 
	    // m=0: 1/2  m>0: 1/sqrt(2)  m<0: 0 (last just skipped in this sum)  
	    (GDouble)(*l_R > 0 ? InvSqrt2 : 1.0 ) *
	    //apply coefficients to the reflectivity basis terms:
	    (GDouble)(l_Rsign==1 ? 1 : *epsilon_R); //v(*epsilon_R) *

	}
      }
      //}
  }
  

  /*printf("DEBUG: perm2=%d -> Iz_b1=%d\tperm3=%d -> Iz_pi=%d\t==>\t%f\t%d %d %d %d %d\n",
	 perm[2],Iz_b1,perm[3],Iz_pi, CB(1, 1, Iz_b1, Iz_pi, mI_X, Iz_b1 + Iz_pi),
	 mIz[2],mIz[3],mIz[4],mIz[5],mIz[6]);
  */

  ThelSum *= N(mL_X) * (GDouble)(mL_X==0 ? 1.0 : (mL_X==1 ? q : pow(q,mL_X))) *
    // to apply polarization fraction weights: 
    (GDouble)sqrt((1.0-pol*mpolFrac)*0.5) * //(1+g) for x-pol, (1-g) for y-pol   
    (pol==1 ? i : COne)*InvSqrt2 * //to account for |eps_g> ~ sqrt(-eps/2)
    CB(1, 1, Iz_b1, Iz_pi, mI_X, Iz_b1 + Iz_pi);


  if(m_ORTHOCHECK) {
    double I=abs(ThelSum);
    printf("ORTHOCHECK %3.1f %3.1f  %3.1f %3.1f  %3.1f %3.1f\t%17.10e ",
	   m_u_rho_1, m_u_rho_3, 
	   m_u_omega_1, m_u_omega_3, m_u_b1_0, m_u_b1_2, I*I);
    printf("%d %e %d %d %d %d %d\n",mpolBeam,(double)mpolFrac,
	   mJ_X,mPar_X,mL_X,mI_X,mepsilon_R);

  }

  return ThelSum;

}

#ifdef GPU_ACCELERATION

void b1piAngAmp::
launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

  const vector< int >& perm = getCurrentPermutation();  
  int Iz_b1 = mIz[perm[2]];
  int Iz_pi = mIz[perm[3]];

  GPUb1piAngAmp_exec(dimGrid, dimBlock, GPU_AMP_ARGS, mpolBeam, mpolFrac,
		     mJ_X, mPar_X, mL_X, mI_X, mepsilon_R, Iz_b1, Iz_pi,
		     m_u_rho_1, m_u_rho_3, m_u_omega_1, m_u_omega_3,
		     m_u_b1_0, m_u_b1_2, mG0_omega, mG0_b1, 
		     /*m_ORTHOCHECK*/ false);
  
}

#endif
