
#include <stdlib.h>

#include <cassert>
#include <iostream>
#include <string>
#include <sstream>

#include "IUAmpTools/AmpParameter.h"
#include "b1piAngAmp.h"
#include "AMPTOOLS_AMPS/barrierFactor.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"
#include "AMPTOOLS_AMPS/breakupMomentum.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/LorentzRotation.h"
#include "CLHEP/Vector/ThreeVector.h"

b1piAngAmp::b1piAngAmp(int polBeam, float polFrac,//const AmpParameter& polFrac,
		       int J_X, int Par_X, int L_X, int I_X, int epsilon_R,
		       int Iz_b1, int Iz_pi,
		       float u_rho_1, float u_rho_3,
		       float u_omega_1, float u_omega_3,
		       float u_b1_0, float u_b1_2, bool orthocheck=false):
  
  Amplitude(),
  mpolBeam( polBeam ),  // beam polarization component (X=0, Y=1)
  mpolFrac( polFrac ),  // fraction of polarization 0=0% 1=100%. 
  mJ_X( J_X ),     // total J of produced resonance
  // parity of produced resonance
  mPar_X( Par_X==1 ? 1 : -1 ), // for convenience let Par_X=0 --> -1 
  // L between bachelor (pi) and isobar (b1)
  mL_X( L_X >= 0 ? L_X : (J_X==0 ? 1 : 0) ),
  mI_X( I_X ),   // isospin of the resonance
  mepsilon_R(epsilon_R==1? 1 : -1), 
  mIz_b1( Iz_b1 ),
  mIz_pi( Iz_pi ),

  m_u_rho_1(u_rho_1),
  m_u_rho_3(u_rho_3),
  m_u_omega_1(u_omega_1),
  m_u_omega_3(u_omega_3),
  m_u_b1_0(u_b1_0),
  m_u_b1_2(u_b1_2),
  m_ORTHOCHECK(orthocheck)

{
  
  assert( ( polBeam == 0 ) || ( polBeam == 1 ) );
  //assert( ( polFrac >= 0 ) && ( polFrac <= 1 ) );
  if(!(( polFrac >= 0 ) && ( polFrac <= 1 ))){
    cout << "ERROR: polFrac set to " << polFrac << endl << "Should be 0.0-1.0" << endl;
    assert(false);
  }

  assert( J_X >= 0 && J_X <=2 );
  assert( abs( (double)Par_X ) <= 1 );
  assert( abs( (double)I_X )   <= 1 );
  //assert( L_X <= J_X );
  assert( (L_X+1 >= J_X && abs(L_X-1) <= J_X) || L_X==-1 );
  assert( abs(epsilon_R)<=1 );
  assert( abs(Iz_b1) <= 1 );
  assert( abs(Iz_pi) <= 1 );

  //registerParameter( mpolFrac );
  
  setDefaultStatus( false );
}

void PrintHEPvector(HepLorentzVector &v){
  printf("(%6.3f, %6.3f, %6.3f; %6.3f m=%6.3f)\n",v.x(),v.y(),v.z(),v.t(),v.m());
}
void PrintArrVector(GDouble *v){
  printf("arr(%6.3f, %6.3f, %6.3f; %6.3f m=%6.3f)\n",v[1],v[2],v[3],v[0],
	 v[0]*v[0]-v[1]*v[1]-v[2]*v[2]-v[3]*v[3]);
}



inline GDouble N(int J){ return sqrt((2*J+1)/(4*M_PI)); }

HepLorentzVector& MoveToRF(HepLorentzVector &parent,
			  HepLorentzVector &daughter)
{
  daughter.rotateZ(-parent.phi());
  daughter.rotateY(-parent.theta());
  daughter.boost(0,0,-parent.rho()/parent.e());    
  return daughter;
}


inline complex <GDouble> BreitWigner(GDouble m0, GDouble Gamma0, int L,
		    HepLorentzVector &P1,HepLorentzVector &P2)
{
  //printf("L passed to BW: %d\n",L);

   HepLorentzVector Ptot=P1+P2;
   GDouble m  = Ptot.m();
   GDouble mass1 = P1.m();
   GDouble mass2 = P2.m();
   

   //printf("m0=%f m=%f\n",m0,m);

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
   //complex<GDouble> bwtop(sqrt(m0*width_coef), 0.0 );
   
   complex<GDouble> bwbottom( ( m0*m0 - m*m ) ,
			      -1.0 * ( m0 * width ) );
  
   return( bwtop / bwbottom );
   //return 1.0;
}


inline float CB(int j1, int j2, int m1, int m2, int J, int M)
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

/*
void ReduceToOne(vector<int> &v,unsigned int n){
  int selval=v[n];

  v.clear();
  v.push_back(selval);

}
*/


float b1piAngAmp::u_rho(int J_rho) const 
{
  return J_rho==1 ? m_u_rho_1 : (J_rho==3 ? m_u_rho_3 : 0);
}
float b1piAngAmp::u_omega(int L_omega) const 
{
  return L_omega==1 ? m_u_omega_1 : (L_omega==3 ? m_u_omega_3 : 0);
}
float b1piAngAmp::u_b1(int L_b1) const 
{
  return L_b1==0 ? m_u_b1_0 : (L_b1==2 ? m_u_b1_2 : 0);
}
/*float b1piAngAmp::v(int epsilon_R) const 
{
assert( abs(epsilon_R)==1 );
return epsilon_R==-1 ? m_v_m : (epsilon_R==+1 ? m_v_p);
}*/



complex< GDouble >
b1piAngAmp::calcAmplitude( GDouble** pKin ) const
{
  
  HepLorentzVector beam  (pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0]); 
  HepLorentzVector recoil(pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0]);
  complex <GDouble> i(0, 1), COne(1, 0),CZero(0,0);  
  GDouble InvSqrt2=1/sqrt(2.0);
  GDouble m0_rho=0.775,G0_rho=0.149;
  GDouble m0_omega=0.783,G0_omega=0.0085;
  //GDouble m0_omega=1.083,G0_omega=0.0085;
  GDouble m0_b1=1.223,G0_b1=0.143;

  //Exprected particle list: 
  // pi- b1(pi+ omega(pi0 "rho"(pi- pi+)))
  //  2      3         4         5   6
                        
  //Acquire the 4-vectors for the current particle permutation
  const vector< int >& perm = getCurrentPermutation();

  //printf("PERM: %d %d %d %d %d\n",perm[6],perm[5],perm[4],perm[3],perm[2]);

  HepLorentzVector rhos_pip(pKin[perm[6]][1], pKin[perm[6]][2],
			    pKin[perm[6]][3], pKin[perm[6]][0]);
  HepLorentzVector rhos_pim(pKin[perm[5]][1], pKin[perm[5]][2],
			    pKin[perm[5]][3], pKin[perm[5]][0]);
  HepLorentzVector rho = rhos_pip + rhos_pim;

  if( rho.m()+0.135 > m0_omega+3*G0_omega) return CZero;

  HepLorentzVector omegas_pi(pKin[perm[4]][1], pKin[perm[4]][2],
			     pKin[perm[4]][3], pKin[perm[4]][0]);
  HepLorentzVector omega = rho + omegas_pi;

  if(fabs(omega.m()-m0_omega) > 3*G0_omega) return CZero; 

  HepLorentzVector b1s_pi(pKin[perm[3]][1], pKin[perm[3]][2],
			  pKin[perm[3]][3], pKin[perm[3]][0]);
  HepLorentzVector b1 = omega + b1s_pi;

  if( fabs(b1.m()-m0_b1) > 3*G0_b1 ||
      b1.m() < (m0_omega - 3*G0_omega) ) return CZero;

  HepLorentzVector Xs_pi(pKin[perm[2]][1], pKin[perm[2]][2],
			  pKin[perm[2]][3], pKin[perm[2]][0]);
  HepLorentzVector X = b1 + Xs_pi;

  GDouble q = breakupMomentum( X.m(), b1.m(), Xs_pi.m() );


  // orientation of production plane in lab
  GDouble alpha = recoil.vect().phi();
  
  //Resonance RF, Godfried-Jackson frame
  HepLorentzRotation XRFboost( -X.boostVector() );

  HepLorentzVector beam_XRF   = XRFboost * beam;
  HepLorentzVector recoil_XRF = XRFboost * recoil;
  
  //Define coordinate system
  Hep3Vector zGJ = beam_XRF.vect().unit();
  Hep3Vector yGJ = zGJ.cross(recoil_XRF.vect()).unit();
  Hep3Vector xGJ = yGJ.cross(zGJ);

  
  HepLorentzVector b1_XRF      = XRFboost * b1;
  HepLorentzVector omega_XRF   = XRFboost * omega;
  HepLorentzVector rho_XRF     = XRFboost * rho;
  HepLorentzVector rhos_pip_XRF= XRFboost * rhos_pip;

  HepLorentzVector omega_b1RF(MoveToRF(b1_XRF, omega_XRF));
  HepLorentzVector rho_omegaRF(MoveToRF(omega_b1RF, 
					MoveToRF(b1_XRF, rho_XRF)));
  HepLorentzVector rhos_pip_rhoRF(MoveToRF(rho_omegaRF,
					    MoveToRF(omega_b1RF, 
						     MoveToRF(b1_XRF,rhos_pip_XRF))));

  Hep3Vector ang_b1( (b1_XRF.vect()).dot(xGJ),
		     (b1_XRF.vect()).dot(yGJ),
		     (b1_XRF.vect()).dot(zGJ) );
  
  //printf("%f %f    %f %f \n",
  // omega.m(), (rho+omegas_pi).m(),
  //	 rho.m(), (rhos_pip+rhos_pim).m());

  //going from x[y]=0[1] notation to reflectivity: -1[+1]
  int m_X,IMLnum=0;


  // SUMMATION GUIDE:
  // notation meant to resemble TeX symbols in derivation
  // exception: pol = \epsilon_\gamma
  // l -> lambda, indicating helicity
  // u_[particle](q.n.) -> amplitude strength coefficient 

  int pol=(mpolBeam==1 ? +1 : -1); // y and x-pol. respectively
  const int* epsilon_R=&mepsilon_R;
  GDouble rho_omegaRF_cosTheta=rho_omegaRF.cosTheta();
  GDouble rho_omegaRF_phi     =rho_omegaRF.phi();
  GDouble rhos_pip_rhoRF_cosTheta=rhos_pip_rhoRF.cosTheta();
  GDouble rhos_pip_rhoRF_phi     =rhos_pip_rhoRF.phi();



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

  if(mJ_X==0) if(mPar_X*pol*(*epsilon_R) ==  -1 ) return ThelSum;
  
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

		    J_rhoDepTerm += u_rho(*J_rho) *
		      l_rhoDepTerm*BreitWigner(m0_rho,G0_rho, *J_rho,
					       rhos_pip,rhos_pim);
		  }
		  
		  L_omegaDepTerm += u_omega(*L_omega)*
		    J_rhoDepTerm*N(*L_omega)*
		    BreitWigner(m0_omega,G0_omega, *L_omega,
				omegas_pi,rho);
		}
		
		l_omegaDepTerm += 
		  L_omegaDepTerm *
		  conj(wignerD(1, *l_b1, *l_omega, omega_b1RF.cosTheta(), 
			       omega_b1RF.phi())) *
		  CB(*L_b1, 1, 0, *l_omega, 1, *l_omega);
	      }
	      
	      L_b1DepTerm += u_b1(*L_b1)*
		l_omegaDepTerm * N(*L_b1)*
		BreitWigner(m0_b1,G0_b1,*L_b1,
			    b1s_pi,omega);
	    }
	    
	    l_b1DepTerm += 
	      L_b1DepTerm * CB(mL_X, 1, 0, *l_b1, mJ_X, *l_b1)*
	      conj(wignerD(mJ_X, m_X, *l_b1, ang_b1.cosTheta(), ang_b1.phi()));
	    
	    
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
  
  
  
  ThelSum *= N(mL_X) * (GDouble)(mL_X==0 ? 1.0 : (mL_X==1 ? q : pow(q,mL_X))) *
    // to apply polarization fraction weights: 
    (GDouble)sqrt((1.0-pol*mpolFrac)*0.5) * //(1+g) for x-pol, (1-g) for y-pol   
    (pol==1 ? i : COne)*InvSqrt2 * //to account for |eps_g> ~ sqrt(-eps/2)
    CB(1, 1, mIz_b1, mIz_pi, mI_X, mIz_b1 + mIz_pi);


  if(m_ORTHOCHECK) {
    double I=abs(ThelSum);
    printf("ORTHOCHECK %3.1f %3.1f  %3.1f %3.1f  %3.1f %3.1f\t%13.10f ", m_u_rho_1, m_u_rho_3, 
	   m_u_omega_1, m_u_omega_3, m_u_b1_0, m_u_b1_2, I*I);
    printf("%d %1.0f %d %d %d %d %d\n",mpolBeam,(double)mpolFrac,mJ_X,mPar_X,mL_X,mI_X,mepsilon_R);

  }

  return ThelSum;

}
    
    
b1piAngAmp*
b1piAngAmp::newAmplitude( const vector< string >& args ) const {
  
  const unsigned int base_arg_num=9;
  assert( args.size() == base_arg_num || args.size() == 15 );
  
  int polBeam = atoi( args[0].c_str() );
  float  polFrac = atof(args[1].c_str());
  //AmpParameter polFrac( args[1] );
  int J_X      = atoi( args[2].c_str() );
  int Par_X    = atoi( args[3].c_str() );
  int L_X      = atoi( args[4].c_str() );
  int I_X      = atoi( args[5].c_str() );
  int epsilon_R= atoi( args[6].c_str() );
  int Iz_b1    = atoi( args[7].c_str() );
  int Iz_pi    = atoi( args[8].c_str() );
  
  bool use_emp = (args.size() == base_arg_num);
  
  float u_rho_1  = use_emp ? sqrt(.9) : atoi( args[9].c_str());
  float u_rho_3  = use_emp ? sqrt(.1) : atoi( args[10].c_str());
  float u_omega_1= use_emp ? sqrt(.9) : atoi( args[11].c_str());
  float u_omega_3= use_emp ? sqrt(.1) : atoi( args[12].c_str());


  const float b1DSratio2 = 0.277*0.277; //from PDG: D/S amp ratio=0.277+/-0.027
  float u_b1_0   = use_emp ? sqrt(1/(1 + b1DSratio2)) 
    : atoi( args[13].c_str());
  float u_b1_2   = use_emp ? sqrt(b1DSratio2/(1 + b1DSratio2)) 
    : atoi( args[14].c_str());

  
  return new b1piAngAmp( polBeam, polFrac, J_X, Par_X, L_X, I_X, epsilon_R,
			 Iz_b1, Iz_pi,
			 u_rho_1, u_rho_3, u_omega_1, u_omega_3, 
			 u_b1_0, u_b1_2, !use_emp);    
  
  
}

b1piAngAmp*
b1piAngAmp::clone() const {
  
  return ( isDefault() ? new b1piAngAmp() : 
	   new b1piAngAmp(mpolBeam,mpolFrac, mJ_X,mPar_X, mL_X, mI_X,mepsilon_R,
			  mIz_b1, mIz_pi,
			  m_u_rho_1, m_u_rho_3, 
			  m_u_omega_1, m_u_omega_3, 
			  m_u_b1_0, m_u_b1_2, m_ORTHOCHECK));
}
