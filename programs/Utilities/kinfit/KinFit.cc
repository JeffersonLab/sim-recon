#include <cctype>
#include "KinFit.h"
using namespace std;
/** @file KinFit.C
 * @brief KinFit class source file.
 */
//_____________________________________________________________________________
/** @class KinFit
 * @brief Utility class used to perform kinematic fitting. 
 * 
 * The interface with 
 * this class is a little primative since it was meant to be called by 
 * ClasEvent. However, it can be used as a stand alone class. 
 *
 * For each event, the KinFit object must be given the tagged photon energy,
 * detected charged particle 4-momenta, target mass and the covariance matrix.
 * All of these can be set using KinFit::SetEvent. The covariance matrix for
 * all types of fits should have the following format:
 *
 * \f$ C = \left(\begin{array}{cccccccc} \sigma^2_{E_{\gamma}} & 0 & 0 & 0 & \ldots & 0 & 0 & 0 \\ 0 & C^{pp}_1 & C^{p\lambda}_1 & C^{p\phi}_1 & \ldots & 0 & 0 & 0 \\ 0 & C^{p\lambda}_1 & C^{\lambda \lambda}_1 & C^{\lambda \phi}_1 & \ldots & 0 & 0 & 0 \\ 0 & C^{p \phi}_1 & C^{\lambda \phi}_1 & C^{\phi \phi}_1 & \ldots & 0 & 0 & 0 \\ \vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\ 0 & 0 & 0 & 0 & \ldots & C^{pp}_k & C^{p\lambda}_k & C^{p\phi}_k \\ 0 & 0 & 0 & 0 & \ldots & C^{p\lambda}_k & C^{\lambda \lambda}_k & C^{\lambda \phi}_k \\ 0 & 0 & 0 & 0 & \ldots &  C^{p \phi}_k & C^{\lambda \phi}_k & C^{\phi \phi}_k \end{array}\right) \f$
 *
 * where \f$ \sigma^2_{E_{\gamma}} \f$ is the tagged photon energy resolution,
 * and \f$ C^{ab}_i \f$ are the <em>i</em>'th particle's covaraince matrix 
 * elements (where the particles are in the same order as they were passed in
 * to set the 4-momentum).
 *
 * Once these quantities are set, then any of the different types of kinematic
 * fits can be run (see Fit and Fit2MissingTagged).
 *
 * The results of the fit can be obtained thru various KinFit member functions.
 * All @a getter member functions return results of the last fit run.
 *
 * <b> Example Usage: </b>
 *
 * \include KinFit.ex
 */
//_____________________________________________________________________________

KinFit::KinFit(){
  /// Default Constructor
  _targetMass = 0.938272;
  _ndf = 0;
  _chi2 = 666.;
  _ephot_in = 0.;
  _ephot_out = 0.;
  _missingMass = -1.;
  _missingParticle = false;
  _extraC = false;
  _invariantMass = -1.;
  _extraC_miss = false;
  _trackingConversion = "EASY";
}
//_____________________________________________________________________________

void KinFit::_Copy(const KinFit &__kfit){
  /// Copy @a kfit data members.
  _pulls = __kfit._pulls;
  _chi2 = __kfit._chi2;
  _ndf = __kfit._ndf;
  _ephot_in = __kfit._ephot_in;
  _ephot_out = __kfit._ephot_out;
  _p4in = __kfit._p4in;
  _p4out = __kfit._p4out;
  _reconstruction = __kfit._reconstruction;
  _targetMass = __kfit._targetMass;
  _cov.ResizeTo(__kfit._cov);
  _cov = __kfit._cov;
  _missingMass = __kfit._missingMass;
  _missingParticle = __kfit._missingParticle;
  _extraC = __kfit._extraC;
  _invariantMass = __kfit._invariantMass;
  _extraC_meas = __kfit._extraC_meas;
  _extraC_miss = __kfit._extraC_miss;
  _trackingConversion = __kfit._trackingConversion;
  for(int i = 0; i < 3; i++) _sigma_missing[i] = __kfit._sigma_missing[i];
}
//_____________________________________________________________________________

void KinFit::_MainFitter(){
  /// Private function used internally to perform kinematic fitting.
  const int numParts = (int)_p4in.size();
  const int dim = 3*numParts + 1; // number of measured quantites
  // check to see that covariance matrix size is consistent with # particles
  if(_cov.GetNrows() != dim || _cov.GetNcols() != dim){
    // wrong size
    std::cout << "Error! <KinFit::_MainFitter> Covariance matrix size is NOT "
      << "correct for current number of particles. For " << numParts
      << " particles, the covariance matrix should be " << dim << " x "
      << dim << ", but the matrix passed in is " << _cov.GetNrows()
      << " x " << _cov.GetNcols() << std::endl;
    abort();
  }

  // Get the number of constraint equations
  int numConstraints = 4;
  if(_extraC) numConstraints = 5;

  // get the degrees of freedom
  _ndf = 4;
  if(_missingParticle) _ndf -= 3;
  if(_extraC) _ndf += 1;

  int i,j, k;
  double alpha[numParts],mass[numParts],lambda[numParts],phi[numParts];
  double p[numParts],erg[numParts],px[numParts],py[numParts],pz[numParts];
  double e_miss=0.,e_inv=0.,px_inv=0.,py_inv=0.,pz_inv=0.;
  double dumt, dumx, dumy, dumz;
  double originalTrackingParameters[numParts][6];
  TVector3 p3_miss;

  /*********** Define matricies needed to perform kinematic fitting **********/

  // Note: for any matrix m, mT is the transpose of m.  

  // y(i) ==> measured particle quantities. i = 0 is the tagged photon energy,
  // i = 3*n+1,3*n+2,3*n+3 are particle n's |p|,lambda,phi respectively.
  TMatrixD yi(dim,1),y(dim,1); // yi will be intial values
  // x(i) ==> missing particle quantities (0,1,2) = (px,py,pz)
  TMatrixD x(3,1);
  // a(i,j) = dc(i)/dx(j) ==> derivatives of constraint eq's w/r to missing
  // particle quantities x.
  TMatrixD a(numConstraints,3),aT(3,numConstraints);
  // b(i,j) = dc(i)/dy(j) ==> derivatives of constraint eq's w/r to measured
  // quantities y.
  TMatrixD b(numConstraints,dim),bT(dim,numConstraints);
  // c(i) is constraint eq. i ==> (0,1,2,3) are (e,px,py,pz) constraints, any
  // extra mass constraints are i = 4 and up.
  TMatrixD c(numConstraints,1);
  // eps(i) = yi(i)-y(i) are the overall changes in the measured quantities
  TMatrixD eps(dim,1),epsT(1,dim);
  // delta(i) is used to update the measured quantities. After an iteration of
  // the fitting process, y -= delta is how y is updated.
  TMatrixD delta(dim,1);
  // xsi(i) same as delta (above) but for x.
  TMatrixD xsi(3,1);
  // gb is a utility matrix used 
  TMatrixD gb(numConstraints,numConstraints);
  // cov_fit is the covariance matrix OUTPUT by the fit.
  TMatrixD cov_fit(dim,dim);

  /***************************** Initial Set Up ******************************/

  // set up vector of initial measured values:
  yi(0,0) = _ephot_in; // y(0,0) is the photon energy
  for(i = 0; i < numParts; i++){
    // each particle has 3 measured quantites, group them together in y
    mass[i] = _p4in[i].M(); // particle's mass
    double *tempPars = momentum2tracking(_p4in[i], _reconstruction[i], _trackingConversion);
    for(int j=0;j<6;j++) originalTrackingParameters[i][j] = tempPars[j];
    for(int j=0;j<3;j++)
    {
      yi(3*i+(j+1),0) = originalTrackingParameters[i][j];
    }
    //alpha[i] = AlphaTrack(_p4in[i]); // this angle doesn't change
  }
  y = yi; // start off y to be yi

  if(_missingParticle){
    // there's a missing particle...set up it's 3-momentum
    // the intial guess will be the total missing 3-momentum
    x(0,0) = 0.; // px
    x(1,0) = 0.; // py
    x(2,0) = _ephot_in; // pz
    for(i = 0; i < numParts; i++){
      x(0,0) -= _p4in[i].Px();
      x(1,0) -= _p4in[i].Py();
      x(2,0) -= _p4in[i].Pz();
    }
  }

  /************************* Begin Iteration Process *************************/
  for(int iter = 0; iter < 10; iter++){ // it should converge by 10 iterations
    // reset all utility matricies
    a.Zero();
    aT.Zero();
    b.Zero();
    bT.Zero();
    c.Zero();
    delta.Zero();
    xsi.Zero();

    // set kinematic quantities for this iteration:
    for(i = 0; i < numParts; i++)
    {
      double currentTrackingParameters[3] = {y(3*i+1,0), y(3*i+2,0), y(3*i+3,0)};
      double *p3temp = tracking2momentum(currentTrackingParameters, originalTrackingParameters[i], _reconstruction[i], _trackingConversion);
      
      px[i] = p3temp[0];
      py[i] = p3temp[1];
      pz[i] = p3temp[2];
      erg[i] = sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i] + mass[i]*mass[i]);

    }
    if(_missingParticle){
      p3_miss.SetXYZ(x(0,0),x(1,0),x(2,0)); // p3 missing
      e_miss = sqrt(p3_miss.Mag2() + pow(_missingMass,2)); // energy missing
    }

    // Set up the constraint eq's for this iteration:
    // c(0) (energy constraint): e_i - e_f = 0
    // c(i) (p3 constraint): p3_f - p3_i = 0
    c(0,0) = y(0,0) + _targetMass; // e_gamma + e_target;
    c(3,0) = -y(0,0); // -e_gamma (pz constraint eq)
    for(i = 0; i < numParts; i++){
      c(0,0) -= erg[i];
      c(1,0) += px[i];
      c(2,0) += py[i];
      c(3,0) += pz[i];
    }
    if(_missingParticle){
      // add the missing particle quantities to e-p constraints
      c(0,0) -= e_miss;
      c(1,0) += p3_miss.X();
      c(2,0) += p3_miss.Y();
      c(3,0) += p3_miss.Z();
    }

#ifdef __DEBUG__    
    if(iter == 0)
      std::cout << std::endl << "<KinFit::_MainFitter> Constraint equation " 
        << "values (e,px,py,pz): " << std::endl;

    std::cout << c(0,0) << " " << c(1,0) << " " << c(2,0) << " " << c(3,0)
      << std::endl;
#endif /* __DEBUG__ */

    if(_extraC){
      // extra mass constraint...make it c(4).
      e_inv = 0.;
      px_inv = 0.;
      py_inv = 0.;
      pz_inv = 0.;
      for(i = 0; i < numParts; i++){
        if(_extraC_meas[i]){
          // this particle IS in the extra mass constraint
          e_inv += erg[i];
          px_inv += px[i];
          py_inv += py[i];
          pz_inv += pz[i];
        }
      }
      if(_extraC_miss){
        // the missing particle is in the extra mass constraint
        e_inv += e_miss;
        px_inv += p3_miss.X();
        py_inv += p3_miss.Y();
        pz_inv += p3_miss.Z();	
      }
      // set c(4) = e^2 - p^2 - m^2 = 0
      c(4,0) = e_inv*e_inv - px_inv*px_inv - py_inv*py_inv - pz_inv*pz_inv 
        - pow(_invariantMass,2);
    }

    // Set up the derivative matricies:  
    if(_missingParticle){
      // set the derivatives matrix a(i,j) = dc(i)/dx(j)
      for(i = 0; i < 3; i++){
        a(0,i) = -x(i,0)/e_miss; // derivatives of energy constraint eq.
        a(i+1,i) = 1.; // derivatives of p3 constraint eq's 
      }
      if(_extraC && _extraC_miss){
        // the missing particle is in the extra mass constraint
        a(4,0) = 2.*(e_inv*(x(0,0)/e_miss) - px_inv);
        a(4,1) = 2.*(e_inv*(x(1,0)/e_miss) - py_inv);
        a(4,2) = 2.*(e_inv*(x(2,0)/e_miss) - pz_inv);
      }
    }
    aT.Transpose(a); // set a transpose

    // set up derivatives matrix b(i,j) = dc(i)/dy(j):
    // 1st get the tagged photon energy derivatives 
    b(0,0) = 1.;
    b(3,0) = -1;    
    for(i = 0; i < numParts; i++)
    {

      // derivatives for particle i
      double currentTrackingParameters[3] = {y(3*i+1,0), y(3*i+2,0), y(3*i+3,0)};
      double **derivs = constraintDerivs(currentTrackingParameters, originalTrackingParameters[i], _reconstruction[i], _trackingConversion);
      for(j=0;j<4;j++)
      {
        for(k=0;k<3;k++)
        {
          b(j,3*i+(k+1)) = derivs[j][k];
        }
      }

      if(_extraC_meas[i]){
        // this particle is in the extra mass constraint
        for(j = 1; j <= 3; j++)
        {
          int ind = 3*i + j;
          b(4,ind) = 2.*(-e_inv*b(0,ind) - px_inv*b(1,ind) - py_inv*b(2,ind) 
              - pz_inv*b(3,ind));
        }
      }	      
    }
    bT.Transpose(b); // set b transpose

    // gb is a utility matrix we'll need alot below
    gb = b * _cov * bT;
    gb.Invert();

    if(_missingParticle){
      // calculate the change to be made to the missing particle quantities
      xsi = ((aT * gb * a).Invert())*(aT * gb * c);
      // update x
      x -= xsi;
      // calculate the changes to be made to the measured quantities
      delta = _cov * bT * gb * (c - a*xsi);
    }
    else delta = _cov * bT * gb * c;

    // update the measured quantities y
    y -= delta;

    eps = yi - y;
    epsT.Transpose(eps);
    if(TMath::Prob((epsT * bT * gb * b * eps)(0,0),_ndf) < 1.e-10) break;

  } /* end iteration process */

  // get the total changes in the measured quantities
  eps = yi - y;
  epsT.Transpose(eps);

  // get the fit covariance matrix
  if(_missingParticle){
    cov_fit = _cov - _cov*bT*gb*b*_cov 
      + _cov*bT*gb*a*((aT*gb*a).Invert())*aT*gb*b*_cov;
  }
  else 
    cov_fit = _cov - _cov*bT*gb*b*_cov;

  // get chi2 and the pulls
  _pulls.resize(dim);
  for(i = 0; i < dim; i++)
    _pulls[i] = -eps(i,0)/sqrt(_cov(i,i) - cov_fit(i,i));

  _chi2 = (epsT * bT * gb * b * eps)(0,0); // the (0,0)...only...element 

  // update output of measured quantities
  _p4out.resize(numParts);
  for(i = 0; i < numParts; i++) _p4out[i].SetPxPyPzE(px[i],py[i],pz[i],erg[i]);
  _ephot_out = y(0,0);

  // set the missing particle covariance matrix
  if(_missingParticle)
    this->_SetMissingParticleErrors((aT * gb * a).Invert(),x);  
}
//_____________________________________________________________________________

void KinFit::_ResetForNewFit(){
  /// Reset the KinFit object for a new fit (private function).
  // 1st make sure vectors are right size for current p4in
  int size = (int)_p4in.size();
  _p4out.resize(size);
  _extraC_meas.resize(size);
  _pulls.resize(3*size + 1);

  // reset bools to default values
  _missingParticle = false;
  _extraC_miss = false;
  for(int i = 0; i < size; i++) _extraC_meas[i] = false;
  for(int i = 0; i < 3; i++) _sigma_missing[i] = 0.;
}
//_____________________________________________________________________________

double KinFit::Fit(const string &__miss,
    const std::vector<bool> &__constrain_meas,
    bool __constrain_miss,double __invariantMass){
  /// Perform a kinematic fit.
  /** @param miss Name of missing particle
   *  @param constrain_meas Which measured particles are in extra constraint
   *  @param constrain_miss Is the missing particle in extra constraint?
   *  @param invariantMass Mass for extra constraint
   *
   *  
   *  Same as Fit(double,const vector<bool>&,bool,double) execept the missing
   *  particle is passed by name instead of by mass (see that function for 
   *  fitting details). See KinFit::NameToMass for info on which names are
   *  recognized.
   */   
  double mass = -1.;
  if(__miss != string()){
    // a missing particle has been specified
    mass = this->NameToMass(__miss);
  }

  this->Fit(mass,__constrain_meas,__constrain_miss,__invariantMass);
  return this->Prob();
}
//_____________________________________________________________________________

double KinFit::Fit(double __missingMass,
    const std::vector<bool> &__constrain_meas,
    bool __constrain_miss,double __invariantMass){
  /// Perform a kinematic fit.
  /** @param missingMass Mass of missing particle
   *  @param constrain_meas Which measured particles are in extra constraint
   *  @param constrain_miss Is the missing particle in extra constraint?
   *  @param invariantMass Mass for extra constraint
   *
   *  If @a invariantMass is @a NOT specified, then there will be no extra mass
   *  constraint (normal usage) and @a constrain_meas and @a constrain_miss 
   *  will be ignored. If no arguments are specified, then a @a 4C fit is run
   *  to nothing missing with no extra constraints.
   *
   *  The fit is run by the private function KinFit::_MainFitter and returns 
   *  the resulting confidence level.
   *
   *  If an extra mass constraint is specified, then @a constrain_meas must
   *  be the same size as the vector of measured 4-momentum used to set up
   *  this event (also, the indexing of @a constrain_meas is the same as that
   *  vector).
   */   
  // reset for new fit
  this->_ResetForNewFit();

  // check for missing particle
  if(__missingMass >= 0.){
    // there is a missing particle and this is its mass
    _missingMass = __missingMass;
    _missingParticle = true;
  }

  // check for extra mass constraint
  if(__invariantMass > 0.){
    _invariantMass = __invariantMass;
    _extraC = true;
    _extraC_miss = __constrain_miss;
    if(__constrain_meas.size() == _extraC_meas.size())
      _extraC_meas = __constrain_meas;
    else{
      std::cout << "Error! <KinFit::Fit> Boolean vector (argument 2) passed to"
        << " this function has incorrect size. The fit will NOT be "
        << "run." << std::endl;
      this->_SetToBadFit();
      return 0.;
    }
  }

  // run the fit
  this->_MainFitter();

  // return the confidence level
  return this->Prob();
}
//_____________________________________________________________________________

double KinFit::NameToMass(const string &__name) {
  /// Function used to translate names to masses for missing particles.
  /** 
   *  This function is used to translate @a name to a mass. It is used by
   *  Fit to change a missing particle name to a missing mass. All names are
   *  lower case and have no white space,-'s or _'s. 
   *
   *  The function first makes all letters lowercase, then removes any white
   *  space, -'s or _'s from the name. It checks the name against the list of
   *  known names. If it doesn't find a match, it outputs a warning message and
   *  returns 0.
   *
   */
  if(__name == string()) return -1.;

  string name = __name;
  // make all letters lower case
  std::transform(name.begin(),name.end(),name.begin(),(int(*)(int))std::tolower);
  // remove any white space, -'s or _'s
  for(string::size_type pos = 0; pos < name.size(); pos++){
    if(name[pos] == ' ' || name[pos] == '-' || name[pos] == '_'){
      name.erase(pos,1);
      pos--;
    }
  }

  // mesons
  if(name == "pi" || name == "pi+" || name == "pi-") return 0.13957;
  if(name == "pi0") return 0.13498;
  if(name == "k" || name == "k+" || name == "k-") return 0.493677;
  if(name == "k0" || name == "kshort" || name == "klong") return 0.497672;
  if(name == "eta") return 0.54730;
  if(name == "rho" || name == "rho+" || name == "rho-" || name == "rho0")
    return 0.7693;
  if(name == "omega") return 0.78257;
  if(name == "eta'" || name == "etaprime") return 0.95778;
  if(name == "phi") return 1.019417;

  // baryons
  if(name == "p" || name == "proton" || name == "antiproton" || name == "pbar"
      || name == "p-" || name == "p+") return 0.93827;
  if(name == "neutron" || name == "n") return 0.93957;

  // if we get here it's not on the list
  std::cout << "Warning! <KinFit::NameToMass> Unknown particle name: " 
    << __name << ". Mass will be returned as 0." << std::endl;
  return 0.;
}
//_____________________________________________________________________________

void KinFit::_SetToBadFit(){
  /// Set fit output to @a bad values.
  /** Sets pulls and \f$ \chi^2 \f$ to 666, sets degrees of freedom and output
   *  kinematic quantites to 0.
   */
  for(int i = 0; i < (int)_pulls.size(); i++) _pulls[i] = 666.;
  _chi2 = 666.;
  _ndf = 0;
  _ephot_out = 0.;
  for(int i = 0; i < (int)_p4out.size(); i++) _p4out[i].SetXYZT(0.,0.,0.,0.);  
}
//_____________________________________________________________________________

void KinFit::_FitToMissingTaggedPhoton(){
  /// Run a fit to no missing particles but the tagged photon is missing

  _ndf = 3;

  const int numParts = (int)_p4in.size();
  const int dim = 3*numParts; // number of measured quantites
  int i,j;
  double alpha[numParts],mass[numParts],lambda[numParts],phi[numParts];
  double p[numParts],erg[numParts],px[numParts],py[numParts],pz[numParts];

  // utility matrix names/notation is the same as in _MainFitter...see there
  // for descriptions. Note that here we'll define a local cov, that will be
  // a copy of the data member _cov, except the tagged photon entries will be
  // removed.
  TMatrixD a(4,1),aT(1,4),b(4,dim),bT(dim,4),c(4,1),eps(dim,1),epsT(1,dim);
  TMatrixD cov(dim,dim),delta(dim,1),xsi(1,1),gb(4,4);
  TMatrixD yi(dim,1),y(dim,1),cov_fit(dim,dim),x(1,1);

  // set up the local copy of the covariance matrix
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      cov(i,j) = _cov(i+1,j+1);
    }
  }

  // set up vector of initial measured values:
  for(i = 0; i < numParts; i++){
    // each particle has 3 measured quantites, group them together in y
    yi(3*i,0) = _p4in[i].P(); // |p|
    yi(3*i+1,0) = LambdaTrack(_p4in[i]); // tracking angle lambda
    yi(3*i+2,0) = PhiTrack(_p4in[i]); // tracking angle phi
    alpha[i] = AlphaTrack(_p4in[i]); // this angle doesn't change
    mass[i] = _p4in[i].M(); // particle's mass
  }
  y = yi; // start off y to be yi


  // make the initial guess for the tagged energy satisfy the energy constraint
  x(0,0) = -_targetMass;
  for(i = 0; i < numParts; i++) x(0,0) += _p4in[i].E();

  /************************* Begin Iteration Process *************************/
  for(int iter = 0; iter < 10; iter++){ // it should converge by 10 iterations
    // reset all utility matricies
    a.Zero();
    aT.Zero();
    b.Zero();
    bT.Zero();
    c.Zero();
    delta.Zero();
    xsi.Zero();

    // set kinematic quantities for this iteration:
    for(i = 0; i < numParts; i++){
      lambda[i] = y(3*i+1,0);
      phi[i] = y(3*i+2,0);
      p[i] = y(3*i,0);
      erg[i] = sqrt(p[i]*p[i] + mass[i]*mass[i]);

      px[i] = p[i]*(cos(lambda[i])*sin(phi[i])*cos(alpha[i]) 
          - sin(lambda[i])*sin(alpha[i]));
      py[i] = p[i]*(cos(lambda[i])*sin(phi[i])*sin(alpha[i]) 
          + sin(lambda[i])*cos(alpha[i]));
      pz[i] = p[i]*cos(lambda[i])*cos(phi[i]);
    }

    // Set up the constraint eq's for this iteration:
    // c(0) (energy constraint): e_i - e_f = 0
    // c(i) (p3 constraint): p3_f - p3_i = 0
    c(0,0) = x(0,0) + _targetMass; // e_gamma + e_target;
    c(3,0) = -x(0,0); // -e_gamma (pz constraint eq)
    for(i = 0; i < numParts; i++){
      c(0,0) -= erg[i];
      c(1,0) += px[i];
      c(2,0) += py[i];
      c(3,0) += pz[i];
    }

#ifdef __DEBUG__    
    if(iter == 0) 
      std::cout << std::endl << "<KinFit::_MainFitter> Constraint equation " 
        << "values (e,px,py,pz): " << std::endl;
    std::cout << c(0,0) << " " << c(1,0) << " " << c(2,0) << " " << c(3,0)
      << std::endl;
#endif /* __DEBUG__ */

    // set up dervative matrix a(i,0) = dc(i)/dx
    a.Zero();
    a(0,0) = 1.;
    a(3,0) = -1.;
    aT.Transpose(a);

    // set up derivatives matrix b(i,j) = dc(i)/dy(j):
    // 1st get the tagged photon energy derivatives 
    for(i = 0; i < numParts; i++){
      // derivatives for particle i
      /* derivatives of constraint eq's w/r to p */
      b(0,3*i) = -p[i]/erg[i];
      b(1,3*i) = cos(lambda[i])*sin(phi[i])*cos(alpha[i]) 
        - sin(lambda[i])*sin(alpha[i]);
      b(2,3*i) = cos(lambda[i])*sin(phi[i])*sin(alpha[i]) 
        + sin(lambda[i])*cos(alpha[i]);
      b(3,3*i) = cos(lambda[i])*cos(phi[i]);
      /* derivatives of constraint eq's w/r to lambda */
      b(0,3*i+1) = 0.;
      b(1,3*i+1) = p[i]*(-sin(lambda[i])*sin(phi[i])*cos(alpha[i])
          - cos(lambda[i])*sin(alpha[i]));
      b(2,3*i+1) = p[i]*(-sin(lambda[i])*sin(phi[i])*sin(alpha[i])
          + cos(lambda[i])*cos(alpha[i]));
      b(3,3*i+1) = p[i]*(-sin(lambda[i])*cos(phi[i]));
      /* derivatives of constraint eq's w/r to phi */
      b(0,3*i+2) = 0.;
      b(1,3*i+2) = p[i]*cos(lambda[i])*cos(phi[i])*cos(alpha[i]);
      b(2,3*i+2) = p[i]*cos(lambda[i])*cos(phi[i])*sin(alpha[i]);
      b(3,3*i+2) = p[i]*cos(lambda[i])*(-sin(phi[i]));

    }
    bT.Transpose(b); // set b transpose

    // gb is a utility matrix we'll need alot below
    gb = b * cov * bT;
    gb.Invert();

    // get correction to x
    xsi = ((aT*gb*a).Invert())*(aT*gb*c);
    x -= xsi;

    // get correction to y
    delta = cov*bT*gb*(c - a*xsi);
    y -= delta;

  } /* end iteration process */

  eps = yi - y;
  epsT.Transpose(eps);

  // get covariance matrix after fitting
  cov_fit = cov - cov*bT*gb*b*cov 
    + cov*bT*gb*a*((aT*gb*a).Invert())*aT*gb*b*cov;

  // get chi2 and the pulls

  // for the pulls, keep them the size they would be if we were using the 
  // measured photon...so pulls[0] will be just set to 0 and be meaningless
  _pulls.resize(dim+1);
  _pulls[0] = 0.;
  for(i = 0; i < dim; i++)
    _pulls[i+1] = -eps(i,0)/sqrt(cov(i,i) - cov_fit(i,i));

  _chi2 = (epsT * bT * gb * b * eps)(0,0); // the (0,0)...only...element 

  // update output of measured quantities
  _p4out.resize(numParts);
  for(i = 0; i < numParts; i++) _p4out[i].SetPxPyPzE(px[i],py[i],pz[i],erg[i]);
  _ephot_out = x(0,0);
}
//_____________________________________________________________________________

void KinFit::_SetMissingParticleErrors(const TMatrixD &__missingCov,
    const TMatrixD &__x){
  /// Calculate the missing particle errors from the last fit

  double p = sqrt(__x(0,0)*__x(0,0) + __x(1,0)*__x(1,0) + __x(2,0)*__x(2,0));
  TLorentzVector p4(__x(0,0),__x(1,0),__x(2,0),sqrt(p*p + pow(_missingMass,2)));

  // kinematic quantities in tracking coordinates
  double phi = PhiTrack(p4);
  double lam = LambdaTrack(p4);
  double alpha = AlphaTrack(p4);

  // missing particle errors in lab coordinates
  double sigma_px_2 = __missingCov(0,0);
  double sigma_py_2 = __missingCov(1,1);
  double sigma_pz_2 = __missingCov(2,2);

  // jacobian elements we need
  double dp_dpx = __x(0,0)/p;
  double dp_dpy = __x(1,0)/p;
  double dp_dpz = __x(2,0)/p;

  double dlam_dpx = (-sin(alpha) - sin(lam)*dp_dpx)/(p*cos(lam));
  double dlam_dpy = (cos(alpha) - sin(lam)*dp_dpy)/(p*cos(lam));
  double dlam_dpz = -tan(lam)*dp_dpz/p;

  double dphi_dpx = (sin(lam)*p*cos(phi)*dlam_dpx - dp_dpx*cos(lam)*cos(phi))
    /(-sin(phi)*p*cos(lam));
  double dphi_dpy = (sin(lam)*p*cos(phi)*dlam_dpy - dp_dpy*cos(lam)*cos(phi))
    /(-sin(phi)*p*cos(lam));
  double dphi_dpz = (1 + sin(lam)*p*cos(phi)*dlam_dpz 
      - dp_dpz*cos(lam)*cos(phi))/(-sin(phi)*p*cos(lam));

  // get error on p
  _sigma_missing[0] = sqrt(sigma_px_2*(dp_dpx*dp_dpx) 
      + sigma_py_2*(dp_dpy*dp_dpy)
      + sigma_pz_2*(dp_dpz*dp_dpz));

  // get error on lambda
  _sigma_missing[1] = sqrt(sigma_px_2*(dlam_dpx*dlam_dpx) 
      + sigma_py_2*(dlam_dpy*dlam_dpy)
      + sigma_pz_2*(dlam_dpz*dlam_dpz));

  // get error on phi
  _sigma_missing[2] = sqrt(sigma_px_2*(dphi_dpx*dphi_dpx) 
      + sigma_py_2*(dphi_dpy*dphi_dpy)
      + sigma_pz_2*(dphi_dpz*dphi_dpz));

}
//_____________________________________________________________________________
