// KinFit class header file. -*- C++ -*-
/** @file kinfit/KinFit.h
 * @brief KinFit class defintion file.
 */
#ifndef _KinFit_H
#define _KinFit_H
// System Headers:
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
// ROOT Headers:
#include "TMatrixD.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
// Local Headers
#include "KinFitUtilities.h"

using namespace std;
//_____________________________________________________________________________

class KinFit {
  
private:
  // Data Members (private):

  /* kinematic fitting statistical quantities */
  std::vector<double> _pulls; ///< Pull quantities of last fit
  double _chi2; ///< \f$ \chi^2 \f$ of last fit
  int _ndf; ///< Number of degrees-of-freedom of last fit

  /* kinematic quantities */
  double _ephot_in; /// < Photon energy (in)
  double _ephot_out; /// < Photon energy (out)
  std::vector<TLorentzVector> _p4in; ///< Particle 4-momenta (in)
  std::vector<TLorentzVector> _p4out; ///< Particle 4-momenta (out)
  std::vector<string> _reconstruction; ///< Particle reconstruction 0 - Drift chamber
                                    ///                          1 - Calorimeter
  double _targetMass; ///< Target mass

  string _trackingConversion; ///< Tracking conversion 

  /* covariance matrix info */
  TMatrixD _cov; ///< Covariance matrix
  double _sigma_missing[3]; ///< Fit errors on missing quantities

  /* missing particle info */
  double _missingMass; ///< Mass of the missing particle
  bool _missingParticle; ///< Is there a missing particle?

  /* extra mass constraint info */
  bool _extraC; ///< Is there an extra mass constraint?
  double _invariantMass; ///< Invariant mass used in extra mass constraint
  std::vector<bool> _extraC_meas; ///< Which measured particles in constraint?
  bool _extraC_miss; ///< Is missing particle used in extra mass constraint?

  // Functions (private):

  // main kinematic fit function
  void _MainFitter();

  // fit to missing tagged photon (3C fit)
  void _FitToMissingTaggedPhoton();

  // get ready for a new fit
  void _ResetForNewFit();

  // copy the KinFit data members
  void _Copy(const KinFit &__kfit);

  // sets output quantities for a bad fit
  void _SetToBadFit();

  // set the missing particle errors
  void _SetMissingParticleErrors(const TMatrixD &__missingCov,
				 const TMatrixD &__x);

public:
  // Create/Copy/Destroy:
  KinFit();

  KinFit(const KinFit &__kfit){ 
    /// Copy Constructor
    this->_Copy(__kfit);
  }

  virtual ~KinFit(){ 
    /// Destructor
    _pulls.clear();
    _p4in.clear();
    _reconstruction.clear();
    _p4out.clear();
    _extraC_meas.clear();
  }

  KinFit& operator=(const KinFit &__kfit){
    /// Assignment operator
    this->_Copy(__kfit);
    return *this;
  }

  // Setters:

  inline void SetCovMatrix(const TMatrixD &__covMat){ 
    /// Set the covariance matrix.
    _cov.ResizeTo(__covMat);
    _cov = __covMat;
  }

  inline void SetP4(const std::vector<TLorentzVector> &__p4){
    /// Set the input 4-momenta.
    _p4in = __p4;
  }
  
  inline void SetReconstruction(const std::vector<string> &__rec){
    /// Set the input 4-momenta.
    _reconstruction = __rec;
  }
  
  inline void SetPhotonEnergy(double __erg){ 
    /// Set the input tagged photon energy.
    _ephot_in = __erg;
  }

  inline void SetTargetMass(double __mass){
    /// Set the target mass (the ctor sets this to be \f$ m_p \f$)
    _targetMass = __mass;
  }

  inline void SetTrackingConversion(string __tc){
    /// Tells the fit how to convert 4-vectors to the tracking parameters.
    _trackingConversion = __tc;
  }

  inline void SetEvent(double __ephot,const std::vector<TLorentzVector> &__p4,const std::vector<string> &__reconstruction,
		       const TMatrixD &__covMat,double __targMass = -1., string __trackingConversion="EASY"){
    /// Set all input quantities.
    /** @param ephot Input tagged photon energy
     *  @param p4 Vector of detected particle 4-momenta
     *  @param covMat Covariance matrix
     *  @param targMass Target mass 
     *
     *  If @a targMass isn't specified, then the target mass is @a NOT set.
     *  The target mass defaults to \f$ m_p \f$ in the ctor.
     *
     */
    if(__targMass > 0.) _targetMass = __targMass;
    _ephot_in = __ephot;
    _p4in = __p4;
    _reconstruction = __reconstruction;
    _cov.ResizeTo(__covMat);
    _cov = __covMat;
    _trackingConversion = __trackingConversion;
  }

  // Getters:

  inline double Chi2() const {
    /// Return \f$ \chi^2 \f$ of the last fit.
    return _chi2;
  }

  inline int Ndf() const {
    /// Return the number of degrees of freedom of the last fit
    return _ndf;
  }

  inline double GetPull(int __n) const {
    /// Returns the @a n pull quantity
    /** Pulls are stored as 
     *	\f$ (E_{\gamma},p_1,\lambda_1,\phi_1,...,p_n,\lambda_n,\phi_n) \f$ 
     */
    return _pulls[__n];
  }

  inline const TLorentzVector& FitP4(int __n) const {
    /// Returns the particle @a n 4-momentum after fitting
    /** Same indexing as vector used to set input 4-momentum */
    return _p4out[__n];
  }

  inline double FitPhotonEnergy() const {
    /// Returns the tagged photon energy after fitting
    return _ephot_out;
  }

  inline double GetMissingError(int __n) const {
    /// Returns the @a n missing particle error (p,lambda,phi)
    return _sigma_missing[__n];
  }

  inline const TMatrixD& GetCovMat() const {
    /// Returns the covariance matrix
    return _cov;
  }

  // Functions:

  double Prob() const {
    /// Returns the confidence level of the last fit
    return TMath::Prob(_chi2,_ndf);
  }

  // converts names to masses for the missing particle
  static double NameToMass(const string &__name);

  // Kinematic fit functions:
  double Fit(const string &__miss,
	     const std::vector<bool> &__constrain_meas = std::vector<bool>(),
	     bool __constrain_miss = false,double __invariantMass = -1.);

  double Fit(double __missingMass = -1.,
	     const std::vector<bool> &__constrain_meas = std::vector<bool>(),
	     bool __constrain_miss = false,double __invariantMass = -1.);

  inline double Fit2MissingTagged(){
    /// Fit to no missing particle but instead a missing tagged photon.
    /** Set the covariance matrix the same way for this fit as for all other
     *  fits. This function will rearrange the matrix itself to account for the
     *  different role of the tagged photon in the fit (from measured to 
     *  missing). Also, the pulls will have the same order and size as for all
     *  other fits...the first term, which is normally the tagged photon pull,
     *  will be set to 0 (since it is meaningless for this fit).
     *
     * <b> Returns: </b> Confidence level of the fit.
     */
    
    // reset for new fit
    this->_ResetForNewFit();
    this->_FitToMissingTaggedPhoton();
    return this->Prob();
  }

};
//_____________________________________________________________________________

#endif /* _KinFit_H */


