// DKinFit class header file. -*- C++ -*-
/** @file DKinFit/DKinFit.h
 * @brief DKinFit class defintion file.
 */
#ifndef _DKinFit_H
#define _DKinFit_H
// System Headers:
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
// ROOT Headers:
#include "TMatrixD.h"
#include "TMath.h"
// Local Headers
#include "DLorentzVector.h"
#include "DVector3.h"
#include "DKinematicData.h"

using namespace std;
//_____________________________________________________________________________

class DKinFit {

  private:
    // Data Members (private):

    int _verbose; ///< Level of verbosity

    /* kinematic fitting statistical quantities */
    std::vector<double> _pulls; ///< Pull quantities of last fit
    double _chi2; ///< \f$ \chi^2 \f$ of last fit
    int _ndf; ///< Number of degrees-of-freedom of last fit
    bool _useTimingInformation; ///< Tells us whether or not to use whatever
                                ///< timing information is available in the final
                                ///< calculation of the confidence level.

    /* kinematic quantities */
    std::vector<DKinematicData> _kDataInitial_in; ///< Initial particle 4-momenta (in)
    std::vector<DKinematicData> _kDataFinal_in; ///< Final particle 4-momenta (in)
    std::vector<DKinematicData> _kDataInitial_out; ///< Initial particle 4-momenta (out)
    std::vector<DKinematicData> _kDataFinal_out; ///< Final particle 4-momenta (out)

    /* covariance matrix info */
    TMatrixD _cov; ///< Covariance matrix
    double _sigma_missing[3]; ///< Fit errors on missing quantities

    /* missing particle info */
    double _missingMass; ///< Mass of the missing particle
    bool _missingParticle; ///< Is there a missing particle?

    /* extra mass constraint info */
    int _extraC; ///< How many extra mass constraints are there?
    double _invariantMass; ///< Invariant mass used in extra mass constraint
    std::vector<bool> _extraC_meas; ///< Which measured particles in constraint?
    std::vector<double> _constraintMasses; ///< Masses for the extra constraint equations
    std::vector< std::vector<int> > _constraintParticles; ///< Particles in the extra constraint equations
                                                          ///< -1 is the missing mass
    std::vector<bool> _extraC_miss; ///< Is missing particle used in extra mass constraint?

    // Functions (private):
    // main kinematic fit function
    void _MainFitter();

    // get ready for a new fit
    void _ResetForNewFit();

    // copy the DKinFit data members
    void _Copy(const DKinFit &__kfit);

    // sets output quantities for a bad fit
    void _SetToBadFit();

    // set the missing particle errors
    void _SetMissingParticleErrors(const TMatrixD &__missingCov, const TMatrixD &__x);

  public:
    // Create/Copy/Destroy:
    DKinFit();

    DKinFit(const DKinFit &__kfit){ 
      /// Copy Constructor
      this->_Copy(__kfit);
    }

    virtual ~DKinFit(){ 
      /// Destructor
      _pulls.clear();
      _kDataInitial_in.clear();
      _kDataFinal_in.clear();
      _kDataInitial_out.clear();
      _kDataFinal_out.clear();
      _extraC_meas.clear();
    }
    
    DKinFit& operator=(const DKinFit &__kfit){
      /// Assignment operator
      this->_Copy(__kfit);
      return *this;
    }

    // Setters:

    // get ready for a new fit
    inline void ResetForNewFit() { this->_ResetForNewFit();}

    inline void SetCovMatrix(const TMatrixD &__covMat){ 
      /// Set the covariance matrix.
      _cov.ResizeTo(__covMat);
      _cov = __covMat;
    }

    inline void SetVerbose(int __verbose){_verbose = __verbose;}

    /// Set the initial tracks
    inline void SetInitial(std::vector<DKinematicData> &__kd)
    { _kDataInitial_in = __kd; }

    /// Set the final tracks
    inline void SetFinal(std::vector<DKinematicData> &__kd)
    { _kDataFinal_in = __kd; }

    /// Get the momenta BEFORE the kinematic fit (initial particles)
    inline std::vector<DKinematicData> GetInitial_in() { return _kDataInitial_in;}
    /// Get the momenta BEFORE the kinematic fit (final state particles)
    inline std::vector<DKinematicData> GetFinal_in()   { return _kDataFinal_in;}

    /// Get the momenta AFTER the kinematic fit (initial particles)
    inline std::vector<DKinematicData> GetInitial_out() { return _kDataInitial_out;}
    /// Get the momenta AFTER the kinematic fit (final particles)
    inline std::vector<DKinematicData> GetFinal_out()   { return _kDataFinal_out;}

    // Getters:

    /// Return \f$ \chi^2 \f$ of the last fit.
    inline double Chi2() const { return _chi2; }

    /// Return the number of degrees of freedom of the last fit
    inline int Ndf() const { return _ndf; }

    /// Returns the @a n pull quantity
    /** Pulls are stored as 
     *	\f$ (p_{0x},p_{0y},p_{0z},p_{1x},...) \f$.
     *  Initial particles first, then final particles.
     */
    inline double GetPull(int __n) const { return _pulls[__n]; }

    //This function actually does use (p,lambda,phi) even though everything everthing else in this class is px,py,pz
    /// Returns the @a n missing particle error (p,lambda,phi)
    inline double GetMissingError(int __n) const { return _sigma_missing[__n]; }

    /// Returns the covariance matrix
    inline const TMatrixD& GetCovMat() const { return _cov; }

    // Functions:
    /// Returns the confidence level of the last fit
    double Prob() const { return TMath::Prob(_chi2,_ndf); }

    // converts names to masses for the missing particle
    static double NameToMass(const string &__name);

    /// Set whether or not there's a missing particle
    void SetMissingParticle(const double __missingMass);

    /// Set each extra mass constraint
    void SetMassConstraint(const double __constraintMass, const vector<int> constrainedParticles);

    /// Use whatever timing information is available
    void UseTimingInformation() { _useTimingInformation = true; }

    // main kinematic fit function
    void Fit();
    void FitTwoGammas(const float __missingMass, const float errmatrixweight);

    // Utility function
    DLorentzVector MissingMomentum(vector<const DKinematicData*> initial, vector<const DKinematicData*> final);

};
//_____________________________________________________________________________

#endif /* _DKinFit_H */


