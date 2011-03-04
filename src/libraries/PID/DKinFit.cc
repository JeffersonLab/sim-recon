#include <stdio.h>
#include <stdlib.h>

#include <cctype>
#include <algorithm>
#include "PID/DKinFit.h"
using namespace std;
/** @file DKinFit.C
 * @brief DKinFit class source file.
 */
//_____________________________________________________________________________
/** @class DKinFit
 * @brief Utility class used to perform kinematic fitting. 
 * 
 * Specify the particles (i.e. DKinematicData) to be fit using SetInitial() and
 * SetFinal() .
 * Each DKinematicData must have a proper covariance matrix set.
 *
 * If needed you also SetMissingParticle() or SetMassConstraint()
 * appropriately
 *
 * Then just Fit() !
 *
 * The results of the fit can be obtained thru various DKinFit member functions.
 * All @a getter member functions return results of the last fit run.
 *
 */
//_____________________________________________________________________________

DKinFit::DKinFit(){
  /// Default Constructor
  _ndf = 0;
  _chi2 = -1.;
  _missingMass = -1.;
  _missingParticle = false;
  _extraC = 0;
  _invariantMass = -1.;
  _extraC_miss.clear();
  _useTimingInformation = false;
  _verbose = 0;
}
//_____________________________________________________________________________

void DKinFit::_Copy(const DKinFit &__kfit){
  /// Copy @a kfit data members.
  _pulls = __kfit._pulls;
  _chi2 = __kfit._chi2;
  _ndf = __kfit._ndf;
  _kDataInitial_in = __kfit._kDataInitial_in;
  _kDataFinal_in = __kfit._kDataFinal_in;
  _kDataInitial_out = __kfit._kDataInitial_out;
  _kDataFinal_out = __kfit._kDataFinal_out;
  _cov.ResizeTo(__kfit._cov);
  _cov = __kfit._cov;
  _missingMass = __kfit._missingMass;
  _missingParticle = __kfit._missingParticle;
  _extraC = __kfit._extraC;
  _invariantMass = __kfit._invariantMass;
  _extraC_meas = __kfit._extraC_meas;
  _extraC_miss = __kfit._extraC_miss;
  _useTimingInformation = __kfit._useTimingInformation;
  _verbose = __kfit._verbose;
  for(int i = 0; i < 3; i++) _sigma_missing[i] = __kfit._sigma_missing[i];
}

//_____________________________________________________________________________
void DKinFit::Fit()
{
  /// Private function used internally to perform kinematic fitting.
  _kDataInitial_out.clear();
  _kDataFinal_out.clear();
  _extraC_miss.clear();
  const int numInitial = (int)_kDataInitial_in.size();
  const int numFinal = (int)_kDataFinal_in.size();
  const int numParts = numInitial + numFinal;
  const int numypar = 3; // number of measured quantites
  const int dim = numypar*numParts; // number of measured quantites
  if(_verbose>0)
  {
    cerr << "numInitial: " << numInitial << endl;
    cerr << "numFinal: " << numFinal << endl;
    cerr << "numParts: " << numParts << endl;
    cerr << "dim: " << dim << endl;
  }
  // Resize the covariance matrix
  _cov.ResizeTo(dim, dim);

  //const int dim = 3*numParts + 1; // number of measured quantites
  // check to see that covariance matrix size is consistent with # particles
  if(_cov.GetNrows() != dim || _cov.GetNcols() != dim)
  {
    // wrong size
    std::cout << "Error! <DKinFit::_MainFitter> Covariance matrix size is NOT "
      << "correct for current number of particles. For " << numParts
      << " particles, the covariance matrix should be " << dim << " x "
      << dim << ", but the matrix passed in is " << _cov.GetNrows()
      << " x " << _cov.GetNcols() << std::endl;
    abort();
  }

  for(int i=0;i<numInitial;i++)
  {
    _kDataInitial_out.push_back(_kDataInitial_in[i]); /// Make sure we're not using the same pointer
    for(int j=0;j<numypar;j++)
    {
      for(int k=0;k<numypar;k++)
      {
        _cov(numypar*i+j,numypar*i+k) = (_kDataInitial_out[i].errorMatrix())(j,k);
      }
    }
  }
  for(int i=0;i<numFinal;i++)
  {
    _kDataFinal_out.push_back(_kDataFinal_in[i]); /// Make sure we're not using the same pointer
    for(int j=0;j<numypar;j++)
    {
      for(int k=0;k<numypar;k++)
      {
        _cov( numypar*i+j + numypar*numInitial, numypar*i+k + numypar*numInitial) = (_kDataFinal_out[i].errorMatrix())(j,k);
      }
    }
  }

  if(_verbose>1)
  {
    cerr << "Total error matrix: " << endl;
    for(int j=0;j<_cov.GetNrows();j++)
    {
      for(int k=0;k<_cov.GetNcols();k++)
      {
        cerr << _cov(j,k) << " ";
      }
      cerr << endl;
    }
  }

  // Get the number of constraint equations
  int numConstraints = 4;
  if(_extraC) numConstraints += _extraC;
  // For Two gammas
  //numConstraints = 4;

  // get the degrees of freedom
  _ndf = 4;
  if(_missingParticle) _ndf -= 3;
  _ndf += _extraC;
  //if(_extraC) _ndf += 1;
  // For Two gammas
  //_ndf = 1;

  if(_verbose>0) cerr << "numConstraints/ndf: " << numConstraints << " " << _ndf << endl;

  int i;
  // The following were changed to vectors to avoid compiler warnings on SunOS 11/12/2007 D.L.
  vector<double> mass(numParts); //double mass[numParts];
  vector<double> energy(numParts); //double energy[numParts];
  double missingEnergy = 999.0;
  vector<int> initialOrFinal(numParts); //int initialOrFinal[numParts];
  // The following 3 lines were commented out to avoid compiler warnings 7/31/07 D.L.
  //double p[numParts],erg[numParts],px[numParts],py[numParts],pz[numParts];
  //double e_miss=0.,e_inv=0.,px_inv=0.,py_inv=0.,pz_inv=0.;
  //double dumt, dumx, dumy, dumz;
  DVector3 p3_miss;
  // The following were changed to vectors to avoid compiler warnings on SunOS 11/12/2007 D.L.
  vector<double> pxMeas(_extraC), pyMeas(_extraC), pzMeas(_extraC), EMeas(_extraC);
  vector<double> pxTot(_extraC), pyTot(_extraC), pzTot(_extraC), ETot(_extraC), massTot(_extraC);
  //double pxMeas[_extraC], pyMeas[_extraC], pzMeas[_extraC], EMeas[_extraC];
  //double pxTot[_extraC], pyTot[_extraC], pzTot[_extraC], ETot[_extraC], massTot[_extraC];

  for(int i=0;i<numParts;i++)
  {
    if(i<numInitial) initialOrFinal[i] = +1;
    else             initialOrFinal[i] = -1;
  }

  /*********** Define matricies needed to perform kinematic fitting **********/

  // Note: for any matrix m, mT is the transpose of m.  
  //
  DLorentzVector missingMomentum;

  // y(i) ==> measured particle quantities.
  // i = 3*n+0,3*n+1,3*n+2 are particle n's p_x,p_y,p_z respectively.
  TMatrixD yi(dim,1),y(dim,1); // yi will be intial values
  // x(i) ==> missing particle quantities (0,1,2) = (px,py,pz)
  TMatrixD x(3,1);
  //TMatrixD x(4,1);
  // a(i,j) = dc(i)/dx(j) ==> derivatives of constraint eq's w/r to missing
  // particle quantities x.
  TMatrixD a(numConstraints,3),aT(3,numConstraints);
  //TMatrixD a(numConstraints,4),aT(4,numConstraints);
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
  //TMatrixD xsi(4,1);
  // gb is a utility matrix used 
  TMatrixD gb(numConstraints,numConstraints);
  // cov_fit is the covariance matrix OUTPUT by the fit.
  TMatrixD cov_fit(dim,dim);

  /***************************** Initial Set Up ******************************/

  // set up vector of initial measured values:
  for(i = 0; i < numInitial; i++)
  {
    mass[i] = _kDataInitial_in[i].mass(); // particle's mass
    yi(numypar*i + 0,0) = _kDataInitial_in[i].px();
    yi(numypar*i + 1,0) = _kDataInitial_in[i].py();
    yi(numypar*i + 2,0) = _kDataInitial_in[i].pz();
    if(_verbose) cerr << "initial masses: " << mass[i] << endl;
  }
  for(i = 0; i < numFinal; i++)
  {
    mass[numInitial + i] = _kDataFinal_in[i].mass(); // particle's mass
    yi(numypar*(numInitial+i) + 0,0) = _kDataFinal_in[i].px();
    yi(numypar*(numInitial+i) + 1,0) = _kDataFinal_in[i].py();
    yi(numypar*(numInitial+i) + 2,0) = _kDataFinal_in[i].pz();
    if(_verbose) cerr << "final masses: " << mass[i+numInitial] << endl;
  }
  y = yi; // start off y to be yi

  // For Two gammas
  if(_missingParticle)
  {
    x(0,0) = 0.0;
    x(1,0) = 0.0;
    x(2,0) = 0.0;
    for(i = 0; i < numInitial; i++)
    {
      x(0,0) += _kDataInitial_in[i].px();
      x(1,0) += _kDataInitial_in[i].py();
      x(2,0) += _kDataInitial_in[i].pz();
    }
    for(i = 0; i < numFinal; i++)
    {
      x(0,0) -= _kDataFinal_in[i].px();
      x(1,0) -= _kDataFinal_in[i].py();
      x(2,0) -= _kDataFinal_in[i].pz();
    }
  }

  /************************* Begin Iteration Process *************************/
  for(int iter = 0; iter < 10; iter++)
  { 
    // it should converge by 10 iterations
    // reset all utility matricies
    a.Zero();
    aT.Zero();
    b.Zero();
    bT.Zero();
    c.Zero();
    delta.Zero();
    xsi.Zero();

    if(_verbose>1)
    {
      cerr << "y: " << endl;
      for(int j=0;j<y.GetNrows();j++)
      {
        for(int k=0;k<y.GetNcols();k++)
        {
          cerr  << j << " " << k << " " << y(j,k) << " ";
        }
        cerr << endl;
      }
    }

    if(_verbose>1)
    {
      cerr << "x: " << endl;
      for(int j=0;j<x.GetNrows();j++)
      {
        for(int k=0;k<x.GetNcols();k++)
        {
          cerr  << j << " " << k << " " << x(j,k) << " ";
        }
        cerr << endl;
      }
    }

    // Set up the constraint eq's for this iteration:
    // c(0) (energy constraint): e_f - e_i = 0
    // c(i) (p3 constraint): p3_f - p3_i = 0
    if(_verbose>1) cerr << "setting up the constraints " << endl;
    for(int j=0;j<4;j++) c(j,0) = 0.; 
    for(int k=0;k<numInitial;k++)
    {
      int offset = k*3;
      energy[k] = sqrt(pow(mass[k],2) + pow(y(offset+0,0),2) + pow(y(offset+1,0),2) + pow(y(offset+2,0),2));
      c(0,0) += y(offset+0,0);
      c(1,0) += y(offset+1,0);
      c(2,0) += y(offset+2,0);
      c(3,0) += energy[k];
    }
    for(int k=0;k<numFinal;k++)
    {
      int offset = (k+numInitial)*3;
      energy[k+numInitial]  = sqrt(pow(mass[k+numInitial],2) + pow(y(offset+0,0),2) + pow(y(offset+1,0),2) + pow(y(offset+2,0),2));
      c(0,0) -= y(offset+0,0);
      c(1,0) -= y(offset+1,0);
      c(2,0) -= y(offset+2,0);
      c(3,0) -= energy[k+numInitial];
    }
    if(_missingParticle)
    {
      missingEnergy = sqrt(pow(_missingMass,2) + pow(x(0,0),2) + pow(x(1,0),2) + pow(x(2,0),2));
      c(0,0) -= x(0,0);
      c(1,0) -= x(1,0);
      c(2,0) -= x(2,0);
      c(3,0) -= missingEnergy;
    }

    /// If there is a extra mass constraints
    for(int j=0;j<_extraC;j++)
    {
      int cindex = 4 + j;
      int yoffset = 3*numInitial;
      int numConstraintParticles = (int)_constraintParticles[j].size();
      pxMeas[j]=pyMeas[j]=pzMeas[j]=EMeas[j]=0.0;
      pxTot[j]=pyTot[j]=pzTot[j]=ETot[j]=0.0;
      massTot[j] = 0.0;
      _extraC_miss.push_back(false); // Assume missing particle is not in constraint
      /// Add in the particles to the constraint
      for(int k=0;k<numConstraintParticles;k++)
      {
        int partIndex = _constraintParticles[j][k];
        int Eindex = numInitial + partIndex;
        /// Use a measured particle
        if(partIndex>=0 && partIndex<numFinal)
        {
          if(_verbose>1) cerr << "Using constraint final particle: " << partIndex << endl;
          pxMeas[j] += y(yoffset + 3*partIndex + 0 ,0);
          pyMeas[j] += y(yoffset + 3*partIndex + 1 ,0);
          pzMeas[j] += y(yoffset + 3*partIndex + 2 ,0);
          EMeas[j] += energy[Eindex];
          if(_verbose>1) cerr << "constraint4vec: " << pxMeas[j] << " " << pyMeas[j] << " " << pzMeas[j] << " " << EMeas[j] << endl;
        }
        else if(partIndex == -1)
        {
          _extraC_miss[j] = true; // Assume missing particle is not in constraint
        }
      }
      pxTot[j] = pxMeas[j];
      pyTot[j] = pyMeas[j];
      pzTot[j] = pzMeas[j];
      ETot[j] = EMeas[j];
      if(_extraC_miss[j])
      {
        pxTot[j] += x(0,0);
        pyTot[j] += x(1,0);
        pzTot[j] += x(2,0);
        ETot[j] += missingEnergy;
      }
      massTot[j] = sqrt(pow(ETot[j],2) - pow(pxTot[j],2) - pow(pyTot[j],2) - pow(pzTot[j],2));
      if(_verbose>1) cerr << "massTot/constraintMass: " << massTot[j] << " " <<  _constraintMasses[j] << endl;
      c(cindex, 0) = massTot[j] - _constraintMasses[j];
    }

    if(_verbose) 
    {
      for(int k=0;k<numParts;k++) cerr << "energy: " << energy[k] << endl;
      cerr << "missing energy: " << missingEnergy << endl;
    }

    // This is for when I have 4 constraints
    //float gammaE0 = sqrt(y(0,0)*y(0,0) + y(1,0)*y(1,0) + y(2,0)*y(2,0));
    //float gammaE1 = sqrt(y(3,0)*y(3,0) + y(4,0)*y(4,0) + y(5,0)*y(5,0));
    //float pi0E = sqrt(x(0,0)*x(0,0) + x(1,0)*x(1,0) + x(2,0)*x(2,0) + _missingMass * _missingMass);
    //c(0,0) = x(0,0) - y(0,0) - y(3,0);
    //c(1,0) = x(1,0) - y(1,0) - y(4,0);
    //c(2,0) = x(2,0) - y(2,0) - y(5,0);
    //c(3,0) = pi0E - gammaE0 - gammaE1;
    //if(_verbose>1) cerr << "set up the constraints " << endl;

    if(_verbose>1)
    {
      cerr << "c: " << endl;
      for(int j=0;j<c.GetNrows();j++)
      {
        for(int k=0;k<c.GetNcols();k++)
        {
          cerr  << j << " " << k << " " << c(j,k) << " ";
        }
        cerr << endl;
      }
    }


    //This is for 4 constraints
    // Derivatives of constraint equations wrt unknown quantities.
    if(_missingParticle) /// NEED TO DO THIS RIGHT!!!!!!!!!!!!! MB
    {
      for(int j=0;j<3;j++)
      {
        a(j,j) = -1.0;                    // derivatives of p3 constraint eq's 
        a(3,j) = -x(j,0)/missingEnergy; // derivatives of energy constraint eq.
      }
      for(int j=0;j<_extraC;j++)
      {
        int aindex = 4 + j;
        if(_extraC_miss[j])
        {
          a(aindex, 0) = (((x(0,0)*ETot[j])/missingEnergy) - pxTot[j])/massTot[j];
          a(aindex, 1) = (((x(1,0)*ETot[j])/missingEnergy) - pyTot[j])/massTot[j];
          a(aindex, 2) = (((x(2,0)*ETot[j])/missingEnergy) - pzTot[j])/massTot[j];
        }
      }
    }

    if(_verbose>1)
    {
      cerr << "a: " << endl;
      for(int j=0;j<a.GetNrows();j++)
      {
        for(int k=0;k<a.GetNcols();k++)
        {
          cerr << a(j,k) << " ";
        }
        cerr << endl;
      }
    }

    aT.Transpose(a); // set a transpose

    if(_verbose>1)
    {
      cerr << "aT: " << endl;
      for(int j=0;j<aT.GetNrows();j++)
      {
        for(int k=0;k<aT.GetNcols();k++)
        {
          cerr << aT(j,k) << " ";
        }
        cerr << endl;
      }
    }

    // Set up derivatives matrix b(i,j) = dc(i)/dy(j):
    //This is for 4 constraints
    for(int k=0;k<numInitial;k++)
    {
      int offset = k*3;
      b(0,offset+0) = 1.0;
      b(1,offset+1) = 1.0;
      b(2,offset+2) = 1.0;
      for(int j=0;j<3;j++) b(3,offset+j) = y(offset+j,0)/energy[k];
    }
    for(int k=0;k<numFinal;k++)
    {
      int offset = (k+numInitial)*3;
      b(0,offset+0) = -1.0;
      b(1,offset+1) = -1.0;
      b(2,offset+2) = -1.0;
      for(int j=0;j<3;j++) b(3,offset+j) = -y(offset+j,0)/energy[k+numInitial];
    }
    for(int j=0;j<_extraC;j++)
    {
      int bindex = 4 + j;
      int yoffset = 3*numInitial;
      int numConstraintParticles = (int)_constraintParticles[j].size();
      /// Add in the particles to the constraint
      for(int k=0;k<numConstraintParticles;k++)
      {
        int partIndex = _constraintParticles[j][k];
        int Eindex = numInitial + partIndex;
        /// Use a measured particle
        if(partIndex>=0 && partIndex<numFinal)
        {
          b(bindex, yoffset + 3*partIndex +0) = (((y(yoffset + 3*partIndex +0,0)*ETot[j])/energy[Eindex]) - pxTot[j])/massTot[j];
          b(bindex, yoffset + 3*partIndex +1) = (((y(yoffset + 3*partIndex +1,0)*ETot[j])/energy[Eindex]) - pyTot[j])/massTot[j];
          b(bindex, yoffset + 3*partIndex +2) = (((y(yoffset + 3*partIndex +2,0)*ETot[j])/energy[Eindex]) - pzTot[j])/massTot[j];
        }
      }
    }

    if(_verbose>1)
    {
      cerr << "b: " << endl;
      for(int j=0;j<b.GetNrows();j++)
      {
        for(int k=0;k<b.GetNcols();k++)
        {
          cerr << b(j,k) << " ";
        }
        cerr << endl;
      }
    }

    bT.Transpose(b); // set b transpose

    if(_verbose>1)
    {
      cerr << "bT: " << endl;
      for(int j=0;j<bT.GetNrows();j++)
      {
        for(int k=0;k<bT.GetNcols();k++)
        {
          cerr << bT(j,k) << " ";
        }
        cerr << endl;
      }
    }

    // gb is a utility matrix we'll need a lot below
    gb = b * _cov * bT;
    if(_verbose>1)
    {
      cerr << "gb: " << endl;
      for(int j=0;j<gb.GetNrows();j++)
      {
        for(int k=0;k<gb.GetNcols();k++)
        {
          cerr << gb(j,k) << " ";
        }
        cerr << endl;
      }
    }

    if( gb.Determinant() == 0 ) break; // for the case when gb is singular

    gb.Invert();
    if(_verbose>1)
    {
      cerr << "gb inverted: " << endl;
      for(int j=0;j<gb.GetNrows();j++)
      {
        for(int k=0;k<gb.GetNcols();k++)
        {
          cerr << gb(j,k) << " ";
        }
        cerr << endl;
      }

    }

    if(_missingParticle)
    {
      if(_verbose>1)
      {
        cerr << "(aT * gb * a): " << endl;
        for(int j=0;j<(aT * gb * a).GetNrows();j++)
        {
          for(int k=0;k<(aT * gb * a).GetNcols();k++)
          {
            cerr << (aT * gb * a)(j,k) << " ";
          }
          cerr << endl;
        }
        cerr << "Before invert: " << endl;
        cerr << "Determinant: " << (aT * gb * a).Determinant() << endl;
      }
      //(aT * gb * a).Invert();
      if(_verbose>1) cerr << "After invert: " << endl;
      // Is this what we should do? 
      if( (aT * gb * a).Determinant() == 0.0 ) break;

      // calculate the change to be made to the missing particle quantities
      xsi = ((aT * gb * a).Invert())*(aT * gb * c);
      // update x
      x -= xsi;
      // calculate the changes to be made to the measured quantities
      delta = _cov * bT * gb * (c - a*xsi);

      if(_verbose>1)
      {
        cerr << "xsi: " << endl;
        for(int j=0;j<xsi.GetNrows();j++)
        {
          for(int k=0;k<xsi.GetNcols();k++)
          {
            cerr << xsi(j,k) << " ";
          }
          cerr << endl;
        }

        cerr << "(c - a*xsi): " << endl;
        for(int j=0;j<(c - a*xsi).GetNrows();j++)
        {
          for(int k=0;k<(c - a*xsi).GetNcols();k++)
          {
            cerr << (c - a*xsi)(j,k) << " ";
          }
          cerr << endl;
        }

      }

    }
    else 
    {
      delta = _cov * bT * gb * c;
    }

    // update the measured quantities y
    y -= delta;

    if(_verbose>1)
    {
      cerr << "delta: " << endl;
      for(int j=0;j<delta.GetNrows();j++)
      {
        for(int k=0;k<delta.GetNcols();k++)
        {
          cerr << delta(j,k) << " ";
        }
        cerr << endl;
      }
    }

    if(_verbose>1)
    {
      cerr << "new y: " << endl;
      for(int j=0;j<y.GetNrows();j++)
      {
        for(int k=0;k<y.GetNcols();k++)
        {
          cerr << y(j,k) << " ";
        }
        cerr << endl;
      }
    }


    eps = yi - y;
    epsT.Transpose(eps);
    if(TMath::Prob((epsT * bT * gb * b * eps)(0,0),_ndf) < 1.e-10) break;
    if(_verbose>1)
    {
      cerr << "eps: " << endl;
      for(int j=0;j<eps.GetNrows();j++)
      {
        for(int k=0;k<eps.GetNcols();k++)
        {
          cerr << eps(j,k) << " ";
        }
        cerr << endl;
      }
    }


  } /* end iteration process */

  // get the total changes in the measured quantities
  eps = yi - y;
  epsT.Transpose(eps);

  // get the fit covariance matrix
  if(_missingParticle)
  {
    if( (aT * gb * a).Determinant() != 0.0 ) {
      cov_fit = _cov - _cov*bT*gb*b*_cov + _cov*bT*gb*a*((aT*gb*a).Invert())*aT*gb*b*_cov;
    }
    else 
    {
    if(_verbose>1) cerr << "Not going to fill the missing cov matrix output because of singularity...." << endl;
    }
  }
  else 
  {
    cov_fit = _cov - _cov*bT*gb*b*_cov;
  }

  DMatrixDSym dumerr(7);
  //dumerr(0,0) = 1.0;
  //cerr << "FILLING ERROR MATRIX" << endl;
  //cerr << "dum row/col: " << dumerr.GetNrows() << " " << dumerr.GetNcols() << endl;
  //cerr << "kin row/col: " << _kDataInitial_out[0].errorMatrix().GetNrows() << " " << _kDataInitial_out[0].errorMatrix().GetNcols() << endl;
  //_kDataInitial_out[0].setErrorMatrix(dumerr);
  //cerr << "FILLED ERROR MATRIX" << endl;

  for(int i=0;i<numInitial;i++)
  {
    dumerr = _kDataInitial_out[i].errorMatrix();
    for(int j=0;j<numypar;j++)
    {
      for(int k=0;k<numypar;k++)
      {
        dumerr(j,k) = cov_fit(numypar*i+j,numypar*i+k);
      }
    }
    _kDataInitial_out[i].setErrorMatrix(dumerr);
  }
  for(int i=0;i<numFinal;i++)
  {
    dumerr = _kDataFinal_out[i].errorMatrix();
    for(int j=0;j<numypar;j++)
    {
      for(int k=0;k<numypar;k++)
      {
        dumerr(j,k) = cov_fit( numypar*i+j + numypar*numInitial, numypar*i+k + numypar*numInitial);
      }
    }
    _kDataFinal_out[i].setErrorMatrix(dumerr);
  }

  // get chi2 and the pulls
  _pulls.resize(dim);
  for(i = 0; i < dim; i++)
  {
    _pulls[i] = -eps(i,0)/sqrt(_cov(i,i) - cov_fit(i,i));
  }

  _chi2 = (epsT * bT * gb * b * eps)(0,0); // the (0,0)...only...element 

  DRandom rnd;
  rnd.SetSeed((unsigned int)(10000*_chi2));


  if(_verbose>0) cerr << "Prob: " << TMath::Prob(_chi2, _ndf) << " " << _chi2 << " " << _ndf << endl;

  // update output of measured quantities
  if(_verbose>1) cerr << "Filling the 4vecs output...." << endl;
  for(i = 0; i < numInitial; i++) 
  {
    int offset = 3*i;
    _kDataInitial_out[i].setMass(mass[i]);
    _kDataInitial_out[i].setMomentum(DVector3(y(offset+0,0), y(offset+1,0),  y(offset+2,0)));
  }
  for(i = 0; i < numFinal; i++) 
  {
    int offset = 3*(i+numInitial);
    _kDataFinal_out[i].setMass(mass[i+numInitial]);
    _kDataFinal_out[i].setMomentum(DVector3(y(offset+0,0), y(offset+1,0),  y(offset+2,0)));
  }
  if(_verbose>1) cerr << "Filled the 4vecs output...." << endl;

  if(_useTimingInformation)
  {
    //_chi2 = 0;
    //_ndf = 0;
    for(i = 0; i < numFinal; i++) 
    {
      if(_verbose>=1) cerr << "detector: " << i << " " << _kDataFinal_in[i].t1_detector() << endl;
      //if(_kDataFinal_in[i].t1_detector() != SYS_NULL)
      //if(_kDataFinal_in[i].t1_detector() == SYS_TOF)
      //if(_kDataFinal_in[i].t1_detector() == SYS_BCAL)
      if(_kDataFinal_in[i].t1_detector() == SYS_BCAL || _kDataFinal_in[i].t1_detector() == SYS_TOF)
      {
        //cerr << "Before: chi2/_ndf/prob: " << _chi2 << " " << _ndf << " " << TMath::Prob(_chi2, _ndf) << endl;
        double beta = _kDataFinal_out[i].lorentzMomentum().Beta();
        double deltaBeta = _kDataFinal_out[i].deltaBeta();
        double measuredBeta_err = _kDataFinal_out[i].measuredBeta_err();
        double px = _kDataFinal_out[i].lorentzMomentum().X();
        double py = _kDataFinal_out[i].lorentzMomentum().Y();
        double pz = _kDataFinal_out[i].lorentzMomentum().Z();
        double pxerr = sqrt(_kDataFinal_out[i].errorMatrix()(0,0));
        double pyerr = sqrt(_kDataFinal_out[i].errorMatrix()(1,1));
        double pzerr = sqrt(_kDataFinal_out[i].errorMatrix()(2,2));
        double X = beta*beta;
        double A = pow(_kDataFinal_out[i].lorentzMomentum().Rho(),2);
        double B = pow(_kDataFinal_out[i].lorentzMomentum().E(),2);
        double errA = sqrt(pow(px*2.0*pxerr,2) + pow(py*2.0*pyerr,2) + pow(pz*2.0*pzerr,2) );
        double errB = errA;
        double errX = X*sqrt( pow(errA/A,2) + pow(errB/B,2) - 2.0*errA*errB/(A*B));
        double calcBeta_err = beta*(1.0/2.0)*( errX/X );
        double timing_chi2 = pow(deltaBeta,2)/(pow(calcBeta_err,2) + pow(measuredBeta_err,2));
        //double timing_chi2 = pow(deltaBeta,2)/(pow(calcBeta_err,2));
        if(_verbose>=1) cerr << "deltaBeta/err/timing_chi2: " << i << " " << deltaBeta << " " << measuredBeta_err << " " << timing_chi2 << endl;
        if(_verbose>=1) cerr << "Beta/measured: " << i << " " << _kDataFinal_out[i].lorentzMomentum().Beta() << "\t";
        if(_verbose>=1) cerr << _kDataFinal_out[i].measuredBeta() << endl;
        if(_verbose>=1) cerr << "t1_detector: " << i << " " << _kDataFinal_in[i].t1_detector() << "\t";
        if(_verbose>=1) cerr << "chi2: " << _chi2 << endl;
        _chi2 += timing_chi2;
        _ndf += 1;
        //cerr << "chi2/_ndf/prob: " << _chi2 << " " << _ndf << " " << TMath::Prob(_chi2, _ndf) << endl;
      }
    }
  }
  // set the missing particle covariance matrix
  if(_verbose>1) cerr << "Filling the missing cov matrix output...." << endl;
  // Not sure about this part
  if(_missingParticle)
  {
  if( (aT * gb * a).Determinant() != 0.0 ) 
  {
    if(_missingParticle)
      this->_SetMissingParticleErrors((aT * gb * a).Invert(),x);  
    if(_verbose>1) cerr << "Filled the missing cov matrix output...." << endl;
  }
  else 
  {
    if(_verbose>1) cerr << "Not going to fill the missing cov matrix output 'cos singular...." << endl;
  }
  }
}
//_____________________________________________________________________________

void DKinFit::FitTwoGammas(const float __missingMass, const float errmatrixweight=1.0)
{
  /// Private function used internally to perform kinematic fitting.
  //const int numParts = (int)_p4in.size();
  //int dum_dim = 0;
  //for(int i=0;i<(int)_kDataInitial_in.size();i++)
  //{
  //dum_numParts += _kDataInitial_in[i].NumQuantitiesForDKinFit();
  //}
  _missingMass = __missingMass;
  _kDataInitial_out.clear();
  _kDataFinal_out.clear();
  //float errmatrixweight = 1.0;
  const int numInitial = (int)_kDataInitial_in.size();
  const int numFinal = (int)_kDataFinal_in.size();
  const int numParts = numInitial + numFinal;
  //const int dim = 4*numParts; // number of measured quantites
  const int numypar = 3; // number of measured quantites
  const int dim = numypar*numParts; // number of measured quantites
  if(_verbose>0)
  {
    cerr << "numInitial: " << numInitial << endl;
    cerr << "numFinal: " << numFinal << endl;
    cerr << "numParts: " << numParts << endl;
    cerr << "dim: " << dim << endl;
  }
  // Resize the covariance matrix
  _cov.ResizeTo(dim, dim);

  //const int dim = 3*numParts + 1; // number of measured quantites
  // check to see that covariance matrix size is consistent with # particles
  if(_cov.GetNrows() != dim || _cov.GetNcols() != dim)
  {
    // wrong size
    std::cout << "Error! <DKinFit::_MainFitter> Covariance matrix size is NOT "
      << "correct for current number of particles. For " << numParts
      << " particles, the covariance matrix should be " << dim << " x "
      << dim << ", but the matrix passed in is " << _cov.GetNrows()
      << " x " << _cov.GetNcols() << std::endl;
    abort();
  }

  //if(_verbose>1)
  {
    for(int i=0;i<(int)_kDataInitial_in.size();i++)
    {
      _kDataInitial_out.push_back(_kDataInitial_in[i]); /// Make sure we're not using the same pointer
      for(int j=0;j<numypar;j++)
      {
        for(int k=0;k<numypar;k++)
        {
          _cov(numypar*i+j,numypar*i+k) = (_kDataInitial_out[i].errorMatrix())(j,k)/errmatrixweight;
        }
      }
    }
    for(int i=0;i<(int)_kDataFinal_in.size();i++)
    {
      _kDataFinal_out.push_back(_kDataFinal_in[i]); /// Make sure we're not using the same pointer
      for(int j=0;j<numypar;j++)
      {
        for(int k=0;k<numypar;k++)
        {
          _cov(numypar*i+j,numypar*i+k) = (_kDataFinal_out[i].errorMatrix())(j,k)/errmatrixweight;
        }
      }
    }
  }

  if(_verbose>1)
  {
    cerr << "Total error matrix: " << endl;
    for(int j=0;j<_cov.GetNrows();j++)
    {
      for(int k=0;k<_cov.GetNcols();k++)
      {
        cerr << _cov(j,k) << " ";
      }
      cerr << endl;
    }
  }


  // The following was commented out to avoid compiler warnings 7/31/07 D.L.
  //bool isEnergyMeasured[numInitial + numFinal];

  // Get the number of constraint equations
  int numConstraints = 4;
  if(_extraC) numConstraints = 5;
  // For Two gammas
  numConstraints = 4;

  // get the degrees of freedom
  _ndf = 4;
  if(_missingParticle) _ndf -= 3;
  if(_extraC) _ndf += 1;
  // For Two gammas
  _ndf = 1;

  if(_verbose>0) cerr << "numConstraints/ndf: " << numConstraints << " " << _ndf << endl;

  int i;
  // The following were changed to vectors to avoid compiler warnings on SunOS 11/12/2007 D.L.
  vector<double> mass(numParts); //double mass[numParts];
  vector<int> initialOrFinal(numParts); //int initialOrFinal[numParts];
  
  // The following 3 lines were commented out to avoid compiler warnings 7/31/07 D.L.
  //double p[numParts],erg[numParts],px[numParts],py[numParts],pz[numParts];
  //double e_miss=0.,e_inv=0.,px_inv=0.,py_inv=0.,pz_inv=0.;
  //double dumt, dumx, dumy, dumz;
  DVector3 p3_miss;

  for(int i=0;i<numParts;i++)
  {
    if(i<numInitial) initialOrFinal[i] = +1;
    else             initialOrFinal[i] = -1;
  }

  /*********** Define matricies needed to perform kinematic fitting **********/

  // Note: for any matrix m, mT is the transpose of m.  
  //
  DLorentzVector missingMomentum;

  // y(i) ==> measured particle quantities.
  // i = 3*n+0,3*n+1,3*n+=2 are particle n's p_x,p_y,p_z respectively.
  TMatrixD yi(dim,1),y(dim,1); // yi will be intial values
  // x(i) ==> missing particle quantities (0,1,2) = (px,py,pz)
  TMatrixD x(3,1);
  //TMatrixD x(4,1);
  // a(i,j) = dc(i)/dx(j) ==> derivatives of constraint eq's w/r to missing
  // particle quantities x.
  TMatrixD a(numConstraints,3),aT(3,numConstraints);
  //TMatrixD a(numConstraints,4),aT(4,numConstraints);
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
  //TMatrixD xsi(4,1);
  // gb is a utility matrix used 
  TMatrixD gb(numConstraints,numConstraints);
  // cov_fit is the covariance matrix OUTPUT by the fit.
  TMatrixD cov_fit(dim,dim);

  /***************************** Initial Set Up ******************************/

  // set up vector of initial measured values:
  for(i = 0; i < numInitial; i++)
  {
    mass[i] = _kDataInitial_in[i].mass(); // particle's mass
    yi(numypar*i + 0,0) = _kDataInitial_in[i].px();
    yi(numypar*i + 1,0) = _kDataInitial_in[i].py();
    yi(numypar*i + 2,0) = _kDataInitial_in[i].pz();
  }
  for(i = 0; i < numFinal; i++)
  {
    mass[numInitial + i] = _kDataFinal_in[i].mass(); // particle's mass
    yi(numypar*(numInitial+i) + 0,0) = _kDataFinal_in[i].px();
    yi(numypar*(numInitial+i) + 1,0) = _kDataFinal_in[i].py();
    yi(numypar*(numInitial+i) + 2,0) = _kDataFinal_in[i].pz();
  }
  y = yi; // start off y to be yi

  // For Two gammas
  //if(_missingParticle)
  {
    //missingMomentum = this->MissingMomentum(_kDataInitial_in, _kDataFinal_in);
    x(0,0) = yi(0,0) + yi(3,0);
    x(1,0) = yi(1,0) + yi(4,0);
    x(2,0) = yi(2,0) + yi(5,0);
  }

  /************************* Begin Iteration Process *************************/
  for(int iter = 0; iter < 10; iter++)
  { 
    // it should converge by 10 iterations
    // reset all utility matricies
    a.Zero();
    aT.Zero();
    b.Zero();
    bT.Zero();
    c.Zero();
    delta.Zero();
    xsi.Zero();

    // Set up the constraint eq's for this iteration:
    // c(0) (energy constraint): e_f - e_i = 0
    // c(i) (p3 constraint): p3_f - p3_i = 0
    for(int j=0;j<numConstraints;j++) c(j,0) = 0.; 

    // This is for when I have 4 constraints
    if(_verbose>1) cerr << "setting up the constraints " << endl;
    float gammaE0 = sqrt(y(0,0)*y(0,0) + y(1,0)*y(1,0) + y(2,0)*y(2,0));
    float gammaE1 = sqrt(y(3,0)*y(3,0) + y(4,0)*y(4,0) + y(5,0)*y(5,0));
    float pi0E = sqrt(x(0,0)*x(0,0) + x(1,0)*x(1,0) + x(2,0)*x(2,0) + __missingMass * __missingMass);
    c(0,0) = x(0,0) - y(0,0) - y(3,0);
    c(1,0) = x(1,0) - y(1,0) - y(4,0);
    c(2,0) = x(2,0) - y(2,0) - y(5,0);
    c(3,0) = pi0E - gammaE0 - gammaE1;
    if(_verbose>1) cerr << "set up the constraints " << endl;

    if(_verbose>1)
    {
      cerr << "y: " << endl;
      for(int j=0;j<y.GetNrows();j++)
      {
        for(int k=0;k<y.GetNcols();k++)
        {
          cerr  << j << " " << k << " " << y(j,k) << " ";
        }
        cerr << endl;
      }
    }

    if(_verbose>1)
    {
      cerr << "x: " << endl;
      for(int j=0;j<x.GetNrows();j++)
      {
        for(int k=0;k<x.GetNcols();k++)
        {
          cerr  << j << " " << k << " " << x(j,k) << " ";
        }
        cerr << endl;
      }
    }

    if(_verbose>1)
    {
      cerr << "c: " << endl;
      for(int j=0;j<c.GetNrows();j++)
      {
        for(int k=0;k<c.GetNcols();k++)
        {
          cerr  << j << " " << k << " " << c(j,k) << " ";
        }
        cerr << endl;
      }
    }


    //This is for 4 constraints
    // Derivatives of constraint equations wrt unknown quantities.
    // In this case the 4-momentum of the particle which decayed to the 
    // two photons.
    a(0,0) = 1.0;
    a(1,1) = 1.0;
    a(2,2) = 1.0;
    for(i = 0; i < 3; i++)
    {
      a(3,i) = x(i,0)/pi0E; 
    }

    if(_verbose>1)
    {
      cerr << "a: " << endl;
      for(int j=0;j<a.GetNrows();j++)
      {
        for(int k=0;k<a.GetNcols();k++)
        {
          cerr << a(j,k) << " ";
        }
        cerr << endl;
      }
    }

    aT.Transpose(a); // set a transpose

    if(_verbose>1)
    {
      cerr << "aT: " << endl;
      for(int j=0;j<aT.GetNrows();j++)
      {
        for(int k=0;k<aT.GetNcols();k++)
        {
          cerr << aT(j,k) << " ";
        }
        cerr << endl;
      }
    }

    // Set up derivatives matrix b(i,j) = dc(i)/dy(j):
    //This is for 4 constraints
    b(0,0) = -1.0;
    b(1,1) = -1.0;
    b(2,2) = -1.0;
    b(0,3) = -1.0;
    b(1,4) = -1.0;
    b(2,5) = -1.0;
    for(i = 0; i < 3; i++)
    {
      b(3,i) =   -y(i,0)/gammaE0; 
      b(3,i+3) = -y(i+3,0)/gammaE1; 
    }

    if(_verbose>1)
    {
      cerr << "b: " << endl;
      for(int j=0;j<b.GetNrows();j++)
      {
        for(int k=0;k<b.GetNcols();k++)
        {
          cerr << b(j,k) << " ";
        }
        cerr << endl;
      }
    }

    bT.Transpose(b); // set b transpose

    if(_verbose>1)
    {
      cerr << "bT: " << endl;
      for(int j=0;j<bT.GetNrows();j++)
      {
        for(int k=0;k<bT.GetNcols();k++)
        {
          cerr << bT(j,k) << " ";
        }
        cerr << endl;
      }
    }

    // gb is a utility matrix we'll need a lot below
    gb = b * _cov * bT;
    if(_verbose>1)
    {
      cerr << "gb: " << endl;
      for(int j=0;j<gb.GetNrows();j++)
      {
        for(int k=0;k<gb.GetNcols();k++)
        {
          cerr << gb(j,k) << " ";
        }
        cerr << endl;
      }
    }

    if( gb.Determinant() == 0 ) break; // for the case when gb is singular

    gb.Invert();
    if(_verbose>1)
    {
      cerr << "gb inverted: " << endl;
      for(int j=0;j<gb.GetNrows();j++)
      {
        for(int k=0;k<gb.GetNcols();k++)
        {
          cerr << gb(j,k) << " ";
        }
        cerr << endl;
      }

      cerr << "(aT * gb * a): " << endl;
      for(int j=0;j<(aT * gb * a).GetNrows();j++)
      {
        for(int k=0;k<(aT * gb * a).GetNcols();k++)
        {
          cerr << (aT * gb * a)(j,k) << " ";
        }
        cerr << endl;
      }
    }

    // Is this what we should do? 
    if( (aT * gb * a).Determinant() == 0.0 ) break;
    if(_verbose>1)
    {
      cerr << "Before invert: " << endl;
      cerr << "Determinant: " << (aT * gb * a).Determinant() << endl;
      (aT * gb * a).Invert();
      cerr << "After invert: " << endl;
    }

    //if(_missingParticle)
    {
      // calculate the change to be made to the missing particle quantities
      xsi = ((aT * gb * a).Invert())*(aT * gb * c);
      // update x
      x -= xsi;
      // calculate the changes to be made to the measured quantities
      delta = _cov * bT * gb * (c - a*xsi);

      if(_verbose>1)
      {
        cerr << "xsi: " << endl;
        for(int j=0;j<xsi.GetNrows();j++)
        {
          for(int k=0;k<xsi.GetNcols();k++)
          {
            cerr << xsi(j,k) << " ";
          }
          cerr << endl;
        }

        cerr << "(c - a*xsi): " << endl;
        for(int j=0;j<(c - a*xsi).GetNrows();j++)
        {
          for(int k=0;k<(c - a*xsi).GetNcols();k++)
          {
            cerr << (c - a*xsi)(j,k) << " ";
          }
          cerr << endl;
        }

      }

    }
    //else delta = _cov * bT * gb * c;

    // update the measured quantities y
    y -= delta;

    if(_verbose>1)
    {
      cerr << "delta: " << endl;
      for(int j=0;j<delta.GetNrows();j++)
      {
        for(int k=0;k<delta.GetNcols();k++)
        {
          cerr << delta(j,k) << " ";
        }
        cerr << endl;
      }
    }

    if(_verbose>1)
    {
      cerr << "new y: " << endl;
      for(int j=0;j<y.GetNrows();j++)
      {
        for(int k=0;k<y.GetNcols();k++)
        {
          cerr << y(j,k) << " ";
        }
        cerr << endl;
      }
    }


    eps = yi - y;
    epsT.Transpose(eps);
    if(TMath::Prob((epsT * bT * gb * b * eps)(0,0),_ndf) < 1.e-10) break;
    if(_verbose>1)
    {
      cerr << "eps: " << endl;
      for(int j=0;j<eps.GetNrows();j++)
      {
        for(int k=0;k<eps.GetNcols();k++)
        {
          cerr << eps(j,k) << " ";
        }
        cerr << endl;
      }
    }


  } /* end iteration process */

  // get the total changes in the measured quantities
  eps = yi - y;
  epsT.Transpose(eps);

  // get the fit covariance matrix
  //if(_missingParticle)
  if( (aT * gb * a).Determinant() != 0.0 )
  {
    cov_fit = _cov - _cov*bT*gb*b*_cov 
      + _cov*bT*gb*a*((aT*gb*a).Invert())*aT*gb*b*_cov;
  }
  else 
  //cov_fit = _cov - _cov*bT*gb*b*_cov;
  {
     if(_verbose>1) cerr << "Not going to fill the missing cov matrix output 'cos is singular...." << endl;
  }

  // get chi2 and the pulls
  _pulls.resize(dim);
  for(i = 0; i < dim; i++)
  {
    _pulls[i] = -eps(i,0)/sqrt(_cov(i,i) - cov_fit(i,i));
  }

  _chi2 = (epsT * bT * gb * b * eps)(0,0); // the (0,0)...only...element 

  if(_verbose>0) cerr << "Prob: " << TMath::Prob(_chi2, _ndf) << " " << _chi2 << " " << _ndf << endl;

  // update output of measured quantities
  if(_verbose>1) cerr << "Filling the 4vecs output...." << endl;
  for(i = 0; i < (int)_kDataFinal_out.size(); i++) 
  {
    _kDataFinal_out[i].setMass(0.0);
    _kDataFinal_out[i].setMomentum(DVector3(y(3*i+0,0), y(3*i+1,0),  y(3*i+2,0)));
  }
  if(_verbose>1) cerr << "Filled the 4vecs output...." << endl;

  // set the missing particle covariance matrix
  if(_verbose>1) cerr << "Filling the missing cov matrix output...." << endl;
  // Not sure about this part
  if( (aT * gb * a).Determinant() != 0.0 ) 
  {
    if(_missingParticle)
      this->_SetMissingParticleErrors((aT * gb * a).Invert(),x);  
    if(_verbose>1) cerr << "Filled the missing cov matrix output...." << endl;
  }
  else 
  {
    if(_verbose>1) cerr << "Not going to fill the missing cov matrix output 'cos singular...." << endl;
  }
}
//_____________________________________________________________________________

//_____________________________________________________________________________

void DKinFit::_ResetForNewFit(){
  /// Reset the DKinFit object for a new fit (private function).
  // 1st make sure vectors are right size for current p4in
  _kDataInitial_in.clear();
  _kDataFinal_in.clear();
  _kDataInitial_out.clear();
  _kDataFinal_out.clear();
  _pulls.clear();

  // reset bools to default values
  _missingParticle = false;
  _extraC_miss.clear();
}
//_____________________________________________________________________________

double DKinFit::NameToMass(const string &__name) {
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
  std::cout << "Warning! <DKinFit::NameToMass> Unknown particle name: " 
    << __name << ". Mass will be returned as 0." << std::endl;
  return 0.;
}
//_____________________________________________________________________________

void DKinFit::_SetToBadFit()
{
  /// Set fit output to @a bad values.
  /** Sets pulls and \f$ \chi^2 \f$ to 666, sets degrees of freedom and output
   *  kinematic quantites to 0.
   */
  for(int i = 0; i < (int)_pulls.size(); i++) _pulls[i] = 666.;
  _chi2 = -1.;
  _ndf = 0;
  for(int i = 0; i < (int)_kDataInitial_out.size(); i++) 
  {
    _kDataInitial_out[i].setMomentum(DVector3(-666.,-666.,-666.));  
    _kDataInitial_out[i].setMass(-666.);
  }
  for(int i = 0; i < (int)_kDataFinal_out.size(); i++) 
  {
    _kDataFinal_out[i].setMomentum(DVector3(-666.,-666.,-666.));  
    _kDataFinal_out[i].setMass(-666.);
  }
}

//_____________________________________________________________________________
void DKinFit::_SetMissingParticleErrors(const TMatrixD &__missingCov, const TMatrixD &__x){
  /// Calculate the missing particle errors from the last fit

  double p = sqrt(__x(0,0)*__x(0,0) + __x(1,0)*__x(1,0) + __x(2,0)*__x(2,0));
  DLorentzVector p4(__x(0,0),__x(1,0),__x(2,0),sqrt(p*p + pow(_missingMass,2)));

  // kinematic quantities in tracking coordinates
  //double phi = PhiTrack(p4);
  //double lam = LambdaTrack(p4);
  //double alpha = AlphaTrack(p4);
  double phi = 1.0;
  double lam = 1.0;
  double alpha = 1.0;

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
/// Set whether or not there's a missing particle
void DKinFit::SetMissingParticle(const double __missingMass)
{
  _missingParticle = true;
  _missingMass = __missingMass;
}

//_____________________________________________________________________________
/// Set each extra mass constraint
void DKinFit::SetMassConstraint(const double __constraintMass, const vector<int> __constrainedParticles)
{
  _extraC++;
  _constraintMasses.push_back(__constraintMass);
  _constraintParticles.push_back(__constrainedParticles); /// Only final state particles for now
  if(_verbose)
  {
    for(int i=0;i<(int)_constraintMasses.size();i++)
    {
      cerr << "Pushing back constraint mass: " << _constraintMasses[i] << endl;
    }
  }
}

//_____________________________________________________________________________
//
DLorentzVector DKinFit::MissingMomentum(vector<const DKinematicData*> initial, vector<const DKinematicData*> final) 
{
  /// Just a function to return missing momentum
  if(_verbose>1) cerr << "Here in missing momentum...." << endl;
  DLorentzVector ret;
  float E = 0;
  float px = 0;
  float py = 0;
  float pz = 0;
  for(int i = 0; i < (int)initial.size(); i++)
  {
    E += initial[i]->energy();
    px += initial[i]->px();
    py += initial[i]->py();
    pz += initial[i]->pz();
  }
  for(int i = 0; i < (int)final.size(); i++)
  {
    E -= final[i]->energy();
    px -= final[i]->px();
    py -= final[i]->py();
    pz -= final[i]->pz();
  }

  ret.SetXYZT(px, py, pz, E);

  return ret;
}
