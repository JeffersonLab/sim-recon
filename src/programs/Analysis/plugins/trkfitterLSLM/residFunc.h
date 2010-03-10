#ifndef _RESIDFUNC_H_
#define _RESIDFUNC_H_

using namespace CLHEP;

class residFunc {
 public:
  virtual ~residFunc();
  virtual void resid(const HepVector *x, void *data, HepVector *f) = 0;
  virtual void deriv(const HepVector *x, void *data, HepMatrix *J) = 0;
  virtual void residAndDeriv(const HepVector *x, void *data, HepVector *f, HepMatrix *J) = 0;
  virtual unsigned int getN() = 0;
  virtual unsigned int getP() = 0;
  inline double getChiSquared(){return chiSquared;};
 protected:
  double chiSquared;
};

extern "C" {
#include <gsl/gsl_vector.h> 
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

  int fGsl(const gsl_vector *x, void *data, gsl_vector *f);
  int dfGsl(const gsl_vector *x, void *data, gsl_matrix *J);
  int fdfGsl(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J);
}

#endif // _RESIDFUNC_H_
