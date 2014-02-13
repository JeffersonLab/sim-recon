#ifndef _CHISQMIN_H_
#define _CHISQMIN_H_

#include "residFunc.h"

class chisqMin {
 public:
  chisqMin(residFunc *f_in, int level = 1);
  ~chisqMin();
  void setStartParams(HepVector &x);
  void getParams(HepVector &x);
  double getChi2();
  void fit();
  HepVector getParams();
  HepMatrix getCovar();
  int debug_level;
  int getP();
  int getN();
  int getIter();
 private:
  const gsl_multifit_fdfsolver_type *T; /* pointer to solver type */
  gsl_multifit_fdfsolver *s; /* pointer to solver */
  gsl_multifit_function_fdf f; /* function structure */
  residFunc *residFuncPtr;
  unsigned int iter; /* fit iteration */
  int status; /* status from fit */
  gsl_matrix *covar; /* covariance matrix */
};

#endif // _CHISQMIN_H_
