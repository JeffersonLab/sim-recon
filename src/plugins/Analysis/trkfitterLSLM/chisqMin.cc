#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>

using namespace std;

#include "chisqMin.h"

#define _DBG__ cerr<<__FILE__<<":"<<__LINE__<<endl;

void print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  printf ("print_state: iter: %3u "
          "|f(x)| = %g\n",
          iter,
          gsl_blas_dnrm2 (s->f));
}

chisqMin::chisqMin(residFunc *resid, int level)
  : debug_level(level), residFuncPtr(resid)
{
  unsigned int p = residFuncPtr->getP();
  unsigned int n = residFuncPtr->getN();
  if (n <= p) {
    cout << "number of points less than or equal to number of parameters\n";
    int error = 5;
    throw error;
  }
  covar = gsl_matrix_alloc(p, p); /* allocate the covariance matrix */
  T = gsl_multifit_fdfsolver_lmsder; /* choose the solver */
  s = gsl_multifit_fdfsolver_alloc (T, (size_t)n, (size_t)p); /* allocate work space,
                                                 get solver pointer */
}

chisqMin::~chisqMin() {
  gsl_multifit_fdfsolver_free(s);
  gsl_matrix_free(covar);
}

void chisqMin::setStartParams(HepVector &xIn) {
  if (debug_level >= 2) {
    cout << "setStartParams called\n";
    cout << "start params: " << xIn << endl;
  }
  unsigned int n = residFuncPtr->getN();
  unsigned int p = residFuncPtr->getP();
  if (debug_level >=2) cout << "number of data points = " << n << " number of parameters = " << p << endl;
  f.f = &fGsl;
  f.df = &dfGsl;
  f.fdf = &fdfGsl;
  f.n = (size_t)n;
  f.p = (size_t)p;
  f.params = NULL; /* pointer to something that gets the data,
		      don't know if this is needed yet */
  gsl_vector *xGslPtr;
  xGslPtr = gsl_vector_alloc((size_t)p);
  for (unsigned int i = 0; i < p; i++) {
    gsl_vector_set(xGslPtr, (size_t)i, xIn(i + 1));
    if (debug_level >=2) cout << i << " " << gsl_vector_get(xGslPtr, (size_t)i) << endl; 
  }
  if (debug_level >=2) cout << "initialize the solver" << endl;
  gsl_multifit_fdfsolver_set(s, &f, xGslPtr); /* initialize the solver */
  gsl_vector_free(xGslPtr);
}

void chisqMin::getParams(HepVector& xOut) {
  unsigned int p = residFuncPtr->getP();
  for (unsigned int i = 0; i < p; i++) {
    xOut(i + 1) = gsl_vector_get (s->x, (size_t)i);
  }
}
double chisqMin::getChi2() {
  double chi = gsl_blas_dnrm2(s->f);
  return chi*chi;
}

void chisqMin::fit() {
  if (debug_level > 1) cout << "fit called" << endl;
  iter = 0;
  if (debug_level > 1) print_state(iter, s);
  do {
    iter++;
    if (debug_level >= 2) cout << "== begin iteration " << iter << " ==\n";
    status = gsl_multifit_fdfsolver_iterate(s); /* go to next point */
    if (debug_level >= 2) {
      printf("chisqMin::fit: iter = %d, status = %d, text = %s\n", iter, status, gsl_strerror (status));
      print_state(iter, s);
    }
    if (status) continue;
    if (debug_level >= 2) cout << "chisqMin::fit: continue not taken, do delta test" << endl;
    status = gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4);
    if (debug_level >= 2) {
      printf("chisqMin::fit: from delta test, status = %d, text = %s\n", status, gsl_strerror (status));
    }
  }
  while (status == GSL_CONTINUE && iter < 100);
}

int chisqMin::getP() {
  return residFuncPtr->getP();
}

int chisqMin::getN() {
  return residFuncPtr->getN();
}

HepMatrix chisqMin::getCovar() {
  int nparams = residFuncPtr->getP();
  // allocate memory for the gsl covariance matrix
  gsl_matrix* gslCovarPtr = gsl_matrix_alloc(nparams, nparams);
  // calculate the covariance matrix
  gsl_multifit_covar(s->J, 0.0, gslCovarPtr);
  HepMatrix covar(nparams, nparams);
  // copy the gsl matrix into the hep matrix
  for (int i = 0; i < nparams; i++) {
    for (int j = 0; j < nparams; j++) {
      covar(i + 1, j + 1) = gsl_matrix_get(gslCovarPtr, i, j);
    }
  }
  gsl_matrix_free(gslCovarPtr); // free the gsl matrix memory
  return covar;
}

int chisqMin::getIter() {
  return iter;
}
