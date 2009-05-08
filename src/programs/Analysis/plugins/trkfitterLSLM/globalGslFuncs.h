#ifndef _GLOBALGSLFUNCS_H_
#define _GLOBALGSLFUNCS_H_

residFunc *residFuncPtr; // global pointer to a fitter class for gsl callback

int fGsl(const gsl_vector *x, void *data, gsl_vector *f) {
  //  cout << "fGsl called\n";
  unsigned int n = residFuncPtr->getN(), p = residFuncPtr->getP(), i;
  HepVector xHep(p), fHep(n);
  for (i = 0; i < p; i++) xHep(i + 1) = gsl_vector_get(x, (size_t)i);
  //  cout << "input params:" << xHep << endl;
  residFuncPtr->resid(&xHep, data, &fHep);
  //cout << "residuals:" << fHep;
  for (i = 0; i < n; i++) gsl_vector_set(f, (size_t)i, fHep(i + 1));
  return GSL_SUCCESS;
}

int dfGsl(const gsl_vector *x, void *data, gsl_matrix *J){
  //  cout << "dfGsl called\n";
  int n = residFuncPtr->getN(), p = residFuncPtr->getP();
  int i, j;
  HepVector xHep(p);
  HepMatrix JHep(n,p);
  for (i = 0; i < p; i++) xHep(i + 1) = gsl_vector_get(x, (size_t)i);
  //cout << "input params:" << xHep << endl;
  residFuncPtr->deriv(&xHep, data, &JHep);
  //cout << "derivs:" << JHep;
  for (i = 0; i < n; i++) {
    for (j = 0; j < p; j++) {
      gsl_matrix_set(J, (size_t)i, (size_t)j, JHep(i + 1, j + 1));
    }
  }
  return GSL_SUCCESS;
}

int fdfGsl(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J){
  //  cout << "fdfGsl called\n";
  fGsl(x, data, f);
  dfGsl(x, data, J);
  return GSL_SUCCESS;
}

#endif // _GLOBALGSLFUNCS_H_
