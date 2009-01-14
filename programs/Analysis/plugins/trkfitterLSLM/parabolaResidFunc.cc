#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/Vector.h>

using namespace std;

#include "residFunc.h"
#include "parabolaResidFunc.h"

double xdata[7] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0};
double ydata[7] = {1.1,4.1,8.9,15.7,24.8,36.2,49.2};

void parabolaResidFunc::resid(const HepVector *x, void *data, HepVector *f){
  cout << "resid called\n";
  cout << "params: " << *x << endl;
  double a = (*x)(1);
  double b = (*x)(2);
  double c = (*x)(3);
  double xthis;
  for (unsigned int i = 0; i < n; i++) {
    xthis = xdata[i];
    (*f)[i] = a*xthis*xthis + b*xthis + c - ydata[i];
  }
};

void parabolaResidFunc::deriv(const HepVector *x, void *data, HepMatrix *J){
  cout << "deriv called\n";
  cout << "params: " << *x << endl;
  double xthis, dfda, dfdb, dfdc;
  for (unsigned int i = 0; i < n; i++) {
    xthis = xdata[i];
    dfda = xthis*xthis;
    dfdb = xthis;
    dfdc = 1.0;
    (*J)((int)i + 1, 1) = dfda;
    (*J)((int)i + 1, 2) = dfdb;
    (*J)((int)i + 1, 3) = dfdc;
  }
};

void parabolaResidFunc::residAndDeriv(const HepVector *x, void *data, HepVector *f,
		   HepMatrix *J){};
