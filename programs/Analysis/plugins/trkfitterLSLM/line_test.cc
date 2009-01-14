#define PI 3.141582654
#include "DLine.h"
#include <CLHEP/Vector/ThreeVector.h>
#include <iostream>
using namespace std;

int main() 
{
  DLine line(0.0, 0.0, 0.0, PI/2.0, 0.0);
  Hep3Vector point(1.0, 1.0, 1.0);
  double doca = line.doca(point);
  cout << "doca = " << doca << endl;
  return 0;
}


