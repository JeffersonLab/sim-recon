int debug_level = 1;

#include <vector>
#include <CLHEP/Matrix/Vector.h>

using namespace std;

#define PI 3.1415926
#include "MyTrajectoryGrkuta.h"

int main() 
{
  double p = 0.5 , theta = PI/6.0, phi = PI/4.0;
  HepVector B(3);
  B(1) = 0.0;
  B(2) = 0.0;
  B(3) = 4.0;
  double xp0 = 0.0;
  double z0 = 0.0;
  double ptinv = 1.0/(p*sin(theta));
#ifdef GRKUTA
  MyTrajectoryGrkuta trajectory(B, debug_level);
#else
  MyTrajectoryBfield trajectory(B, debug_level);
#endif
  trajectory.swim(xp0, z0, theta, phi, ptinv);
  trajectory.print();
  return 0;
}
