// $Id$

#include <iostream>
#include <iomanip>
using namespace std;

#include "DTrackFit.h"



//-----------------------
// Dist
//-----------------------
double DTrackFit::Dist(TVector3 *v)
{
	/// Returns the shortest distance to the track from the point
	/// passed.
	
	// The way this attempts to do this is the following:
	// First, we have to get a point on the track close enough to the hit
	// that we can calculate the distance assuming the track follows a
	// perfect helix. Since we want the minimal distance, we write the
	// distance of the hit to the helix in terms of an angle delta_phi which
	// represents the difference of the phi of the hit itself (measured using
	// a corrdinate system whose origin is on the center of the helix) and
	// the phi of the point on the helix. Differentiating this distance wrt
	// delta_phi and setting equal to zero allows us to solve for delta_phi
	// (sort of). The result is actually a transendental equation:
	//
	//  sin(delta_phi) = k*delta_phi
	//
	// where k = -(dz/dphi)^2/(r*Ro) which comes from the helix definition.
	// We solve this by expanding sin(delta_phi) to 3 terms and reducing it
	// to a quadratic equation. The expansion is good to within 0.5% for
	// delta_phi less than pi/2. This actually should be very good though
	// in most cases since the angular difference from the hit to the track
	// should be very small given the large dz/dphi.
	
	double r, Ro, dzdphi;
	
	double k = -(dzdphi*dzdphi)/r/Ro;

	// coefficients of quadratic equation
	double a = 1.0/120.0;
	double b = -1.0/6.0;
	double c = 1.0 - k;
	double delta_phi_m = (-b - sqrt(b*b - 4.0*a*c))/2.0/a;
	double delta_phi_p = (-b + sqrt(b*b - 4.0*a*c))/2.0/a;
	double delta_phi = fabs(delta_phi_p)<fabs(delta_phi_m) ? delta_phi_p:delta_phi_m;
	
	// calculate distance
	double delta_z = dzdphi*delta_phi;
	double dist = sqrt(r*r + Ro*Ro - 2.0*r*Ro*cos(delta_phi) + delta_z*delta_z);

	return dist;
}
	
