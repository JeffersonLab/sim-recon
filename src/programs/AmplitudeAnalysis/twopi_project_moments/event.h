#ifndef EVENT_H__
#define EVENT_H__

#include <map>

#include "wave.h"

extern double massLow, massHigh;

class event {
public:
  float mass;
  float tPrime;
  float theta;
  float phi;

  // Lookup table.  This contains the Ylm for 0<=l<=4, 0<=m<=l.  They
  // are arranged densely.  It would be smarter to actually know about
  // the waveset in use, then one could only store the waves actually
  // in use, but this is good enough for the numbers of events I'm
  // dealing with.
  float ampls[6];

  event() { mass = tPrime = theta = phi = 0; memset(ampls, 0, sizeof(ampls)); }
  event(double th, double ph) { mass = tPrime = 0; theta = th; phi = ph; initAmpls(); }
  event(double mass_, double tPrime_, double th, double ph)
  { mass = mass_; tPrime = tPrime_; theta = th; phi = ph; initAmpls(); }
  event(const event& o) { mass = o.mass; tPrime = o.tPrime; theta = o.theta; phi = o.phi; memcpy(ampls, o.ampls, sizeof(ampls)); }
    
  void initAmpls();
  double decayAmplitude(int reflectivity, const wave& w) const
  { return decayAmplitude(reflectivity, w.l, w.m); };
  double decayAmplitude(int reflectivity, int l, int m) const;
  double MCweight(int reflectivity, const wave& w1, const wave& w2) const
  { return MCweight(reflectivity, w1.l, w1.m, w2.l, w2.m); }
  double MCweight(int reflectivity, int l1, int m1, int l2, int m2) const;
  std::complex<double> momentWeight(int L, int M) const;

  bool accepted() const {
    return (this->mass >= massLow && this->mass < massHigh
	    //&& this->tPrime > 0.1 && this->tPrime < 0.3
	    //&& this->tPrime > 0.1
	    //&& this->tPrime > 0.1
	    //&& fabs(cos(theta)) < 0.5
            && 1);
  }
};

#endif
