// -----------------------------------------
// DDIRCTruthBarHit.h
// created on: 05.04.2017
// initial athor: r.dzhygadlo at gsi.de
// -----------------------------------------

#ifndef DDIRCTRUTHBARHIT_H_
#define DDIRCTRUTHBARHIT_H_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DDIRCTruthBarHit: public JObject {

public:
  JOBJECT_PUBLIC (DDIRCTruthBarHit);

  float x, y, z;     //  coordinate where ch. track hits the radiator
  float px, py, pz;  //  components of the track momentum
  float t;	     // time
  float E;	     // energy

  int pdg;           // PDG of the particle
  int bar;           // index of the bar
  int track;         // index of the MC track

  void toStrings(vector<pair<string, string> >&items) const {
    AddString(items, "x", "%1.3f", x);
    AddString(items, "y", "%1.3f", y);
    AddString(items, "z", "%1.3f", z);
    AddString(items, "px", "%1.3f", px);
    AddString(items, "py", "%1.3f", py);
    AddString(items, "pz", "%1.3f", pz);
    AddString(items, "t", "%1.3f", t);
    AddString(items, "E", "%1.3f", E);
    AddString(items, "pdg", "%d", pdg);
    AddString(items, "bar", "%d", bar);
    AddString(items, "track", "%d", track);
  }
};

#endif /* DDIRCTRUTHBARHIT_H_ */
