// -----------------------------------------
// DDIRCTruthMcpHit.h
// created on: 05.04.2017
// initial athor: r.dzhygadlo at gsi.de
// -----------------------------------------

#ifndef DDIRCMCPHIT_H_
#define DDIRCMCPHIT_H_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DDIRCTruthMcpHit: public JObject {

public:
  JOBJECT_PUBLIC (DDIRCTruthMcpHit);

  float x, y, z;   // hit position
  float t;	   // detection time
  float E;	   // poton energy
  int   ch;        // MCP channel of the hit
  int   key_bar;   // key of the corresponding bar hit
  
  void toStrings(vector<pair<string, string> >&items) const {
    AddString(items, "x", "%1.3f", x);
    AddString(items, "y", "%1.3f", y);
    AddString(items, "z", "%1.3f", z);
    AddString(items, "t", "%1.3f", t);
    AddString(items, "E", "%1.3f", E);
    AddString(items, "ch", "%d", ch);
    AddString(items, "key_bar", "%d", key_bar);
  }
};

#endif /* DDIRCMCPHIT_H_ */
