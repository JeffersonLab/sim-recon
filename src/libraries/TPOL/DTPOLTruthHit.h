#ifndef _DTPOLTruthHit_
#define _DTPOLTruthHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DTPOLTruthHit:public JObject{
public:
    JOBJECT_PUBLIC(DTPOLTruthHit);

    float dEdx;
    bool primary;
    int track;
    int itrack;
    int ptype;
    float r;
    float phi;
    float z;
    float t;
    int sector;

    void toStrings(vector<pair<string,string> > &items)const{
        AddString(items, "track", "%d", track);
        AddString(items, "itrack", "%d", itrack);
        AddString(items, "primary", "%d", primary);
        AddString(items, "ptype", "%d", ptype);
        AddString(items, "dEdx(MeV/cm)", "%1.3f", dEdx*1.0E3);
        AddString(items, "t", "%3.2f", t);
        AddString(items, "r", "%3.1f", r);
        AddString(items, "phi", "%1.3f", phi);
        AddString(items, "z", "%3.1f", z);
        AddString(items, "sector", "%d", sector);
    }
};

#endif // _DTPOLTruthHit_

