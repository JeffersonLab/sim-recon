#ifndef _DTPOLHit_
#define _DTPOLHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DTPOLHit:public jana::JObject{
public:
    JOBJECT_PUBLIC(DTPOLHit);

    int sector;    // sector number 1-32
    double phi;
    int ring;      // ring number 1-24
    double theta;
    double pulse_peak;
    double integral;
    double dE;     // Energy loss in keV
    double t;

    void toStrings(vector<pair<string,string> > &items)const{
        AddString(items, "sector", "%d", sector);
        AddString(items, "phi", "%3.3f", phi);
        AddString(items, "ring", "%d", ring);
        AddString(items, "theta", "%3.3f", theta);
        AddString(items, "pulse_peak", "%3.3f", pulse_peak);
        AddString(items, "integral", "%3.3f", integral);
        AddString(items, "dE", "%3.3f", dE);
        AddString(items, "t", "%3.3f", t);
    }
};

#endif // _DTPOLHit_
