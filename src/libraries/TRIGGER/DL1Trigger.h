#ifndef _DL1Trigger_
#define _DL1Trigger_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DL1Trigger:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DL1Trigger);
  
  DL1Trigger():timestamp(0),trig_mask(0),fp_trig_mask(0),livetime(0),busytime(0),sync_evt(0){}

  
  int event_type;  /* 1 - trigger
		      2 - SYNC and trig */
  
  
  uint64_t timestamp; 
  uint32_t trig_mask;
  uint32_t fp_trig_mask;
  
  
  uint32_t livetime;
  uint32_t busytime;
  
  uint32_t sync_evt;
  
  typedef struct {
    int line;
    int type;
    int fcal;
    int bcal;
    float en_thr;
    int nhit;
    unsigned int pattern;
    int prescale;
  } gtp_par;
  
  vector<gtp_par> gtp_conf; 
  

  // the second argument to AddString is printf style format
  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items,    "timestamp",     "%ld",  timestamp );     
    AddString(items,    "event_type",    "%d",   event_type );
    AddString(items,    "trig_mask",     "%d",   trig_mask );
    AddString(items,    "fp_trig_mask",  "%d",   fp_trig_mask );    
  }
  
};

#endif // _DL1Trigger_

