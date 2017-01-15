#ifndef _DL1Trigger_
#define _DL1Trigger_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DL1Trigger:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DL1Trigger);
  
  DL1Trigger():event_type(0),timestamp(0),trig_mask(0),fp_trig_mask(0),nsync(0),trig_number(0),live(0),busy(0),live_inst(0),unix_time(0){}

  
  int event_type;  // 0 - unknown , 1 - trigger   2 - SYNC and trig
  
  
  uint64_t timestamp; 
  uint32_t trig_mask;
  uint32_t fp_trig_mask;
  

  uint32_t nsync;
  uint32_t trig_number;
  uint32_t live;
  uint32_t busy;
  uint32_t live_inst;
  uint32_t unix_time;

  vector<uint32_t> gtp_sc;
  vector<uint32_t> fp_sc;
  vector<uint32_t> gtp_rate;
  vector<uint32_t> fp_rate;

 
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
    AddString(items,  "timestamp",     "%ld",  timestamp );     
    AddString(items,  "event_type",    "%d",   event_type );
    AddString(items,  "trig_mask",     "0x%08x",   trig_mask );
    AddString(items,  "fp_trig_mask",  "0x%08x",   fp_trig_mask );    

    AddString(items,  "nsync"       , "%d" , nsync); 
    AddString(items,  "trig_number" , "%d" , trig_number); 
    AddString(items,  "live"        , "%d" , live); 
    AddString(items,  "busy"        , "%d" , busy); 
    AddString(items,  "live_inst"   , "%d" , live_inst); 
    AddString(items,  "unix_time"   , "%d" , unix_time); 
		  
    AddString(items, "gtp_sc"    ,   "%d" ,   gtp_sc.size());
    AddString(items, "fp_sc"     ,   "%d" ,   fp_sc.size());    
    AddString(items, "gtp_rate"   ,  "%d" ,   gtp_rate.size());	    
    AddString(items, "fp_rate"    ,  "%d" ,   fp_rate.size());

  }
  
};

#endif // _DL1Trigger_

