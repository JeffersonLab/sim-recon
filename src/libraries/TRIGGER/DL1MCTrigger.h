#ifndef _DL1MCTrigger_
#define _DL1MCTrigger_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DL1MCTrigger:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DL1MCTrigger);
  
 DL1MCTrigger():trig_mask(0),fcal_en(0),fcal_adc(0),fcal_adc_en(0),fcal_gtp(0),fcal_gtp_en(0),
    bcal_en(0),bcal_adc(0),bcal_adc_en(0),bcal_gtp(0),bcal_gtp_en(0) {}
 
  
  uint32_t trig_mask;
  
  float   fcal_en;
  int     fcal_adc;
  float   fcal_adc_en;
  int     fcal_gtp;
  float   fcal_gtp_en;

  float   bcal_en;
  int     bcal_adc;
  float   bcal_adc_en;
  int     bcal_gtp;
  float   bcal_gtp_en;

  int trig_time[32];


  // the second argument to AddString is printf style format
  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items,  "trig_mask",       "0x%08x", trig_mask );
    
    AddString(items,  "FCAL E(GeV)",     "%6.3f",  fcal_en);
    AddString(items,  "FCAL ADC E(cnt)", "%d",     fcal_adc);
    AddString(items,  "FCAL ADC E(GeV)", "%6.3f",  fcal_adc_en);
    AddString(items,  "FCAL GTP E(cnt)", "%d",     fcal_gtp);
    AddString(items,  "FCAL GTP E(GeV)", "%6.3f",  fcal_gtp_en);
    
    AddString(items,  "BCAL E(GeV)",     "%6.3f",  bcal_en);
    AddString(items,  "BCAL ADC E(cnt)", "%d",     bcal_adc);
    AddString(items,  "BCAL ADC E(GeV)", "%6.3f",  bcal_adc_en);
    AddString(items,  "BCAL GTP E(cnt)", "%d",     bcal_gtp);
    AddString(items,  "BCAL GTP E(GeV)", "%6.3f",  bcal_gtp_en);

    AddString(items,  "Trig Time (samp)", "%d",    trig_time[0]);
    

  }
  
};

#endif // _DL1MCTrigger_

