// $Id$
//
//    File: DTAGHTDCDigiHit.h
// Created: Tue Aug  6 13:02:22 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DTAGHTDCDigiHit_
#define _DTAGHTDCDigiHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DTAGHTDCDigiHit: public jana::JObject {
   public:
      JOBJECT_PUBLIC(DTAGHTDCDigiHit);
      
      // Add data members here. For example:
      int counter_id;  ///< counter id 1-274

	  uint32_t time;
      
      // This method is used primarily for pretty printing
      // the second argument to AddString is printf style format
      void toStrings(vector<pair<string,string> > &items)const{
         AddString(items, "counter_id", "%d", counter_id);
         AddString(items, "time", "%d", time);
      }
      
};

#endif // _DTAGHTDCDigiHit_

