// $Id$
//
//    File: DTAGMTDCDigiHit.h
// Created: Tue Aug  6 13:02:22 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DTAGMTDCDigiHit_
#define _DTAGMTDCDigiHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DTAGMTDCDigiHit: public jana::JObject {
   public:
      JOBJECT_PUBLIC(DTAGMTDCDigiHit);
      
      // Add data members here. For example:
      int row;         ///< row number 1-5
      int column;      ///< column number 1-102

	  uint32_t time;
      
      // This method is used primarily for pretty printing
      // the second argument to AddString is printf style format
      void toStrings(vector<pair<string,string> > &items)const{
         AddString(items, "row", "%d", row);
         AddString(items, "column", "%d", column);
         AddString(items, "time", "%d", time);
      }
      
};

#endif // _DTAGMTDCDigiHit_
