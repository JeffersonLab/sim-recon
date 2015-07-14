// DDANAEVIODOMTree.h

// E.Wolin, 19-Mar-2010


#ifndef _DDANAEVIODOMTree_
#define _DDANAEVIODOMTree_


#include "JANA/JObject.h"
#include "evioUtil.hxx"


using namespace jana;
using namespace evio;


//------------------------------------------------------------------------------------


class DDANAEVIODOMTree : public JObject {
  
 public:
  JOBJECT_PUBLIC(DDANAEVIODOMTree);
  
  DDANAEVIODOMTree(uint16_t tag, uint8_t num) : tree(tag,num) {}; 
  ~DDANAEVIODOMTree() {};
   

  // EVIO DOM tree holds event
  evioDOMTree tree;


 private:
  DDANAEVIODOMTree();
};



#endif // _DDANAEVIODOMTree_


//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
