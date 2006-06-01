#include "TObject.h"

class FCALHitCopy : public TObject {
 public:
	Float_t x;
	Float_t y;
	Float_t E;
	Float_t t;

  FCALHitCopy(){};

  ~FCALHitCopy(){};

  ClassDef(FCALHitCopy,1)
};

ClassImp(FCALHitCopy)
