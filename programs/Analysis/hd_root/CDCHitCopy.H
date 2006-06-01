#include "TObject.h"

class CDCHitCopy : public TObject {
 public:
  Float_t  radius;
  Float_t  phim;
  Float_t  dE;
  Float_t  t;
  Float_t  x;
  Float_t  y;  

  CDCHitCopy(){};

  ~CDCHitCopy(){};

  ClassDef(CDCHitCopy,1)
};

ClassImp(CDCHitCopy)
