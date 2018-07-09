#include "DrcEvent.h"

ClassImp(DrcEvent)

// // -----   Default constructor   -------------------------------------------
DrcEvent::DrcEvent(): fId(-1),fType(0),fPdg(0),fParent(0),
  fTime(-1),fHitSize(0),fTest1(0),fTest2(0)
{ 
}

void DrcEvent::AddHit(DrcHit hit){
  fHitArray.push_back(hit);
  fHitSize++;
}
