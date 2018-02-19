// -----------------------------------------
// DrcEvent.h
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef DrcEvent_h
#define DrcEvent_h 1

#include <vector>

#include "TObject.h"
#include "TVector3.h"

#include "DrcHit.h"

class DrcEvent: public TObject  {

public:

  DrcEvent(); 	//the default constructor
  ~DrcEvent() {}; 

  void AddHit(DrcHit hit);
  DrcHit GetHit(Int_t ind) { return fHitArray[ind]; }

  // Accessors 
  Int_t GetId() const { return fId; }
  Int_t GetType() const { return fType; }
  Double_t GetTime() const { return fTime; }

  Int_t GetPdg()          const { return fPdg; }
  Int_t GetParent()       const { return fParent; }
  TVector3 GetMomentum()  const { return fMomentum; }
  TVector3 GetPosition()  const { return fPosition; }
  Int_t GetHitSize()      const { return fHitSize; }
  Double_t GetTest1()     const { return fTest1; }
  Double_t GetTest2()     const { return fTest2; }
  
  // Mutators
  void SetId(Int_t val)        { fId=val; }
  void SetType(Int_t val)        { fType=val; }
  void SetTime(Double_t val)      { fTime=val; }

  void SetPdg(Int_t val) { fPdg = val; }
  void SetParent(Int_t val) { fParent = val; }
  void SetMomentum(TVector3 val) { fMomentum = val; }
  void SetPosition(TVector3 val) { fPosition = val; }
  void SetTest1(Double_t val) { fTest1 = val; }
  void SetTest2(Double_t val) { fTest2 = val; }

private: 
  Int_t fId;
  Int_t fType;
  Int_t fPdg;
  Int_t fParent;
  Double_t fTime;
  
  Int_t fHitSize;
  std::vector<DrcHit> fHitArray;

  TVector3 fMomentum;
  TVector3 fPosition;
  Double_t fTest1;
  Double_t fTest2;
  
  ClassDef(DrcEvent, 1);
};

#endif

// #if defined(__ROOTCLING__)
// #pragma link C++ class DrcEvent+;
// #endif
