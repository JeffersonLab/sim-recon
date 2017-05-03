// -----------------------------------------
// DrcHit.h
// created on: 07.04.2017
// initial author: r.dzhygadlo at gsi.de
// -----------------------------------------

#ifndef DrcHit_h
#define DrcHit_h 1

#include <vector>

#include "TObject.h"
#include "TVector3.h"

class DrcHit : public TObject {

public:   
 
  //Constructor
  DrcHit();

  ~DrcHit(){};
 
  // Accessors 
  Int_t GetType()        { return fType; }
  Int_t GetNreflectionsInPrizm()  { return fNreflectionsInPrizm; }
  Double_t GetPathInPrizm()  { return fPathInPrizm; }
  TVector3 GetMomentum()     { return fMomentum; }
  TVector3 GetPosition()     { return fPosition; }
  Double_t GetCherenkovMC()  { return fCherenkovMC;}
  
  Int_t GetPmtId()       { return fPmtId; }
  Int_t GetPixelId()     { return fPixelId; }
  Int_t GetChannel() { return fChannel;}
  Double_t GetLeadTime() { return fLeadTime; } 
    
  // Mutators
  void SetType(Int_t val)    { fType = val; }
  void SetNreflectionsInPrizm(Int_t val)  { fNreflectionsInPrizm = val; }
  void SetPathInPrizm(Double_t val) { fPathInPrizm = val; }
  void SetMomentum(TVector3 val)    { fMomentum = val; }
  void SetPosition(TVector3 val)    { fPosition = val; }
  void SetCherenkovMC(Double_t val) { fCherenkovMC = val; }
  
  void SetPmtId(Int_t val)   { fPmtId = val; }
  void SetPixelId(Int_t val) { fPixelId = val; }
  void SetChannel(Int_t val) { fChannel=val; }
  void SetLeadTime(Double_t val) { fLeadTime=val; } 

protected:

  Int_t fType;
  Int_t fNreflectionsInPrizm;
  Double_t fPathInPrizm;
  TVector3 fMomentum;
  TVector3 fPosition;
  Double_t fCherenkovMC;
  
  Int_t fPmtId;
  Int_t fPixelId;
  Int_t fChannel;
  Double_t fLeadTime;    

  ClassDef(DrcHit,1)
};

#endif

// #if defined(__ROOTCLING__)
// #pragma link C++ class DrcHit+;
// #endif
