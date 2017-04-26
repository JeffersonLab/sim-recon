// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "DrcEvent.h"
#include "DrcHit.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_DrcHit(void *p = 0);
   static void *newArray_DrcHit(Long_t size, void *p);
   static void delete_DrcHit(void *p);
   static void deleteArray_DrcHit(void *p);
   static void destruct_DrcHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DrcHit*)
   {
      ::DrcHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::DrcHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("DrcHit", ::DrcHit::Class_Version(), "DrcHit.h", 15,
                  typeid(::DrcHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::DrcHit::Dictionary, isa_proxy, 4,
                  sizeof(::DrcHit) );
      instance.SetNew(&new_DrcHit);
      instance.SetNewArray(&newArray_DrcHit);
      instance.SetDelete(&delete_DrcHit);
      instance.SetDeleteArray(&deleteArray_DrcHit);
      instance.SetDestructor(&destruct_DrcHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DrcHit*)
   {
      return GenerateInitInstanceLocal((::DrcHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::DrcHit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_DrcEvent(void *p = 0);
   static void *newArray_DrcEvent(Long_t size, void *p);
   static void delete_DrcEvent(void *p);
   static void deleteArray_DrcEvent(void *p);
   static void destruct_DrcEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DrcEvent*)
   {
      ::DrcEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::DrcEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("DrcEvent", ::DrcEvent::Class_Version(), "DrcEvent.h", 17,
                  typeid(::DrcEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::DrcEvent::Dictionary, isa_proxy, 4,
                  sizeof(::DrcEvent) );
      instance.SetNew(&new_DrcEvent);
      instance.SetNewArray(&newArray_DrcEvent);
      instance.SetDelete(&delete_DrcEvent);
      instance.SetDeleteArray(&deleteArray_DrcEvent);
      instance.SetDestructor(&destruct_DrcEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DrcEvent*)
   {
      return GenerateInitInstanceLocal((::DrcEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::DrcEvent*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr DrcHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *DrcHit::Class_Name()
{
   return "DrcHit";
}

//______________________________________________________________________________
const char *DrcHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DrcHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int DrcHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DrcHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *DrcHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DrcHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *DrcHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DrcHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr DrcEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *DrcEvent::Class_Name()
{
   return "DrcEvent";
}

//______________________________________________________________________________
const char *DrcEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DrcEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int DrcEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DrcEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *DrcEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DrcEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *DrcEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DrcEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void DrcHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class DrcHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(DrcHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(DrcHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_DrcHit(void *p) {
      return  p ? new(p) ::DrcHit : new ::DrcHit;
   }
   static void *newArray_DrcHit(Long_t nElements, void *p) {
      return p ? new(p) ::DrcHit[nElements] : new ::DrcHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_DrcHit(void *p) {
      delete ((::DrcHit*)p);
   }
   static void deleteArray_DrcHit(void *p) {
      delete [] ((::DrcHit*)p);
   }
   static void destruct_DrcHit(void *p) {
      typedef ::DrcHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::DrcHit

//______________________________________________________________________________
void DrcEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class DrcEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(DrcEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(DrcEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_DrcEvent(void *p) {
      return  p ? new(p) ::DrcEvent : new ::DrcEvent;
   }
   static void *newArray_DrcEvent(Long_t nElements, void *p) {
      return p ? new(p) ::DrcEvent[nElements] : new ::DrcEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_DrcEvent(void *p) {
      delete ((::DrcEvent*)p);
   }
   static void deleteArray_DrcEvent(void *p) {
      delete [] ((::DrcEvent*)p);
   }
   static void destruct_DrcEvent(void *p) {
      typedef ::DrcEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::DrcEvent

namespace ROOT {
   static TClass *vectorlEDrcHitgR_Dictionary();
   static void vectorlEDrcHitgR_TClassManip(TClass*);
   static void *new_vectorlEDrcHitgR(void *p = 0);
   static void *newArray_vectorlEDrcHitgR(Long_t size, void *p);
   static void delete_vectorlEDrcHitgR(void *p);
   static void deleteArray_vectorlEDrcHitgR(void *p);
   static void destruct_vectorlEDrcHitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<DrcHit>*)
   {
      vector<DrcHit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<DrcHit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<DrcHit>", -2, "vector", 214,
                  typeid(vector<DrcHit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEDrcHitgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<DrcHit>) );
      instance.SetNew(&new_vectorlEDrcHitgR);
      instance.SetNewArray(&newArray_vectorlEDrcHitgR);
      instance.SetDelete(&delete_vectorlEDrcHitgR);
      instance.SetDeleteArray(&deleteArray_vectorlEDrcHitgR);
      instance.SetDestructor(&destruct_vectorlEDrcHitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<DrcHit> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<DrcHit>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEDrcHitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<DrcHit>*)0x0)->GetClass();
      vectorlEDrcHitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEDrcHitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEDrcHitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<DrcHit> : new vector<DrcHit>;
   }
   static void *newArray_vectorlEDrcHitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<DrcHit>[nElements] : new vector<DrcHit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEDrcHitgR(void *p) {
      delete ((vector<DrcHit>*)p);
   }
   static void deleteArray_vectorlEDrcHitgR(void *p) {
      delete [] ((vector<DrcHit>*)p);
   }
   static void destruct_vectorlEDrcHitgR(void *p) {
      typedef vector<DrcHit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<DrcHit>

namespace {
  void TriggerDictionaryInitialization_Dict_Impl() {
    static const char* headers[] = {
"DrcEvent.h",
"DrcHit.h",
0
    };
    static const char* includePaths[] = {
"/data.local1/install/root-6.08.00/include",
"/data.local1/dirc/halldgluex/sim-recon/master/src/plugins/Analysis/pid_dirc/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$DrcEvent.h")))  DrcHit;
class __attribute__((annotate("$clingAutoload$DrcEvent.h")))  DrcEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "DrcEvent.h"
#include "DrcHit.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"DrcEvent", payloadCode, "@",
"DrcHit", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_Dict() {
  TriggerDictionaryInitialization_Dict_Impl();
}
