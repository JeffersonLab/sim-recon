#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <iostream>
using namespace std;

XERCES_CPP_NAMESPACE_USE

#include "XString.hpp"
#include "XParsers.hpp"

class makeTargetTable
{
public :
   makeTargetTable(const int capacity);
   ~makeTargetTable();

   int add(DOMElement* const targetEl,
           DOMElement* const parentEl = 0);
   DOMElement* lookup(const DOMElement* const targetEl,
                      const XString& type, const XString& name);

private :
   int fTableLen;
   int fTableSize;
   DOMElement** fTargetEl;
   DOMElement** fParentEl;

   int addElement(DOMElement* const targetEl);
   int addComposite(DOMElement* const targetEl);
   int addModel(DOMElement* const targetEl);

   void set_name(DOMElement* const targetEl);
   float set_a(DOMElement* const targetEl);
   float set_z(DOMElement* const targetEl);
   float set_density(DOMElement* const targetEl);
   float set_radlen(DOMElement* const targetEl);
   float set_collen(DOMElement* const targetEl);
   float set_abslen(DOMElement* const targetEl);
   float set_dedx(DOMElement* const targetEl);

   void dump(const DOMElement* const targetEl);
   void dump(DOMElement* const targetEl);
};
