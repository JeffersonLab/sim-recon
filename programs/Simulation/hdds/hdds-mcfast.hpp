#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <iostream>
using namespace std;

XERCES_CPP_NAMESPACE_USE

/* a simple error handler deriviative to install on parser */

class MyOwnErrorHandler : public ErrorHandler
{
public:
   MyOwnErrorHandler();
   ~MyOwnErrorHandler();

   bool getSawErrors() const;

/* Implementation of the SAX ErrorHandler interface */
   void warning(const SAXParseException& e);
   void error(const SAXParseException& e);
   void fatalError(const SAXParseException& e);
   void resetErrors();

private :
   MyOwnErrorHandler(const MyOwnErrorHandler&);
   void operator=(const MyOwnErrorHandler&);

   bool    fSawErrors;     // flag to record that an error occurred
};

class MyDOMErrorHandler : public DOMErrorHandler
{
public:
   MyDOMErrorHandler();
   ~MyDOMErrorHandler();

   bool getSawErrors() const;
   bool handleError(const DOMError& domError);
   void resetErrors();

private :
    MyDOMErrorHandler(const MyDOMErrorHandler&);
    void operator=(const MyDOMErrorHandler&);
    bool fSawErrors;
};

inline bool MyDOMErrorHandler::getSawErrors() const
{
       return fSawErrors;
}

/*  a simple class for translation between XMLCh strings and local coding */

class XString
{
public :
   XString(void);
   XString(const XMLCh* const x);
   XString(const char* const s);
   XString(const XString& X);
   ~XString();

   const char* localForm() const;
   const XMLCh* unicodeForm() const;
   bool equals(const XString& X) const;
   bool equals(const char* const s) const;
   bool equals(const XMLCh* const x) const;
   int stringLen() const;
   bool operator==(const int len) const;
   bool operator!=(const int len) const;
   XString& operator=(const XString& X);
   XString& operator+=(const XString& X);

private :
    XMLCh* fUnicodeForm;	// string in XMLCh coding, eg. Unicode
    char* fLocalForm;		// string in local coding, eg. ASCII
};

/* class for holding a parameter list */

class MCfastParameterList
{
public:
   MCfastParameterList(const int capacity);
   ~MCfastParameterList();

   void inherits(MCfastParameterList& anc);
   void append(const DOMElement* par);
   const DOMElement* get(char* name, char* type,
                          int number=0, char* unit=0) const;

private:
   MCfastParameterList* fAncestor;	// inheritance pointer, if any
   const DOMElement** fListEl;		// array of DOM_Element contents
   int fListLen;			// current parameter count
   int fListCap;			// allocation size
};


inline ostream& operator<<(ostream& target, const XString& toDump)
{
   target << toDump.localForm();
   return target;
}

inline bool MyOwnErrorHandler::getSawErrors() const
{
   return fSawErrors;
}

DOMDocument* parseInputDocument(const char* file);
DOMDocument* buildDOMDocument(const char* file);
