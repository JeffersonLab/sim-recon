#include <sax/ErrorHandler.hpp>
#include <iostream.h>

class SAXParseException;

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


/*  a simple class for translation between XMLCh data to local coding */

class StrX
{
public :
   StrX(const XMLCh* const toTranscode);
   ~StrX();

   const char* localForm() const;

private :
    char* fLocalForm;		// string in local coding, eg. ASCII
};


#define ELEMENT_NODE 1
#define ATTRIBUTE_NODE 2


inline bool MyOwnErrorHandler::getSawErrors() const
{
   return fSawErrors;
}

inline StrX::StrX(const XMLCh* const toTranscode)
{
   fLocalForm = XMLString::transcode(toTranscode);
}

inline StrX::~StrX()
{
   delete [] fLocalForm;
}

inline const char* StrX::localForm() const
{
   return fLocalForm;
}

inline ostream& operator<<(ostream& target, const StrX& toDump)
{
   target << toDump.localForm();
   return target;
}

