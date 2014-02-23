/*
 * XParsers: service classes to support parsing of xml domuments
 *           using standard DOM parsing tools
 *
 * Class implementation
 * September 21, 2003
 * Richard Jones
 *
 * Bundled with the error handler classes are two utility functions:
 *
 * 1. parseInputDocument - parser implemented using the old-style
 *    XercesDOMParser interface based on the example code in
 *    $XERCESCROOT/samples/DOMPrint
 *
 * 2. buildDOMDocument - parser implemented using the w3c standard
 *    DOMBuilder interface based on the example code in
 *    $XERCESCROOT/samples/DOMCount
 *
 * Implementation Notes:
 * ---------------------
 * To prevent memory leaks, each of these parsers only retains a single
 * document in memory at a time.  The next call will destroy the DOM
 * tree created on the previous call and return the resources to the
 * pool.  To prevent this behavior, call the parser with the argument
 * perm=true, in which case the resulting DOMDocument will persist for
 * the rest of the lifetime of the program.
 *
 *
 * Modification Notes:
 * --------------------
 * 11/7/2012 DL
 *   Added EntityResolver class to keep track of all of the XML files
 *   pulled in by the parser so an md5 checksum could be performed.
 *   results are written to a FORTRAN function called "md5geom" so the
 *   checksum can be accessed programatically.
 *
 * 6/12/2012  DL
 *   Xerces 3 has done away with the DOMBuilder API, yet retains
 *   the DOMParser. It seems the code using the routines in this file
 *   looked to the pre-processor variable OLD_STYLE_XERCES_PARSER to
 *   decide whether to call parseInputDocument() or buildDOMDocument().
 *   The former being called if the variable was defined implying
 *   the former was likely to be deprecated. The simplest change that
 *   could be made to get this working with XERCES 3 was to turn the
 *   buildDOMDocument() routine into a wrapper for the parseInputDocument()
 *   routine. This is done below.
 *
 */

#include <fstream>
using namespace std;

#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>

#include "XParsers.hpp"
#include "XString.hpp"
#include "md5.h"

std::string last_md5_checksum = "";

/*
 * FIX_XERCES_getElementById_BUG does a store/load cycle at parsing time
 * to fully instantiate entity references on the document tree.
 * See xerces-c++ bug 12800 at http://nagoya.apache.org
 */
#define FIX_XERCES_getElementById_BUG true

#define X(str) XString(str).unicode_str()
#define S(str) str.c_str()

xercesc::DOMDocument* parseInputDocument(const XString& xmlFile, bool keep)
{
   static xercesc::XercesDOMParser* scratchParser=0;
   xercesc::XercesDOMParser* parser;
   if (keep)
   {
      parser = new xercesc::XercesDOMParser;
   }
   else if (scratchParser == 0)
   {
      parser = scratchParser = new xercesc::XercesDOMParser;
   }
   else
   {
      parser = scratchParser;
   }
   
   MyEntityResolver myEntityResolver(xmlFile);
   
   parser->setValidationScheme(xercesc::XercesDOMParser::Val_Auto);
   parser->setCreateEntityReferenceNodes(false);
   parser->setValidationSchemaFullChecking(true);
   parser->setDoNamespaces(true);
   parser->setDoSchema(true);
   parser->setEntityResolver(&myEntityResolver);

   MyOwnErrorHandler errorHandler;
   parser->setErrorHandler(&errorHandler);

   try
   {
      parser->parse(xmlFile.c_str());
	  myEntityResolver.GetMD5_checksum();
   }
   catch (const xercesc::XMLException& toCatch)
   {
      std::cerr
           << "\nparseInputDocument: Error during parsing: '" << xmlFile
	   << "'\n" << "Exception message is:  \n"
           << toCatch.getMessage() << "\n" << std::endl;
      return 0;
   }
   catch (const xercesc::DOMException& toCatch)
   {
      std::cerr
           << "\nXParsers: Error during parsing: '" << xmlFile << "'\n"
           << "Exception message is:  \n"
           << toCatch.msg << "\n" << std::endl;
      xercesc::XMLPlatformUtils::Terminate();
      return 0;
   }
   catch (...)
   {
      std::cerr
           << "\nparseInputDocument: Unexpected exception during parsing: '"
           << xmlFile << "'\n";
      xercesc::XMLPlatformUtils::Terminate();
      return 0;
   }

   if (errorHandler.getSawErrors())
   {
      std::cerr << "\nErrors occured, no output available\n" << std::endl;
      return 0;
   }

   return parser->getDocument();
}

xercesc::DOMDocument* buildDOMDocument(const XString& xmlFile, bool keep)
{
return parseInputDocument(xmlFile, keep);
#if 0 // below no longer works in XERCES 3
	
   xercesc::DOMImplementation *impl =
         xercesc:: DOMImplementationRegistry::getDOMImplementation(X("LS"));
   static xercesc::DOMBuilder* scratchBuilder=0;
   xercesc::DOMBuilder* builder;
   if (keep)
   {
      builder = ((xercesc::DOMImplementationLS*)impl)->createDOMBuilder(
                  xercesc::DOMImplementationLS::MODE_SYNCHRONOUS, 0);
   }
   else if (scratchBuilder == 0)
   {
      builder = scratchBuilder = ((xercesc::DOMImplementationLS*)impl)->
	        createDOMBuilder(xercesc::DOMImplementationLS::MODE_SYNCHRONOUS,
                 0);
   }
   else
   {
      builder = scratchBuilder;
   }
   XString tmpFileS = ".tmp-"+xmlFile.basename();

   builder->setFeature(xercesc::XMLUni::fgDOMValidation, true);
   builder->setFeature(xercesc::XMLUni::fgDOMNamespaces, true);
   builder->setFeature(xercesc::XMLUni::fgDOMDatatypeNormalization, true);
   builder->setFeature(xercesc::XMLUni::fgDOMEntities, false);
   builder->setFeature(xercesc::XMLUni::fgXercesSchemaFullChecking, true);
   builder->setFeature(xercesc::XMLUni::fgXercesSchema, true);

   MyDOMErrorHandler errHandler;
   builder->setErrorHandler(&errHandler);

   xercesc::DOMDocument* doc = 0;

   try {
      builder->resetDocumentPool();
      doc = builder->parseURI(xmlFile.c_str());
#if defined FIX_XERCES_getElementById_BUG
      xercesc::DOMWriter* writer = ((xercesc::DOMImplementationLS*)impl)->
                                    createDOMWriter();
      xercesc::LocalFileFormatTarget* lfft =
                     new xercesc::LocalFileFormatTarget(X(tmpFileS));
      writer->writeNode(lfft,*(doc->getDocumentElement()));
      delete lfft;
      delete writer;
      builder->resetDocumentPool();
      doc = builder->parseURI(X(tmpFileS));
#endif
   }
   catch (const xercesc::XMLException& toCatch) {
      std::cout << "Exception message is: \n" << toCatch.getMessage() << "\n";
      return 0;
   }
   catch (const xercesc::DOMException& toCatch) {
      std::cout << "Exception message is: \n" << toCatch.msg << "\n";
      return 0;
   }
   catch (...) {
      std::cout << "Unexpected Exception \n" ;
      return 0;
   }

   if (errHandler.getSawErrors())
   {
      std::cerr << "\nErrors occured, no output available\n" << std::endl;
      return 0;
   }

   return doc;
#endif // 0
}

MyOwnErrorHandler::MyOwnErrorHandler() : 
   fSawErrors(false)
{
}

MyOwnErrorHandler::~MyOwnErrorHandler()
{
}

// Overrides of the SAX ErrorHandler interface

void MyOwnErrorHandler::error(const xercesc::SAXParseException& e)
{
   fSawErrors = true;
   XString systemId(e.getSystemId());
   XString message(e.getMessage());
   std::cerr
        << "\nparseInputDocument: Error at file " << S(systemId)
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << S(message) << std::endl;
}

void MyOwnErrorHandler::fatalError(const xercesc::SAXParseException& e)
{
   fSawErrors = true;
   XString systemId(e.getSystemId());
   XString message(e.getMessage());
   std::cerr
        << "\nparseInputDocument: Fatal Error at file " << S(systemId)
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << S(message) << std::endl;
}

void MyOwnErrorHandler::warning(const xercesc::SAXParseException& e)
{
   XString systemId(e.getSystemId());
   XString message(e.getMessage());
   std::cerr
        << "\nparseInputDocument: Warning at file " << S(systemId)
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << S(message) << std::endl;
}

void MyOwnErrorHandler::resetErrors()
{
}

MyDOMErrorHandler::MyDOMErrorHandler() :

    fSawErrors(false)
{
}

MyDOMErrorHandler::~MyDOMErrorHandler()
{
}

//  MyDOMHandlers: Overrides of the DOM ErrorHandler interface

bool MyDOMErrorHandler::handleError(const xercesc::DOMError& domError)
{
   fSawErrors = true;
   if (domError.getSeverity() == xercesc::DOMError::DOM_SEVERITY_WARNING)
      std::cerr << "\nWarning at file ";
   else if (domError.getSeverity() == xercesc::DOMError::DOM_SEVERITY_ERROR)
      std::cerr << "\nError at file ";
   else
      std::cerr << "\nFatal Error at file ";

   std::cerr
        << XString(domError.getLocation()->getURI()).c_str()
        << ", line " << domError.getLocation()->getLineNumber()
        << ", char " << domError.getLocation()->getColumnNumber()
        << "\n  Message: " << XString(domError.getMessage()).c_str()
       	<< std::endl;

   return true;
}

void MyDOMErrorHandler::resetErrors()
{
   fSawErrors = false;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

//----------------------------------
// MyEntityResolver (constructor)
//----------------------------------
MyEntityResolver::MyEntityResolver(const XString& xmlFile)
{
	xml_filenames.push_back(xmlFile);
	
	string fname = xmlFile;
	size_t pos = fname.find_last_of('/');
	if(pos != string::npos){
		path = fname.substr(0,pos) + "/";
	}
}

//----------------------------------
// MyEntityResolver (destructor)
//----------------------------------
MyEntityResolver::~MyEntityResolver()
{

}

//----------------------------------
// resolveEntity
//----------------------------------
xercesc::InputSource* MyEntityResolver::resolveEntity(const XMLCh* const publicId, const XMLCh* const systemId)
{
	/// This method gets called from the xerces parser each time it
	/// opens a file (except for the top-level file). For each of these,
	/// record the name of the file being opened, then just return NULL
	/// to have xerces handle opening the file in the normal way.

	// Do some backflips to get strings into std::string format
	std::string my_publicId = "";
	std::string my_systemId = "";
	if(publicId){
		char *my_publicId_ptr = xercesc::XMLString::transcode(publicId);
		my_publicId = my_publicId_ptr;
		xercesc::XMLString::release(&my_publicId_ptr);
	}
	if(systemId){
		char *my_systemId_ptr = xercesc::XMLString::transcode(systemId);
		my_systemId = my_systemId_ptr;
		xercesc::XMLString::release(&my_systemId_ptr);
	}
	//std::cerr<<"publicId="<<my_publicId<<"  systemId="<<my_systemId<<std::endl;

	// The systemId seems to be the one we want
	xml_filenames.push_back(path + my_systemId);

	return NULL; // have xerces handle this using its defaults
}

//----------------------------------
// GetXMLFilenames
//----------------------------------
std::vector<std::string> MyEntityResolver::GetXMLFilenames(void)
{
	return xml_filenames;
}

//----------------------------------
// GetMD5_checksum
//----------------------------------
std::string MyEntityResolver::GetMD5_checksum(void)
{
	/// This will calculate an MD5 checksum using all of the files currently
	/// in the list of XML files. To do this, it opens each file and reads it
	/// in, in its entirety, updating the checksum as it goes. The checksum is
	/// returned as a hexadecimal string.

	md5_state_t pms;
	md5_init(&pms);
	for(unsigned int i=0; i<xml_filenames.size(); i++){

		//std::cerr<<".... Adding file to MD5 checksum : " << xml_filenames[i] << std::endl;
	
		ifstream ifs(xml_filenames[i].c_str());
		if(!ifs.is_open())continue;

		// get length of file:
		ifs.seekg (0, ios::end);
		unsigned int length = ifs.tellg();
		ifs.seekg (0, ios::beg);

		// allocate memory:
		char *buff = new char [length];

		// read data as a block:
		ifs.read (buff,length);
		ifs.close();

		md5_append(&pms, (const md5_byte_t *)buff, length);

		delete[] buff;

		//std::cerr<<".... Adding file to MD5 checksum : " << xml_filenames[i] << "  (size=" << length << ")" << std::endl;
	}
	
	md5_byte_t digest[16];
	md5_finish(&pms, digest);
	
	char hex_output[16*2 + 1];
	for(int di = 0; di < 16; ++di) sprintf(hex_output + di * 2, "%02x", digest[di]);

	return last_md5_checksum = hex_output;
}

