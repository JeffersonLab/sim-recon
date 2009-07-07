/*  HDDS Browser Classes
 *
 *  Author: richard.t.jones@uconn.edu
 *
 *  Original version - Richard Jones, June 3, 2008.
 *
 */

#ifndef SAW_HDDSBROWSER_DEF
#define SAW_HDDSBROWSER_DEF true

#include <vector>

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/XMLStringTokenizer.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/sax/ErrorHandler.hpp>

using namespace xercesc;

#include "XString.hpp"
#include "XParsers.hpp"
#include "hddsCommon.hpp"

class hddsBrowser
{
 /* The hddsBrowser class is a general utility class for looking up
  * geometry information in the hdds geometry tree.
  */
 public:
   hddsBrowser(const XString xmlFile);	// constructor from xml document
   std::vector<Refsys>* find(const XString volume, const Refsys *ref = 0,
                             const DOMElement *contEl = 0);
                                   	// look up a volume in the geometry
 private:
   DOMDocument *fGeomDoc;		// the geometry tree document
};

#endif
