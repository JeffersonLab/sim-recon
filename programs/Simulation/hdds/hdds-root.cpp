/*
 *  hdds-root :   an interface utility that reads in a HDDS document
 *		   (Hall D Detector Specification) and writes out a
 *		   ROOT macro to instantiate the geometry within ROOT.
 *
 *  Original version - Edward Brash, November 1 2003.
 *  Based on hdds-geant by Richard Jones, May 19 2001.
 *
 *  Notes:
 *  ------
 * 1. The HDDS specification is an xml document, as described by HDDS.dtd.
 * 2. Access by hdds-root to the xml source is through the industry-
 *    standard DOM-1 interface.
 * 3. The code has been tested with the xerces-c DOM implementation from
 *    Apache, and is intended to be used with the xerces-c library.
 * 4. Output is sent to standard out through the ordinary c++ i/o library.
 * 5. As a by-product of using the DOM parser to access the xml source,
 *    hdds-root verifies the source against the dtd before translating it.
 *    Therefore it may also be used as a validator of the xml specification
 *    (see the -v option).
 *
 *  Implementation:
 *  ---------------
 * Most of the translation was straight-forward, but there were a couple
 * of tricky cases where decisions had to be made.  I think that these
 * choices should work out in most cases.  If not, further tuning of the
 * algorithm will be necessary.  Here are the tricky cases.
 *
 * 1. When to use divisions instead of placing multiple copies.
 * 
 *  Most of the time when a <mpos...> command appears it can be translated
 *  into a division of the mother volume in Geant.  This is good to do
 *  because it makes both more compact description and is more efficient
 *  at tracking time.  The difficulty here is that there is no easy way
 *  to check if the contents fit entirely inside the division.  This is
 *  not a problem in the HDDS geometry description because the <mpos...>
 *  command is only for positioning, and makes no statement about what
 *  slice of the mother it occupies.  But it is a problem for Geant because
 *  contents of divisions have to fit inside the division.  This can
 *  happen at any time, but it most frequently happens when the object
 *  is being rotated before placement, as in the case of stereo layers
 *  in a drift chamber.  My solution is to make a strict set of rules
 *  that are required before hdds-root will create a division in response
 *  to a <mpos...> command, and do individual placement by default.
 *  The rules for creation of divisions are as follows:
 *      (a) the <composition> command must have a solid container, either
 *          via the envelope="..." attribute or by itself being a division.
 *      (b) the <mpos...> command must be alone inside its <composition>
 *      (c) the kind of shape of the container must be compatible with the
 *          <mpos...> type, eg. <mposPhi> would work if its container is
 *          a "tubs" (or division theroef) but not if it is a "box".
 *      (d) for <mposPhi> the impliedRot attribute must be "true".
 *      (e) the rot="..." attribute must be zeros or missing.
 *  The last condition is not logically necessary for it to work, but it
 *  avoids the problems that often occur with rotated placements failing
 *  to fit inside the division.  To bypass this limitation and place
 *  rotated volumes inside divisions, one can simply create a new volume
 *  as a <composition> into which the content is placed with rotation,
 *  and then place the new volume with the <mpos...> command without rot.
 *
 * 2. How to recognize which media contain magnetic fields.
 *
 *  There is no provision in the hdds geometry model for magnetic field
 *  information.  Ultimately that is something that will be stored as a map
 *  somewhere in a database.  Geant needs to distinguish between 4 cases:
 *       (0) no magnetic field
 *       (1) general case of inhomogenous field (Runge-Kutta)
 *       (2) quasi-homogenous field with map (helical segments)
 *       (3) uniform field (helices along local z-axis)
 *  My solution is as follows.  By default, I assume case 0.  For all
 *  contents of a composition named "fieldVolume" I assign case 2.  For
 *  all contents of a composition named "*Magnet*" I assign case 3.
 *  For the magnitude of the field, I simply store a constant field value
 *  (kG) for each case, and rely on the GUFLD user routine in Geant3 to
 *  handle the actual field values at tracking time.
 *
 * 3. What to do about stackX/stackY/stackZ tags
 *
 *  In the case of Boolean tags (union/intersection/subtraction) the choice
 *  was easy: no support in Geant3 for these cases.  For the stacks it is
 *  possible to construct such structures in Geant3 but it is complicated
 *  by the fact that the stacks do not give information about the kind or
 *  size of volume that should be used for the mother.  Since stacks can
 *  be implemented easily using compositions anyway, I decided not to include
 *  support for them in hdds-root.
 */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE true
#endif

/*
 * FIX_XERCES_getElementById_BUG does a store/load cycle at parsing time
 * to fully instantiate entity references on the document tree.
 * See xerces-c++ bug 12800 at http://nagoya.apache.org
 */
#define FIX_XERCES_getElementById_BUG true


#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/XMLStringTokenizer.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/dom/DOM.hpp>

#include "hdds-root.hpp"

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

#include <fstream>

#define X(XString) XString.unicodeForm()
#define S(XString) XString.localForm()

#ifdef _Tru64
#undef basename
#define basename _RC_basename
static char * basename(const char *f)
{ const char *base;
                                                                                
  for(base = f; *f; f++)
  { if (*f == '/')
      base = f+1;
  }
                                                                                
  return (char *)base;
}
#endif

double fieldStrength[] =
{
   0.0,	  // zero field regions
   0.0,   // inhomogenous field regions (unused)
   22.4,  // mapped field regions: solenoid (kG, approximate)
   2.0    // uniform field regions: sweep magnets (kG)
};

int first_volume_placement = 0;

void usage()
{
    cerr << "\nUsage:\n"
            "    hdds-root [-v] {HDDS file}\n\n"
            "Options:\n"
            "    -v   validate only\n"
         << endl;
}

class Refsys
{
 public:
   DOMElement* fMother;		// current mother volume element
   int fMagField;		// flag indicating magnetic field
   double fOrigin[3];		// x,y,z coordinate of volume origin (cm)
   double fPhiOffset;		// azimuthal angle of volume origin (deg)
   double fRmatrix[3][3];	// rotation matrix (daughter -> mother)
   int fRotation;		// unique Rmatrix flag

   static char* fIdentifierList; // list of identifier strings (space-sep)
   static int fVolumes;		// total number of volumes to far

   struct VolIdent
   {
      XString fieldS;	
      int value;
      int step;
   };
   struct VolIdent fIdentifier[999];	// identifier tag list 
   int fIdentifiers;			// length of identifier tag list

   Refsys();				// empty constructor
   Refsys(const Refsys& src);		// copy constructor
   Refsys& operator=(Refsys& src);	// copy operator
   Refsys& reset();			// reset origin, Rmatrix
   Refsys& reset(const Refsys& ref);	// reset origin, Rmatrix to ref
   Refsys& shift(const double vector[3]); // translate origin
   Refsys& shift(const Refsys& ref);	  // copy origin from ref
   Refsys& shift(const Refsys& ref,
                 const double vector[3]); // translate origin in ref frame
   Refsys& rotate(const double omega[3]); // rotate by vector omega (radians)
   Refsys& rotate(const Refsys& ref);	  // copy Rmatrix from ref
   Refsys& rotate(const Refsys& ref,
                  const double omega[3]); // rotate by omega in ref frame

   int createMaterial(DOMElement* el);	// generate code for materials
   double getMaterialA(DOMElement* el);	// return A for a material
   double getMaterialZ(DOMElement* el);	// return Z for a material
   double getMixtureA(DOMElement* el);	// return A_eff for a mixture
   double getMixtureZ(DOMElement* el);	// return Z_eff for a mixture
   int createSolid(DOMElement* el);	// generate code for solids
   int createRotation();		// generate code for rotations
   int createDivision(char* divStr,
                      int ncopy,
                      int iaxis,
                      double start,
                      double step);	// generate code for divisions
   int createVolume(DOMElement* el);	// generate code for placement

private:
   static int fRotations;	// non-trivial rotations defined so far
};

void rootMacroHeader()
{
  cout << "void hddsroot()" << endl
       << "{" << endl
       << "//" << endl
       << "//  This file has been generated automatically via the " << endl
       << "//  utility hdds-root directly from main_HDDS.xml " << endl
       << "//   (see ROOT class TGeoManager for an example of use) " << endl
       << "//" << endl
       << "gSystem->Load(\"libGeom\");" << endl
       << "TGeoRotation *rot;" << endl
       << "TGeoNode *Node, *Node1;" << endl
       << " " << endl
       << "TGeoManager *detector = new TGeoManager(\"hddsroot\",\"hddsroot.C\");" << endl
       << " " << endl
       << " " << endl
       << "//-----------List of Materials and Mixtures--------------" << endl
       << " " << endl;
}

void rootMacroTrailer()
{
  cout << "gGeoManager->CloseGeometry();" << endl
       << "Double_t *origin = new Double_t[3];" << endl
       << "origin[0] = 450; origin[1] = -50; origin[2] = -200;" << endl
       << "TGeoBBox *clip = new TGeoBBox(\"CLIP\",300,300,300,origin);" << endl
       << "gGeoManager->SetClippingShape(clip);" << endl
       << "gGeoManager->DefaultColors();" << endl
       << "gGeoManager->SetVisLevel(9);" << endl
       << "HALL->Raytrace();" << endl
       << "}" << endl;
}

int main(int argC, char* argV[])
{
   try
   {
      XMLPlatformUtils::Initialize();
   }
   catch (const XMLException& toCatch)
   {
      XString message(toCatch.getMessage());
      cerr << "hdds-root: Error during initialization! :\n"
           << S(message) << endl;
      return 1;
   }

   if (argC < 2)
   {
      usage();
      return 1;
   }
   else if ((argC == 2) && (strcmp(argV[1], "-?") == 0))
   {
      usage();
      return 2;
   }

   const char* xmlFile = 0;
   bool rootMacroOutput = true;
   int argInd;
   for (argInd = 1; argInd < argC; argInd++)
   {
      if (argV[argInd][0] != '-')
         break;

      if (strcmp(argV[argInd], "-v") == 0)
         rootMacroOutput = false;
      else
         cerr << "Unknown option \'" << argV[argInd]
              << "\', ignoring it\n" << endl;
   }

   if (argInd != argC - 1)
   {
      usage();
      return 1;
   }
   xmlFile = argV[argInd];

#if defined OLD_STYLE_XERCES_PARSER
   DOMDocument* document = parseInputDocument(xmlFile);
#else
   DOMDocument* document = buildDOMDocument(xmlFile,false);
#endif
   if (document == 0)
   {
      cerr << "hdds-root : Error parsing HDDS document, "
           << "cannot continue" << endl;
      return 1;
   }

   DOMNode* docEl;
   try {
      docEl = document->getDocumentElement();
   }
   catch (DOMException& e) {
      XString msgS(e.msg);
      cerr << "Woops " << S(msgS) << endl;
      return 1;
   }

   XString everythingS("everything");
   DOMElement* rootEl = document->getElementById(X(everythingS));
   if (rootEl == 0)
   {
      cerr << "hdds-root : Error scanning HDDS document, " << endl
           << "  no element named \"eveything\" found" << endl;
      return 1;
   }

   if (rootMacroOutput)
   {
      Refsys mrs;
      rootMacroHeader();
      mrs.createVolume(rootEl);
      rootMacroTrailer();
   }

   XMLPlatformUtils::Terminate();
   return 0;
}

#if 0
/* Parser implemented using the old-style XercesDOMParser interface
 * based on the example code in $XERCESCROOT/samples/DOMPrint
 */
DOMDocument* parseInputDocument(const char* xmlFile)
{
   XercesDOMParser* parser = new XercesDOMParser;
   parser->setValidationScheme(XercesDOMParser::Val_Auto);
   parser->setCreateEntityReferenceNodes(false);
   parser->setValidationSchemaFullChecking(true);
   parser->setDoNamespaces(true);
   parser->setDoSchema(true);

   MyOwnErrorHandler errorHandler;
   parser->setErrorHandler(&errorHandler);

   try
   {
      parser->parse(xmlFile);
   }
   catch (const XMLException& toCatch)
   {
      XString message(toCatch.getMessage());
      cerr << "\nhdds-root: Error during parsing: '" << xmlFile << "'\n"
           << "Exception message is:  \n"
           << S(message) << "\n" << endl;
      return 0;
   }
   catch (const DOMException& toCatch)
   {
      XString message(toCatch.msg);
      cerr << "\nhdds-root: Error during parsing: '" << xmlFile << "'\n"
           << "Exception message is:  \n"
           << S(message) << "\n" << endl;
      XMLPlatformUtils::Terminate();
      return 0;
   }
   catch (...)
   {
      cerr << "\nhdds-root: Unexpected exception during parsing: '"
           << xmlFile << "'\n";
      XMLPlatformUtils::Terminate();
      return 0;
   }

   if (errorHandler.getSawErrors())
   {
      cerr << "\nErrors occured, no output available\n" << endl;
      return 0;
   }

   return parser->getDocument();
}

/* Parser implemented using the w3c standard DOMBuilder interface
 * based on the example code in $XERCESCROOT/samples/DOMCount
 */
DOMDocument* buildDOMDocument(const char* xmlFile)
{
   XString lsS("LS");
   DOMImplementation *impl = DOMImplementationRegistry::getDOMImplementation(X(lsS));
   DOMBuilder* parser = ((DOMImplementationLS*)impl)->createDOMBuilder(DOMImplementationLS::MODE_SYNCHRONOUS, 0);
   DOMWriter* writer = ((DOMImplementationLS*)impl)->createDOMWriter();
   XString tmpFileS(".tmp-");
   XString suffix(basename(xmlFile));
   tmpFileS += suffix;

   parser->setFeature(XMLUni::fgDOMValidation, true);
   parser->setFeature(XMLUni::fgDOMNamespaces, true);
   parser->setFeature(XMLUni::fgDOMDatatypeNormalization, true);
   parser->setFeature(XMLUni::fgDOMEntities, false);
   parser->setFeature(XMLUni::fgXercesSchemaFullChecking, true);
   parser->setFeature(XMLUni::fgXercesSchema, true);

   MyDOMErrorHandler* errHandler = new MyDOMErrorHandler();
   parser->setErrorHandler(errHandler);

   DOMDocument* doc = 0;

   try {
      parser->resetDocumentPool();
      doc = parser->parseURI(xmlFile);
#if defined FIX_XERCES_getElementById_BUG
      LocalFileFormatTarget* lfft = new LocalFileFormatTarget(X(tmpFileS));
      writer->writeNode(lfft,*(doc->getDocumentElement()));
      delete lfft;
      parser->resetDocumentPool();
      doc = parser->parseURI(X(tmpFileS));
#endif
   }
   catch (const XMLException& toCatch) {
      XString message(toCatch.getMessage());
      cout << "Exception message is: \n"
           << S(message) << "\n";
      return 0;
   }
   catch (const DOMException& toCatch) {
      XString message(toCatch.msg);
      cout << "Exception message is: \n"
           << S(message) << "\n";
      return 0;
   }
   catch (...) {
      cout << "Unexpected Exception \n" ;
      return 0;
   }

   if (errHandler->getSawErrors())
   {
      cerr << "\nErrors occured, no output available\n" << endl;
      return 0;
   }

   return doc;
}


MyOwnErrorHandler::MyOwnErrorHandler() : 
   fSawErrors(false)
{
}

MyOwnErrorHandler::~MyOwnErrorHandler()
{
}

// Overrides of the SAX ErrorHandler interface

void MyOwnErrorHandler::error(const SAXParseException& e)
{
   fSawErrors = true;
   XString systemId(e.getSystemId());
   XString message(e.getMessage());
   cerr << "\nhdds-root: Error at file " << S(systemId)
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << S(message) << endl;
}

void MyOwnErrorHandler::fatalError(const SAXParseException& e)
{
   fSawErrors = true;
   XString systemId(e.getSystemId());
   XString message(e.getMessage());
   cerr << "\nhdds-root: Fatal Error at file " << S(systemId)
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << S(message) << endl;
}

void MyOwnErrorHandler::warning(const SAXParseException& e)
{
   XString systemId(e.getSystemId());
   XString message(e.getMessage());
   cerr << "\nhdds-root: Warning at file " << S(systemId)
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << S(message) << endl;
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

bool MyDOMErrorHandler::handleError(const DOMError& domError)
{
   fSawErrors = true;
   if (domError.getSeverity() == DOMError::DOM_SEVERITY_WARNING)
      cerr << "\nWarning at file ";
   else if (domError.getSeverity() == DOMError::DOM_SEVERITY_ERROR)
       cerr << "\nError at file ";
   else
       cerr << "\nFatal Error at file ";

   cerr << XString(domError.getLocation()->getURI()).localForm()
        << ", line " << domError.getLocation()->getLineNumber()
        << ", char " << domError.getLocation()->getColumnNumber()
        << "\n  Message: " << XString(domError.getMessage()).localForm()
       	<< endl;

   return true;
}

void MyDOMErrorHandler::resetErrors()
{
   fSawErrors = false;
}

#endif

int Refsys::fVolumes = 0;
int Refsys::fRotations = 0;
char* Refsys::fIdentifierList = 0;

Refsys::Refsys()			// empty constructor
{
   fMother = 0;
   fMagField = 0;
   fPhiOffset = 0;
   fIdentifiers = 0;
   this->reset();
}

Refsys::Refsys(const Refsys& src)	// copy constructor
{
   reset(src);
   fMother = src.fMother;
   fMagField = src.fMagField;
   fPhiOffset = src.fPhiOffset;
   fIdentifiers = src.fIdentifiers;
   for (int i = 0; i < fIdentifiers; i++)
   {
      fIdentifier[i] = src.fIdentifier[i];
   }
} 

Refsys& Refsys::operator=(Refsys& src)	// copy operator (deep sematics)
{
   Refsys* dst = new Refsys(src);
   return *dst;
}

Refsys& Refsys::reset()			// reset Origin, Rmatrix to null
{
   fOrigin[0] = fOrigin[1] = fOrigin[2] = 0;
   fRmatrix[0][0] = fRmatrix[1][1] = fRmatrix[2][2] = 1;
   fRmatrix[0][1] = fRmatrix[1][0] = fRmatrix[1][2] =
   fRmatrix[0][2] = fRmatrix[2][0] = fRmatrix[2][1] = 0;
   fRotation = 0;
   return *this;
}

Refsys& Refsys::reset(const Refsys& ref) // reset Origin, Rmatrix to ref
{
   for (int i = 0; i < 3; i++)
   {
      fOrigin[i] = ref.fOrigin[i];
      fRmatrix[i][0] = ref.fRmatrix[i][0];
      fRmatrix[i][1] = ref.fRmatrix[i][1];
      fRmatrix[i][2] = ref.fRmatrix[i][2];
   }
   fRotation = ref.fRotation;
   return *this;
}

Refsys& Refsys::shift(const double vector[3])     // translate origin
{
   for (int i = 0; i < 3; i++)
   {
      fOrigin[i] += fRmatrix[i][0] * vector[0] +
                    fRmatrix[i][1] * vector[1] +
                    fRmatrix[i][2] * vector[2];
   }
   return *this;
}

Refsys& Refsys::shift(const Refsys& ref)      // copy origin from ref
{
   fOrigin[0] = ref.fOrigin[0];
   fOrigin[1] = ref.fOrigin[1];
   fOrigin[2] = ref.fOrigin[2];
   return *this;
}

Refsys& Refsys::shift(const Refsys& ref,
                      const double vector[3]) // translate origin in ref frame
{
   Refsys myRef(ref);
   myRef.shift(vector);
   return shift(myRef);
}

Refsys& Refsys::rotate(const double omega[3]) // rotate by vector omega (rad)
{
   if ( (omega[0] != 0.0) || (omega[1] != 0.0) || (omega[2] != 0.0) )
   {
      double cosx = cos(omega[0]);
      double sinx = sin(omega[0]);
      double cosy = cos(omega[1]);
      double siny = sin(omega[1]);
      double cosz = cos(omega[2]);
      double sinz = sin(omega[2]);

      for (int i = 0; i < 3; i++)
      {
         double x[3];
         double xx[3];
         x[0] = fRmatrix[i][0] * cosz + fRmatrix[i][1] * sinz;
         x[1] = fRmatrix[i][1] * cosz - fRmatrix[i][0] * sinz;
         x[2] = fRmatrix[i][2];
         xx[0] = x[0] * cosy - x[2] * siny;
         xx[1] = x[1];
         xx[2] = x[2] * cosy + x[0] * siny;
         fRmatrix[i][0] = xx[0];
         fRmatrix[i][1] = xx[1] * cosx + xx[2] * sinx;
         fRmatrix[i][2] = xx[2] * cosx - xx[1] * sinx;
      }

      fRotation = -1;
   }

   return *this;
}

Refsys& Refsys::rotate(const Refsys& ref)      // copy Rmatrix from ref
{
   for (int i = 0; i < 3; i++)
   {
      fRmatrix[i][0] = ref.fRmatrix[i][0];
      fRmatrix[i][1] = ref.fRmatrix[i][1];
      fRmatrix[i][2] = ref.fRmatrix[i][2];
   }
   fRotation = ref.fRotation;
   return *this;
}

Refsys& Refsys::rotate(const Refsys& ref,
                       const double omega[3])  // rotate by omega in ref frame
{
   Refsys myRef(ref);
   myRef.rotate(omega);
   return rotate(myRef);
}

double Refsys::getMaterialA(DOMElement* el)
{
  // XString tagS(el->getTagName());
  // XString nameAttS("name");
  // XString matS(el->getAttribute(X(nameAttS)));

   XString aAttS("a");
   XString aS(el->getAttribute(X(aAttS)));
   double a = atof(S(aS));

   return a;
}

double Refsys::getMaterialZ(DOMElement* el)
{
  // XString tagS(el->getTagName());
  // XString nameAttS("name");
  // XString matS(el->getAttribute(X(nameAttS)));

   XString zAttS("z");
   XString zS(el->getAttribute(X(zAttS)));
   double z = atof(S(zS));

   return z;
}

double Refsys::getMixtureA(DOMElement* el)
{
   XString tagS(el->getTagName());
   XString nameAttS("name");
   XString matS(el->getAttribute(X(nameAttS)));

   int nList = 0;
   int iList[999];
   double wList[999];
   double aList[999];
   double zList[999];
   XString addmaterialTagS("addmaterial");
   DOMDocument* document = el->getOwnerDocument();
   DOMNodeList* compList = el->getElementsByTagName(X(addmaterialTagS));
   //cout << "In getMixtureA: " << compList->getLength() << " <- value of compList->getLength()" << endl;
   for (int ic = 0; ic < compList->getLength(); ic++)
   {
      DOMNode* node = compList->item(ic);
      DOMElement* elem = (DOMElement*) node;
      XString compAttS("material");
      XString compS(elem->getAttribute(X(compAttS)));
      DOMElement* compEl = document->getElementById(X(compS));
      XString imateAttS("Geant3imate");
      XString cS(compEl->getAttribute(X(imateAttS)));

      aList[nList] = getMaterialA(compEl);
      zList[nList] = getMaterialZ(compEl);
      if ((aList[nList] == 0.0) || (zList[nList] == 0.0))
      {
	aList[nList] = getMixtureA(compEl);
	zList[nList] = getMixtureZ(compEl);
      }

      XString natomsTagS("natoms");
      XString frmassTagS("fractionmass");
      DOMNodeList* atomList = elem->getElementsByTagName(X(natomsTagS));
      DOMNodeList* fracList = elem->getElementsByTagName(X(frmassTagS));
      if (atomList->getLength() == 1)
      {
         DOMNode* node = atomList->item(0);
         DOMElement* elem = (DOMElement*) node;
         XString nAttS("n");
         XString nS(elem->getAttribute(X(nAttS)));
         wList[nList] = atoi(S(nS))*aList[nList];
      }
      else if (fracList->getLength() == 1)
      {
         DOMNode* node = fracList->item(0);
         DOMElement* elem = (DOMElement*) node;
         XString fractionAttS("fraction");
         XString fS(elem->getAttribute(X(fractionAttS)));
         wList[nList] = atof(S(fS));
      }
      else
      {
         cerr << "hdds-root error: " << S(tagS) << " " << S(matS)
              << " is missing some proportion data." << endl;
         exit(1);
      }
      //            cout << "In 2nd MIXTURE LOOP: " << nList << ", " << iList[nList] << ", " << wList[nList]
      //	 << ", " << aList[nList] << ", " << zList[nList] << endl;
      nList++;
   }

   double numer=0.0;
   double denom=0.0;
   for (int im = 0; im < nList; im++)
     { 
       numer = numer + aList[im]*wList[im];
       denom = denom + wList[im];
     }

   numer=numer/denom;

   return numer;

}

double Refsys::getMixtureZ(DOMElement* el)
{

   XString tagS(el->getTagName());
   XString nameAttS("name");
   XString matS(el->getAttribute(X(nameAttS)));

   int nList = 0;
   int iList[999];
   double wList[999];
   double aList[999];
   double zList[999];
   XString addmaterialTagS("addmaterial");
   DOMDocument* document = el->getOwnerDocument();
   DOMNodeList* compList = el->getElementsByTagName(X(addmaterialTagS));
   //cout << "In getMixtureA: " << compList->getLength() << " <- value of compList->getLength()" << endl;
   for (int ic = 0; ic < compList->getLength(); ic++)
   {
      DOMNode* node = compList->item(ic);
      DOMElement* elem = (DOMElement*) node;
      XString compAttS("material");
      XString compS(elem->getAttribute(X(compAttS)));
      DOMElement* compEl = document->getElementById(X(compS));
      XString imateAttS("Geant3imate");
      XString cS(compEl->getAttribute(X(imateAttS)));

      aList[nList] = getMaterialA(compEl);
      zList[nList] = getMaterialZ(compEl);
      if ((aList[nList] == 0.0) || (zList[nList] == 0.0))
      {
	aList[nList] = getMixtureA(compEl);
	zList[nList] = getMixtureZ(compEl);
      }

      XString natomsTagS("natoms");
      XString frmassTagS("fractionmass");
      DOMNodeList* atomList = elem->getElementsByTagName(X(natomsTagS));
      DOMNodeList* fracList = elem->getElementsByTagName(X(frmassTagS));
      if (atomList->getLength() == 1)
      {
         DOMNode* node = atomList->item(0);
         DOMElement* elem = (DOMElement*) node;
         XString nAttS("n");
         XString nS(elem->getAttribute(X(nAttS)));
         wList[nList] = atoi(S(nS))*aList[nList];
      }
      else if (fracList->getLength() == 1)
      {
         DOMNode* node = fracList->item(0);
         DOMElement* elem = (DOMElement*) node;
         XString fractionAttS("fraction");
         XString fS(elem->getAttribute(X(fractionAttS)));
         wList[nList] = atof(S(fS));
      }
      else
      {
         cerr << "hdds-root error: " << S(tagS) << " " << S(matS)
              << " is missing some proportion data." << endl;
         exit(1);
      }
      //            cout << "In 2nd MIXTURE LOOP: " << nList << ", " << iList[nList] << ", " << wList[nList]
      //	 << ", " << aList[nList] << ", " << zList[nList] << endl;
      nList++;
   }

   double numer=0.0;
   double denom=0.0;
   for (int im = 0; im < nList; im++)
     { 
       numer = numer + zList[im]*wList[im];
       denom = denom + wList[im];
     }

   numer=numer/denom;

   return numer;

}

int Refsys::createMaterial(DOMElement* el)
{
   static int imateCount = 0;
   int imate = ++imateCount;
   char imateStr[30];
   sprintf(imateStr, "%d", imate);
   XString imateAttS("Geant3imate");
   XString imateS(imateStr);
   el->setAttribute(X(imateAttS),X(imateS));

   XString tagS(el->getTagName());
   XString nameAttS("name");
   XString matS(el->getAttribute(X(nameAttS)));

   XString aAttS("a");
   XString aS(el->getAttribute(X(aAttS)));
   double a = atof(S(aS));

   XString zAttS("z");
   XString zS(el->getAttribute(X(zAttS)));
   double z = atof(S(zS));

   //cout << "At the start of createMaterial: " << a << ", " << z << endl;

   double dens = -1;
   double radl = -1;
   double absl = -1;
   XString realAttS("real");
   DOMNodeList* paramList = el->getElementsByTagName(X(realAttS));
   for (int ip = 0; ip < paramList->getLength(); ip++)
   {
      DOMNode* node = paramList->item(ip);
      DOMElement* elem = (DOMElement*) node;
      XString paramAttS("name");
      XString valueAttS("value");
      XString paramS(elem->getAttribute(X(paramAttS)));
      XString valueS(elem->getAttribute(X(valueAttS)));
      if (paramS.equals("density"))
      {
         dens = atof(S(valueS));
      }
      else if (paramS.equals("radlen"))
      {
         radl = atof(S(valueS));
      }
      else if (paramS.equals("abslen"))
      {
         absl = atof(S(valueS));
      }
   }

   int nList = 0;
   int iList[999];
   double wList[999];
   double aList[999];
   double zList[999];
   DOMDocument* document = el->getOwnerDocument();
   XString addmaterialTagS("addmaterial");
   DOMNodeList* compList = el->getElementsByTagName(X(addmaterialTagS));
   //cout << "Here I am !!! " << compList->getLength() << " <- value of compList->getLength()" << endl;
   for (int ic = 0; ic < compList->getLength(); ic++)
   {
      DOMNode* node = compList->item(ic);
      DOMElement* elem = (DOMElement*) node;
      XString compAttS("material");
      XString compS(elem->getAttribute(X(compAttS)));
      DOMElement* compEl = document->getElementById(X(compS));
      XString imateAttS("Geant3imate");
      XString cS(compEl->getAttribute(X(imateAttS)));
      if (cS != 0)
      {
	iList[nList] = atoi(S(cS));
      }
      else
      {
	iList[nList] = createMaterial(compEl);
      }

      aList[nList] = getMaterialA(compEl);
      zList[nList] = getMaterialZ(compEl);
      if ((aList[nList] == 0.0) || (zList[nList] == 0.0))
      {
	aList[nList] = getMixtureA(compEl);
	zList[nList] = getMixtureZ(compEl);
      }

      //cout << "In top level creatematerial loop: " << nList 
      //   << ", " << iList[nList] << ", " << aList[nList] << ", " 
      //   << zList[nList] << ", " << wList[nList] << endl;

      XString natomsTagS("natoms");
      XString frmassTagS("fractionmass");
      DOMNodeList* atomList = elem->getElementsByTagName(X(natomsTagS));
      DOMNodeList* fracList = elem->getElementsByTagName(X(frmassTagS));
      if (atomList->getLength() == 1)
      {
         XString typeS(compEl->getTagName());
         if (! typeS.equals("element"))
         {
            cerr << "hdds-root: error processing composite " << S(matS) << endl
                 << "natoms can only be specified for elements." << endl;
            exit(1);
         }
         else if (dens < 0)
         {
            cerr << "hdds-root error: " << S(tagS) << " " << S(matS)
                 << " is missing a density specification." << endl;
            exit(1);
         }
         DOMNode* node = atomList->item(0);
         DOMElement* elem = (DOMElement*) node;
         XString nAttS("n");
         XString nS(elem->getAttribute(X(nAttS)));
         wList[nList] = atoi(S(nS))*aList[nList];
      }
      else if (fracList->getLength() == 1)
      {
         DOMNode* node = fracList->item(0);
         DOMElement* elem = (DOMElement*) node;
         XString fractionAttS("fraction");
         XString fS(elem->getAttribute(X(fractionAttS)));
         wList[nList] = atof(S(fS));
	 XString aAttS("a");
	 XString aS(elem->getAttribute(X(aAttS)));

         double rho = 0;
         XString realTagS("real");
         DOMNodeList* propList = compEl->getElementsByTagName(X(realTagS));
         for (int ip = 0; ip < propList->getLength(); ip++)
         {
            DOMNode* pnode = propList->item(ip);
            DOMElement* pelem = (DOMElement*) pnode;
            XString pAttS("name");
            XString vAttS("value");
            XString pS(pelem->getAttribute(X(pAttS)));
            XString vS(pelem->getAttribute(X(vAttS)));
            if (pS.equals("density"))
            {
               rho = atof(S(vS));
            }
         }
         assert (rho > 0);
         if (dens < 0)
         {
            dens -= wList[nList] / rho;
         }
      }
      else
      {
         cerr << "hdds-root error: " << S(tagS) << " " << S(matS)
              << " is missing some proportion data." << endl;
         exit(1);
      }
      //      cout << "In MIXTURE LOOP: " << nList << ", " << iList[nList] << ", " << wList[nList]
      //<< ", " << dens << endl;
      nList++;
   }

   if (tagS.equals("element"))
   {
      if (dens < 0)
      {
         cerr << "hdds-root error: " << S(tagS) << " " << S(matS)
              << " is missing a density specification." << endl;
         exit(1);
      }
      if ((radl == -1) || (absl == -1))
      {
	//         cout << endl
        //      << "      imate = " << imate << endl
        //      << "      namate = \'" << S(matS) << "\'" << endl
        //      << "      a = " << a << endl
        //      << "      z = " << z << endl
        //      << "      dens = " << dens << endl
        //      << "      nlmat = 1" << endl
        //      << "      wmat(1) = 1" << endl
        //      << "      call gsmixt(imate,namate,a,z,dens,nlmat,wmat)"
        //      << endl;
	// cout << endl
	   cout << "TGeoMixture *mat" << imate << "= new TGeoMixture(\"" << S(matS) 
	      << "\",1," << dens << ");" << endl;
         cout << "mat" << imate << "->SetUniqueID(   " << imate << ");" << endl;
         cout << "mat" << imate << "->DefineElement(0,"<< a << "," << z << ",1.000);" << endl;
      }
      else
      {
	//cout << endl
	//     << "      imate = " << imate << endl
	//     << "      chnama = \'" << S(matS) << "\'" << endl
	//     << "      a = " << a << endl
	//     << "      z = " << z << endl
	//     << "      dens = " << dens << endl
	//     << "      radl = " << radl / (dens + 1e-30) << endl
	//     << "      absl = " << absl / (dens + 1e-30) << endl
	//     << "      nwbuf = 0" << endl
	//     << "      call gsmate(imate,chnama,a,z,dens,radl,absl,ubuf,nwbuf)"
	//     << endl;
	// cout << endl
	   cout << "TGeoMaterial *mat" << imate << "= new TGeoMaterial(\"" << S(matS) 
	      << "\"," << a << "," << z << "," << dens << ");" << endl;
         cout << "mat" << imate << "->SetUniqueID(   " << imate << ");" << endl;
      }
   }
   else
   {
      double rho = dens;
      if (rho < 0)
      {
         rho = -1/(dens + 0.999999);
      }
      //cout << endl
      //     << "      imate = " << imate << endl
      //    << "      namate = \'" << S(matS) << "\'" << endl;
      double wsum = 0.0;
      for (int im = 0; im < nList; im++)
	{
	  wsum = wsum + wList[im];
	}
      for (int im = 0; im < nList; im++)
      {
        // cout << "      wmat(" << im + 1 << ") = " << wList[im] << endl
	//    << "      call gfmate(" << iList[im] << ",chnama,"
	//    << "amat(" << im + 1 << "),zmat(" << im + 1 << "),"
	//    << "dens,radl,absl,ubuf,nwbuf)" << endl;
      }
      //cout << "      dens = " << rho << endl
      //   << "      nlmat = " << ((wList[0] < 1) ? nList : -nList) << endl
      //   << "      call gsmixt(imate,namate,amat,zmat,dens,nlmat,wmat)"
      //   << endl;
      //cout << endl
	cout << "TGeoMixture *mat" << imate << "= new TGeoMixture(\"" << S(matS) 
	  //	   << "\"," << ((wList[0] < 1) ? nList : -nList) << "," << rho << ");" << endl;
	  	   << "\"," << nList << "," << rho << ");" << endl;
      cout << "mat" << imate << "->SetUniqueID(   " << imate << ");" << endl;

      for (int im = 0; im < nList; im++)
      {
	//cout << "Writing mixture list: " << im << ", " << iList[im] << ", " << aList[im]<< ", " << zList[im] << endl; 
	cout << "mat" << imate << "->DefineElement(" << im << ","<< aList[im] << "," << zList[im] << "," << (wList[im]/wsum) << ");" << endl;
      }
   }
   
   if (dens < 0)
   {
      char densStr[30];
      sprintf(densStr, "%f", -1/(dens + 0.999999));
      XString realTagS("real");
      XString nameAttS("name");
      XString valueAttS("value");
      XString nameS("density");
      XString valueS(densStr);
      DOMElement* densEl = document->createElement(X(realTagS));
      densEl->setAttribute(X(nameAttS),X(nameS));
      densEl->setAttribute(X(valueAttS),X(valueS));
      el->appendChild(densEl);
   }

   return imate;
}

void getConversions(DOMElement* el, double& tocm, double& todeg)
{
   XString unitAttS("unit_length");
   XString unitS(el->getAttribute(X(unitAttS)));
   if (unitS.equals("mm"))
   {
      tocm = 0.1;
   }
   else if (unitS.equals("cm"))
   {
      tocm = 1.0;
   }
   else if (unitS.equals("m"))
   {
      tocm = 100.;
   }
   else
   {
      XString tagS(el->getTagName());
      cerr << "hdds-root error: unknown length unit " << S(unitS)
           << " on tag " << S(tagS) << endl;
      exit(1);
   }

   XString unaAttS("unit_angle");
   unitS = el->getAttribute(X(unaAttS));
   if (unitS.equals("deg"))
   {
      todeg = 1.0;
   }
   else if (unitS.equals("mrad"))
   {
      todeg = 0.180/M_PI;
   }
   else
   {
      XString tagS(el->getTagName());
      cerr << "hdds-root error: unknown angle unit " << S(unitS)
           << " on volume " << S(tagS) << endl;
      exit(1);
   }
}

int Refsys::createSolid(DOMElement* el)
{
   XString nameTagS("name");
   XString matAttS("material");
   XString nameS(el->getAttribute(X(nameTagS)));
   XString matS(el->getAttribute(X(matAttS)));

   int imate;
   DOMDocument* document = el->getOwnerDocument();
   DOMElement* matEl = document->getElementById(X(matS));
   XString imateAttS("Geant3imate");
   XString imateS(matEl->getAttribute(X(imateAttS)));
   if (imateS != 0)
   {
      imate = atoi(S(imateS));
   }
   else
   {
      imate = createMaterial(matEl);
   }
   
   static int itmedCount = 0;
   int itmed = ++itmedCount;
   XString sensiAttS("sensitive");
   XString sensiS(el->getAttribute(X(sensiAttS)));
   //cout << endl
   //    << "      itmed = " << itmed << endl
   //   << "      natmed = \'" << S(nameS) << " " << S(matS) << "\'" << endl
   //   << "      nmat = " << imate << endl
   //   << "      isvol = " << (sensiS.equals("true") ? 1 : 0) << endl
   //   << "      ifield = " << ((fMagField == 0) ? 0 : 2) << endl
   //   << "      fieldm = " << fieldStrength[fMagField] << endl
   //   << "      tmaxfd = " << ((fMagField == 0) ? 0 : 1) << endl
   //   << "      stemax = 0" << endl
   //   << "      deemax = 0" << endl
   //   << "      epsil = 1e-3" << endl
   //   << "      stmin = 0" << endl
   //   << "      nwbuf = 0" << endl
   //   << "      call gstmed(itmed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,"
   //   << endl
   //   << "     +            stemax,deemax,epsil,stmin,ubuf,nwbuf)" << endl;
   cout << "TGeoMedium *med" << itmed << " = new TGeoMedium(\"" << S(nameS) << " " << S(matS) << "\"," 
	<< itmed << "," << imate << "," << (sensiS.equals("true") ? 1 : 0) << "," 
	<< ((fMagField == 0) ? 0 : 2) << "," << fieldStrength[fMagField] << "," << ((fMagField == 0) ? 0 : 1) 
	<< ",-1,-1,0.1000000E-02,-1);"
	<< endl;

   double tocm, todeg;
   getConversions(el, tocm, todeg);

   double par[99];
   int npar = 0;
   XString shapeS(el->getTagName());
   if (shapeS.equals("box"))
   {
      shapeS = "BOX ";
      double xl, yl, zl;
      XString xyzAttS("X_Y_Z");
      XString xyzS(el->getAttribute(X(xyzAttS)));
      sscanf(S(xyzS), "%lf %lf %lf", &xl, &yl, &zl);

      npar = 3;
      par[0] = xl/2 * tocm;
      par[1] = yl/2 * tocm;
      par[2] = zl/2 * tocm;

      cout << "TGeoVolume *" << S(nameS) << "= gGeoManager->MakeBox(\"" << S(nameS) << "\",med" << itmed << "," 
	   << par[0] << "," << par[1] << "," << par[2] << ");" << endl;
   }
   else if (shapeS.equals("tubs"))
   {
      shapeS = "TUBS";
      double ri, ro, zl, phi0, dphi;
      XString riozAttS("Rio_Z");
      XString riozS(el->getAttribute(X(riozAttS)));
      sscanf(S(riozS), "%lf %lf %lf", &ri, &ro, &zl);
      XString profAttS("profile");
      XString profS(el->getAttribute(X(profAttS)));
      sscanf(S(profS), "%lf %lf", &phi0, &dphi);

      npar = 5;
      par[0] = ri * tocm;
      par[1] = ro * tocm;
      par[2] = zl/2 * tocm;
      par[3] = phi0 * todeg;
      par[4] = (phi0 + dphi) * todeg;
      if (dphi == 360)
	{
	  shapeS = "TUBE";
	  npar = 3;
	  
	  cout << "TGeoVolume *" << S(nameS) << "= gGeoManager->MakeTube(\"" 
	       << S(nameS) << "\",med" << itmed << "," 
	       << par[0] << "," << par[1] << "," << par[2] << ");" << endl;
	}
      else
	{
	  cout << "TGeoVolume *" << S(nameS) << "= gGeoManager->MakeTubs(\"" << S(nameS) 
	       << "\",med" << itmed << "," << par[0] << "," << par[1] << "," << par[2] 
	       << "," << par[3] << "," << par[4] << ");" << endl;
	}
   }
   else if (shapeS.equals("trd"))
   {
      shapeS = "TRAP";
      double xm, ym, xp, yp, zl;
      XString xyzAttS("Xmp_Ymp_Z");
      XString xyzS(el->getAttribute(X(xyzAttS)));
      sscanf(S(xyzS), "%lf %lf %lf %lf %lf", &xm, &xp, &ym, &yp, &zl);
      double alph_xz, alph_yz;
      XString incAttS("inclination");
      XString incS(el->getAttribute(X(incAttS)));
      sscanf(S(incS), "%lf %lf", &alph_xz, &alph_yz);

      npar = 11;
      double x = tan(alph_xz * M_PI/180);
      double y = tan(alph_yz * M_PI/180);
      double r = sqrt(x*x + y*y);
      par[0] = zl/2 * tocm;
      par[1] = atan2(r,1) * 180/M_PI;
      par[2] = atan2(y,x) * 180/M_PI;
      par[3] = ym/2 * tocm;
      par[4] = xm/2 * tocm;
      par[5] = xm/2 * tocm;
      par[6] = 0;
      par[7] = yp/2 * tocm;
      par[8] = xp/2 * tocm;
      par[9] = xp/2 * tocm;
      par[10] = 0;
      
      cout << "TGeoVolume *" << S(nameS) << "= gGeoManager->MakeTrap(\"" << S(nameS) 
	   << "\",med" << itmed << "," << par[0] << "," << par[1] << "," << par[2] 
	   << "," << par[3] << "," << par[4] << "," << par[5] << "," << par[6] 
	   << "," << par[7] << "," << par[8] << "," << par[9] << "," << par[10] << ");" << endl;
   }
   else if (shapeS.equals("pcon"))
   {
      shapeS = "PCON";
      double phi0, dphi;
      XString profAttS("profile");
      XString profS(el->getAttribute(X(profAttS)));
      sscanf(S(profS), "%lf %lf", &phi0, &dphi);
      XString planeTagS("polyplane");
      DOMNodeList* planeList = el->getElementsByTagName(X(planeTagS));

      npar = 3;
      par[0] = phi0 * todeg;
      par[1] = dphi * todeg;
      par[2] = planeList->getLength();
      for (int p = 0; p < planeList->getLength(); p++)
      {
         double ri, ro, zl;
         DOMNode* node = planeList->item(p);
         DOMElement* elem = (DOMElement*) node;
         XString riozAttS("Rio_Z");
         XString riozS(elem->getAttribute(X(riozAttS)));
         sscanf(S(riozS), "%lf %lf %lf", &ri, &ro, &zl);
         par[npar++] = zl * tocm;
         par[npar++] = ri * tocm;
         par[npar++] = ro * tocm;
      }

      cout << "TGeoVolume *" << S(nameS) << "= gGeoManager->MakePcon(\"" << S(nameS) 
	   << "\",med" << itmed << "," << par[0] << "," << par[1] << "," << par[2] << ");" << endl;
      for (int mycounter=0; mycounter < par[2]; mycounter++)
	{
	  cout << "  ((TGeoPcon*)" << S(nameS) << "->GetShape())->DefineSection(" << mycounter 
	       << "," << par[3+3*mycounter] << "," << par[4+3*mycounter] << "," 
	       << par[5+3*mycounter] << ");" << endl;
	}
   }
   else if (shapeS.equals("cons"))
   {
      shapeS = "CONS";
      double rim, rip, rom, rop, zl;
      XString riozAttS("Rio1_Rio2_Z");
      XString riozS(el->getAttribute(X(riozAttS)));
      sscanf(S(riozS), "%lf %lf %lf %lf %lf", &rim, &rom, &rip, &rop, &zl);
      double phi0, dphi;
      XString profAttS("profile");
      XString profS(el->getAttribute(X(profAttS)));
      sscanf(S(profS), "%lf %lf", &phi0, &dphi);

      npar = 7;
      par[0] = zl/2 * tocm;
      par[1] = rim * tocm;
      par[2] = rom * tocm;
      par[3] = rip * tocm;
      par[4] = rop * tocm;
      par[5] = phi0 * todeg;
      par[6] = (phi0 + dphi) * todeg;
      if (dphi == 360)
      {
         shapeS = "CONE";
         npar = 5;

	 cout << "TGeoVolume *" << S(nameS) << "= gGeoManager->MakeCone(\"" 
	      << S(nameS) << "\",med" << itmed << "," 
	      << par[0] << "," << par[1] << "," << par[2] 
	      << par[3] << "," << par[4] << ");" << endl;
      }
      else
	{
	  cout << "TGeoVolume *" << S(nameS) << "= gGeoManager->MakeCons(\"" << S(nameS) 
	       << "\",med" << itmed << "," << par[0] << "," << par[1] << "," << par[2] 
	       << "," << par[3] << "," << par[4] << "," << par[5] << "," << par[6] << ");" << endl;
	}
   }
   else
   {
      cerr << "hdds-root error: volume " << S(nameS)
           << " should be one of the valid shapes, not " << S(shapeS) << endl;
      exit(1);
   }

   if (nameS.stringLen() > 4)
   {
      cerr << "hdds-root error: volume name " << S(nameS)
           << " should be no more than 4 characters long." << endl;
      exit(1);
   }

   //   cout << endl
   //   << "      chname = \'" << S(nameS) << "\'" << endl
   //   << "      chshap = \'" << S(shapeS) << "\'" << endl
   //   << "      nmed = " << itmed << endl
   //   << "      npar = " << npar << endl;
   //for (int ipar = 0; ipar < npar; ipar++)
   //{
   //   cout << "      par(" << ipar + 1 << ") = " << par[ipar] << endl;
   //}
   //cout << "      call gsvolu(chname,chshap,nmed,par,npar,ivolu)" << endl;

   char ivoluStr[10];
   int ivolu = ++Refsys::fVolumes;
   sprintf(ivoluStr, "%d", ivolu);
   XString ivoluAttS("Geant3ivolu");
   XString icopyAttS("Geant3icopy");
   XString ivoluS(ivoluStr);
   XString icopyS("0");
   el->setAttribute(X(ivoluAttS),X(ivoluS));  
   el->setAttribute(X(icopyAttS),X(icopyS));  

/* consistency check #1: require Geant's volume index to match mine
 * 
 * This is required if the getX() lookup functions are going to work.
 * I count volumes in the order I define them, starting from 1.  If
 * Geant does the same thing then this error should never occur.
 */
//   cout << "      if (ivolu.ne." << ivolu << ")"
   //      << " stop \'consistency check #1 failed\'" << endl;

   return Refsys::fVolumes;
}

int Refsys::createRotation()
{
   if (fRotation < 0)
   {
      fRotation = ++fRotations;
   }
   else
   {
      return fRotation;
   }

   //   cout << endl
   //   << "      irot = " << fRotation << endl;

   double theta[3], phi[3];
   for (int i = 0; i < 3; i++)
   {
      double r = sqrt(fRmatrix[0][i] * fRmatrix[0][i]
                    + fRmatrix[1][i] * fRmatrix[1][i]);
      theta[i] = atan2(r, fRmatrix[2][i]) * 180/M_PI;
      phi[i] = atan2(fRmatrix[1][i], fRmatrix[0][i]) * 180/M_PI;
      // cout << "      theta" << i + 1 << " = " << theta[i] << endl
      //   << "      phi" << i + 1 << " = " << phi[i] << endl;
   }

   //cout << "      "
   //    << "call gsrotm(irot,theta1,phi1,theta2,phi2,theta3,phi3)"
   //     << endl;
   cout << "TGeoRotation *rot" << fRotation <<" = new TGeoRotation(\"rot" << fRotation << "\","
	<< theta[0] << "," << phi[0] << "," << theta[1] << "," << phi[1] << "," 
	<< theta[2] << "," << phi[2] << ");" << endl;

   return fRotation;
}

int Refsys::createDivision(char* divStr,
                           int ncopy,
                           int iaxis,
                           double start,
                           double step)
{
   divStr[0] = toupper(divStr[0]);
   divStr[1] = toupper(divStr[1]);
   divStr[2] = toupper(divStr[2]);
   divStr[3] = toupper(divStr[3]);

   assert (fMother != 0);

   XString nameAttS("name");
   XString motherS(fMother->getAttribute(X(nameAttS)));

   //cout << endl
   //   << "      chname = \'" << divStr << "\'" << endl
   //   << "      chmoth = \'" << S(motherS) << "\'" << endl
   //   << "      ndiv = " << ncopy << endl
   //   << "      iaxis = " << iaxis << endl
   //   << "      step = " << step << endl
   //   << "      c0 = " << start << endl
   //   << "      numed = 0" << endl
   //   << "      ndvmax = 0" << endl
   //   << "      call gsdvx(chname,chmoth,ndiv,iaxis,step,c0,numed,ndvmax)"
   //   << endl;
   cout << "TGeoVolume *" << divStr << "= "<< S(motherS) << "->Divide(\"" << divStr << "\","
	<< iaxis << "," << ncopy << "," << start << "," << step << ");" << endl;


   char attStr[30];
   DOMDocument* document = fMother->getOwnerDocument();
   XString divTagS("Geant3division");
   DOMElement* divEl = document->createElement(X(divTagS));
   XString divS(divStr);
   XString voluAttS("volume");
   divEl->setAttribute(X(nameAttS),X(divS));
   divEl->setAttribute(X(voluAttS),X(motherS));
   XString ivoluAttS("Geant3ivolu");
   sprintf(attStr, "%d", ++Refsys::fVolumes);
   XString ivoluS(attStr);
   divEl->setAttribute(X(ivoluAttS),X(ivoluS));  
   XString icopyAttS("Geant3icopy");
   sprintf(attStr, "%d", ncopy);
   XString icopyS(attStr);
   divEl->setAttribute(X(icopyAttS),X(icopyS));  
   fMother->appendChild(divEl);
   fMother = divEl;

   for (int id = 0; id < fIdentifiers; id++)
   {
      XString fieldS(fIdentifier[id].fieldS);
      int value = fIdentifier[id].value;
      int step = fIdentifier[id].step;
      XString idlistS;
      for (int ic = 0; ic < ncopy; ic++)
      {
         char str[30];
         sprintf(str, "%d ", value);
         idlistS += str;
         value += step;
      }
      divEl->setAttribute(X(fieldS),X(idlistS));
   }
   fIdentifiers = 0;
   return ncopy;
}

int Refsys::createVolume(DOMElement* el)
{
   int icopy = 0;

   Refsys myRef(*this);
   XString tagS(el->getTagName());
   XString nameAttS("name");
   XString nameS(el->getAttribute(X(nameAttS)));
   if (nameS.equals("fieldVolume"))
   {
      myRef.fMagField = 2;
   }
   else if (strstr(S(nameS),"Magnet"))
   {
      myRef.fMagField = 3;
   }

   DOMElement* env = 0;
   DOMDocument* document = el->getOwnerDocument();
   XString envAttS("envelope");
   XString envS(el->getAttribute(X(envAttS)));
   if (envS != 0)
   {
      env = document->getElementById(X(envS));
      XString contAttS("contains");
      XString containS(env->getAttribute(X(contAttS)));
      if (containS.equals(nameS))
      {
         return myRef.createVolume(env);
      }
      else if (containS != 0)
      {
         cerr << "hdds-root error: re-use of shape " << S(envS)
              << " is not allowed." << endl;
         exit(1);
      }
      env->setAttribute(X(contAttS),X(nameS));
      icopy = myRef.createVolume(env);
      myRef.fIdentifiers = 0;
      myRef.fMother = env;
      myRef.reset();
   }

   if (tagS.equals("intersection") ||
       tagS.equals("subtraction") ||
       tagS.equals("union"))
   {
      cerr << "hdds-root error: boolean " << S(tagS)
           << " operator is not supported." << endl;
      exit(1);
   }
   else if (tagS.equals("composition"))
   {
      DOMNode* cont;
      int nSiblings = 0;
      for (cont = el->getFirstChild(); 
           cont != 0;
           cont = cont->getNextSibling())
      {
         if (cont->getNodeType() == DOMNode::ELEMENT_NODE)
         {
            ++nSiblings;
         }
      }

      for (cont = el->getFirstChild(); 
           cont != 0;
           cont = cont->getNextSibling())
      {
         if (cont->getNodeType() != DOMNode::ELEMENT_NODE)
         {
            continue;
         }
         DOMElement* contEl = (DOMElement*) cont;
         XString comdS(contEl->getTagName());
         XString voluAttS("volume");
         XString targS(contEl->getAttribute(X(voluAttS)));
         DOMElement* targEl = document->getElementById(X(targS));

         Refsys drs(myRef);
         double origin[3], angle[3];
         XString rotAttS("rot");
         XString rotS(contEl->getAttribute(X(rotAttS)));
         sscanf(S(rotS), "%lf %lf %lf", &angle[0], &angle[1], &angle[2]);
         double tocm, todeg;
         getConversions(contEl, tocm, todeg);
         double torad = todeg * M_PI/180;
         angle[0] *= torad;
         angle[1] *= torad;
         angle[2] *= torad;
         bool noRotation = (angle[0] == 0) &&
                           (angle[1] == 0) &&
                           (angle[2] == 0) ;

         DOMNode* ident;
         for (ident = cont->getFirstChild(); 
              ident != 0;
              ident = ident->getNextSibling())
         {
            if (ident->getNodeType() != DOMNode::ELEMENT_NODE)
            {
               continue;
            }
            int id = drs.fIdentifiers++;
            DOMElement* identEl = (DOMElement*) ident;
            XString fieldAttS("field");
            XString valueAttS("value");
            XString stepAttS("step");
            XString fieldS(identEl->getAttribute(X(fieldAttS)));
            XString valueS(identEl->getAttribute(X(valueAttS)));
            XString stepS(identEl->getAttribute(X(stepAttS)));
            drs.fIdentifier[id].fieldS = fieldS;
            drs.fIdentifier[id].value = atoi(S(valueS));
            drs.fIdentifier[id].step = atoi(S(stepS));
            if (fIdentifierList == 0)
            {
               fIdentifierList = new char[32];
               fIdentifierList[0] = 0;
               strncpy(fIdentifierList, S(fieldS), 30);
            }
            else if (strstr(fIdentifierList, S(fieldS)) == 0)
            {
               int len = strlen(fIdentifierList);
               char* newList = new char[len + 32];
               strcpy(newList, fIdentifierList);
               strcat(newList, " ");
               strncat(newList, S(fieldS), 30);
               delete [] fIdentifierList;
               fIdentifierList = newList;
            }
         }

         if (comdS.equals("posXYZ"))
         {
            XString xyzAttS("X_Y_Z");
            XString xyzS(contEl->getAttribute(X(xyzAttS)));
            sscanf(S(xyzS), "%lf %lf %lf", &origin[0], &origin[1], &origin[2]);
            origin[0] *= tocm;
            origin[1] *= tocm;
            origin[2] *= tocm;
            drs.shift(origin);
            drs.rotate(angle);
            drs.createVolume(targEl);
         }
         else if (comdS.equals("posRPhiZ"))
         {
            double r, phi, z;
            XString rphizAttS("R_Phi_Z");
            XString rphizS(contEl->getAttribute(X(rphizAttS)));
            sscanf(S(rphizS), "%lf %lf %lf", &r, &phi, &z);
            double s;
            XString sAttS("S");
            XString sS(contEl->getAttribute(X(sAttS)));
            s = atof(S(sS));
            phi *= torad;
            r *= tocm;
            z *= tocm;
            s *= tocm;
            origin[0] = r * cos(phi) - s * sin(phi);
            origin[1] = r * sin(phi) + s * cos(phi);
            origin[2] = z;
            XString irotAttS("impliedRot");
            XString implrotS(contEl->getAttribute(X(irotAttS)));
            if (implrotS.equals("true") && (phi != 0))
            {
               angle[2] += phi;
            }
            drs.shift(origin);
            drs.rotate(angle);
            drs.createVolume(targEl);
         }
         else if (comdS.equals("mposPhi"))
         {
            XString ncopyAttS("ncopy");
            XString ncopyS(contEl->getAttribute(X(ncopyAttS)));
            int ncopy = atoi(S(ncopyS));
            if (ncopy <= 0)
            {
               cerr << "hdds-root error: volume " << S(nameS)
                    << " is positioned with " << ncopy << " copies!" << endl;
               exit(1);
            }

            double phi0, dphi;
            XString phi0AttS("Phi0");
            XString phi0S(contEl->getAttribute(X(phi0AttS)));
            phi0 = atof(S(phi0S)) * torad;
            XString dphiAttS("dPhi");
            XString dphiS(contEl->getAttribute(X(dphiAttS)));
            if (dphiS != 0)
            {
               dphi = atof(S(dphiS)) * torad;
            }
            else
            {
               dphi = 2 * M_PI / ncopy;
            }

            double r, s, z;
            XString rzAttS("R_Z");
            XString rzS(contEl->getAttribute(X(rzAttS)));
            sscanf(S(rzS), "%lf %lf", &r, &z);
            XString sAttS("S");
            XString sS(contEl->getAttribute(X(sAttS)));
            s = atof(S(sS));
            r *= tocm;
            z *= tocm;
            s *= tocm;

            XString containerS;
            if (env != 0)
            {
               containerS = env->getTagName();
            }
            else
            {
               XString divAttS("divides");
               containerS = el->getAttribute(X(divAttS));
            }
            XString irotAttS("impliedRot");
            XString implrotS(contEl->getAttribute(X(irotAttS)));
            if (noRotation && (nSiblings == 1) &&
                (containerS.equals("pcon") ||
                 containerS.equals("cons") ||
                 containerS.equals("tubs")) &&
                 implrotS.equals("true"))
            {
               static int phiDivisions = 0xd00;
               char* divStr = new char[5];
               sprintf(divStr, "s%3.3x", ++phiDivisions);
               phi0 *= 180/M_PI;
               dphi *= 180/M_PI;
               XString envAttS("envelope");
               XString targEnvS(targEl->getAttribute(X(envAttS)));
               DOMElement* targEnv;
               if (targEnvS != 0)
               {
                  targEnv = document->getElementById(X(targEnvS));
               }
               else
               {
                  targEnv = targEl;
               }
               XString profAttS("profile");
               XString profS(targEnv->getAttribute(X(profAttS)));
               if ((r == 0) && (profS != 0))
               {
                  double phi1, dphi1;
                  sscanf(S(profS), "%lf %lf", &phi1, &dphi1);
                  getConversions(targEnv, tocm, todeg);
                  double phiOffset = (phi1 + dphi1/2) * todeg;
		  drs.fPhiOffset = phiOffset * (M_PI/180);
                  phi0 += phiOffset;
               }
               int iaxis = 2;
               drs.createDivision(divStr, ncopy, iaxis, phi0 - dphi/2, dphi);
               XString divAttS("divides");
               targEl->setAttribute(X(divAttS),X(containerS));
               origin[0] = r;
               origin[1] = s;
               origin[2] = z;
               drs.reset();
               drs.shift(origin);
               drs.createVolume(targEl);
            }
            else
            {
               Refsys drs0(drs);
               drs.rotate(angle);
               for (int inst = 0; inst < ncopy; inst++)
               {
                  double phi = phi0 + inst * dphi;
                  origin[0] = r * cos(phi) - s * sin(phi);
                  origin[1] = r * sin(phi) + s * cos(phi);
                  origin[2] = z;
                  drs.shift(drs0, origin);
                  if (implrotS.equals("true"))
                  {
                     angle[2] += ((inst == 0) ? phi0 : dphi);
                     drs.rotate(drs0, angle);
                  }
                  drs.createVolume(targEl);
               }
            }
         }
         else if (comdS.equals("mposR"))
         {
            XString ncopyAttS("ncopy");
            XString ncopyS(contEl->getAttribute(X(ncopyAttS)));
            int ncopy = atoi(S(ncopyS));
            if (ncopy <= 0)
            {
               cerr << "hdds-root error: volume " << S(nameS)
                    << " is positioned with " << ncopy << " copies!" << endl;
               exit(1);
            }

            double r0, dr;
            XString r0AttS("R0");
            XString r0S(contEl->getAttribute(X(r0AttS)));
            r0 = atof(S(r0S)) * tocm;
            XString drAttS("dR");
            XString drS(contEl->getAttribute(X(drAttS)));
            dr = atof(S(drS)) * tocm;

            double phi, z, s;
            XString zphiAttS("Z_Phi");
            XString zphiS(contEl->getAttribute(X(zphiAttS)));
            sscanf(S(zphiS), "%lf %lf", &z, &phi);
            XString sAttS("S");
            XString sS(contEl->getAttribute(X(sAttS)));
            s = atof(S(sS));
            phi *= torad;
            z *= tocm;
            s *= tocm;

            XString containerS;
            if (env != 0)
            {
               containerS = env->getTagName();
            }
            else
            {
               XString divAttS("divides");
               containerS = el->getAttribute(X(divAttS));
            }
            if (noRotation && (nSiblings == 1) &&
                (containerS.equals("pcon") ||
                 containerS.equals("cons") ||
                 containerS.equals("tubs")))
            {
               static int rDivisions = 0xd00;
               char* divStr = new char[5];
               sprintf(divStr, "r%3.3x", ++rDivisions);
               int iaxis = 1;
               drs.createDivision(divStr, ncopy, iaxis, r0 - dr/2, dr);
               XString divAttS("divides");
               targEl->setAttribute(X(divAttS),X(containerS));
               origin[0] = r0 * cos(phi) - s * sin(phi);
               origin[1] = r0 * sin(phi) + s * cos(phi);
               origin[2] = z;
               drs.reset();
               drs.shift(origin);
               drs.rotate(angle);
               drs.createVolume(targEl);
            }
            else
            {
               Refsys drs0(drs);
               drs.rotate(angle);
               for (int inst = 0; inst < ncopy; inst++)
               {
                  double r = r0 + inst * dr;
                  origin[0] = r * cos(phi) - s * sin(phi);
                  origin[1] = r * sin(phi) + s * cos(phi);
                  origin[2] = z;
                  drs.shift(drs0, origin);
                  drs.createVolume(targEl);
               }
            }
         }
         else if (comdS.equals("mposX"))
         {
            XString ncopyAttS("ncopy");
            XString ncopyS(contEl->getAttribute(X(ncopyAttS)));
            int ncopy = atoi(S(ncopyS));
            if (ncopy <= 0)
            {
               cerr << "hdds-root error: volume " << S(nameS)
                    << " is positioned with " << ncopy << " copies!" << endl;
               exit(1);
            }

            double x0, dx;
            XString x0AttS("X0");
            XString x0S(contEl->getAttribute(X(x0AttS)));
            x0 = atof(S(x0S)) * tocm;
            XString dxAttS("dX");
            XString dxS(contEl->getAttribute(X(dxAttS)));
            dx = atof(S(dxS)) * tocm;

            double y, z, s;
            XString yzAttS("Y_Z");
            XString yzS(contEl->getAttribute(X(yzAttS)));
            sscanf(S(yzS), "%lf %lf", &y, &z);
            XString sAttS("S");
            XString sS(contEl->getAttribute(X(sAttS)));
            s = atof(S(sS));
            y *= tocm;
            z *= tocm;
            s *= tocm;

            XString containerS;
            if (env != 0)
            {
               containerS = env->getTagName();
            }
            else
            {
               XString divAttS("divides");
               containerS = el->getAttribute(X(divAttS));
            }
            if (noRotation && (nSiblings == 1) && 
                containerS.equals("box"))
            {
               static int xDivisions = 0xd00;
               char* divStr = new char[5];
               sprintf(divStr, "x%3.3x", ++xDivisions);
               int iaxis = 1;
               drs.createDivision(divStr, ncopy, iaxis, x0 - dx/2, dx);
               XString divAttS("divides");
               targEl->setAttribute(X(divAttS),X(containerS));
               origin[0] = 0;
               origin[1] = y + s;
               origin[2] = z;
               drs.reset();
               drs.shift(origin);
               drs.rotate(angle);
               drs.createVolume(targEl);
            }
            else
            {
               Refsys drs0(drs);
               drs.rotate(angle);
               for (int inst = 0; inst < ncopy; inst++)
               {
                  double x = x0 + inst * dx;
                  origin[0] = x;
                  origin[1] = y + s;
                  origin[2] = z;
                  drs.shift(drs0, origin);
                  drs.createVolume(targEl);
               }
            }
         }
         else if (comdS.equals("mposY"))
         {
            XString ncopyAttS("ncopy");
            XString ncopyS(contEl->getAttribute(X(ncopyAttS)));
            int ncopy = atoi(S(ncopyS));
            if (ncopy <= 0)
            {
               cerr << "hdds-root error: volume " << S(nameS)
                    << " is positioned with " << ncopy << " copies!" << endl;
               exit(1);
            }

            double y0, dy;
            XString y0AttS("Y0");
            XString y0S(contEl->getAttribute(X(y0AttS)));
            y0 = atof(S(y0S)) * tocm;
            XString dyAttS("dY");
            XString dyS(contEl->getAttribute(X(dyAttS)));
            dy = atof(S(dyS)) * tocm;

            double x, z, s;
            XString zxAttS("Z_X");
            XString zxS(contEl->getAttribute(X(zxAttS)));
            sscanf(S(zxS), "%lf %lf", &z, &x);
            XString sAttS("S");
            XString sS(contEl->getAttribute(X(sAttS)));
            s = atof(S(sS));
            x *= tocm;
            z *= tocm;
            s *= tocm;

            XString containerS;
            if (env != 0)
            {
               containerS = env->getTagName();
            }
            else
            {
               XString divAttS("divides");
               containerS = el->getAttribute(X(divAttS));
            }
            if (noRotation && (nSiblings == 1) && 
                containerS.equals("box"))
            {
               static int yDivisions = 0xd00;
               char* divStr = new char[5];
               sprintf(divStr, "y%3.3x", ++yDivisions);
               int iaxis = 2;
               drs.createDivision(divStr, ncopy, iaxis, y0 - dy/2, dy);
               XString divAttS("divides");
               targEl->setAttribute(X(divAttS),X(containerS));
               origin[0] = x + s;
               origin[1] = 0;
               origin[2] = z;
               drs.reset();
               drs.shift(origin);
               drs.rotate(angle);
               drs.createVolume(targEl);
            }
            else
            {
               Refsys drs0(drs);
               drs.rotate(angle);
               for (int inst = 0; inst < ncopy; inst++)
               {
                  double y = y0 + inst * dy;
                  double phi = atan2(y,x);
                  origin[0] = x;
                  origin[1] = y;
                  origin[2] = z;
                  drs.shift(drs0, origin);
                  drs.createVolume(targEl);
               }
            }
         }
         else if (comdS.equals("mposZ"))
         {
            XString ncopyAttS("ncopy");
            XString ncopyS(contEl->getAttribute(X(ncopyAttS)));
            int ncopy = atoi(S(ncopyS));
            if (ncopy <= 0)
            {
               cerr << "hdds-root error: volume " << S(nameS)
                    << " is positioned with " << ncopy << " copies!" << endl;
               exit(1);
            }

            double z0, dz;
            XString z0AttS("Z0");
            XString z0S(contEl->getAttribute(X(z0AttS)));
            z0 = atof(S(z0S)) * tocm;
            XString dzAttS("dZ");
            XString dzS(contEl->getAttribute(X(dzAttS)));
            dz = atof(S(dzS)) * tocm;

            double x, y, s;
            XString xyAttS("X_Y");
            XString xyS(contEl->getAttribute(X(xyAttS)));
            if (xyS.stringLen() > 0)
            {
               sscanf(S(xyS), "%lf %lf", &x, &y);
            }
            else
            {
               double r, phi;
               XString rphiAttS("R_Phi");
               XString rphiS(contEl->getAttribute(X(rphiAttS)));
               sscanf(S(rphiS), "%lf %lf", &r, &phi);
	       phi *= torad;
               x = r * cos(phi);
               y = r * sin(phi);
            }
            XString sAttS("S");
            XString sS(contEl->getAttribute(X(sAttS)));
            s = atof(S(sS));
            x *= tocm;
            y *= tocm;
            s *= tocm;

            XString containerS;
            if (env != 0)
            {
               containerS = env->getTagName();
            }
            else
            {
               XString divAttS("divides");
               containerS = el->getAttribute(X(divAttS));
            }
            if (noRotation && (nSiblings == 1) &&
                (containerS != 0))
            {
               static int zDivisions = 0xd00;
               char* divStr = new char[5];
               sprintf(divStr, "z%3.3x", ++zDivisions);
               int iaxis = 3;
               drs.createDivision(divStr, ncopy, iaxis, z0 - dz/2, dz);
               XString divAttS("divides");
               targEl->setAttribute(X(divAttS),X(containerS));
               double phi = atan2(y,x);
               origin[0] = x - s * sin(phi);
               origin[1] = y + s * cos(phi);
               origin[2] = 0;
               drs.reset();
               drs.shift(origin);
               drs.rotate(angle);
               drs.createVolume(targEl);
            }
            else
            {
               Refsys drs0(drs);
               drs.rotate(angle);
               double phi = atan2(y,x);
               origin[0] = x - s * sin(phi);
               origin[1] = y + s * cos(phi);
               for (int inst = 0; inst < ncopy; inst++)
               {
                  double z = z0 + inst * dz;
                  origin[2] = z;
                  drs.shift(drs0, origin);
                  drs.createVolume(targEl);
               }
            }
         }
         else
         {
            cerr << "hdds-root error: composition of volume " << S(nameS)
                 << " contains unknown tag " << S(comdS) << endl;
            exit(1);
         }
      }
   }
   else if (tagS.equals("stackX") || 
            tagS.equals("stackY") ||
            tagS.equals("stackZ"))
   {
      cerr << "hdds-root error: stacks are not supported." << endl
           << "Use compositions instead." << endl;
      exit(1);
   }
   else
   {
      XString icopyAttS("Geant3icopy");
      XString icopyS(el->getAttribute(X(icopyAttS)));
      if (icopyS != 0)
      {
         icopy = atoi(S(icopyS));
      }
      else
      {
         XString profAttS("profile");
         XString profS(el->getAttribute(X(profAttS)));
         if (profS != 0)
         {
            double phi0, dphi;
	    double tocm, todeg;
            sscanf(S(profS), "%lf %lf", &phi0, &dphi);
            getConversions(el, tocm, todeg);
	    if ( (myRef.fOrigin[0] == 0) && (myRef.fOrigin[1] == 0) )
	    {
               phi0 -= myRef.fPhiOffset*(180/M_PI)/todeg;
	    }
            char pStr[80];
            sprintf(pStr, "%lf %lf", phi0, dphi);
            XString pS(pStr);
            el->setAttribute(X(profAttS),X(pS));
	 }
         myRef.createSolid(el);
         icopy = 0;
      }

      if (myRef.fMother != 0)      {
         XString nameAttS("name");
         XString motherS(myRef.fMother->getAttribute(X(nameAttS)));
         int irot = myRef.createRotation();
	 //         cout << endl
	 //   << "      chname = \'" << S(nameS) << "\'" << endl
	 //   << "      nr = " << ++icopy << endl
	 //   << "      chmoth = \'" << S(motherS) << "\'" << endl
	 //   << "      x = " << myRef.fOrigin[0] << endl
	 //   << "      y = " << myRef.fOrigin[1] << endl
	 //   << "      z = " << myRef.fOrigin[2] << endl
	 //   << "      irot = " << irot << endl
	 //   << "      chonly = \'ONLY\'" << endl
	 //   << "      call gspos(chname,nr,chmoth,x,y,z,irot,chonly)"
	 //   << endl;

	 icopy++;
	 if (first_volume_placement == 0) 
	   {
	     cout << "gGeoManager->SetTopVolume(" << S(motherS) << ");" << endl;
	     first_volume_placement = 1;
	   }
	 if (irot == 0) 
	   {
	     if ( (myRef.fOrigin[0] == 0) && (myRef.fOrigin[1] == 0) && (myRef.fOrigin[2] == 0))
	       {
		 cout << S(motherS) << "->AddNode(" << S(nameS) << "," << icopy << ",gGeoIdentity);" << endl;
	       } 
	     else  { 
		 cout << S(motherS) <<"->AddNode(" << S(nameS) << "," << icopy << ",new TGeoTranslation(" 
		      <<  myRef.fOrigin[0] << "," 
		      <<  myRef.fOrigin[1] << "," 
		      <<  myRef.fOrigin[2] << "));" << endl;
	       }
	   }
	 else
	   {
	     cout << S(motherS) <<"->AddNode(" << S(nameS)<< "," << icopy << ",new TGeoCombiTrans(" 
		  <<  myRef.fOrigin[0] << "," 
		  <<  myRef.fOrigin[1] << "," 
		  <<  myRef.fOrigin[2] << "," << "rot" << irot <<"));" << endl;
	   }
      }

      char icopyStr[30];
      sprintf(icopyStr, "%d", icopy);
      XString copyS(icopyStr);
      XString copyAttS("Geant3icopy");
      el->setAttribute(X(copyAttS),X(copyS));
      for (int id = 0; id < myRef.fIdentifiers; id++)
      {
         XString fieldS(myRef.fIdentifier[id].fieldS);
         XString idlistS(el->getAttribute(X(fieldS)));
         XString spaceS(" ");
         XMLStringTokenizer picker(X(idlistS),X(spaceS));
         int count = icopy;
         XString idS;
         for (idS = picker.nextToken(); idS != 0; idS = picker.nextToken())
         {
            count--;
         }
         XString mylistS(idlistS);
         for ( ; count > 1; --count)
         {
            mylistS += "0 ";
         }
         char str[30];
         sprintf(str, "%d ", myRef.fIdentifier[id].value);
         mylistS += str;
         el->setAttribute(X(fieldS),X(mylistS));
         myRef.fIdentifier[id].value += myRef.fIdentifier[id].step;
      }
   }
   return icopy;
}

#if 0

XString::XString(void)
{
   fUnicodeForm = XMLString::transcode("");
   fLocalForm = XMLString::replicate("");
}

XString::XString(const XMLCh* const x)
{
   if (x) {
      fUnicodeForm = XMLString::replicate(x);
      fLocalForm = XMLString::transcode(x);
   }
   else {
      fUnicodeForm = XMLString::transcode("");
      fLocalForm = XMLString::replicate("");
   }
}

XString::XString(const char* const s)
{
   if (s) {
      fUnicodeForm = XMLString::transcode(s);
      fLocalForm = XMLString::replicate(s);
   }
   else {
      fUnicodeForm = XMLString::transcode("");
      fLocalForm = XMLString::replicate("");
   }
}

XString::XString(const XString& X)
{
   if (X.fUnicodeForm) {
      fUnicodeForm = XMLString::replicate(X.fUnicodeForm);
      fLocalForm = XMLString::transcode(X.fUnicodeForm);
   }
   else {
      fUnicodeForm = XMLString::transcode("");
      fLocalForm = XMLString::replicate("");
   }
}

XString::~XString()
{
   XMLString::release(&fUnicodeForm);
   XMLString::release(&fLocalForm);
}

const char* XString::localForm() const
{
   return fLocalForm;
}

const XMLCh* XString::unicodeForm() const
{
   return fUnicodeForm;
}

bool XString::equals(const XString& X) const
{
   return XMLString::equals(fUnicodeForm,X.fUnicodeForm);
}

bool XString::equals(const char* const s) const
{
   return XMLString::equals(fLocalForm,s);
}

bool XString::equals(const XMLCh* const x) const
{
   return XMLString::equals(fUnicodeForm,x);
}

int XString::stringLen() const
{
   return XMLString::stringLen(fUnicodeForm);
}

bool XString::operator==(const int len) const
{
   return (stringLen() == len);
}

bool XString::operator!=(const int len) const
{
   return (stringLen() != len);
}

XString& XString::operator=(const XString& X)
{
   XMLString::release(&fUnicodeForm);
   XMLString::release(&fLocalForm);
   fUnicodeForm = XMLString::replicate(X.fUnicodeForm);
   fLocalForm = XMLString::transcode(X.fUnicodeForm);
   return *this;
}

XString& XString::operator+=(const XString& X)
{
   int len = stringLen() + X.stringLen();
   XMLCh* sum = new XMLCh[len+1];
   XMLString::copyString(sum,fUnicodeForm);
   XMLString::catString(sum,X.fUnicodeForm);
   XMLString::release(&fUnicodeForm);
   XMLString::release(&fLocalForm);
   fUnicodeForm = XMLString::replicate(sum);
   fLocalForm = XMLString::transcode(sum);
   delete [] sum;
   return *this;
}
#endif
