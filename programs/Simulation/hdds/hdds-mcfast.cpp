/*
 *  hdds-mcfast :  an interface utility that reads in a HDDS document
 *		   (Hall D Detector Specification) and writes out a
 *		   geometry description suitable for input to the
 *		   mcfast Monte Carlo detector simulation.
 *
 *  Original version - Richard Jones, May 16 2001.
 *
 *  Notes:
 *  ------
 * 1. The HDDS specification is an xml document, as described by HDDS.dtd.
 * 2. Access by hdds-mcfast to the xml source is through the industry-
 *    standard DOM-1 interface.
 * 3. The code has been tested with the xerces-c DOM implementation from
 *    Apache, and is intended to be used with the xerces-c library.
 * 4. Output is sent to standard out through the ordinary c++ i/o library.
 * 5. Within the HDDS document are references to mcfast db files that list
 *    the variables required on each output line.  These are looked for
 *    starting from the current working directory, and are typically given
 *    as "db/*".  Thus hdds-mcfast must generally be invoked after setting
 *    the current working directory to the one where "db" is located.
 * 6. As a by-product of using the DOM parser to access the xml source,
 *    hdds-mcfast verifies the source against the dtd before translating it.
 *    Therefore it may also be used as a validator of the xml specification
 *    (see the -v option).
 */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax/SAXException.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/dom/DOM.hpp>

#include "hdds-mcfast.hpp"

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

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

char* templateFileList = 0;

void usage()
{
    cerr << "\nUsage:\n"
            "    hdds-mcfast [-v] {HDDS file}\n\n"
            "Options:\n"
            "    -v   validate only\n"
         << endl;
}

void processTemplateFile(const char* fname, MCfastParameterList& pars)
{
   ifstream dbFile(fname);
   if (dbFile == 0)
   {
      cerr << "hdds-mcfast: Error opening input file " << fname << endl;
      exit(2);
   }

   while (! dbFile.eof())
   {
      char line[250];
      dbFile.getline(line,250);
      if ((strlen(line) == 0) || (line[0] == '!'))
      {
         continue;
      }
      char* token = strtok(line," ");
      if (strcasecmp(token,"end") == 0)
      {
         break;
      }
      else if (strcasecmp(token,"include") == 0)
      {
         processTemplateFile(strtok(0," "),pars);
      }
      else if (strcasecmp(token,"template") == 0)
      {
         continue;
      }
      else if (strcasecmp(token,"int") == 0)
      {
         char* var = strtok(0,", (");
         char* arg = strtok(0,")");
         int   dim = (arg == 0) ? 0 : atoi(arg);
         const DOMElement* el;
         XString valueS;
         if (dim == 0)
         {
            el = pars.get(var,"int");
	    XString valAttS("value");
            valueS = el->getAttribute(X(valAttS));
         }
         else
         {
            el = pars.get(var,"int_array",dim);
	    XString valAttS("values");
            valueS = el->getAttribute(X(valAttS));
         }
         if (el)
         {
            cout << " " << S(valueS);
         }
         else
         {
            cerr << "hdds-mcfast: Parameter \"" << var << "\""
                 << " required in template file " << fname << endl
                 << "is missing from HDDS" << endl;
            exit(3);
         }
      }
      else if (strcasecmp(token,"real") == 0)
      {
         char* var = strtok(0,", (");
         char* arg = strtok(0,")");
         int   dim = (arg == 0) ? 0 : atoi(arg);
         const DOMElement* el;
         XString valueS;
         if (dim == 0)
         {
            el = pars.get(var,"real");
	    XString valAttS("value");
            valueS = el->getAttribute(X(valAttS));
         }
         else
         {
            el = pars.get(var,"real_array",dim);
	    XString valAttS("values");
            valueS = el->getAttribute(X(valAttS));
         }
         if (el)
         {
            cout << " " << S(valueS);
         }
         else
         {
            cerr << "hdds-mcfast: Parameter \"" << var << "\""
                 << " required in template file " << fname << endl
                 << "is missing from HDDS" << endl;
            exit(3);
         }
      }
      else if (strcasecmp(token,"char") == 0)
      {
         char* var = strtok(0,", (");
         char* arg = strtok(0,")");
         int   dim = (arg == 0) ? 0 : atoi(arg);
         const DOMElement* el;
         if (dim == 0)
         {
            if (! (el = pars.get(var,"string")) )
            {
               cerr << "hdds-mcfast: Parameter \"" << var << "\""
                    << " required in template file " << fname << endl
                    << "is missing from HDDS" << endl;
               exit(3);
            }
	    XString valAttS("value");
            XString valueS(el->getAttribute(X(valAttS)));
            cout << " \"" << S(valueS) << "\"";
         }
         else
         {
            if (! (el = pars.get(var,"string_vector")) )
            {
               cerr << "hdds-mcfast: Parameter \"" << var << "\""
                    << " required in template file " << fname << endl
                    << "is missing from HDDS" << endl;
               exit(3);
            }
            int vcount = 0;
            DOMNode* vect;
            for ( vect = el->getFirstChild(); 
                  vect != 0;
                  vect = vect->getNextSibling() )
            {
               if (vect->getNodeType() != DOMNode::ELEMENT_NODE) continue;
               DOMElement* vectEl = (DOMElement*) vect;
               XString tagS(vectEl->getTagName());
               if (tagS.equals("string_data"))
               {
	          XString valAttS("value");
                  XString valueS(vectEl->getAttribute(X(valAttS)));
                  cout << " \"" << S(valueS) << "\"";
                  vcount++;
               }
            }
            if (vcount != dim)
            {
	       XString nameAttS("name");
               XString vnameS(el->getAttribute(X(nameAttS)));
               cerr << "hdds-mcfast: String vector "
                    << S(vnameS) << " has too few elements, "
                    << dim << " required." << endl;
               exit(3);
            }
         }
      }
      else if ((strcasecmp(token,"parent") == 0) ||
               (strcasecmp(token,"child") == 0))
      {
         continue;
      }
      else if (strcasecmp(token,"material") == 0)
      {
         char* var = strtok(0,", (");
         char* arg = strtok(0,")");
         int   dim = (arg == 0) ? 0 : atoi(arg);
         const DOMElement* el;
         if (dim == 0)
         {
            if (! (el = pars.get(var,"reference")) )
            {
               cerr << "hdds-mcfast: Parameter \"" << var << "\""
                    << " required in template file " << fname << endl
                    << "is missing from HDDS" << endl;
               exit(3);
            }
	    XString valAttS("value");
            XString valueS(el->getAttribute(X(valAttS)));
            cout << " \"" << S(valueS) << "\"";
         }
         else
         {
            if (! (el = pars.get(var,"reference_vector")) )
            {
               cerr << "hdds-mcfast: Parameter \"" << var << "\""
                    << " required in template file " << fname << endl
                    << "is missing from HDDS" << endl;
               exit(3);
            }
            DOMNode* vect;
            int vcount = 0;
            for ( vect = el->getFirstChild();
                  vect != 0;
                  vect = vect->getNextSibling() )
            {
               if (vect->getNodeType() != DOMNode::ELEMENT_NODE) continue;
               DOMElement* vectEl = (DOMElement*) vect;
               XString tagS(vectEl->getTagName());
               if (tagS.equals("reference_data"))
               {
	          XString valAttS("value");
                  XString valueS(vectEl->getAttribute(X(valAttS)));
                  cout << " \"" << S(valueS) << "\"";
                  vcount++;
               }
            }
            if (vcount < dim)
            {
	       XString nameAttS("name");
               XString rnameS(el->getAttribute(X(nameAttS)));
               cerr << "hdds-mcfast: Reference vector "
                    << S(rnameS) << " has too few elements, "
                    << dim << " required." << endl;
               exit(3);
            }
         }
      }
      else
      {
         cerr << "hdds-mcfast: Template file " << fname
              << " contains unknown parameter type " << token << endl;
         exit(3);
      }
   }
}

void writeMCfastRecord(const DOMElement* el, MCfastParameterList* pars)
{
   XString parAttS("parameters");
   XString parBlockS(el->getAttribute(X(parAttS)));
   DOMElement* parEl = el->getOwnerDocument()->getElementById(X(parBlockS));
   DOMNode* parNode;
   int parCount;
   if (parEl != 0)
   {
      parNode = parEl->getFirstChild();
      XString wildS("*");
      parCount = parEl->getElementsByTagName(X(wildS))->getLength();
   }
   else
   {
      parNode = 0;
      parCount = 0;
   }
   MCfastParameterList basePars(parCount);
   basePars.inherits(*pars);
   for ( ; parNode != 0; parNode = parNode->getNextSibling() )
   {
      if (parNode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
      DOMElement* thisEl = (DOMElement*) parNode;
      basePars.append(thisEl);
   }

   XString wildS("*");
   parCount = el->getElementsByTagName(X(wildS))->getLength();
   MCfastParameterList myPars(parCount);
   myPars.inherits(basePars);
   for ( parNode = el->getFirstChild();
         parNode != 0;
         parNode = parNode->getNextSibling() )
   {
      if (parNode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
      DOMElement* thisEl = (DOMElement*) parNode;
      XString tagS(thisEl->getTagName());
      if (tagS.equals("mcfast"))
      {
         writeMCfastRecord(thisEl, &myPars);
      }
      else
      {
         myPars.append(thisEl);
      }
   }

   XString tempAttS("template");
   XString templS(el->getAttribute(X(tempAttS)));
   if (templateFileList == 0)
   {
      templateFileList = new char[99999];	// big enough to forget
      templateFileList[0] = 0;
   }
   if (strstr(templateFileList,S(templS)) == 0)
   {
      strcat(templateFileList,S(templS));
      cout << "include " << S(templS) << endl;
   }
     
   XString modAttS("model");
   XString modelS(el->getAttribute(X(modAttS)));
   cout << "make " << S(modelS);
   processTemplateFile(S(templS),myPars);
   cout << endl;
}

void writeMaterialRecord(DOMElement* el)
{
   MCfastParameterList parList(3);

   XString stringS("string");
   DOMElement* nameEl = el->getOwnerDocument()->createElement(X(stringS));
   XString nameAttS("name");
   XString valueAttS("value");
   nameEl->setAttribute(X(nameAttS),X(nameAttS));
   nameEl->setAttribute(X(valueAttS),el->getAttribute(X(nameAttS)));
   parList.append(nameEl);

   XString realS("real");
   DOMElement* aEl = el->getOwnerDocument()->createElement(X(realS));
   XString aS("a");
   aEl->setAttribute(X(nameAttS),X(aS));
   aEl->setAttribute(X(valueAttS),el->getAttribute(X(aS)));
   parList.append(aEl);

   DOMElement* zEl = el->getOwnerDocument()->createElement(X(realS));
   XString zS("z");
   zEl->setAttribute(X(nameAttS),X(zS));
   zEl->setAttribute(X(valueAttS),el->getAttribute(X(zS)));
   parList.append(zEl);

   XString modelAttS("model");
   XString materialS("Material");
   el->setAttribute(X(modelAttS),X(materialS));
   XString templAttS("template");
   XString dbS("db/materials.db");
   el->setAttribute(X(templAttS),X(dbS));
   writeMCfastRecord(el,&parList);
}

void writeMixtureRecord(DOMElement* el)
{
   XString realS("real");
   DOMNodeList* propList = el->getElementsByTagName(X(realS));
   int propCount = propList->getLength();
   if (propCount > 4)
   {
      writeMaterialRecord(el);
      return;
   }

   MCfastParameterList parList(4);

   XString stringS("string");
   DOMElement* nameEl = el->getOwnerDocument()->createElement(X(stringS));
   XString nameAttS("name");
   XString valueAttS("value");
   nameEl->setAttribute(X(nameAttS),X(nameAttS));
   nameEl->setAttribute(X(valueAttS),el->getAttribute(X(nameAttS)));
   parList.append(nameEl);

   XString addmatS("addmaterial");
   DOMNodeList* matList = el->getElementsByTagName(X(addmatS));
   int matCount = matList->getLength();
   assert (matCount > 0);
   float fVol[matCount], fSum = 0.;
   XString refvecS("reference_vector");
   DOMElement* matVecEl =
               el->getOwnerDocument()->createElement(X(refvecS));
   DOMElement* dataEl[5];
   int m;
   for (m = 0; m < matCount; m++)
   {
      DOMNode* mat = matList->item(m);
      DOMElement* matEl = (DOMElement*) mat;
      XString materialS("material");
      XString refIdS(matEl->getAttribute(X(materialS)));
      XString refdataS("reference_data");
      dataEl[m] = el->getOwnerDocument()->createElement(X(refdataS));
      XString valueAttS("value");
      dataEl[m]->setAttribute(X(valueAttS),X(refIdS));
      matVecEl->appendChild(dataEl[m]);

      DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refIdS));
      XString aAttS("a");
      XString aS(refEl->getAttribute(X(aAttS)));
      float a = atof(S(aS));
      XString zAttS("z");
      XString zS(refEl->getAttribute(X(zAttS)));
      float z = atof(S(zS));

      DOMNodeList* parList = refEl->getElementsByTagName(X(realS));
      int parCount = parList->getLength();
      float density = 1.0e30;
      for (int p = 0; p < parCount; p++)
      {
         DOMNode* par = parList->item(p);
         DOMElement* parEl = (DOMElement*) par;
	 XString nameAttS("name");
         XString parName(parEl->getAttribute(X(nameAttS)));
         if (parName.equals("density"))
         {
	    XString valueAttS("value");
            XString valS(parEl->getAttribute(X(valueAttS)));
            density = atof(S(valS));
         }
      }

      XString wildS("*");
      DOMNodeList* specList = matEl->getElementsByTagName(X(wildS));
      assert (specList->getLength() > 0);
      DOMNode* spec = specList->item(0);
      DOMElement* specEl = (DOMElement*) spec;
      XString specTypeS(specEl->getTagName());
      float fMass = 0;
      if (specTypeS.equals("natoms"))
      {
	 XString nAttS("n");
         XString nS(specEl->getAttribute(X(nAttS)));
         int n = atoi(S(nS));
         fMass = n*a;
      }
      else if (specTypeS.equals("fractionmass"))
      {
	 XString fracAttS("fraction");
         XString fS(specEl->getAttribute(X(fracAttS)));
         fMass = atof(S(fS));
      }
      fVol[m] = fMass/density;
      fSum += fVol[m];
   }

   for (; m < 5; m++)
   {
      XString tagS("reference_data");
      XString valAttS("value");
      XString minuS("-");
      dataEl[m] = el->getOwnerDocument()->createElement(X(tagS));
      dataEl[m]->setAttribute(X(valAttS),X(minuS));
      matVecEl->appendChild(dataEl[m]);
      fVol[m] = 0;
   }

   char nmatStr[10];
   sprintf(nmatStr,"%d",matCount);
   XString nmatS(nmatStr);
   XString nmatAttS("nmat");
   XString intS("int");
   DOMElement* nmatEl = el->getOwnerDocument()->createElement(X(intS));
   nmatEl->setAttribute(X(nameAttS),X(nmatAttS));
   nmatEl->setAttribute(X(valueAttS),X(nmatS));
   parList.append(nmatEl);

   XString matnameS("matnames");
   matVecEl->setAttribute(X(nameAttS),X(matnameS));
   parList.append(matVecEl);

   char fVolStr[300];
   sprintf(fVolStr,"%f %f %f %f %f", fVol[0]/fSum, fVol[1]/fSum,
           fVol[2]/fSum, fVol[3]/fSum, fVol[4]/fSum);
   XString fvolS(fVolStr);
   XString rarrayS("real_array");
   DOMElement* fVolEl = el->getOwnerDocument()->createElement(X(rarrayS));
   XString propS("prop");
   fVolEl->setAttribute(X(nameAttS),X(propS));
   fVolEl->setAttribute(X(valueAttS),X(fvolS));
   parList.append(fVolEl);

   XString modelAttS("model");
   XString templAttS("template");
   XString mixtureS("Mixture");
   XString dbS("db/mixtures.db");
   el->setAttribute(X(modelAttS),X(mixtureS));
   el->setAttribute(X(templAttS),X(dbS));
   writeMCfastRecord(el,&parList);
}

void doMaterials(DOMElement* rootEl)
{
   XString elemS("element");
   DOMNodeList* matList = rootEl->getElementsByTagName(X(elemS));
   int matCount = matList->getLength();
   for (int m = 0; m < matCount; m++)
   {
      DOMNode* mat = matList->item(m);
      DOMElement* matEl = (DOMElement*) mat;
      XString realS("real");
      int nRealPars = matEl->getElementsByTagName(X(realS))->getLength();
      if (nRealPars > 4)
      {
         writeMaterialRecord(matEl);
      }
   }

   XString compoS("composite");
   DOMNodeList* compList = rootEl->getElementsByTagName(X(compoS));
   int compCount = compList->getLength();
   for (int c = 0; c < compCount; c++)
   {
      DOMNode* comp = compList->item(c);
      DOMElement* compEl = (DOMElement*) comp;
      XString realS("real");
      int nRealPars = compEl->getElementsByTagName(X(realS))->getLength();
      writeMixtureRecord(compEl);
   }
}

int main(int argC, char* argV[])
{
   try
   {
      XMLPlatformUtils::Initialize();
   }
   catch (const XMLException& toCatch)
   {
      XString msgS(toCatch.getMessage());
      cerr << "hdds-mcfast: Error during initialization! :\n"
           << S(msgS) << endl;
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

   const char*  xmlFile = 0;
   bool mcfastOutput = true;
   int argInd;
   for (argInd = 1; argInd < argC; argInd++)
   {
      if (argV[argInd][0] != '-')
         break;

      if (strcmp(argV[argInd],"-v") == 0)
         mcfastOutput = false;
      else
         cerr << "Unknown option '" << argV[argInd]
              << "', ignoring it\n" << endl;
   }

   if (argInd != argC - 1)
   {
      usage();
      return 1;
   }
   xmlFile = argV[argInd];

#if defined OLD_STYLE_XERCES_PARSER
   DOMDocument* doc = parseInputDocument(xmlFile);
#else
   DOMDocument* doc = buildDOMDocument(xmlFile);
#endif
   if (doc == 0)
   {
      cerr << "hdds-geant : Error parsing HDDS document, "
           << "cannot continue" << endl;
      return 1;
   }

   if (! mcfastOutput)
   {
      return 0;
   }

   cout << "database mcfast 0000" << endl;

   DOMElement* rootEl = doc->getDocumentElement();
   DOMNode* sect;
   for ( sect = rootEl->getLastChild();
         sect != 0;
         sect = sect->getPreviousSibling() )
   {
      if (sect->getNodeType() != DOMNode::ELEMENT_NODE) continue;
      DOMElement* sectEl = (DOMElement*) sect;
      XString stagS = sectEl->getTagName();
      if (stagS.equals("section"))
      {
         DOMNode* cont;
         for ( cont = sect->getFirstChild();
               cont != 0;
               cont = cont->getNextSibling() )
         {
            if (cont->getNodeType() != DOMNode::ELEMENT_NODE) continue;
            DOMElement* contEl = (DOMElement*) cont;
            XString tagS = contEl->getTagName();
            if (tagS.equals("mcfast"))
            {
               static bool firstTime = true;
               writeMCfastRecord(contEl, 0);
               if (firstTime)
               {
                  doMaterials(rootEl);
                  firstTime = false;
               }
            }
         }
      }
   }

   XMLPlatformUtils::Terminate();
   return 0;
}

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
      cerr << "\nhdds-geant: Error during parsing: '" << xmlFile << "'\n"
           << "Exception message is:  \n"
           << S(message) << "\n" << endl;
      return 0;
   }
   catch (const DOMException& toCatch)
   {
      XString message(toCatch.msg);
      cerr << "\nhdds-geant: Error during parsing: '" << xmlFile << "'\n"
           << "Exception message is:  \n"
           << S(message) << "\n" << endl;
      XMLPlatformUtils::Terminate();
      return 0;
   }
   catch (...)
   {
      cerr << "\nhdds-geant: Unexpected exception during parsing: '"
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

   XString badAttS("CDSI");
   DOMElement* targ = doc->getElementById(X(badAttS));
   return doc;
}

// Overrides of the SAX ErrorHandler interface

void MyOwnErrorHandler::error(const SAXParseException& e)
{
   fSawErrors = true;
   XString systemId(e.getSystemId());
   XString message(e.getMessage());
   cerr << "\nhdds-geant: Error at file " << S(systemId)
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << S(message) << endl;
}

void MyOwnErrorHandler::fatalError(const SAXParseException& e)
{
   fSawErrors = true;
   XString systemId(e.getSystemId());
   XString message(e.getMessage());
   cerr << "\nhdds-geant: Fatal Error at file " << S(systemId)
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << S(message) << endl;
}

void MyOwnErrorHandler::warning(const SAXParseException& e)
{
   XString systemId(e.getSystemId());
   XString message(e.getMessage());
   cerr << "\nhdds-geant: Warning at file " << S(systemId)
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
MyOwnErrorHandler::MyOwnErrorHandler() : 
   fSawErrors(false)
{
}

MyOwnErrorHandler::~MyOwnErrorHandler()
{
}


MCfastParameterList::MCfastParameterList(const int capacity)
{
   fListCap = (capacity > 0) ? capacity : 1;
   fListEl = new const DOMElement* [fListCap];
   fAncestor = 0;
   fListLen = 0;
}

MCfastParameterList::~MCfastParameterList()
{
   delete [] fListEl;
}

void MCfastParameterList::inherits(MCfastParameterList& anc)
{
   fAncestor = &anc;
}

void MCfastParameterList::append(const DOMElement* par)
{
   DOMElement* myCopy = (DOMElement*) par->cloneNode(true);
   if (fListLen < fListCap)
   {
      fListEl[fListLen++] = myCopy;
   }
   else
   {
      cerr << "hdds-mcfast: Error in MCfastParameterList::append" << endl
           << "   list overflow, insufficient capacity" << endl;
      exit(1);
   }
}

const DOMElement* MCfastParameterList::get(char* name, char* type,
                                            int number, char* unit) const
{
   int p = fListLen;
   while (--p >= 0)
   {
      const DOMElement* el = fListEl[p];
      XString nameAttS("name");
      XString namS(el->getAttribute(X(nameAttS)));
      XString typS(el->getTagName());
      XString unitAttS("unit");
      XString uniS(el->getAttribute(X(unitAttS)));
      bool foundIt = ((name == 0) || namS.equals(name)) &&
                     ((type == 0) || typS.equals(type)) &&
                     ((unit == 0) || uniS.equals(unit)) ;
      if (foundIt)
      {
         break;
      }
   }
   
   if (p < 0)
   {
      if (fAncestor)
      {
         return fAncestor->get(name, type, number, unit);
      }
      else
      {
         return 0;
      }
   }
   else
   {
      return fListEl[p];
   } 
}

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
