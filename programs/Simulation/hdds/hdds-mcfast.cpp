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

#include <util/PlatformUtils.hpp>
#include <sax/SAXException.hpp>
#include <sax/SAXParseException.hpp>
#include <parsers/DOMParser.hpp>
#include <dom/DOM_DOMException.hpp>

#include "hdds-mcfast.hpp"

#include <assert.h>
#include <fstream.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

char* templateFileList = 0;

void usage()
{
    cerr << "\nUsage:\n"
            "    hdds-mcfast [-v] {HDDS file}\n\n"
            "Options:\n"
            "    -v   validate only\n"
         << endl;
}

void processTemplateFile(char* fname, McfastParameterList& pars)
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
         const DOM_Element* el;
         char* value;
         if (dim == 0)
         {
            el = pars.get(var,"int");
            value = el->getAttribute("value").transcode();
         }
         else
         {
            el = pars.get(var,"int_array",dim);
            value = el->getAttribute("values").transcode();
         }
         if (el)
         {
            cout << " " << value;
         }
         else
         {
            cerr << "hdds-mcfast: Parameter \"" << var << "\""
                 << " required in template file " << fname << endl
                 << "is missing from HDDS" << endl;
            exit(3);
         }
         delete [] value;
      }
      else if (strcasecmp(token,"real") == 0)
      {
         char* var = strtok(0,", (");
         char* arg = strtok(0,")");
         int   dim = (arg == 0) ? 0 : atoi(arg);
         const DOM_Element* el;
         char* value;
         if (dim == 0)
         {
            el = pars.get(var,"real");
            value = el->getAttribute("value").transcode();
         }
         else
         {
            el = pars.get(var,"real_array",dim);
            value = el->getAttribute("values").transcode();
         }
         if (el)
         {
            cout << " " << value;
         }
         else
         {
            cerr << "hdds-mcfast: Parameter \"" << var << "\""
                 << " required in template file " << fname << endl
                 << "is missing from HDDS" << endl;
            exit(3);
         }
         delete [] value;
      }
      else if (strcasecmp(token,"char") == 0)
      {
         char* var = strtok(0,", (");
         char* arg = strtok(0,")");
         int   dim = (arg == 0) ? 0 : atoi(arg);
         const DOM_Element* el;
         if (dim == 0)
         {
            if (! (el = pars.get(var,"string")) )
            {
               cerr << "hdds-mcfast: Parameter \"" << var << "\""
                    << " required in template file " << fname << endl
                    << "is missing from HDDS" << endl;
               exit(3);
            }
            char* value = el->getAttribute("value").transcode();
            cout << " \"" << value << "\"";
            delete [] value;
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
            DOM_Node vect;
            for ( vect = el->getFirstChild(); 
                  vect != 0;
                  vect = vect.getNextSibling() )
            {
               if (vect.getNodeType() != ELEMENT_NODE) continue;
               DOM_Element vectEl = (DOM_Element&) vect;
               char* tag = vectEl.getTagName().transcode();
               if (strcmp(tag,"string_data") == 0)
               {
                  char* value = vectEl.getAttribute("value").transcode();
                  cout << " \"" << value << "\"";
                  delete [] value;
                  vcount++;
               }
            }
            if (vcount != dim)
            {
               char* vname = el->getAttribute("name").transcode();
               cerr << "hdds-mcfast: String vector "
                    << vname << " has too few elements, "
                    << dim << " required." << endl;
               delete [] vname;
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
         const DOM_Element* el;
         if (dim == 0)
         {
            if (! (el = pars.get(var,"reference")) )
            {
               cerr << "hdds-mcfast: Parameter \"" << var << "\""
                    << " required in template file " << fname << endl
                    << "is missing from HDDS" << endl;
               exit(3);
            }
            char* value = el->getAttribute("value").transcode();
            cout << " \"" << value << "\"";
            delete [] value;
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
            DOM_Node vect;
            int vcount = 0;
            for ( vect = el->getFirstChild();
                  vect != 0;
                  vect = vect.getNextSibling() )
            {
               if (vect.getNodeType() != ELEMENT_NODE) continue;
               DOM_Element vectEl = (DOM_Element&) vect;
               char* tag = vectEl.getTagName().transcode();
               if (strcmp(tag,"reference_data") == 0)
               {
                  char* value = vectEl.getAttribute("value").transcode();
                  cout << " \"" << value << "\"";
                  delete [] value;
                  vcount++;
               }
            }
            if (vcount < dim)
            {
               char* rname = el->getAttribute("name").transcode();
               cerr << "hdds-mcfast: Reference vector "
                    << rname << " has too few elements, "
                    << dim << " required." << endl;
               delete [] rname;
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

void writeMcfastRecord(const DOM_Element& el, McfastParameterList* pars)
{
   DOMString parBlock = el.getAttribute("parameters");
   DOM_Element parEl = el.getOwnerDocument().getElementById(parBlock);
   DOM_Node parNode;
   int parCount = 0;
   if (parEl != 0)
   {
      parNode = parEl.getFirstChild();
      parCount = parEl.getElementsByTagName("*").getLength();
   }
   McfastParameterList basePars(parCount);
   basePars.inherits(*pars);
   for ( ; parNode != 0; parNode = parNode.getNextSibling() )
   {
      if (parNode.getNodeType() != ELEMENT_NODE) continue;
      DOM_Element* thisEl = (DOM_Element*) &parNode;
      basePars.append(*thisEl);
   }

   parCount = el.getElementsByTagName("*").getLength();
   McfastParameterList myPars(parCount);
   myPars.inherits(basePars);
   for ( parNode = el.getFirstChild();
         parNode != 0;
         parNode = parNode.getNextSibling() )
   {
      if (parNode.getNodeType() != ELEMENT_NODE) continue;
      DOM_Element* thisEl = (DOM_Element*) &parNode;
      char* tag = thisEl->getTagName().transcode();
      if (strcmp(tag,"mcfast") == 0)
      {
         writeMcfastRecord(*thisEl, &myPars);
      }
      else
      {
         myPars.append(*thisEl);
      }
      delete [] tag;
   }

   char* templ = el.getAttribute("template").transcode();
   if (templateFileList == 0)
   {
      templateFileList = new char[99999];	// big enough to forget
      templateFileList[0] = 0;
   }
   if (strstr(templateFileList,templ) == 0)
   {
      strcat(templateFileList,templ);
      cout << "include " << templ << endl;
   }
     
   char* model = el.getAttribute("model").transcode();
   cout << "make " << model;
   processTemplateFile(templ, myPars);
   cout << endl;
   delete [] model;
   delete [] templ;
}

void writeMaterialRecord(DOM_Element& el)
{
   McfastParameterList parList(3);

   DOM_Element nameEl = el.getOwnerDocument().createElement("string");
   nameEl.setAttribute("name","name");
   nameEl.setAttribute("value",el.getAttribute("name"));
   parList.append(nameEl);

   DOM_Element aEl = el.getOwnerDocument().createElement("real");
   aEl.setAttribute("name","a");
   aEl.setAttribute("value",el.getAttribute("a"));
   parList.append(aEl);

   DOM_Element zEl = el.getOwnerDocument().createElement("real");
   zEl.setAttribute("name","z");
   zEl.setAttribute("value",el.getAttribute("z"));
   parList.append(zEl);

   el.setAttribute("model","Material");
   el.setAttribute("template","db/materials.db");
   writeMcfastRecord(el, &parList);
}

void writeMixtureRecord(DOM_Element& el)
{
   DOM_NodeList propList = el.getElementsByTagName("real");
   int propCount = propList.getLength();
   if (propCount > 4)
   {
      writeMaterialRecord(el);
      return;
   }

   McfastParameterList parList(4);

   DOM_Element nameEl = el.getOwnerDocument().createElement("string");
   nameEl.setAttribute("name","name");
   nameEl.setAttribute("value",el.getAttribute("name"));
   parList.append(nameEl);

   DOM_NodeList matList = el.getElementsByTagName("addmaterial");
   int matCount = matList.getLength();
   assert (matCount > 0);
   float fVol[matCount], fSum = 0.;
   DOM_Element matVecEl =
               el.getOwnerDocument().createElement("reference_vector");
   DOM_Element dataEl[5];
   int m;
   for (m = 0; m < matCount; m++)
   {
      DOM_Node mat = matList.item(m);
      DOM_Element matEl = (DOM_Element&) mat;
      DOMString refID = matEl.getAttribute("material");
      dataEl[m] = el.getOwnerDocument().createElement("reference_data");
      dataEl[m].setAttribute("value",refID);
      matVecEl.appendChild(dataEl[m]);

      DOM_Element refEl = matEl.getOwnerDocument().getElementById(refID);
      char* aStr = refEl.getAttribute("a").transcode();
      float a = atof(aStr);
      delete [] aStr;
      char* zStr = refEl.getAttribute("z").transcode();
      float z = atof(zStr);
      delete [] zStr;

      DOM_NodeList parList = refEl.getElementsByTagName("real");
      int parCount = parList.getLength();
      float density = 1.0e30;
      for (int p = 0; p < parCount; p++)
      {
         DOM_Node par = parList.item(p);
         DOM_Element parEl = (DOM_Element&) par;
         char* parName = parEl.getAttribute("name").transcode();
         if (strcmp(parName,"density") == 0)
         {
            char* valStr = parEl.getAttribute("value").transcode();
            density = atof(valStr);
            delete [] valStr;
         }
         delete [] parName;
      }

      DOM_NodeList specList = matEl.getElementsByTagName("*");
      assert (specList.getLength() > 0);
      DOM_Node spec = specList.item(0);
      DOM_Element specEl = (DOM_Element&) spec;
      char* specType = specEl.getTagName().transcode();
      float fMass = 0;
      if (strcmp(specType,"natoms") == 0)
      {
         char* nStr = specEl.getAttribute("n").transcode();
         int n = atoi(nStr);
         delete [] nStr;
         fMass = n*a;
      }
      else if (strcmp(specType,"fractionmass") == 0)
      {
         char* fStr = specEl.getAttribute("fraction").transcode();
         fMass = atof(fStr);
         delete [] fStr;
      }
      delete [] specType;
      fVol[m] = fMass/density;
      fSum += fVol[m];
   }

   for (; m < 5; m++)
   {
      dataEl[m] = el.getOwnerDocument().createElement("reference_data");
      dataEl[m].setAttribute("value","-");
      matVecEl.appendChild(dataEl[m]);
      fVol[m] = 0;
   }

   char nmatStr[10];
   sprintf(nmatStr,"%d",matCount);
   DOM_Element nmatEl = el.getOwnerDocument().createElement("int");
   nmatEl.setAttribute("name","nmat");
   nmatEl.setAttribute("value",nmatStr);
   parList.append(nmatEl);

   matVecEl.setAttribute("name","matnames");
   parList.append(matVecEl);

   char fVolStr[300];
   sprintf(fVolStr,"%f %f %f %f %f", fVol[0]/fSum, fVol[1]/fSum,
           fVol[2]/fSum, fVol[3]/fSum, fVol[4]/fSum);
   DOM_Element fVolEl = el.getOwnerDocument().createElement("real_array");
   fVolEl.setAttribute("name","prop");
   fVolEl.setAttribute("value",fVolStr);
   parList.append(fVolEl);

   el.setAttribute("model","Mixture");
   el.setAttribute("template","db/mixtures.db");
   writeMcfastRecord(el, &parList);
}

void doMaterials(DOM_Element& rootEl)
{
   DOM_NodeList matList = rootEl.getElementsByTagName("element");
   int matCount = matList.getLength();
   for (int m = 0; m < matCount; m++)
   {
      DOM_Node mat = matList.item(m);
      DOM_Element matEl = (DOM_Element&) mat;
      int nRealPars = matEl.getElementsByTagName("real").getLength();
      if (nRealPars > 4)
      {
         writeMaterialRecord(matEl);
      }
   }

   DOM_NodeList compList = rootEl.getElementsByTagName("composite");
   int compCount = compList.getLength();
   for (int c = 0; c < compCount; c++)
   {
      DOM_Node comp = compList.item(c);
      DOM_Element compEl = (DOM_Element&) comp;
      int nRealPars = compEl.getElementsByTagName("real").getLength();
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
      cerr << "hdds-mcfast: Error during initialization! :\n"
           << StrX(toCatch.getMessage()) << endl;
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

   DOMParser parser;
   parser.setValidationScheme(DOMParser::Val_Auto);
   parser.setCreateEntityReferenceNodes(false);
   parser.setDoNamespaces(false);

   MyOwnErrorHandler errorHandler;
   parser.setErrorHandler(&errorHandler);

   try
   {
      parser.parse(xmlFile);
   }
   catch (const XMLException& toCatch)
   {
      cerr << "\nhdds-mcfast: Error during parsing: '" << xmlFile << "'\n"
           << "Exception message is:  \n"
           << StrX(toCatch.getMessage()) << "\n" << endl;
      return -1;
   }
   catch (const DOM_DOMException& toCatch)
   {
      cerr << "\nhdds-mcfast: Error during parsing: '" << xmlFile << "'\n"
           << "Exception message is:  \n"
           << toCatch.msg.transcode() << "\n" << endl;
      XMLPlatformUtils::Terminate();
      return 4;
   }
   catch (...)
   {
      cerr << "\nhdds-mcfast: Unexpected exception during parsing: '"
           << xmlFile << "'\n";
      XMLPlatformUtils::Terminate();
      return 4;
   }

   if (errorHandler.getSawErrors())
   {
      cerr << "\nErrors occured, no output available\n" << endl;
   }

   DOM_Document doc = parser.getDocument();
   if (! mcfastOutput)
   {
      return 0;
   }

   cout << "database mcfast 0000" << endl;

   DOM_Element rootEl = doc.getDocumentElement();
   DOM_Node sect;
   for ( sect = rootEl.getLastChild();
         sect != 0;
         sect = sect.getPreviousSibling() )
   {
      if (sect.getNodeType() != ELEMENT_NODE) continue;
      DOM_Element sectEl = (DOM_Element&) sect;
      char* stag = sectEl.getTagName().transcode();
      if (strcmp(stag,"section") == 0)
      {
         DOM_Node cont;
         for ( cont = sect.getFirstChild();
               cont != 0;
               cont = cont.getNextSibling() )
         {
            if (cont.getNodeType() != ELEMENT_NODE) continue;
            DOM_Element contEl = (DOM_Element&) cont;
            char* tag = contEl.getTagName().transcode();
            if (strcmp(tag,"mcfast") == 0)
            {
               static bool firstTime = true;
               writeMcfastRecord(contEl, 0);
               if (firstTime)
               {
                  doMaterials(rootEl);
                  firstTime = false;
               }
            }
            delete [] tag;
         }
      }
      delete [] stag;
   }

   XMLPlatformUtils::Terminate();
   return 0;
}


MyOwnErrorHandler::MyOwnErrorHandler() : 
   fSawErrors(false)
{
}

MyOwnErrorHandler::~MyOwnErrorHandler()
{
}

/* Overrides of the SAX ErrorHandler interface */

void MyOwnErrorHandler::error(const SAXParseException& e)
{
   fSawErrors = true;
   cerr << "\nhdds-mcfast: Error at file " << StrX(e.getSystemId())
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << StrX(e.getMessage()) << endl;
}

void MyOwnErrorHandler::fatalError(const SAXParseException& e)
{
   fSawErrors = true;
   cerr << "\nhdds-mcfast: Fatal Error at file " << StrX(e.getSystemId())
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << StrX(e.getMessage()) << endl;
}

void MyOwnErrorHandler::warning(const SAXParseException& e)
{
   cerr << "\nhdds-mcfast: Warning at file " << StrX(e.getSystemId())
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << StrX(e.getMessage()) << endl;
}

void MyOwnErrorHandler::resetErrors()
{
}

McfastParameterList::McfastParameterList(const int capacity)
{
   fListCap = (capacity > 0) ? capacity : 1;
   fListEl = new const DOM_Element* [fListCap];
   fAncestor = 0;
   fListLen = 0;
}

McfastParameterList::~McfastParameterList()
{
   for (int i=0; i<fListLen; i++)
   {
      delete fListEl[i];
   }
   delete [] fListEl;
}

void McfastParameterList::inherits(McfastParameterList& anc)
{
   fAncestor = &anc;
}

void McfastParameterList::append(const DOM_Element& par)
{
   DOM_Element* myCopy = new DOM_Element(par);
   if (fListLen < fListCap)
   {
      fListEl[fListLen++] = myCopy;
   }
   else
   {
      cerr << "hdds-mcfast: Error in McfastParameterList::append" << endl
           << "   list overflow, insufficient capacity" << endl;
      exit(1);
   }
}

const DOM_Element* McfastParameterList::get(char* name, char* type,
                                            int number=0, char* unit=0) const
{
   int p = fListLen;
   while (--p >= 0)
   {
      const DOM_Element* el = fListEl[p];
      char* nam = el->getAttribute("name").transcode();
      char* typ = el->getTagName().transcode();
      char* uni = el->getAttribute("unit").transcode();
      bool foundIt = ((name == 0) || (strcmp(name,nam) == 0)) &&
                     ((type == 0) || (strcmp(type,typ) == 0)) &&
                     ((unit == 0) || (strcmp(unit,uni) == 0)) ;
      delete [] nam;
      delete [] typ;
      delete [] uni;
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
