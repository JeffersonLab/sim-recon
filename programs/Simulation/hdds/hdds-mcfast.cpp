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
#include <sstream>

#define X(XString) XString.unicodeForm()
#define S(XString) XString.localForm()

#ifdef _Tru64
#undef basename
#define basename _RC_basename
static char* basename(const char *f)
{
   const char* base;
   for (base = f; *f; f++)
   {
      if (*f == '/')
      base = f+1;
   }
   return (char *)base;
}
#endif

struct _modelTableEntry
{
   char* model;
   int refcount;
   int maxcount;
   ostringstream db;
} modelTable[999];
int modelTableLen=0;

makeTargetTable targetTable(999);


void usage()
{
    cerr << "\nUsage:\n"
            "    hdds-mcfast [-v] {HDDS file}\n\n"
            "Options:\n"
            "    -v   validate only\n"
         << endl;
}

void processTemplateFile(const DOMElement* const targetEl,
                         const char* const fname)
{
   XString* processingTarget=0;
   int processedTemplates=0;

   const XString modelAttS("model");
   const XString modelS = targetEl->getAttribute(X(modelAttS));
   int model;
   for (model=0; model < modelTableLen; model++)
   {
      if (modelS.equals(modelTable[model].model))
         break;
   }
   if (model == modelTableLen)
   {
      char* str = new char[strlen(S(modelS))+1];
      strcpy(str,S(modelS));
      modelTable[model].model = str;
      modelTable[model].refcount = 0;
      modelTable[model].maxcount = 0;
      modelTableLen++;
   }
   if (modelTable[model].refcount == 0)
   {
      modelTable[model].db << "include " << fname << endl;
   }

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
      char token[250];
      istringstream sline(line);
      sline >> token;
      if (strcasecmp(token,"end") == 0)
      {
         if (processingTarget)
         {
            delete processingTarget;
            processingTarget = 0;
         }
         else
         {
            cerr << "hdds-mcfast: end statement without matching template"
                 << " in input file " << fname << endl;
            exit(2);
         }
         modelTable[model].db << endl;
      }
      else if (strcasecmp(token,"include") == 0)
      {
         char templFile[250];
         sline >> templFile;
         processTemplateFile(targetEl,templFile);
      }
      else if (strcasecmp(token,"template") == 0)
      {
         char tgt[250];
         sline >> tgt;
         processingTarget = new XString(strtok(tgt,", ("));
         const char* maxcount = strtok(0,")");
         int dim = (maxcount == 0) ? 0 : atoi(maxcount);
         if (dim == 0)
         {
            cerr << "hdds-mcfast: bad format in template statement"
                 << " in input file " << fname << endl;
            exit(2);
         }
         else if (processingTarget && processingTarget->equals(modelS))
         {
            if (++modelTable[model].refcount > dim)
            { 
               cerr << "hdds-mcfast: number of objects of type " << S(modelS)
                    << " exceeds maximum of " << dim << endl;
               cerr << "defined in input file " << fname << endl;
               cerr << "Increase the array size in the template statement"
                    << " found in the above file" << endl;
               cerr << "and try again." << endl;
               exit(2);
            }
            modelTable[model].maxcount = dim;
            modelTable[model].db << "make " << S(modelS);
            ++processedTemplates;
         }
      }
      else if (strcasecmp(token,"make") == 0)
      {
         char tgt[250];
         sline >> tgt;
         int m;
         for (m=0; m < modelTableLen; m++)
         {
            if (strcmp(modelTable[m].model,tgt) == 0)
               break;
         }
         if (m == modelTableLen)
         { 
            cerr << "hdds-mcfast: error in template file " << fname << endl;
            cerr << "Statement template " << S(modelS)
                 << " must appear before first make instance." << endl;
            exit(2);
         }
         else if (modelTable[m].refcount > modelTable[m].maxcount)
         {
            cerr << "hdds-mcfast: number of objects of type " << tgt
                 << " overflows table." << endl;
            cerr << "Increase the array size in the template file"
                 << " and try again." << endl;
            exit(2);
         }
      }
      else if (!(processingTarget && processingTarget->equals(modelS)))
      {
         continue;
      }
      else if (strcasecmp(token,"int") == 0)
      {
         char tgt[250];
         sline >> tgt;
         const char* var = strtok(tgt,", (");
         const char* arg = strtok(0,")");
         int dim = (arg == 0) ? 0 : atoi(arg);
         const DOMElement* el;
         XString valueS;
         if (dim == 0)
         {
            const XString typeS("int");
            const XString varS(var);
            if (el = targetTable.lookup(targetEl,typeS,varS))
            {
	       const XString valAttS("value");
               valueS = el->getAttribute(X(valAttS));
            }
         }
         else
         {
            const XString typeS("int_array");
            const XString varS(var);
            if (el = targetTable.lookup(targetEl,typeS,varS))
            {
	       const XString valAttS("values");
               valueS = el->getAttribute(X(valAttS));
            }
         }
         if (valueS == 0)
         {
            cerr << "hdds-mcfast: Parameter \"" << var << "\""
                 << " required in template file " << fname << endl
                 << "is missing from HDDS" << endl;
            exit(3);
         }
         char* intstr = strtok((char*)S(valueS)," ");
         modelTable[model].db << " " << intstr;
         for (int i=1; i < dim; i++)
         {
            if (intstr = strtok(0," "))
            {
               modelTable[model].db << " " << intstr;
            }
            else
            {
               modelTable[model].db << " 0";
            }
         }
      }
      else if (strcasecmp(token,"real") == 0)
      {
         char tgt[250];
         sline >> tgt;
         const char* var = strtok(tgt,", (");
         const char* arg = strtok(0,")");
         int dim = (arg == 0) ? 0 : atoi(arg);
         const DOMElement* el;
         XString valueS;
         if (dim == 0)
         {
            const XString typeS("real");
            const XString varS(var);
            if (el = targetTable.lookup(targetEl,typeS,varS))
            {
	       const XString valAttS("value");
               valueS = el->getAttribute(X(valAttS));
            }
         }
         else
         {
            const XString typeS("real_array");
            const XString varS(var);
            if (el = targetTable.lookup(targetEl,typeS,varS))
            {
	       const XString valAttS("values");
               valueS = el->getAttribute(X(valAttS));
            }
         }
         if (valueS == 0)
         {
            cerr << "hdds-mcfast: Parameter \"" << var << "\""
                 << " required in template file " << fname << endl
                 << "is missing from HDDS" << endl;
            exit(3);
         }
         const char* fltstr = strtok((char*)S(valueS)," ");
         modelTable[model].db << " " << fltstr;
         for (int i=1; i < dim; i++)
         {
            if (fltstr = strtok(0," "))
            {
               modelTable[model].db << " " << fltstr;
            }
            else
            {
               modelTable[model].db << " 0";
            }
         }
      }
      else if (strcasecmp(token,"char") == 0)
      {
         char tgt[250];
         sline >> tgt;
         const char* var = strtok(tgt,", (");
         const char* arg = strtok(0,")");
         int dim = (arg == 0) ? 0 : atoi(arg);
         const DOMElement* el;
         XString valueS;
         if (dim == 0)
         {
            const XString typeS("string");
            const XString varS(var);
            if ((el = targetTable.lookup(targetEl,typeS,varS)) == 0)
            {
               cerr << "hdds-mcfast: Parameter \"" << var << "\""
                    << " required in template file " << fname << endl
                    << "is missing from HDDS" << endl;
               exit(3);
            }
	    const XString valAttS("value");
            valueS = el->getAttribute(X(valAttS));
            modelTable[model].db << " \"" << S(valueS) << "\"";
         }
         else
         {
            const XString typeS("string_vector");
            const XString varS(var);
            if ((el = targetTable.lookup(targetEl,typeS,varS)) == 0)
            {
               cerr << "hdds-mcfast: Parameter \"" << var << "\""
                    << " required in template file " << fname << endl
                    << "is missing from HDDS" << endl;
               exit(3);
            }
            int vcount = 0;
            const DOMNode* vect;
            for ( vect = el->getFirstChild(); 
                  vect != 0;
                  vect = vect->getNextSibling() )
            {
               if (vect->getNodeType() != DOMNode::ELEMENT_NODE)
                  continue;
               const DOMElement* vectEl = (DOMElement*) vect;
               const XString tagS(vectEl->getTagName());
               if (tagS.equals("string_data"))
               {
	          const XString valAttS("value");
                  const XString valueS(vectEl->getAttribute(X(valAttS)));
                  modelTable[model].db << " \"" << S(valueS) << "\"";
                  vcount++;
               }
            }
            for ( ; vcount < dim; vcount++)
            {
               modelTable[model].db << " \"\"";
            }
            if (vcount != dim)
            {
	       const XString nameAttS("name");
               const XString vnameS(el->getAttribute(X(nameAttS)));
               cerr << "hdds-mcfast: mcfast array size of " << dim
                    << " for variable " << S(vnameS) << " is too small"
                    << " to hold " << vcount << " elements" << endl;
               cerr << "Please increase array size in " << fname << endl;
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
         char tgt[250];
         sline >> tgt;
         const char* var = strtok(tgt,", (");
         const char* arg = strtok(0,")");
         int dim = (arg == 0) ? 0 : atoi(arg);
         const DOMElement* el;
         if (dim == 0)
         {
            const XString typeS("reference");
            const XString varS(var);
            if ((el = targetTable.lookup(targetEl,typeS,varS)) == 0)
            {
               cerr << "hdds-mcfast: Parameter \"" << var << "\""
                    << " required in template file " << fname << endl
                    << "is missing from HDDS" << endl;
               exit(3);
            }
	    const XString valAttS("value");
            const XString valueS(el->getAttribute(X(valAttS)));
            modelTable[model].db << " \"" << S(valueS) << "\"";
         }
         else
         {
            const XString typeS("reference_vector");
            const XString varS(var);
            if ((el = targetTable.lookup(targetEl,typeS,varS)) == 0)
            {
               cerr << "hdds-mcfast: Parameter \"" << var << "\""
                    << " required in template file " << fname << endl
                    << "is missing from HDDS" << endl;
               exit(3);
            }
            const DOMNode* vect;
            int vcount = 0;
            for ( vect = el->getFirstChild();
                  vect != 0;
                  vect = vect->getNextSibling() )
            {
               if (vect->getNodeType() != DOMNode::ELEMENT_NODE)
                  continue;
               const DOMElement* vectEl = (DOMElement*) vect;
               const XString tagS(vectEl->getTagName());
               if (tagS.equals("reference_data"))
               {
	          const XString valAttS("value");
                  const XString valueS(vectEl->getAttribute(X(valAttS)));
                  modelTable[model].db << " \"" << S(valueS) << "\"";
                  vcount++;
               }
            }
            for ( ; vcount < dim; vcount++)
            {
               modelTable[model].db << " \"-\"";
            }
            if (vcount != dim)
            {
	       const XString nameAttS("name");
               const XString vnameS(el->getAttribute(X(nameAttS)));
               cerr << "hdds-mcfast: mcfast array size of " << dim
                    << " for variable " << S(vnameS) << " is too small"
                    << " to hold " << vcount << " elements" << endl;
               cerr << "Please increase array size in " << fname << endl;
               exit(3);
            }
         }
      }
      else if (strstr(token,"!") == token)
      {
         continue;
      }
      else
      {
         cerr << "hdds-mcfast: Template file " << fname
              << " contains unknown parameter type " << token << endl;
         exit(3);
      }
   }
   if (processedTemplates == 0)
   {
      cerr << "hdds-mcfast: template for " << S(modelS)
           << " not found in input file " << fname << endl;
      exit(2);
   }
   else if (processingTarget)
   {
      cerr << "hdds-mcfast: template statement without matching end"
           << " in input file " << fname << endl;
      exit(2);
   }
}

void makedb(DOMElement* el)
{
   DOMNode* cont;
   for (cont = el->getLastChild();
        cont != 0;
        cont = cont->getPreviousSibling())
   {
      if (cont->getNodeType() != DOMNode::ELEMENT_NODE)
         continue;
      DOMElement* contEl = (DOMElement*)cont;
      const XString tagS = contEl->getTagName();
      if (tagS.equals("mcfast"))
      {
         targetTable.add(contEl);
      }
      else
      {
         makedb(contEl);
      }
   }
}

void printdb()
{
   char line[999];
   for (int m=0; m < modelTableLen; m++)
   {
      istringstream idb(modelTable[m].db.str());
      while (! idb.eof())
      {
         idb.getline(line,999);
         cout << line << endl;
      }
   }
   cout << "end" << endl;
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
   DOMDocument* doc = buildDOMDocument(xmlFile,false);
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

   int model=modelTableLen++;
   char hdrStr[] = "database";             // header line must come first
   modelTable[model].model = hdrStr;
   modelTable[model].refcount = 1;
   modelTable[model].maxcount = 1;
   modelTable[model].db << "database mcfast 0000" << endl;

   model=modelTableLen++;
   char detStr[] = "detector";             // detector is first declaration
   modelTable[model].model = detStr;
   modelTable[model].refcount = 0;
   modelTable[model].maxcount = 0;

   model=modelTableLen++;
   char matStr[] = "Material";             // force materials next
   modelTable[model].model = matStr;
   modelTable[model].refcount = 0;
   modelTable[model].maxcount = 0;

   model=modelTableLen++;
   char mixStr[] = "Mixture";              // then mixtures
   modelTable[model].model = mixStr;
   modelTable[model].refcount = 0;
   modelTable[model].maxcount = 0;

   DOMElement* rootEl = doc->getDocumentElement();
   makedb(rootEl);

   model=modelTableLen++;
   char hitsStr[] = "hitsontrack";         // last comes histontrack
   modelTable[model].model = hitsStr;
   modelTable[model].refcount = 1;
   modelTable[model].maxcount = 1;
   modelTable[model].db << "include db/hitsontrack.db" << endl
                        << "make HitsOnTrack 4 0 0" << endl;
   printdb();

   XMLPlatformUtils::Terminate();
   return 0;
}


makeTargetTable::makeTargetTable(const int capacity)
{
   if (capacity > 0) {
      fTargetEl = new DOMElement*[capacity];
      fParentEl = new DOMElement*[capacity];
      fTableSize = capacity;
   }
   else {
      fTableSize = 0;
   }
   fTableLen = 0;
}

makeTargetTable::~makeTargetTable()
{
   delete [] fTargetEl;
   delete [] fParentEl;
}

int makeTargetTable::add(DOMElement* const targetEl,
                         DOMElement* const parentEl)
{
   if (targetEl == 0)
   {
      cerr << "makeTargetTable::add - called with null pointer" << endl;
      return 0;
   }
   int i;
   for (i=0; i < fTableLen; i++)
   {
      if (fTargetEl[i] == targetEl)
         break;
   }
   if (i >= fTableSize) {
      cerr << "makeTargetTable::add - internal table overflow, "
           << "please increase table size." << endl;
      return 0;
   }
   else if (i == fTableLen)
   {
      fTargetEl[i] = targetEl;
      fParentEl[i] = parentEl;
      ++fTableLen;
      const XString tagS = targetEl->getTagName();
      if (tagS.equals("mcfast"))
      {
         addModel(targetEl);
      }
      else if (tagS.equals("element"))
      {
         addElement(targetEl);
      }
      else if (tagS.equals("composite"))
      {
         addComposite(targetEl);
      }
      else
      {
         cerr << "makeTargetTable::add - called with unexpected tag "
              << S(tagS) << endl;
         return 0;
      }
   }
}

int makeTargetTable::addElement(DOMElement* const targetEl)
{
   set_name(targetEl);
   set_a(targetEl);
   set_z(targetEl);
   const XString modelAttS("model");
   const XString materialS("Material");
   targetEl->setAttribute(X(modelAttS),X(materialS));
   const XString templAttS("template");
   const XString dbS("db/materials.db");
   targetEl->setAttribute(X(templAttS),X(dbS));
   processTemplateFile(targetEl,S(dbS));
   return 0;
}

int makeTargetTable::addComposite(DOMElement* const targetEl)
{
   const XString densityS("density");
   const XString radlenS("radlen");
   const XString collenS("collen");
   const XString abslenS("abslen");
   const XString dedxS("dedx");
   const XString realS("real");
   if (lookup(targetEl,realS,densityS) ||
       lookup(targetEl,realS,radlenS) ||
       lookup(targetEl,realS,collenS) ||
       lookup(targetEl,realS,abslenS) ||
       lookup(targetEl,realS,dedxS))
   {
      set_density(targetEl);
      set_radlen(targetEl);
      set_abslen(targetEl);
      set_collen(targetEl);
      set_dedx(targetEl);
      return addElement(targetEl);
   }
   else
   {
      set_name(targetEl);
      set_a(targetEl);
      set_z(targetEl);
      set_density(targetEl);
   }

   const XString addmatS("addmaterial");
   DOMNodeList* matList = targetEl->getElementsByTagName(X(addmatS));
   const int matCount = matList->getLength();
   char nmatStr[10];
   sprintf(nmatStr,"%d",matCount);
   const XString nmatS(nmatStr);
   const XString nmatAttS("nmat");
   const XString intS("int");
   DOMElement* nmatEl = targetEl->getOwnerDocument()->createElement(X(intS));
   const XString nameAttS("name");
   const XString valueAttS("value");
   nmatEl->setAttribute(X(nameAttS),X(nmatAttS));
   nmatEl->setAttribute(X(valueAttS),X(nmatS));
   targetEl->appendChild(nmatEl);

   float fVol[matCount];
   DOMElement* dataEl[matCount];
   const XString refvecS("reference_vector");
   DOMElement* matVecEl = targetEl->getOwnerDocument()->createElement(X(refvecS));
   const XString matnameS("matnames");
   matVecEl->setAttribute(X(nameAttS),X(matnameS));
   for (int m=0; m < matCount; m++)
   {
      DOMElement* matEl = (DOMElement*)matList->item(m);
      const XString materialS("material");
      const XString refIdS(matEl->getAttribute(X(materialS)));
      const XString refdataS("reference_data");
      dataEl[m] = targetEl->getOwnerDocument()->createElement(X(refdataS));
      const XString valueAttS("value");
      dataEl[m]->setAttribute(X(valueAttS),X(refIdS));
      matVecEl->appendChild(dataEl[m]);
      DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refIdS));
      add(refEl);
      const XString fractionS("fractionmass");
      DOMNodeList* specList = matEl->getElementsByTagName(X(fractionS));
      if (specList->getLength() == 0)
      {
         const XString nameAttS("name");
         const XString nameS = matEl->getAttribute(X(nameAttS));
         cerr << "makeTargetTable::addComposite - "
              << "composite " << S(nameS) << " is missing information."
              << endl;
         exit(2);
      }
      DOMElement* specEl = (DOMElement*)specList->item(0);
      const XString fracAttS("fraction");
      const XString fS(specEl->getAttribute(X(fracAttS)));
      float fMass_m = atof(S(fS));
      float density_m = set_density(refEl);
      float density = set_density(targetEl);
      fVol[m] = fMass_m*(density/density_m);
   }
   targetEl->appendChild(matVecEl);

   const XString rarrayS("real_array");
   DOMElement* fVolEl = targetEl->getOwnerDocument()->createElement(X(rarrayS));
   const XString propS("prop");
   fVolEl->setAttribute(X(nameAttS),X(propS));
   char str[20];
   sprintf(str,"%f",fVol[0]);
   XString fVolS(str);
   for (int m=1; m < matCount; m++)
   {
      sprintf(str," %f",fVol[m]);
      const XString xS(str);
      fVolS += xS;
   }
   const XString valuesAttS("values");
   fVolEl->setAttribute(X(valuesAttS),X(fVolS));
   targetEl->appendChild(fVolEl);

   const XString modelAttS("model");
   const XString templAttS("template");
   const XString mixtureS("Mixture");
   const XString dbS("db/mixtures.db");
   targetEl->setAttribute(X(modelAttS),X(mixtureS));
   targetEl->setAttribute(X(templAttS),X(dbS));
   processTemplateFile(targetEl,S(dbS));
   return 0;
}

int makeTargetTable::addModel(DOMElement* const targetEl)
{
   const XString paramAttS("parameters");
   const XString parS = targetEl->getAttribute(X(paramAttS));
   if (parS != 0)
   {
      DOMElement* parEl = targetEl->getOwnerDocument()->getElementById(X(parS));
      for (DOMNode* var = parEl->getFirstChild(); 
           var != 0;
           var = var->getNextSibling() )
      {
         if (var->getNodeType() != DOMNode::ELEMENT_NODE)
            continue;
         DOMElement* varEl = (DOMElement*)var;
         const XString refS("reference");
         const XString varS = varEl->getTagName();
         if (varS.equals(refS))
         {
            const XString valueAttrS("value");
            const XString valueS = varEl->getAttribute(X(valueAttrS));
            DOMElement* matEl = varEl->getOwnerDocument()->getElementById(X(valueS));
            add(matEl);
         }
      }
   }
   for (DOMNode* var = targetEl->getFirstChild(); 
        var != 0;
        var = var->getNextSibling() )
   {
      if (var->getNodeType() != DOMNode::ELEMENT_NODE)
         continue;
      DOMElement* varEl = (DOMElement*)var;
      const XString varS = varEl->getTagName();
      if (varS.equals("reference"))
      {
         const XString valueAttrS("value");
         const XString valueS = varEl->getAttribute(X(valueAttrS));
         DOMElement* matEl = varEl->getOwnerDocument()->getElementById(X(valueS));
         add(matEl);
      }
      else if (varS.equals("mcfast"))
      {
         add(varEl,targetEl);
      }
   }
   const XString templAttS("template");
   const XString templS = targetEl->getAttribute(X(templAttS));
   processTemplateFile(targetEl,S(templS));
   return 0;
}

void makeTargetTable::set_name(DOMElement* const targetEl)
{
   const XString nameAttS("name");
   const XString valueAttS("value");
   const XString stringS("string");
   DOMElement* nameEl = lookup(targetEl,stringS,nameAttS);
   if (nameEl == 0)
   {
      const XString nameS = targetEl->getAttribute(X(nameAttS));
      nameEl = targetEl->getOwnerDocument()->createElement(X(stringS));
      nameEl->setAttribute(X(nameAttS),X(nameAttS));
      nameEl->setAttribute(X(valueAttS),X(nameS));
      targetEl->appendChild(nameEl);
   }
}

float makeTargetTable::set_a(DOMElement* const targetEl)
{
   const XString nameAttS("name");
   const XString valueAttS("value");
   const XString realS("real");
   const XString aAttS("a");
   XString aS = targetEl->getAttribute(X(aAttS));
   DOMElement* aEl = lookup(targetEl,realS,aAttS);
   if (aEl == 0)
   {
      if (aS == 0)
      {
         const XString addmatS("addmaterial");
         DOMNodeList* matList = targetEl->getElementsByTagName(X(addmatS));
         int matCount = matList->getLength();
         float a_sum=0;
         float wgt_sum=0;
         for (int m=0; m < matCount; m++)
         {
            DOMElement* matEl = (DOMElement*)matList->item(m);
            const XString materialS("material");
            const XString refS = matEl->getAttribute(X(materialS));
            DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refS));
            float a_m = set_a(refEl);
            float wgt_m=0;
            const XString natomsS("natoms");
            DOMNodeList* wgtList = matEl->getElementsByTagName(X(natomsS));
            if (wgtList->getLength())
            {
               DOMElement* natomsEl = (DOMElement*)wgtList->item(0);
               const XString nAttS("n");
               const XString nS = natomsEl->getAttribute(X(nAttS));
               int n = atoi(S(nS));
               wgt_m = n*a_m;
            }
            const XString fractionS("fractionmass");
            wgtList = matEl->getElementsByTagName(X(fractionS));
            if (wgtList->getLength())
            {
               DOMElement* fractionEl = (DOMElement*)wgtList->item(0);
               const XString fractionAttS("fraction");
               const XString fracS = fractionEl->getAttribute(X(fractionAttS));
               wgt_m = atof(S(fracS));
            }
            a_sum += a_m*wgt_m;
            wgt_sum += wgt_m;
         }
         char aStr[20];
         sprintf(aStr,"%g",a_sum/wgt_sum);
         const XString aStrS(aStr);
         aS = aStrS;
      }
      DOMElement* aEl = targetEl->getOwnerDocument()->createElement(X(realS));
      aEl->setAttribute(X(nameAttS),X(aAttS));
      aEl->setAttribute(X(valueAttS),X(aS));
      targetEl->appendChild(aEl);
   }
   return atof(S(aS));
}

float makeTargetTable::set_z(DOMElement* const targetEl)
{
   const XString nameAttS("name");
   const XString valueAttS("value");
   const XString realS("real");
   const XString zAttS("z");
   XString zS = targetEl->getAttribute(X(zAttS));
   DOMElement* zEl = lookup(targetEl,realS,zAttS);
   if (zEl == 0)
   {
      if (zS == 0)
      {
         const XString addmatS("addmaterial");
         DOMNodeList* matList = targetEl->getElementsByTagName(X(addmatS));
         int matCount = matList->getLength();
         float z_sum=0;
         float wgt_sum=0;
         for (int m=0; m < matCount; m++)
         {
            DOMElement* matEl = (DOMElement*)matList->item(m);
            const XString materialS("material");
            const XString refS = matEl->getAttribute(X(materialS));
            DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refS));
            float z_m = set_z(refEl);
            float wgt_m=0;
            const XString natomsS("natoms");
            DOMNodeList* wgtList = matEl->getElementsByTagName(X(natomsS));
            if (wgtList->getLength())
            {
               DOMElement* natomsEl = (DOMElement*)wgtList->item(0);
               const XString nAttS("n");
               const XString nS = natomsEl->getAttribute(X(nAttS));
               int n = atoi(S(nS));
               const XString materialS("material");
               const XString refIdS = matEl->getAttribute(X(materialS));
               DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refIdS));
               wgt_m = n*set_a(refEl);
            }
            const XString fractionS("fractionmass");
            wgtList = matEl->getElementsByTagName(X(fractionS));
            if (wgtList->getLength())
            {
               DOMElement* fractionEl = (DOMElement*)wgtList->item(0);
               const XString fractionAttS("fraction");
               const XString fracS = fractionEl->getAttribute(X(fractionAttS));
               wgt_m = atof(S(fracS));
            }
            z_sum += z_m*wgt_m;
            wgt_sum += wgt_m;
         }
         char zStr[20];
         sprintf(zStr,"%g",z_sum/wgt_sum);
         const XString zStrS(zStr);
         zS = zStrS;
      }
      DOMElement* zEl = targetEl->getOwnerDocument()->createElement(X(realS));
      zEl->setAttribute(X(nameAttS),X(zAttS));
      zEl->setAttribute(X(valueAttS),X(zS));
      targetEl->appendChild(zEl);
   }
   return atof(S(zS));
}

float makeTargetTable::set_density(DOMElement* const targetEl)
{
   const XString nameAttS("name");
   const XString valueAttS("value");
   const XString realS("real");
   const XString densityS("density");
   DOMElement* densityEl = lookup(targetEl,realS,densityS);
   if (densityEl == 0)
   {
      const XString addmatS("addmaterial");
      DOMNodeList* matList = targetEl->getElementsByTagName(X(addmatS));
      int matCount = matList->getLength();
      float ohr_sum=0;
      float wgt_sum=0;
      for (int m=0; m < matCount; m++)
      {
         DOMElement* matEl = (DOMElement*)matList->item(m);
         const XString materialS("material");
         const XString refS = matEl->getAttribute(X(materialS));
         DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refS));
         float rho_m = set_density(refEl);
         float wgt_m=0;
         const XString natomsS("natoms");
         DOMNodeList* wgtList = matEl->getElementsByTagName(X(natomsS));
         if (wgtList->getLength())
         {
            const XString nameS = matEl->getAttribute(X(nameAttS));
            cerr << "makeTargetTable::set_density - "
                 << "no automatic density calculation for atomic mixture "
                 << S(nameS) << endl;
            exit(2);
         }
         const XString fractionS("fractionmass");
         wgtList = matEl->getElementsByTagName(X(fractionS));
         if (wgtList->getLength())
         {
            DOMElement* fractionEl = (DOMElement*)wgtList->item(0);
            const XString fractionAttS("fraction");
            const XString fracS = fractionEl->getAttribute(X(fractionAttS));
            wgt_m = atof(S(fracS));
         }
         ohr_sum += wgt_m/rho_m;
         wgt_sum += wgt_m;
      }
      char rhoStr[20];
      sprintf(rhoStr,"%g",wgt_sum/ohr_sum);
      const XString rhoS(rhoStr);
      densityEl = targetEl->getOwnerDocument()->createElement(X(realS));
      densityEl->setAttribute(X(nameAttS),X(densityS));
      densityEl->setAttribute(X(valueAttS),X(rhoS));
      targetEl->appendChild(densityEl);
   }
   const XString resultS = densityEl->getAttribute(X(valueAttS));
   return atof(S(resultS));
}

float makeTargetTable::set_radlen(DOMElement* const targetEl)
{
   const XString nameAttS("name");
   const XString valueAttS("value");
   const XString realS("real");
   const XString radlenS("radlen");
   DOMElement* radlenEl = lookup(targetEl,realS,radlenS);
   if (radlenEl == 0)
   {
      const XString addmatS("addmaterial");
      DOMNodeList* matList = targetEl->getElementsByTagName(X(addmatS));
      int matCount = matList->getLength();
      float rho = set_density(targetEl);
      float adbmal_sum=0;
      float wgt_sum=0;
      for (int m=0; m < matCount; m++)
      {
         DOMElement* matEl = (DOMElement*)matList->item(m);
         const XString materialS("material");
         const XString refS = matEl->getAttribute(X(materialS));
         DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refS));
         float lambda_m = set_radlen(refEl);
         float rho_m = set_density(refEl);
         float wgt_m=0;
         const XString natomsS("natoms");
         DOMNodeList* wgtList = matEl->getElementsByTagName(X(natomsS));
         if (wgtList->getLength())
         {
            DOMElement* natomsEl = (DOMElement*)wgtList->item(0);
            const XString nAttS("n");
            const XString nS = natomsEl->getAttribute(X(nAttS));
            int n = atoi(S(nS));
            wgt_m = n*set_a(refEl);
         }
         const XString fractionS("fractionmass");
         wgtList = matEl->getElementsByTagName(X(fractionS));
         if (wgtList->getLength())
         {
            DOMElement* fractionEl = (DOMElement*)wgtList->item(0);
            const XString fractionAttS("fraction");
            const XString fracS = fractionEl->getAttribute(X(fractionAttS));
            wgt_m = atof(S(fracS));
         }
         adbmal_sum += (wgt_m/lambda_m)*(rho/rho_m);
         wgt_sum += wgt_m;
      }
      char lambdaStr[20];
      sprintf(lambdaStr,"%g",wgt_sum/adbmal_sum);
      const XString lambdaS(lambdaStr);
      radlenEl = targetEl->getOwnerDocument()->createElement(X(realS));
      radlenEl->setAttribute(X(nameAttS),X(radlenS));
      radlenEl->setAttribute(X(valueAttS),X(lambdaS));
      targetEl->appendChild(radlenEl);
   }
   const XString resultS = radlenEl->getAttribute(X(valueAttS));
   return atof(S(resultS));
}

float makeTargetTable::set_collen(DOMElement* const targetEl)
{
   const XString nameAttS("name");
   const XString valueAttS("value");
   const XString realS("real");
   const XString collenS("collen");
   DOMElement* collenEl = lookup(targetEl,realS,collenS);
   if (collenEl == 0)
   {
      const XString addmatS("addmaterial");
      DOMNodeList* matList = targetEl->getElementsByTagName(X(addmatS));
      int matCount = matList->getLength();
      float rho = set_density(targetEl);
      float adbmal_sum=0;
      float wgt_sum=0;
      for (int m=0; m < matCount; m++)
      {
         DOMElement* matEl = (DOMElement*)matList->item(m);
         const XString materialS("material");
         const XString refS = matEl->getAttribute(X(materialS));
         DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refS));
         float lambda_m = set_collen(refEl);
         float rho_m = set_density(refEl);
         float wgt_m=0;
         const XString natomsS("natoms");
         DOMNodeList* wgtList = matEl->getElementsByTagName(X(natomsS));
         if (wgtList->getLength())
         {
            DOMElement* natomsEl = (DOMElement*)wgtList->item(0);
            const XString nAttS("n");
            const XString nS = natomsEl->getAttribute(X(nAttS));
            int n = atoi(S(nS));
            wgt_m = n*set_a(refEl);
         }
         const XString fractionS("fractionmass");
         wgtList = matEl->getElementsByTagName(X(fractionS));
         if (wgtList->getLength())
         {
            DOMElement* fractionEl = (DOMElement*)wgtList->item(0);
            const XString fractionAttS("fraction");
            const XString fracS = fractionEl->getAttribute(X(fractionAttS));
            wgt_m = atof(S(fracS));
         }
         adbmal_sum += (wgt_m/lambda_m)*(rho/rho_m);
         wgt_sum += wgt_m;
      }
      char lambdaStr[20];
      sprintf(lambdaStr,"%g",wgt_sum/adbmal_sum);
      const XString lambdaS(lambdaStr);
      collenEl = targetEl->getOwnerDocument()->createElement(X(realS));
      collenEl->setAttribute(X(nameAttS),X(collenS));
      collenEl->setAttribute(X(valueAttS),X(lambdaS));
      targetEl->appendChild(collenEl);
   }
   const XString resultS = collenEl->getAttribute(X(valueAttS));
   return atof(S(resultS));
}

float makeTargetTable::set_abslen(DOMElement* const targetEl)
{
   const XString nameAttS("name");
   const XString valueAttS("value");
   const XString realS("real");
   const XString abslenS("abslen");
   DOMElement* abslenEl = lookup(targetEl,realS,abslenS);
   if (abslenEl == 0)
   {
      const XString addmatS("addmaterial");
      DOMNodeList* matList = targetEl->getElementsByTagName(X(addmatS));
      int matCount = matList->getLength();
      float rho = set_density(targetEl);
      float adbmal_sum=0;
      float wgt_sum=0;
      for (int m=0; m < matCount; m++)
      {
         DOMElement* matEl = (DOMElement*)matList->item(m);
         const XString materialS("material");
         const XString refS = matEl->getAttribute(X(materialS));
         DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refS));
         float lambda_m = set_abslen(refEl);
         float rho_m = set_density(refEl);
         float wgt_m=0;
         const XString natomsS("natoms");
         DOMNodeList* wgtList = matEl->getElementsByTagName(X(natomsS));
         if (wgtList->getLength())
         {
            DOMElement* natomsEl = (DOMElement*)wgtList->item(0);
            const XString nAttS("n");
            const XString nS = natomsEl->getAttribute(X(nAttS));
            int n = atoi(S(nS));
            wgt_m = n*set_a(refEl);
         }
         const XString fractionS("fractionmass");
         wgtList = matEl->getElementsByTagName(X(fractionS));
         if (wgtList->getLength())
         {
            DOMElement* fractionEl = (DOMElement*)wgtList->item(0);
            const XString fractionAttS("fraction");
            const XString fracS = fractionEl->getAttribute(X(fractionAttS));
            wgt_m = atof(S(fracS));
         }
         adbmal_sum += (wgt_m/lambda_m)*(rho/rho_m);
         wgt_sum += wgt_m;
      }
      char lambdaStr[20];
      sprintf(lambdaStr,"%g",wgt_sum/adbmal_sum);
      const XString lambdaS(lambdaStr);
      abslenEl = targetEl->getOwnerDocument()->createElement(X(realS));
      abslenEl->setAttribute(X(nameAttS),X(abslenS));
      abslenEl->setAttribute(X(valueAttS),X(lambdaS));
      targetEl->appendChild(abslenEl);
   }
   const XString resultS = abslenEl->getAttribute(X(valueAttS));
   return atof(S(resultS));
}

float makeTargetTable::set_dedx(DOMElement* const targetEl)
{
   const XString nameAttS("name");
   const XString valueAttS("value");
   const XString realS("real");
   const XString dedxS("dedx");
   DOMElement* dedxEl = lookup(targetEl,realS,dedxS);
   if (dedxEl == 0)
   {
      const XString addmatS("addmaterial");
      DOMNodeList* matList = targetEl->getElementsByTagName(X(addmatS));
      int matCount = matList->getLength();
      float rho = set_density(targetEl);
      float dedx_sum=0;
      float wgt_sum=0;
      for (int m=0; m < matCount; m++)
      {
         DOMElement* matEl = (DOMElement*)matList->item(m);
         const XString materialS("material");
         const XString refS = matEl->getAttribute(X(materialS));
         DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refS));
         float dedx_m = set_dedx(refEl);
         float rho_m = set_density(refEl);
         float wgt_m=0;
         const XString natomsS("natoms");
         DOMNodeList* wgtList = matEl->getElementsByTagName(X(natomsS));
         if (wgtList->getLength())
         {
            DOMElement* natomsEl = (DOMElement*)wgtList->item(0);
            const XString nAttS("n");
            const XString nS = natomsEl->getAttribute(X(nAttS));
            int n = atoi(S(nS));
            wgt_m = n*set_a(refEl);
         }
         const XString fractionS("fractionmass");
         wgtList = matEl->getElementsByTagName(X(fractionS));
         if (wgtList->getLength())
         {
            DOMElement* fractionEl = (DOMElement*)wgtList->item(0);
            const XString fractionAttS("fraction");
            const XString fracS = fractionEl->getAttribute(X(fractionAttS));
            wgt_m = atof(S(fracS));
         }
         dedx_sum += (wgt_m*dedx_m)*(rho/rho_m);
         wgt_sum += wgt_m;
      }
      char dedxStr[20];
      sprintf(dedxStr,"%g",dedx_sum/wgt_sum);
      const XString dedxStrS(dedxStr);
      dedxEl = targetEl->getOwnerDocument()->createElement(X(realS));
      dedxEl->setAttribute(X(nameAttS),X(dedxS));
      dedxEl->setAttribute(X(valueAttS),X(dedxStrS));
      targetEl->appendChild(dedxEl);
   }
   const XString resultS = dedxEl->getAttribute(X(valueAttS));
   return atof(S(resultS));
}

DOMElement* makeTargetTable::lookup(const DOMElement* const targetEl,
                                    const XString& type, const XString& name)
{
   if (targetEl == 0)
      return 0;

   const XString nameAttS("name");
   DOMNodeList* varList = targetEl->getElementsByTagName(X(type));
   int varCount = varList->getLength();
   for (int v=0; v < varCount; v++)
   {
      DOMElement* varEl = (DOMElement*)varList->item(v);
      const XString nameS = varEl->getAttribute(X(nameAttS));
      if (nameS.equals(name))
         return varEl;
   }

   const XString paramAttS("parameters");
   const XString paramS = targetEl->getAttribute(X(paramAttS));
   if (paramS != 0)
   {
      DOMElement* paramEl = targetEl->getOwnerDocument()->getElementById(X(paramS));
      varList = paramEl->getElementsByTagName(X(type));
      varCount = varList->getLength();
      for (int v=0; v < varCount; v++)
      {
         DOMElement* varEl = (DOMElement*)varList->item(v);
         const XString nameS = varEl->getAttribute(X(nameAttS));
         if (nameS.equals(name))
            return varEl;
      }
   }

   for (int i=0; i < fTableLen; i++)
   {
      if (fTargetEl[i] == targetEl)
      {
         return lookup(fParentEl[i],type,name);
      }
   }
   return 0;
}

void makeTargetTable::dump(const DOMElement* const targetEl)
{
   static int indent = 0;
   for (int n = 0; n < indent; n++)
   {
      cout << "  ";
   }
   
   XString tagS(targetEl->getTagName());
   cout << "<" << S(tagS);
   DOMNamedNodeMap* attrList = targetEl->getAttributes();
   int attrListLength = attrList->getLength();
   for (int a = 0; a < attrListLength; a++)
   {
      DOMNode* node = attrList->item(a);
      XString nameS(node->getNodeName());
      XString valueS(node->getNodeValue());
      cout << " " << S(nameS) << "=\"" << S(valueS) << "\"";
   }

   DOMNodeList* contList = targetEl->getChildNodes();
   int contListLength = contList->getLength();
   if (contListLength > 0)
   {
      cout << ">" << endl;
      indent++;
      for (int c = 0; c < contListLength; c++)
      {
         DOMNode* node = contList->item(c);
         if (node->getNodeType() == DOMNode::ELEMENT_NODE)
         {
            DOMElement* contEl = (DOMElement*) node;
            dump(contEl);
         }
      }
      indent--;
      for (int n = 0; n < indent; n++)
      {
         cout << "  ";
      }
      cout << "</" << S(tagS) << ">" << endl;
   }
   else
   {
      cout << " />" << endl;
   }
}

void makeTargetTable::dump(DOMElement* const targetEl)
{
   return dump((const DOMElement* const) targetEl);
}
