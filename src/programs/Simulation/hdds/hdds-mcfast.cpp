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

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax/SAXException.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/sax/ErrorHandler.hpp>

#include "XString.hpp"
#include "XParsers.hpp"

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>

using namespace std;

XERCES_CPP_NAMESPACE_USE

#define X(XString) XString.unicode_str()
#define S(XString) XString.c_str()


class makeTargetTable
{
public :
   int add(DOMElement* const targetEl,
           DOMElement* const parentEl = 0);
   DOMElement* lookup(const DOMElement* const targetEl,
                      XString& type, const XString& name);

private :
   std::map<const DOMElement*,DOMElement*> fParentElList;

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

#ifdef BASENAME_IN_LIBGEN
#include <libgen.h>
#elif defined BASENAME_USE_BUILTIN
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

class modelTableEntry
{
 public:
   int refcount;
   int maxcount;
   ostringstream db;

   modelTableEntry() : refcount(0), maxcount(0) {};
   modelTableEntry(const modelTableEntry& src)
   {
      refcount = src.refcount;
      maxcount = src.maxcount;
      db.str(src.db.str());
   }
   modelTableEntry& operator=(const modelTableEntry& src)
   {
      refcount = src.refcount;
      maxcount = src.maxcount;
      db.str(src.db.str());
   }
};

std::map<const std::string,modelTableEntry> modelTable;

makeTargetTable targetTable;


void usage()
{
    cerr << "\nUsage:\n"
            "    hdds-mcfast [-v] {HDDS file}\n\n"
            "Options:\n"
            "    -v   validate only\n"
         << endl;
}

class DbMaker
{
/* This class contains the utility functions that read the
 * MCFast materials and geometry db files and generate the
 * output structure based on that information.
 */ 
 public:
   static void processTemplateFile(const DOMElement* const targetEl,
                                            const char* const fname);
   static void makedb(DOMElement* el);
   static void printdb();
};

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

   modelTableEntry entry;

   entry.refcount = 1;
   entry.maxcount = 1;
   entry.db.str("");
   entry.db << "database mcfast 0000" << endl;
   modelTable["database"] = entry;         // header line must come first

   entry.refcount = 0;
   entry.maxcount = 0;
   entry.db.str("");
   modelTable["detector"] = entry;         // detector is first declaration

   entry.refcount = 0;
   entry.maxcount = 0;
   entry.db.str("");
   modelTable["Material"] = entry;         // force materials next

   entry.refcount = 0;
   entry.maxcount = 0;
   entry.db.str("");
   modelTable["Mixture"] = entry;          // then mixtures

   DOMElement* rootEl = doc->getDocumentElement();
   DbMaker::makedb(rootEl);

   entry.refcount = 1;
   entry.maxcount = 1;
   entry.db.str("");
   entry.db << "include db/hitsontrack.db" << endl
            << "make HitsOnTrack 4 0 0" << endl;
   modelTable["hitsontrack"] = entry;      // last comes histontrack

   DbMaker::printdb();

   XMLPlatformUtils::Terminate();
   return 0;
}


int makeTargetTable::add(DOMElement* const targetEl,
                         DOMElement* const parentEl)
{
   if (targetEl == 0)
   {
      cerr << "makeTargetTable::add - called with null pointer" << endl;
      return 0;
   }

   std::map<const DOMElement*,DOMElement*>::iterator iter;
   iter = fParentElList.find(targetEl);
   if (iter == fParentElList.end())
   {
      fParentElList[targetEl] = parentEl;
      XString tagS = targetEl->getTagName();
      if (tagS == "mcfast")
      {
         addModel(targetEl);
      }
      else if (tagS == "element")
      {
         addElement(targetEl);
      }
      else if (tagS == "composite")
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
   XString modelAttS("model");
   XString materialS("Material");
   targetEl->setAttribute(X(modelAttS),X(materialS));
   XString templAttS("template");
   XString dbS("db/materials.db");
   targetEl->setAttribute(X(templAttS),X(dbS));
   DbMaker::processTemplateFile(targetEl,S(dbS));
   return 0;
}

int makeTargetTable::addComposite(DOMElement* const targetEl)
{
   XString densityS("density");
   XString radlenS("radlen");
   XString collenS("collen");
   XString abslenS("abslen");
   XString dedxS("dedx");
   XString realS("real");
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

   XString addmatS("addmaterial");
   DOMNodeList* matList = targetEl->getElementsByTagName(X(addmatS));
   const int matCount = matList->getLength();
   char nmatStr[10];
   sprintf(nmatStr,"%d",matCount);
   XString nmatS(nmatStr);
   XString nmatAttS("nmat");
   XString intS("int");
   DOMElement* nmatEl = targetEl->getOwnerDocument()->createElement(X(intS));
   XString nameAttS("name");
   XString valueAttS("value");
   nmatEl->setAttribute(X(nameAttS),X(nmatAttS));
   nmatEl->setAttribute(X(valueAttS),X(nmatS));
   targetEl->appendChild(nmatEl);

   std::vector<float> fVol;
   std::vector<DOMElement*> dataEl;
   XString refvecS("reference_vector");
   DOMElement* matVecEl = targetEl->getOwnerDocument()->createElement(X(refvecS));
   XString matnameS("matnames");
   matVecEl->setAttribute(X(nameAttS),X(matnameS));
   for (int m=0; m < matCount; m++)
   {
      DOMElement* matEl = (DOMElement*)matList->item(m);
      XString materialS("material");
      XString refIdS(matEl->getAttribute(X(materialS)));
      XString refdataS("reference_data");
      dataEl[m] = targetEl->getOwnerDocument()->createElement(X(refdataS));
      XString valueAttS("value");
      dataEl[m]->setAttribute(X(valueAttS),X(refIdS));
      matVecEl->appendChild(dataEl[m]);
      DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refIdS));
      add(refEl);
      XString fractionS("fractionmass");
      DOMNodeList* specList = matEl->getElementsByTagName(X(fractionS));
      if (specList->getLength() == 0)
      {
         XString nameAttS("name");
         XString nameS = matEl->getAttribute(X(nameAttS));
         cerr << "makeTargetTable::addComposite - "
              << "composite " << S(nameS) << " is missing information."
              << endl;
         exit(2);
      }
      DOMElement* specEl = (DOMElement*)specList->item(0);
      XString fracAttS("fraction");
      XString fS(specEl->getAttribute(X(fracAttS)));
      float fMass_m = atof(S(fS));
      float density_m = set_density(refEl);
      float density = set_density(targetEl);
      fVol[m] = fMass_m*(density/density_m);
   }
   targetEl->appendChild(matVecEl);

   XString rarrayS("real_array");
   DOMElement* fVolEl = targetEl->getOwnerDocument()->createElement(X(rarrayS));
   XString propS("prop");
   fVolEl->setAttribute(X(nameAttS),X(propS));
   char str[20];
   sprintf(str,"%f",fVol[0]);
   XString fVolS(str);
   for (int m=1; m < matCount; m++)
   {
      sprintf(str," %f",fVol[m]);
      XString xS(str);
      fVolS += xS;
   }
   XString valuesAttS("values");
   fVolEl->setAttribute(X(valuesAttS),X(fVolS));
   targetEl->appendChild(fVolEl);

   XString modelAttS("model");
   XString templAttS("template");
   XString mixtureS("Mixture");
   XString dbS("db/mixtures.db");
   targetEl->setAttribute(X(modelAttS),X(mixtureS));
   targetEl->setAttribute(X(templAttS),X(dbS));
   DbMaker::processTemplateFile(targetEl,S(dbS));

   return 0;
}

int makeTargetTable::addModel(DOMElement* const targetEl)
{
   XString paramAttS("parameters");
   XString parS = targetEl->getAttribute(X(paramAttS));
   if (parS.size() != 0)
   {
      DOMElement* parEl = targetEl->getOwnerDocument()->getElementById(X(parS));
      XString refS("reference");
      DOMNodeList* matList = parEl->getElementsByTagName(X(refS));
      int matCount = matList->getLength();
      for (int m=0; m < matCount; m++)
      {
         DOMElement* ptrEl = (DOMElement*)matList->item(m);
         XString valueAttrS("value");
         XString valueS = ptrEl->getAttribute(X(valueAttrS));
         DOMElement* matEl = ptrEl->getOwnerDocument()->getElementById(X(valueS));
         add(matEl);
      }
      XString refdatS("reference_data");
      matList = parEl->getElementsByTagName(X(refdatS));
      matCount = matList->getLength();
      for (int m=0; m < matCount; m++)
      {
         DOMElement* ptrEl = (DOMElement*)matList->item(m);
         XString valueAttrS("value");
         XString valueS = ptrEl->getAttribute(X(valueAttrS));
         DOMElement* matEl = ptrEl->getOwnerDocument()->getElementById(X(valueS));
         add(matEl);
      }
   }
   XString refS("reference");
   DOMNodeList* matList = targetEl->getElementsByTagName(X(refS));
   int matCount = matList->getLength();
   for (int m=0; m < matCount; m++)
   {
      DOMElement* ptrEl = (DOMElement*)matList->item(m);
      XString valueAttrS("value");
      XString valueS = ptrEl->getAttribute(X(valueAttrS));
      DOMElement* matEl = ptrEl->getOwnerDocument()->getElementById(X(valueS));
      add(matEl);
   }
   XString refdatS("reference_data");
   matList = targetEl->getElementsByTagName(X(refdatS));
   matCount = matList->getLength();
   for (int m=0; m < matCount; m++)
   {
      DOMElement* ptrEl = (DOMElement*)matList->item(m);
      XString valueAttrS("value");
      XString valueS = ptrEl->getAttribute(X(valueAttrS));
      DOMElement* matEl = ptrEl->getOwnerDocument()->getElementById(X(valueS));
      add(matEl);
   }
   for (DOMNode* var = targetEl->getFirstChild(); 
        var != 0;
        var = var->getNextSibling() )
   {
      if (var->getNodeType() != DOMNode::ELEMENT_NODE)
         continue;
      DOMElement* varEl = (DOMElement*)var;
      XString varS = varEl->getTagName();
      if (varS == "mcfast")
      {
         add(varEl,targetEl);
      }
   }
   XString templAttS("template");
   XString templS = targetEl->getAttribute(X(templAttS));
   DbMaker::processTemplateFile(targetEl,S(templS));
   return 0;
}

void makeTargetTable::set_name(DOMElement* const targetEl)
{
   XString nameAttS("name");
   XString valueAttS("value");
   XString stringS("string");
   DOMElement* nameEl = lookup(targetEl,stringS,nameAttS);
   if (nameEl == 0)
   {
      XString nameS = targetEl->getAttribute(X(nameAttS));
      nameEl = targetEl->getOwnerDocument()->createElement(X(stringS));
      nameEl->setAttribute(X(nameAttS),X(nameAttS));
      nameEl->setAttribute(X(valueAttS),X(nameS));
      targetEl->appendChild(nameEl);
   }
}

float makeTargetTable::set_a(DOMElement* const targetEl)
{
   XString nameAttS("name");
   XString valueAttS("value");
   XString realS("real");
   XString aAttS("a");
   XString aS = targetEl->getAttribute(X(aAttS));
   DOMElement* aEl = lookup(targetEl,realS,aAttS);
   if (aEl == 0)
   {
      if (aS.size() == 0)
      {
         XString addmatS("addmaterial");
         DOMNodeList* matList = targetEl->getElementsByTagName(X(addmatS));
         int matCount = matList->getLength();
         float a_sum=0;
         float wgt_sum=0;
         for (int m=0; m < matCount; m++)
         {
            DOMElement* matEl = (DOMElement*)matList->item(m);
            XString materialS("material");
            XString refS = matEl->getAttribute(X(materialS));
            DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refS));
            float a_m = set_a(refEl);
            float wgt_m=0;
            XString natomsS("natoms");
            DOMNodeList* wgtList = matEl->getElementsByTagName(X(natomsS));
            if (wgtList->getLength())
            {
               DOMElement* natomsEl = (DOMElement*)wgtList->item(0);
               XString nAttS("n");
               XString nS = natomsEl->getAttribute(X(nAttS));
               int n = atoi(S(nS));
               wgt_m = n*a_m;
            }
            XString fractionS("fractionmass");
            wgtList = matEl->getElementsByTagName(X(fractionS));
            if (wgtList->getLength())
            {
               DOMElement* fractionEl = (DOMElement*)wgtList->item(0);
               XString fractionAttS("fraction");
               XString fracS = fractionEl->getAttribute(X(fractionAttS));
               wgt_m = atof(S(fracS));
            }
            a_sum += a_m*wgt_m;
            wgt_sum += wgt_m;
         }
         char aStr[20];
         sprintf(aStr,"%g",a_sum/wgt_sum);
         XString aStrS(aStr);
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
   XString nameAttS("name");
   XString valueAttS("value");
   XString realS("real");
   XString zAttS("z");
   XString zS = targetEl->getAttribute(X(zAttS));
   DOMElement* zEl = lookup(targetEl,realS,zAttS);
   if (zEl == 0)
   {
      if (zS.size() == 0)
      {
         XString addmatS("addmaterial");
         DOMNodeList* matList = targetEl->getElementsByTagName(X(addmatS));
         int matCount = matList->getLength();
         float z_sum=0;
         float wgt_sum=0;
         for (int m=0; m < matCount; m++)
         {
            DOMElement* matEl = (DOMElement*)matList->item(m);
            XString materialS("material");
            XString refS = matEl->getAttribute(X(materialS));
            DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refS));
            float z_m = set_z(refEl);
            float wgt_m=0;
            XString natomsS("natoms");
            DOMNodeList* wgtList = matEl->getElementsByTagName(X(natomsS));
            if (wgtList->getLength())
            {
               DOMElement* natomsEl = (DOMElement*)wgtList->item(0);
               XString nAttS("n");
               XString nS = natomsEl->getAttribute(X(nAttS));
               int n = atoi(S(nS));
               XString materialS("material");
               XString refIdS = matEl->getAttribute(X(materialS));
               DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refIdS));
               wgt_m = n*set_a(refEl);
            }
            XString fractionS("fractionmass");
            wgtList = matEl->getElementsByTagName(X(fractionS));
            if (wgtList->getLength())
            {
               DOMElement* fractionEl = (DOMElement*)wgtList->item(0);
               XString fractionAttS("fraction");
               XString fracS = fractionEl->getAttribute(X(fractionAttS));
               wgt_m = atof(S(fracS));
            }
            z_sum += z_m*wgt_m;
            wgt_sum += wgt_m;
         }
         char zStr[20];
         sprintf(zStr,"%g",z_sum/wgt_sum);
         XString zStrS(zStr);
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
   XString nameAttS("name");
   XString valueAttS("value");
   XString realS("real");
   XString densityS("density");
   XString unitAttS("unit");
   DOMElement* densityEl = lookup(targetEl,realS,densityS);
   if (densityEl == 0)
   {
      XString addmatS("addmaterial");
      DOMNodeList* matList = targetEl->getElementsByTagName(X(addmatS));
      int matCount = matList->getLength();
      float ohr_sum=0;
      float wgt_sum=0;
      for (int m=0; m < matCount; m++)
      {
         DOMElement* matEl = (DOMElement*)matList->item(m);
         XString materialS("material");
         XString refS = matEl->getAttribute(X(materialS));
         DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refS));
         float rho_m = set_density(refEl);
         float wgt_m=0;
         XString natomsS("natoms");
         DOMNodeList* wgtList = matEl->getElementsByTagName(X(natomsS));
         if (wgtList->getLength())
         {
            XString nameS = matEl->getAttribute(X(nameAttS));
            cerr << "makeTargetTable::set_density - "
                 << "no automatic density calculation for atomic mixture "
                 << S(nameS) << endl;
            exit(2);
         }
         XString fractionS("fractionmass");
         wgtList = matEl->getElementsByTagName(X(fractionS));
         if (wgtList->getLength())
         {
            DOMElement* fractionEl = (DOMElement*)wgtList->item(0);
            XString fractionAttS("fraction");
            XString fracS = fractionEl->getAttribute(X(fractionAttS));
            wgt_m = atof(S(fracS));
         }
         ohr_sum += wgt_m/rho_m;
         wgt_sum += wgt_m;
      }
      char rhoStr[20];
      sprintf(rhoStr,"%g",wgt_sum/ohr_sum);
      XString rhoS(rhoStr);
      densityEl = targetEl->getOwnerDocument()->createElement(X(realS));
      densityEl->setAttribute(X(nameAttS),X(densityS));
      densityEl->setAttribute(X(valueAttS),X(rhoS));
      XString unitS("g/cm^3");
      densityEl->setAttribute(X(unitAttS),X(unitS));
      targetEl->appendChild(densityEl);
   }
   XString resultS = densityEl->getAttribute(X(valueAttS));
   XString unitS = densityEl->getAttribute(X(unitAttS));
   if (unitS != "g/cm^3")
   {
      cerr << "makeTargetTable::set_density - "
           << "unsupported density units " << S(unitS)
           << " found in materials specification" << endl;
      exit(2);
   }
   return atof(S(resultS));
}

float makeTargetTable::set_radlen(DOMElement* const targetEl)
{
   XString nameAttS("name");
   XString valueAttS("value");
   XString realS("real");
   XString unitAttS("unit");
   XString radlenS("radlen");
   DOMElement* radlenEl = lookup(targetEl,realS,radlenS);
   if (radlenEl == 0)
   {
      XString addmatS("addmaterial");
      DOMNodeList* matList = targetEl->getElementsByTagName(X(addmatS));
      int matCount = matList->getLength();
      float adbmal_sum=0;
      float wgt_sum=0;
      for (int m=0; m < matCount; m++)
      {
         DOMElement* matEl = (DOMElement*)matList->item(m);
         XString materialS("material");
         XString refS = matEl->getAttribute(X(materialS));
         DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refS));
         float lambda_m = set_radlen(refEl);
         float wgt_m=0;
         XString natomsS("natoms");
         DOMNodeList* wgtList = matEl->getElementsByTagName(X(natomsS));
         if (wgtList->getLength())
         {
            DOMElement* natomsEl = (DOMElement*)wgtList->item(0);
            XString nAttS("n");
            XString nS = natomsEl->getAttribute(X(nAttS));
            int n = atoi(S(nS));
            wgt_m = n*set_a(refEl);
         }
         XString fractionS("fractionmass");
         wgtList = matEl->getElementsByTagName(X(fractionS));
         if (wgtList->getLength())
         {
            DOMElement* fractionEl = (DOMElement*)wgtList->item(0);
            XString fractionAttS("fraction");
            XString fracS = fractionEl->getAttribute(X(fractionAttS));
            wgt_m = atof(S(fracS));
         }
         adbmal_sum += wgt_m/lambda_m;
         wgt_sum += wgt_m;
      }
      char lambdaStr[20];
      sprintf(lambdaStr,"%g",wgt_sum/adbmal_sum);
      XString lambdaS(lambdaStr);
      radlenEl = targetEl->getOwnerDocument()->createElement(X(realS));
      radlenEl->setAttribute(X(nameAttS),X(radlenS));
      radlenEl->setAttribute(X(valueAttS),X(lambdaS));
      XString unitS("g/cm^2");
      radlenEl->setAttribute(X(unitAttS),X(unitS));
      targetEl->appendChild(radlenEl);
   }
   XString resultS = radlenEl->getAttribute(X(valueAttS));
   XString unitS = radlenEl->getAttribute(X(unitAttS));
   if (unitS != "g/cm^2")
   {
      cerr << "makeTargetTable::set_radlen - "
           << "unsupported radlen units " << S(unitS)
           << " found in materials specification" << endl;
      exit(2);
   }
   return atof(S(resultS));
}

float makeTargetTable::set_collen(DOMElement* const targetEl)
{
   XString nameAttS("name");
   XString valueAttS("value");
   XString realS("real");
   XString collenS("collen");
   XString unitAttS("unit");
   DOMElement* collenEl = lookup(targetEl,realS,collenS);
   if (collenEl == 0)
   {
      XString addmatS("addmaterial");
      DOMNodeList* matList = targetEl->getElementsByTagName(X(addmatS));
      int matCount = matList->getLength();
      float adbmal_sum=0;
      float wgt_sum=0;
      for (int m=0; m < matCount; m++)
      {
         DOMElement* matEl = (DOMElement*)matList->item(m);
         XString materialS("material");
         XString refS = matEl->getAttribute(X(materialS));
         DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refS));
         float lambda_m = set_collen(refEl);
         float wgt_m=0;
         XString natomsS("natoms");
         DOMNodeList* wgtList = matEl->getElementsByTagName(X(natomsS));
         if (wgtList->getLength())
         {
            DOMElement* natomsEl = (DOMElement*)wgtList->item(0);
            XString nAttS("n");
            XString nS = natomsEl->getAttribute(X(nAttS));
            int n = atoi(S(nS));
            wgt_m = n*set_a(refEl);
         }
         XString fractionS("fractionmass");
         wgtList = matEl->getElementsByTagName(X(fractionS));
         if (wgtList->getLength())
         {
            DOMElement* fractionEl = (DOMElement*)wgtList->item(0);
            XString fractionAttS("fraction");
            XString fracS = fractionEl->getAttribute(X(fractionAttS));
            wgt_m = atof(S(fracS));
         }
         adbmal_sum += (wgt_m/lambda_m);
         wgt_sum += wgt_m;
      }
      char lambdaStr[20];
      sprintf(lambdaStr,"%g",wgt_sum/adbmal_sum);
      XString lambdaS(lambdaStr);
      collenEl = targetEl->getOwnerDocument()->createElement(X(realS));
      collenEl->setAttribute(X(nameAttS),X(collenS));
      collenEl->setAttribute(X(valueAttS),X(lambdaS));
      XString unitS("g/cm^2");
      collenEl->setAttribute(X(unitAttS),X(unitS));
      targetEl->appendChild(collenEl);
   }
   XString resultS = collenEl->getAttribute(X(valueAttS));
   XString unitS = collenEl->getAttribute(X(unitAttS));
   if (unitS != "g/cm^2")
   {
      cerr << "makeTargetTable::set_collen - "
           << "unsupported collen units " << S(unitS)
           << " found in materials specification" << endl;
      exit(2);
   }
   return atof(S(resultS));
}

float makeTargetTable::set_abslen(DOMElement* const targetEl)
{
   XString nameAttS("name");
   XString valueAttS("value");
   XString realS("real");
   XString abslenS("abslen");
   XString unitAttS("unit");
   DOMElement* abslenEl = lookup(targetEl,realS,abslenS);
   if (abslenEl == 0)
   {
      XString addmatS("addmaterial");
      DOMNodeList* matList = targetEl->getElementsByTagName(X(addmatS));
      int matCount = matList->getLength();
      float adbmal_sum=0;
      float wgt_sum=0;
      for (int m=0; m < matCount; m++)
      {
         DOMElement* matEl = (DOMElement*)matList->item(m);
         XString materialS("material");
         XString refS = matEl->getAttribute(X(materialS));
         DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refS));
         float lambda_m = set_abslen(refEl);
         float wgt_m=0;
         XString natomsS("natoms");
         DOMNodeList* wgtList = matEl->getElementsByTagName(X(natomsS));
         if (wgtList->getLength())
         {
            DOMElement* natomsEl = (DOMElement*)wgtList->item(0);
            XString nAttS("n");
            XString nS = natomsEl->getAttribute(X(nAttS));
            int n = atoi(S(nS));
            wgt_m = n*set_a(refEl);
         }
         XString fractionS("fractionmass");
         wgtList = matEl->getElementsByTagName(X(fractionS));
         if (wgtList->getLength())
         {
            DOMElement* fractionEl = (DOMElement*)wgtList->item(0);
            XString fractionAttS("fraction");
            XString fracS = fractionEl->getAttribute(X(fractionAttS));
            wgt_m = atof(S(fracS));
         }
         adbmal_sum += wgt_m/lambda_m;
         wgt_sum += wgt_m;
      }
      char lambdaStr[20];
      sprintf(lambdaStr,"%g",wgt_sum/adbmal_sum);
      XString lambdaS(lambdaStr);
      abslenEl = targetEl->getOwnerDocument()->createElement(X(realS));
      abslenEl->setAttribute(X(nameAttS),X(abslenS));
      abslenEl->setAttribute(X(valueAttS),X(lambdaS));
      XString unitS("g/cm^2");
      abslenEl->setAttribute(X(unitAttS),X(unitS));
      targetEl->appendChild(abslenEl);
   }
   XString resultS = abslenEl->getAttribute(X(valueAttS));
   XString unitS = abslenEl->getAttribute(X(unitAttS));
   if (unitS != "g/cm^2")
   {
      cerr << "makeTargetTable::set_abslen - "
           << "unsupported abslen units " << S(unitS)
           << " found in materials specification" << endl;
      exit(2);
   }
   return atof(S(resultS));
}

float makeTargetTable::set_dedx(DOMElement* const targetEl)
{
   XString nameAttS("name");
   XString valueAttS("value");
   XString realS("real");
   XString dedxS("dedx");
   XString unitAttS("unit");
   DOMElement* dedxEl = lookup(targetEl,realS,dedxS);
   if (dedxEl == 0)
   {
      XString addmatS("addmaterial");
      DOMNodeList* matList = targetEl->getElementsByTagName(X(addmatS));
      int matCount = matList->getLength();
      float dedx_sum=0;
      float wgt_sum=0;
      for (int m=0; m < matCount; m++)
      {
         DOMElement* matEl = (DOMElement*)matList->item(m);
         XString materialS("material");
         XString refS = matEl->getAttribute(X(materialS));
         DOMElement* refEl = matEl->getOwnerDocument()->getElementById(X(refS));
         float dedx_m = set_dedx(refEl);
         float wgt_m=0;
         XString natomsS("natoms");
         DOMNodeList* wgtList = matEl->getElementsByTagName(X(natomsS));
         if (wgtList->getLength())
         {
            DOMElement* natomsEl = (DOMElement*)wgtList->item(0);
            XString nAttS("n");
            XString nS = natomsEl->getAttribute(X(nAttS));
            int n = atoi(S(nS));
            wgt_m = n*set_a(refEl);
         }
         XString fractionS("fractionmass");
         wgtList = matEl->getElementsByTagName(X(fractionS));
         if (wgtList->getLength())
         {
            DOMElement* fractionEl = (DOMElement*)wgtList->item(0);
            XString fractionAttS("fraction");
            XString fracS = fractionEl->getAttribute(X(fractionAttS));
            wgt_m = atof(S(fracS));
         }
         dedx_sum += wgt_m*dedx_m;
         wgt_sum += wgt_m;
      }
      char dedxStr[20];
      sprintf(dedxStr,"%g",dedx_sum/wgt_sum);
      XString dedxStrS(dedxStr);
      dedxEl = targetEl->getOwnerDocument()->createElement(X(realS));
      dedxEl->setAttribute(X(nameAttS),X(dedxS));
      dedxEl->setAttribute(X(valueAttS),X(dedxStrS));
      XString unitS("MeV/g/cm^2");
      dedxEl->setAttribute(X(unitAttS),X(unitS));
      targetEl->appendChild(dedxEl);
   }
   XString resultS = dedxEl->getAttribute(X(valueAttS));
   XString unitS = dedxEl->getAttribute(X(unitAttS));
   if (unitS != "MeV/g/cm^2")
   {
      cerr << "makeTargetTable::set_dedx - "
           << "unsupported dedx units " << S(unitS)
           << " found in materials specification" << endl;
      exit(2);
   }
   return atof(S(resultS));
}

DOMElement* makeTargetTable::lookup(const DOMElement* const targetEl,
                                    XString& type, const XString& name)
{
   if (targetEl == 0)
      return 0;

   XString nameAttS("name");
   DOMNodeList* varList = targetEl->getElementsByTagName(X(type));
   int varCount = varList->getLength();
   for (int v=0; v < varCount; v++)
   {
      DOMElement* varEl = (DOMElement*)varList->item(v);
      XString nameS = varEl->getAttribute(X(nameAttS));
      if (nameS == name)
         return varEl;
   }

   XString paramAttS("parameters");
   XString paramS = targetEl->getAttribute(X(paramAttS));
   if (paramS.size() != 0)
   {
      DOMElement* paramEl = targetEl->getOwnerDocument()->getElementById(X(paramS));
      varList = paramEl->getElementsByTagName(X(type));
      varCount = varList->getLength();
      for (int v=0; v < varCount; v++)
      {
         DOMElement* varEl = (DOMElement*)varList->item(v);
         XString nameS = varEl->getAttribute(X(nameAttS));
         if (nameS == name)
            return varEl;
      }
   }

   std::map<const DOMElement*,DOMElement*>::iterator iter;
   iter = fParentElList.find(targetEl);
   if (iter != fParentElList.end())
   {
      return lookup(iter->second,type,name);
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

void DbMaker::processTemplateFile(const DOMElement* const targetEl,
                                  const char* const fname)
{
   XString* processingTarget=0;
   int processedTemplates=0;

   XString modelAttS("model");
   XString modelS = targetEl->getAttribute(X(modelAttS));

   std::map<const std::string,modelTableEntry>::iterator iter;
   iter = modelTable.find(modelS);
   if (iter == modelTable.end())
   {
      modelTableEntry entry;
      entry.refcount = 0;
      entry.maxcount = 0;
      modelTable[modelS] = entry;
   }
   if (modelTable[modelS].refcount == 0)
   {
      modelTable[modelS].db << "include " << fname << endl;
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
         modelTable[modelS].db << endl;
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
         else if (processingTarget && *processingTarget == modelS)
         {
            if (++modelTable[modelS].refcount > dim)
            { 
               cerr << "hdds-mcfast: number of objects of type " << S(modelS)
                    << " exceeds maximum of " << dim << endl;
               cerr << "defined in input file " << fname << endl;
               cerr << "Increase the array size in the template statement"
                    << " found in the above file" << endl;
               cerr << "and try again." << endl;
               exit(2);
            }
            modelTable[modelS].maxcount = dim;
            modelTable[modelS].db << "make " << S(modelS);
            ++processedTemplates;
         }
      }
      else if (strcasecmp(token,"make") == 0)
      {
         std::map<const std::string,modelTableEntry>::iterator iter;
         iter = modelTable.find(sline.str());
         if (iter == modelTable.end())
         { 
            cerr << "hdds-mcfast: error in template file " << fname << endl;
            cerr << "Statement template " << S(modelS)
                 << " must appear before first make instance." << endl;
            exit(2);
         }
         else if (iter->second.refcount > iter->second.maxcount)
         {
            cerr << "hdds-mcfast: number of objects of type " << sline.str()
                 << " overflows table." << endl;
            cerr << "Increase the array size in the template file"
                 << " and try again." << endl;
            exit(2);
         }
      }
      else if (!(processingTarget && *processingTarget == modelS))
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
            XString typeS("int");
            XString varS(var);
            if (el = targetTable.lookup(targetEl,typeS,varS))
            {
	       XString valAttS("value");
               valueS = el->getAttribute(X(valAttS));
            }
         }
         else
         {
            XString typeS("int_array");
            XString varS(var);
            if (el = targetTable.lookup(targetEl,typeS,varS))
            {
	       XString valAttS("values");
               valueS = el->getAttribute(X(valAttS));
            }
         }
         if (valueS.size() == 0)
         {
            cerr << "hdds-mcfast: Parameter \"" << var << "\""
                 << " required in template file " << fname << endl
                 << "is missing from HDDS" << endl;
            exit(3);
         }
         char* intstr = strtok((char*)S(valueS)," ");
         modelTable[modelS].db << " " << intstr;
         for (int i=1; i < dim; i++)
         {
            if (intstr = strtok(0," "))
            {
               modelTable[modelS].db << " " << intstr;
            }
            else
            {
               modelTable[modelS].db << " 0";
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
            XString typeS("real");
            XString varS(var);
            if (el = targetTable.lookup(targetEl,typeS,varS))
            {
	       XString valAttS("value");
               valueS = el->getAttribute(X(valAttS));
            }
         }
         else
         {
            XString typeS("real_array");
            XString varS(var);
            if (el = targetTable.lookup(targetEl,typeS,varS))
            {
	       XString valAttS("values");
               valueS = el->getAttribute(X(valAttS));
            }
         }
         if (valueS.size() == 0)
         {
            cerr << "hdds-mcfast: Parameter \"" << var << "\""
                 << " required in template file " << fname << endl
                 << "is missing from HDDS" << endl;
            exit(3);
         }
	 XString unitAttS("unit");
         XString unitS = el->getAttribute(X(unitAttS));
         float rho=1e-30;
         DOMElement* densEl;
         XString realS("real");
         XString densAttS("density");
         if (densEl = targetTable.lookup(targetEl,realS,densAttS))
         {
	    XString uS = densEl->getAttribute(X(unitAttS));
            if (uS != "g/cm^3")
            {
               cerr << "hdds-mcfast: density parameter found with"
                    << " incorrect units " << S(uS) << endl;
               exit(3);
            }
	    XString valAttS("value");
            XString rhoS = densEl->getAttribute(X(valAttS));
            if (rhoS.size() != 0)
            {
               rho=atof(S(rhoS))+1e-30;
            }
         }
         float fconvert=
         /* standard unit for length in MCfast is cm */
            (unitS == "m")? 100
          : (unitS == "cm")? 1
          : (unitS == "mm")? 0.1
         /* standard unit for angles in MCfast is radians */
          : (unitS == "rad")? 1
          : (unitS == "mrad")? 0.001
          : (unitS == "deg")? M_PI/180
         /* standard unit for energy in MCfast is GeV */
          : (unitS == "MeV")? 0.001
          : (unitS == "GeV")? 1
         /* standard unit for density in MCfast is g/cm^3 */
          : (unitS == "g/cm^3")? 1
         /* standard unit for interaction lengths in MCfast is cm */
          : (unitS == "g/cm^2")? 1/rho
         /* standard unit for dE/dx in MCfast is GeV/cm */
          : (unitS == "MeV/g/cm^2")? 0.001*rho
          : (unitS == "MeV/cm")? 0.001
          : (unitS == "GeV/cm")? 1
         /* standard unit for magnetic field in MCfast is Tesla */
          : (unitS == "Tesla")? 1
          : (unitS == "T")? 1
          : (unitS == "kG")? 0.1
          : (unitS == "kGs")? 0.1
          : (unitS == "G")? 0.0001
          : (unitS == "Gs")? 0.0001
         /* fractions are always unit-normalized */
          : (unitS == "percent")? 0.01
          : (unitS == "none" || unitS.size() == 0)? 1
         /* any other dimensions are not valid at this point */
          : 0;
         if (fconvert == 0)
         {
            cerr << "hdds-mcfast: Parameter \"" << var << "\""
                 << " has unrecognized units " << S(unitS) << endl;
            exit(3);
         }
         const char* fltstr = strtok((char*)S(valueS)," ");
         float flt=atof(fltstr);
         modelTable[modelS].db << " " << showpoint << atof(fltstr)*fconvert;
         for (int i=1; i < dim; i++)
         {
            if (fltstr = strtok(0," "))
            {
               modelTable[modelS].db << " " << atof(fltstr)*fconvert;
            }
            else
            {
               modelTable[modelS].db << " 0";
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
            XString typeS("string");
            XString varS(var);
            if ((el = targetTable.lookup(targetEl,typeS,varS)) == 0)
            {
               cerr << "hdds-mcfast: Parameter \"" << var << "\""
                    << " required in template file " << fname << endl
                    << "is missing from HDDS" << endl;
               exit(3);
            }
	    XString valAttS("value");
            valueS = el->getAttribute(X(valAttS));
            modelTable[modelS].db << " \"" << S(valueS) << "\"";
         }
         else
         {
            XString typeS("string_vector");
            XString varS(var);
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
               XString tagS(vectEl->getTagName());
               if (tagS == "string_data")
               {
	          XString valAttS("value");
                  XString valueS(vectEl->getAttribute(X(valAttS)));
                  modelTable[modelS].db << " \"" << S(valueS) << "\"";
                  vcount++;
               }
            }
            for ( ; vcount < dim; vcount++)
            {
               modelTable[modelS].db << " \"\"";
            }
            if (vcount != dim)
            {
	       XString nameAttS("name");
               XString vnameS(el->getAttribute(X(nameAttS)));
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
            XString typeS("reference");
            XString varS(var);
            if ((el = targetTable.lookup(targetEl,typeS,varS)) == 0)
            {
               cerr << "hdds-mcfast: Parameter \"" << var << "\""
                    << " required in template file " << fname << endl
                    << "is missing from HDDS" << endl;
               exit(3);
            }
	    XString valAttS("value");
            XString valueS(el->getAttribute(X(valAttS)));
            modelTable[modelS].db << " \"" << S(valueS) << "\"";
         }
         else
         {
            XString typeS("reference_vector");
            XString varS(var);
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
               XString tagS(vectEl->getTagName());
               if (tagS == "reference_data")
               {
	          XString valAttS("value");
                  XString valueS(vectEl->getAttribute(X(valAttS)));
                  modelTable[modelS].db << " \"" << S(valueS) << "\"";
                  vcount++;
               }
            }
            for ( ; vcount < dim; vcount++)
            {
               modelTable[modelS].db << " \"-\"";
            }
            if (vcount != dim)
            {
	       XString nameAttS("name");
               XString vnameS(el->getAttribute(X(nameAttS)));
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

void DbMaker::makedb(DOMElement* el)
{
   DOMNode* cont;
   for (cont = el->getLastChild();
        cont != 0;
        cont = cont->getPreviousSibling())
   {
      if (cont->getNodeType() != DOMNode::ELEMENT_NODE)
         continue;
      DOMElement* contEl = (DOMElement*)cont;
      XString tagS = contEl->getTagName();
      if (tagS == "mcfast")
      {
         targetTable.add(contEl);
      }
      else
      {
         makedb(contEl);
      }
   }
}

void DbMaker::printdb()
{
   char line[999];
   std::map<const std::string,modelTableEntry>::iterator iter;
   for (iter = modelTable.begin();
        iter != modelTable.end();
        ++iter)
   {
      istringstream idb(iter->second.db.str());
      while (! idb.eof())
      {
         idb.getline(line,999);
         cout << line << endl;
      }
   }
   cout << "end" << endl;
}
