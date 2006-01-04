/*
 *  hddm-cpp :	tool that reads in a HDDM document (Hall D Data Model)
 *		and writes a c++ header file that embodies the model in
 *		c++ classes.  It also generates input/output member
 *		functions to translate the model between the memory
 *		representation and a default binary representation that
 *		is suitable for passing over a pipe or storing on disk.
 *
 *  Original version - Richard Jones, May 25 2001.
 *
 *
 *  Programmer's Notes:
 *  -------------------
 * This translator has yet to be implemented, but it should be simple
 * to do.  Just take the hddm-c.cpp code and make the following changes.
 *
 * 1. In the header file replace the make_x_yyy function declarations
 *    with class definitions named x_yyy, using the signature of the
 *    make_x_yyy for the constructor, and the default destructor.
 *
 * 2. Implement the data members of the c structures as public data 
 *    members of the parent class, for efficient access.
 *
 * 3. The pointers become pointers to class instance.
 *
 * 4. The i/o functions are implemented as member functions of the
 *    by top-level class x_HDDM.
 *
 * 5. The pop-stack helper functions are implemented in an auxilliary
 *    class popstack.
 *
 */


#define MAX_POPLIST_LENGTH 99

#include "XString.hpp"
#include "XParsers.hpp"
#include <assert.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>

#include <fstream>

#define X(XString) XString.unicode_str()
#define S(XString) XString.c_str()

using namespace xercesc;

char* hFilename = 0;
std::ofstream hFile;
std::ofstream cFile;

const char* classPrefix;
int tagListLength = 0;
DOMElement* tagList[100000];
bool verifyOnly = false;
char containerList[65536]="";
char constructorCalls[65536]="";
char destructorCalls[65536]="";

void usage();
char* plural(const char* str);
char* simpleStructType(const char* tag);
char* listStructType(const char* tag);
void checkConsistency(DOMElement* el, int t);
void writeHeader(DOMElement* el);
int constructGroup(DOMElement* el);
void constructMakeFuncs();
void constructUnpackers();
void constructReadFunc(DOMElement* topEl);
void constructPackers();
void constructFlushFunc(DOMElement* el);
void writeMatcher();
void constructOpenFunc(DOMElement* el);
void constructInitFunc(DOMElement* el);
void constructCloseFunc(DOMElement* el);
void constructDocument(DOMElement* el);

void usage()
{
   std::cerr
        << "\nUsage:\n"
        << "    hddm-cpp [-v | -o <filename>] {HDDM file}\n\n"
        << "Options:\n"
        <<  "    -v			validate only\n"
        <<  "    -o <filename>	write to <filename>.h"
        << std::endl;
}

/* Generate the plural form of a noun */

char* plural(const char* str)
{
   int len = strlen(str);
   char* p = new char [len+10];
   strcpy(p,str);
   if ((len > 3) && (strcmp(&p[len-3],"tum")  == 0))
   {
      strcpy(&p[len-3],"ta");
   }
   else if ((len > 2) && (strcmp(&p[len-2],"ex")  == 0))
   {
      strcpy(&p[len-2],"ices");
   }
   else if ((len > 2) && (strcmp(&p[len-2],"sh")  == 0))
   {
      strcpy(&p[len-2],"shes");
   }
   else if ((len > 1) && (strcmp(&p[len-1],"s")  == 0))
   {
      strcpy(&p[len-1],"ses");
   }
   else if (len > 1)
   {
      strcat(p,"s");
   }
   return p;
}

/* Map from tag name to name of the corresponding c++-class
 * for the case of simple tags (those that do not repeat)
 */
char* simpleStructType(const char* tag)
{
   int len = strlen(tag) + strlen(classPrefix);
   char* p = new char [len + 10];
   char* q = new char [len];
   strcpy(q,tag);
   q[0] = toupper(q[0]);
   sprintf(p,"%s_%s_t",classPrefix,q);
   delete [] q;
   return p;
}

/* Map from tag name to name of the corresponding c++-class
 * for the case of list tags (those that may repeat)
 */
char* listStructType(const char* tag)
{
   int len = strlen(tag) + strlen(classPrefix);
   char* p = new char [len + 10];
   char* tags = plural(tag);
   tags[0] = toupper(tags[0]);
   sprintf(p,"%s_%s_t",classPrefix,tags);
   delete [] tags;
   return p;
}

/* Verify that the tag group under this element does not collide
 * with existing tag group t, otherwise exit with fatal error
 */
void checkConsistency(DOMElement* el, int t)
{
   XString tagS(el->getTagName());
   DOMNamedNodeMap* oldAttr = tagList[t]->getAttributes();
   DOMNamedNodeMap* newAttr = el->getAttributes();
   int listLength = oldAttr->getLength();
   for (int n = 0; n < listLength; n++)
   {
      XString nameS(oldAttr->item(n)->getNodeName());
      XString oldS(tagList[t]->getAttribute(X(nameS)));
      XString newS(el->getAttribute(X(nameS)));
      if (nameS == "minOccurs")
      {
         continue;
      }
      else if (nameS == "maxOccurs")
      {
         int maxold = (oldS == "unbounded")? 9999 : atoi(S(oldS));
         int maxnew = (newS == "unbounded")? 9999 : atoi(S(newS));
	 if (maxold*maxnew <= maxold)
         {
            std::cerr
                 << "hddm-cpp error: inconsistent maxOccurs usage by tag "
                 << "\"" << S(tagS) << "\" in xml document." << std::endl;
            exit(1);
         }
      }
      else if (newS != oldS)
      {
         std::cerr
              << "hddm-cpp error: inconsistent usage of attribute "
              << "\"" << S(nameS) << "\" in tag "
              << "\"" << S(tagS) << "\" in xml document." << std::endl;
         exit(1);
      }
   }
   listLength = newAttr->getLength();
   for (int n = 0; n < listLength; n++)
   {
      XString nameS(newAttr->item(n)->getNodeName());
      XString oldS(tagList[t]->getAttribute(X(nameS)));
      XString newS(el->getAttribute(X(nameS)));
      if (nameS == "minOccurs")
      {
         continue;
      }
      else if (nameS == "maxOccurs")
      {
         int maxold = (oldS == "unbounded")? 9999 : atoi(S(oldS));
         int maxnew = (newS == "unbounded")? 9999 : atoi(S(newS));
	 if (maxold*maxnew <= maxnew)
         {
            std::cerr
                 << "hddm-cpp error: inconsistent maxOccurs usage by tag "
                 << "\"" << S(tagS) << "\" in xml document." << std::endl;
            exit(1);
         }
      }
      else if (newS != oldS)
      {
         std::cerr
              << "hddm-cpp error: inconsistent usage of attribute "
              << "\"" << S(nameS) << "\" in tag "
              << "\"" << S(tagS) << "\" in xml document." << std::endl;
         exit(1);
      }
   }
   DOMNodeList* oldList = tagList[t]->getChildNodes();
   DOMNodeList* newList = el->getChildNodes();
   listLength = oldList->getLength();
   if (newList->getLength() != listLength)
   {
      std::cerr
           << "hddm-cpp error: inconsistent usage of tag "
           << "\"" << S(tagS) << "\" in xml document." << std::endl;
   exit(1);
   }
   for (int n = 0; n < listLength; n++)
   {
      DOMNode* cont = oldList->item(n);
      XString nameS(cont->getNodeName());
      short type = cont->getNodeType();
      if (type == DOMNode::ELEMENT_NODE)
      {
         DOMNodeList* contList = el->getElementsByTagName(X(nameS));
         if (contList->getLength() != 1)
         {
             std::cerr
                  << "hddm-cpp error: inconsistent usage of tag "
                  << "\"" << S(tagS) << "\" in xml document." << std::endl;
             exit(1);
         }
      }
   }
}

/* Write declaration of c-structure for this tag to c-header file */

void writeHeader(DOMElement* el)
{
   XString tagS(el->getTagName());
   char* ctypeDef = simpleStructType(S(tagS));

   XString repAttS("maxOccurs");
   XString repS(el->getAttribute(X(repAttS)));
   int rep = (repS == "unbounded")? 9999 : atoi(S(repS));
   if (rep > 1)
   {
      char* ctypeRef = listStructType(S(tagS));
		char cDef[256];
		strcpy(cDef,ctypeDef);
		cDef[strlen(cDef)-2] = 0;
		char cRef[256];
		strcpy(cRef,&ctypeRef[2]);
		cRef[strlen(cRef)-2] = 0;
		char cpptypeRef[256];
		strcpy(cpptypeRef,ctypeRef);
		cpptypeRef[strlen(cpptypeRef)-1] = 'c';
		hFile << "//------------- "<<cpptypeRef<<" --------------"<<std::endl
            << std::endl << "class "<<cpptypeRef<<":public DContainer"<<std::endl
		      << "{" << std::endl
				<< "	public:"<<std::endl
            << "		"<<cpptypeRef<<"(void)"
				<< ":DContainer((void**)&"<<cDef<<", sizeof("<<ctypeDef<<"), \""<<cDef<<"\"){}" << std::endl
            << "   	" << ctypeDef << " *"<<cDef<<";" << std::endl
            << "};"<< std::endl
				<< "//-------------------------------------------"<<std::endl;
		sprintf(containerList,"%s\n	%s *%s;",containerList,cpptypeRef, cRef);
		sprintf(constructorCalls,"%s\n	hddm->%s 	= new %s();", constructorCalls, cRef, cpptypeRef);
		sprintf(destructorCalls,"%s\n	delete hddm->%s;", destructorCalls, cRef);
      delete [] ctypeRef;
   }

   //hFile << "#endif /* " << ctypeDef << " */" 			<< std::endl;
}

/* Generate c-structure declarations for this tag and its descendants;
 * this function calls itself recursively
 */
int constructGroup(DOMElement* el)
{
   XString tagS(el->getTagName());
   int t;
   for (t = 0; t < tagListLength; t++)
   {
      XString targS(tagList[t]->getTagName());
      if (tagS == targS)
      {
         checkConsistency(el,t);
         return t;
      }
   }

   tagList[t] = el;
   tagListLength++;

   DOMNodeList* contList = el->getChildNodes();
   int contLength = contList->getLength();
   for (int c = 0; c < contLength; c++)
   {
      DOMNode* cont = contList->item(c);
      short type = cont->getNodeType();
      if (type == DOMNode::ELEMENT_NODE)
      {
         DOMElement* contEl = (DOMElement*) cont;
         constructGroup(contEl);
      }
   }

   writeHeader(el);
   return t;
}

/* Generate the xml template in normal form and store in a string */

void constructDocument(DOMElement* el)
{
   static int indent = 0;
   cFile << "\"";
   for (int n = 0; n < indent; n++)
   {
      cFile << "  ";
   }
   
   XString tagS(el->getTagName());
   cFile << "<" << S(tagS);
   DOMNamedNodeMap* attrList = el->getAttributes();
   int attrListLength = attrList->getLength();
   for (int a = 0; a < attrListLength; a++)
   {
      DOMNode* node = attrList->item(a);
      XString nameS(node->getNodeName());
      XString valueS(node->getNodeValue());
      cFile << " " << S(nameS) << "=\\\"" << S(valueS) << "\\\"";
   }

   DOMNodeList* contList = el->getChildNodes();
   int contListLength = contList->getLength();
   if (contListLength > 0)
   {
      cFile << ">\\n\"" << std::endl;
      indent++;
      for (int c = 0; c < contListLength; c++)
      {
         DOMNode* node = contList->item(c);
         if (node->getNodeType() == DOMNode::ELEMENT_NODE)
         {
            DOMElement* contEl = (DOMElement*) node;
            constructDocument(contEl);
         }
      }
      indent--;
      cFile << "\"";
      for (int n = 0; n < indent; n++)
      {
         cFile << "  ";
      }
      cFile << "</" << S(tagS) << ">\\n\"" << std::endl;
   }
   else
   {
      cFile << " />\\n\"" << std::endl;
   }
}

int main(int argC, char* argV[])
{
   try
   {
      XMLPlatformUtils::Initialize();
   }
   catch (const XMLException* toCatch)
   {
      XString msg(toCatch->getMessage());
      std::cerr
           << "hddm-cpp: Error during initialization! :\n"
           << S(msg) << std::endl;
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
   int argInd;
   for (argInd = 1; argInd < argC; argInd++)
   {
      if (argV[argInd][0] != '-')
      {
         break;
      }
      if (strcmp(argV[argInd],"-v") == 0)
      {
         verifyOnly = true;
      }
      else if (strcmp(argV[argInd],"-o") == 0)
      {
         hFilename = argV[++argInd];
      }
      else
      {
         std::cerr
              << "Unknown option \'" << argV[argInd]
              << "\', ignoring it\n" << std::endl;
      }
   }

   if (argInd != argC - 1)
   {
      usage();
      return 1;
   }
   xmlFile = argV[argInd];

#if defined OLD_STYLE_XERCES_PARSER
   DOMDocument* document = parseInputDocument(xmlFile,false);
#else
   DOMDocument* document = buildDOMDocument(xmlFile,false);
#endif
   if (document == 0)
   {
      std::cerr
           << "hddm-cpp : Error parsing HDDM document, "
           << "cannot continue" << std::endl;
      return 1;
   }

   DOMElement* rootEl = document->getDocumentElement();
   XString rootS(rootEl->getTagName());
   if (rootS != "HDDM")
   {
      std::cerr
           << "hddm-cpp error: root element of input document is "
           << "\"" << S(rootS) << "\", expected \"HDDM\""
           << std::endl;
      return 1;
   }

   XString classAttS("class");
   XString classS(rootEl->getAttribute(X(classAttS)));
   classPrefix = S(classS);

   char hname[510];
	char hppname[510];
   if (verifyOnly)
   {
      sprintf(hname,"/dev/null");
   }
   else if (hFilename)
   {
      sprintf(hname,"%s.h",hFilename);
   }
   else
   {
      sprintf(hname,"hddm_%s.h",classPrefix);
   }
	sprintf(hppname,"%spp",hname);

   hFile.open(hppname);
   if (! hFile.is_open())
   {
      std::cerr
           << "hddm-cpp error: unable to open output file "
           << hppname << std::endl;
      return 1;
   }

   char cname[510];
	sprintf(cname,"hddm_containers_%s.cc",classPrefix);
   cFile.open(cname);
   if (! cFile.is_open())
   {
      std::cerr
           << "hddm-cpp error: unable to open output file "
           << cname << std::endl;
      return 1;
   }

   hFile << "/*"						<< std::endl
	 << " * " << hppname << " - DO NOT EDIT THIS FILE"	<< std::endl
	 << " *"						<< std::endl
	 << " * This file was generated automatically by hddm-cpp"
	 << " from the file"					<< std::endl
    << " * " << xmlFile					<< std::endl
    << " * This header file defines the c++ structures that"
	 << " hold the data"					<< std::endl
	 << " * described in the data model"
    << " (from " << xmlFile << "). "			<< std::endl
	 << " *"						<< std::endl
	 << " * The hddm data model tool set was written by"	<< std::endl
	 << " * Richard Jones, University of Connecticut."	<< std::endl
	 << " *"						<< std::endl
	 << " *"						<< std::endl
	 << " * The C++ container system was written by"	<< std::endl
	 << " * David Lawrence, Jefferson Lab."	<< std::endl
	 << " *"						<< std::endl
	 << " * For more information see the following web site"<< std::endl
	 << " *"						<< std::endl
	 << " * http://zeus.phys.uconn.edu/halld/datamodel/doc"	<< std::endl
	 << " *"						<< std::endl
	 << " */"						<< std::endl
	 							<< std::endl;

   cFile	<< "/*"						<< std::endl
			<< " * " << cname << " - DO NOT EDIT THIS FILE"	<< std::endl
			<< " *"						<< std::endl
			<< " * This file was generated automatically by hddm-cpp"
			<< " from the file"					<< std::endl
			<< " * " << xmlFile					<< std::endl
			<< " *"						<< std::endl
			<< " * The hddm data model tool set was written by"	<< std::endl
			<< " * Richard Jones, University of Connecticut."	<< std::endl
			<< " *"						<< std::endl
			<< " *"						<< std::endl
			<< " * The C++ container system was written by"	<< std::endl
			<< " * David Lawrence, Jefferson Lab."	<< std::endl
			<< " *"						<< std::endl
			<< " * For more information see the following web site"<< std::endl
			<< " *"						<< std::endl
			<< " * http://zeus.phys.uconn.edu/halld/datamodel/doc"	<< std::endl
			<< " */"						<< std::endl
			<< std::endl;

   hFile << "#include \""<<hname<<"\"" 		<< std::endl
			<< "#include \"DContainer.h\"" 		<< std::endl
			<< std::endl
			<< "#ifndef _HDDM_HPP_"<<std::endl
			<< "#define _HDDM_HPP_"<<std::endl;

   constructGroup(rootEl);
	
	hFile	<< std::endl
			<< std::endl
			<< "//----------------------------------------------------------------------------"<<std::endl
			<< "//------------------------------- hddm_containers_t -------------------------------"<<std::endl
			<< "//----------------------------------------------------------------------------"<<std::endl
			<< "typedef struct{"<<std::endl
			<< std::endl
			<< "	/// This struct should consist ONLY of class pointers derived from DContainer"<<std::endl
			<< "	"<<containerList<<std::endl
			<< "}"<<classPrefix<<"_hddm_containers_t;"<<std::endl
			<< std::endl
			<< "// in HDDM/hddm_containers.cc"<<std::endl
			<< "derror_t init_hddm_containers("<<classPrefix<<"_hddm_containers_t *hddm);"<<std::endl
			<< "derror_t delete_hddm_containers("<<classPrefix<<"_hddm_containers_t *hddm);"<<std::endl
			<< std::endl
			<< "#endif // _HDDM_HPP_"<<std::endl
			<< std::endl;

   cFile << "#include \"" << hppname << "\"" 			<< std::endl
			<< "//----------------------"<<std::endl
			<< "// init_hddm_containers_t"<<std::endl
			<< "//----------------------"<<std::endl
			<< "derror_t init_hddm_containers("<<classPrefix<<"_hddm_containers_t *hddm)"<<std::endl
			<< "{"<<std::endl
			<< "	/// Call constructors for all DContainer derived classes"<<std::endl
			<< constructorCalls<<std::endl
			<< "}"<<std::endl
			<< std::endl
			<< "//----------------------"<<std::endl
			<< "// delete_hddm_containers_t"<<std::endl
			<< "//----------------------"<<std::endl
			<< "derror_t delete_hddm_containers("<<classPrefix<<"_hddm_containers_t *hddm)"<<std::endl
			<< "{"<<std::endl
			<< "	/// Call destructors for all DContainer derived classes"<<std::endl
			<< destructorCalls<<std::endl
			<< "}"<<std::endl
			<< std::endl;

   XMLPlatformUtils::Terminate();
   return 0;
}
