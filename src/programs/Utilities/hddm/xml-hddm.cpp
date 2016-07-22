/*
 *  xml-hddm :	tool that reads in a plain-text xml data document
 *              and writes it out as a hddm (Hall D Data Model)
 *              following the template provided.
 *
 *  Version 1.2 - Richard Jones, December 2005.
 *  - Updated code to use STL strings and vectors instead of old c-style
 *    pre-allocated arrays and strXXX functions.
 *  - Moved functions into classes grouped by function for better clarity.
 *  - Introduced the XStream class library instead of the direct interface
 *    to the rpc/xdr c-library function.  This also gives access to a nice
 *    integrated set of compression/decompression streambuf classes.
 *
 *  Version 1.1 - Richard Jones, September 2003.
 *  - Updated code to work with the new DOM-2 implementation Xerces-c
 *    from apache.org.  Significant changes have taken place in the API
 *    since DOM-1.
 *  - Added support for new types "long" (int64), "string" (char arrays of
 *    arbitrary length), and "anyURI" (special case of string).
 *  - Switched from native encoding to the use of the XDR library to make
 *    hddm files machine-independent.
 *
 *  Original version - Richard Jones, October 1 2001.
 *
 *
 *  Programmer's Notes:
 *  -------------------
 * 1. The output from xml-hddm is a valid hddm data stream.
 *
 * 2. Two inputs are required: the input xml data document, and a xml
 *    template describing the structure of the data in the document.
 *    Both documents must be well-formed.  In addition, the data document
 *    must conform to the hddm specification document.  Only if both of
 *    these requirements are met is it possible to insure that the data
 *    document can be expressed in a hddm stream.
 *
 * 3. The code has been tested with the xerces-c DOM implementation from
 *    Apache, and is intended to be used with the xerces-c library.
 *
 * 4. Output is sent by default to stdout and can be changed with the
 *    -o option.
 */

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>

#include "XParsers.hpp"
#include "XString.hpp"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <rpc/rpc.h>
#include <unistd.h>

#include <fstream>
#include <string>
#include <sstream>
#include <xstream/xdr.h>

#include "particleType.h"

using namespace xercesc;

#define X(str) XString(str).unicode_str()
#define S(str) str.c_str()

int explicit_repeat_count = 1;

class HDDMmaker
{
 public:
   std::ostream* ofs; 

   HDDMmaker() {};
   ~HDDMmaker() {};

   void constructDocument(DOMElement* el);
   void outputStream(DOMElement* thisEl, DOMElement* modelEl,
                     std::ostream& ofs);
};

void usage()
{
   std::cerr
        << "\nUsage:\n"
        << "    xml-hddm [-o <filename>] -t <template> [input file]\n\n"
        << "Options:\n"
        <<  "    -t <template>    read template from <template>\n"
        <<  "    -o <filename>    write to hddm file <filename>"
        << std::endl;
}


int main(int argC, char* argV[])
{
   std::ifstream* ifs;
   XString templFilename;
   XString outFilename;

   try
   {
      XMLPlatformUtils::Initialize();
   }
   catch (const XMLException* toCatch)
   {
      std::cerr
           << "hddm-xml: Error during initialization! :\n"
           << toCatch->getMessage() << std::endl;
      return 1;
   }

   int argInd;
   for (argInd = 1; argInd < argC; argInd++)
   {
      if (argV[argInd][0] != '-')
      {
         break;
      }
      else if (strcmp(argV[argInd],"-t") == 0)
      {
         templFilename = argV[++argInd];
      }
      else if (strcmp(argV[argInd],"-o") == 0)
      {
         outFilename = argV[++argInd];
      }
      else
      {
         usage();
         return 1;
      }
   }

   if (templFilename.size() == 0)
   {
      usage();
      exit(1);
   }

   HDDMmaker builder;
   if (outFilename.size())
   {
      XString fname(outFilename + ".hddm");
      builder.ofs = new std::ofstream(fname.c_str());
   }
   else
   {
      builder.ofs = &std::cout;
   }

#if defined OLD_STYLE_XERCES_PARSER
   DOMDocument* document = parseInputDocument(templFilename.c_str(),true);
#else
   DOMDocument* document = buildDOMDocument(templFilename.c_str(),true);
#endif
   if (document == 0)
   {
      std::cerr
           << "xml-hddm : Error parsing template HDDM document, "
           << "cannot continue" << std::endl;
      return 1;
   }

   DOMElement* rootEl = document->getDocumentElement();
   XString rootS(rootEl->getTagName());
   if (rootS != "HDDM")
   {
      std::cerr
           << "xml-hddm error: root element of input document is "
           << "\"" << S(rootS) << "\", expected \"HDDM\""
           << std::endl;
      exit(1);
   }

   builder.constructDocument(rootEl);

   XString xmlFile;
   if (argInd == argC)
   {
      ifs = new std::ifstream(0);
   }
   else if (argInd == argC - 1)
   {
      xmlFile = XString(argV[argInd]);
      ifs = new std::ifstream(xmlFile.c_str());
   }
   else
   {
      usage();
      return 1;
   }

   if (! ifs->is_open())
   {
      std::cerr
           << "xml-hddm: Error opening input stream "
           << xmlFile << std::endl;
      exit(1);
   }

   XString xmlHeader;
   XString docHeader;
   XString line;
   if (std::getline(*ifs,line))
   {
      if (line.substr(0,5) != "<?xml")
      {
         std::cerr
              << "xml-hddm: Error reading input stream "
              << xmlFile << std::endl;
         std::cerr
              << "Input file does not appear to be an xml document!"
              << std::endl;
         exit(1);
      }
   }
   else
   {
      std::cerr
           << "xml-hddm: Error reading from input stream " 
           << xmlFile << std::endl;
      exit(1);
   }
   xmlHeader = line;
   if (std::getline(*ifs,line) && line.substr(0,5) != "<HDDM")
   {
      std::cerr
           << "xml-hddm: Input document tag is not HDDM!"
           << std::endl;
      exit(1);
   }
   docHeader = line;

   xstream::xdr::ostream ofx(*builder.ofs);
   std::stringstream tmpFileStr;
   tmpFileStr << "tmp" << getpid();
   while (getline(*ifs,line))
   {
      if (line.size() > 500000)
      {
         std::cerr
              << "xml-hddm: line too long in input document" << std::endl;
         exit(1);
      }

      XString text(line);

      std::ofstream ofs(tmpFileStr.str().c_str());
      if (! ofs.is_open())
      {
         std::cerr
              << "xml-hddm: Error opening temp file "
              << tmpFileStr.str() << std::endl;
         exit(2);
      }
      ofs << xmlHeader << std::endl;
      ofs << docHeader << std::endl;

      while (text.size())
      {
         XString::size_type start = text.find_first_of("<");
         if (start == XString::npos)
         {
            break;
         }
         else if (text.substr(start,2) == "</")
         {
            break;
         }
         XString::size_type end = text.find_first_of('>');
         while (end == XString::npos)
         {
            if (line.size() > 400000)
            {
               std::cerr
                    << "xml-hddm: tag too long in input document" << std::endl;
               exit(1);
            }
            else
            {
               std::getline(*ifs,line);
               text += line;
            }
            end = text.find_first_of('>');
         }
         if (text.substr(end-1,2) == "/>")
         {
            ofs << text.substr(0,end+1) << std::endl;
            text.erase(0,end+1);
         }
         else
         {
            XString endTag;
            endTag = "</" + text.substr(start+1,
                            text.find_first_of(" \t",start)-start-1);
            while (text.find(endTag) == XString::npos)
            {
               ofs << text << std::endl;
               std::getline(*ifs,text);
            }
            ofs << text << std::endl;
         }
         ofs << "</HDDM>" << std::endl;
         ofs.close();

#if defined OLD_STYLE_XERCES_PARSER
         document = parseInputDocument(tmpFile.str().c_str(),false);
#else
         document = buildDOMDocument(tmpFileStr.str().c_str(),false);
#endif
         if (document == 0)
         {
            std::cerr
                 << "xml-hddm : Error parsing HDDM document, "
                 << "cannot continue" << std::endl;
            delete ifs;
            return 1;
         }
         unlink(tmpFileStr.str().c_str());

         DOMElement* thisEl = document->getDocumentElement();

         std::ostringstream ofsbuf;
         xstream::xdr::ostream ofxbuf(ofsbuf);
         builder.outputStream(thisEl,rootEl,ofsbuf);
         int size = (int)ofsbuf.tellp();
         ofx << ((size > 0)? size : 0);
         if (size > 0)
         {
            *builder.ofs << ofsbuf.str();
         }
      }
   }

   if (builder.ofs != &std::cout) {
      ((std::ofstream*)builder.ofs)->close();
   }
   unlink(tmpFileStr.str().c_str());
   XMLPlatformUtils::Terminate();
   return 0;
}

/* Generate the xml document template in normal form */

void HDDMmaker::constructDocument(DOMElement* el)
{
   static int indent = 0;
   for (int n = 0; n < indent; n++)
   {
      *ofs << "  ";
   }
   
   XString tagS(el->getTagName());
   *ofs << "<" << tagS;
   DOMNamedNodeMap* attrList = el->getAttributes();
   int attrListLength = attrList->getLength();
   for (int a = 0; a < attrListLength; a++)
   {
      DOMNode* node = attrList->item(a);
      XString nameS(node->getNodeName());
      XString valueS(node->getNodeValue());
      *ofs << " " << nameS << "=\"" << valueS << "\"";
   }

   DOMNodeList* contList = el->getChildNodes();
   int contListLength = contList->getLength();
   if (contListLength > 0)
   {
      *ofs << ">" << std::endl;
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
      for (int n = 0; n < indent; n++)
      {
         *ofs << "  ";
      }
      *ofs << "</" << tagS << ">" << std::endl;
   }
   else
   {
      *ofs << " />" << std::endl;
   }
}

/* Generate the output binary stream according the HDDM template */

void HDDMmaker::outputStream(DOMElement* thisEl, DOMElement* modelEl,
                             std::ostream& ofs)
{
   XString modelS(modelEl->getTagName());
   XString thisS(thisEl->getTagName());

   DOMNamedNodeMap* modelAttrList = modelEl->getAttributes();
   int modelAttrListLength = modelAttrList->getLength();
   DOMNamedNodeMap* thisAttrList = thisEl->getAttributes();
   int thisAttrListLength = thisAttrList->getLength();
   XString minS(modelEl->getAttribute(X("minOccurs")));
   XString maxS(modelEl->getAttribute(X("maxOccurs")));
   int expectAttrCount = modelAttrList->getLength()
                         - (minS == ""? 0 : 1)
                         - (maxS == ""? 0 : 1);
   if (thisAttrListLength != expectAttrCount)
   {
      std::cerr
           << "xml-hddm: Inconsistency in input xml document" << std::endl
           << "tag " << S(thisS) << " in input document with "
           << thisAttrListLength << " attributes " << std::endl
           << "matched to tag " << S(modelS) << " in hddm template "
           << "with " << expectAttrCount << " attributes." << std::endl;
      exit(1);
   }

   xstream::xdr::ostream ofx(ofs);
   for (int a = 0; a < modelAttrListLength; a++)
   {
      XString attrS(modelAttrList->item(a)->getNodeName());
      XString typeS(modelAttrList->item(a)->getNodeValue());
      XString valueS(thisEl->getAttribute(X(attrS)));
      if (attrS == "maxOccurs" || attrS == "minOccurs")
      {
         continue;
      }
      else if (valueS == "" and typeS != "string")
      {
         std::cerr
              << "xml-hddm: Inconsistency in input xml document" << std::endl
              << "tag " << S(thisS) << " in input document is missing "
              << "attribute " << S(attrS) << std::endl;
         exit(1);
      }
      std::stringstream valueStr(valueS);
      if (typeS == "int")
      {
         int32_t val;
         valueStr >> val;
         ofx << val;
      }
      if (typeS == "long")
      {
         int64_t val;
         valueStr >> val;
         ofx << val;
      }
      else if (typeS == "float")
      {
         float val;
         if (valueS == "nan") {
            val = NAN;
         }
         else if (valueS == "-nan") {
            val = -NAN;
         }
         else if (valueS == "inf") {
            val = INFINITY;
         }
         else if (valueS == "-inf") {
            val = -INFINITY;
         }
         else {
            valueStr >> val;
         }
         ofx << val;
      }
      else if (typeS == "double")
      {
         double val;
         if (valueS == "nan") {
            val = NAN;
         }
         else if (valueS == "-nan") {
            val = -NAN;
         }
         else if (valueS == "inf") {
            val = INFINITY;
         }
         else if (valueS == "-inf") {
            val = -INFINITY;
         }
         else {
            valueStr >> val;
         }
         ofx << val;
      }
      else if (typeS == "boolean")
      {
         int val;
         valueStr >> val;
         ofx << val;
      }
      else if (typeS == "Particle_t")
      {
         int32_t val;
         for (val = 0; val < 99; val++)
         {
            if (valueS == ParticleType((Particle_t)val))
            {
               break;
            }
         }
         ofx << val;
      }
      else if (typeS == "string" || typeS == "anyURI")
      {
         ofx << valueS;
      }
      else
      {
         // other types are treated as comments
      }
   }

   DOMNodeList* thisList = thisEl->getChildNodes();
   int thisListLength = thisList->getLength();
   DOMNodeList* modelList = modelEl->getChildNodes();
   int modelListLength = modelList->getLength();
   for (int m = 0; m < modelListLength; m++)
   {
      DOMNode* mode = modelList->item(m);
      short type = mode->getNodeType();
      if (type == DOMNode::ELEMENT_NODE)
      {
         DOMElement* model = (DOMElement*) mode;
         XString modelS(model->getTagName());
         /*
         XString reqS(model->getAttribute(X("minOccurs")));
	 int req = (reqS == "unbounded")? INT_MAX : 
                   (reqS == "")? 1 :
                   atoi(S(reqS));
         */
         XString repS(model->getAttribute(X("maxOccurs")));
	 int rep = (repS == "unbounded")? INT_MAX :
                   (repS == "")? 1 :
                   atoi(S(repS));
         int repCount=0;
         std::stringstream ofsbuf;
         xstream::xdr::ostream ofxbuf(ofsbuf);
         for (int i = 0; i < thisListLength; i++)
         {
            DOMNode* instance = thisList->item(i);
            short type = instance->getNodeType();
            if (type == DOMNode::ELEMENT_NODE)
            {
               DOMElement* instanceEl = (DOMElement*) instance;
               XString nameS(instanceEl->getTagName());
               if (nameS == modelS)
               {
                  outputStream(instanceEl,model,ofsbuf);
                  if (repCount++ && (rep == 1))
                  {
                     std::cerr
                          << "xml-hddm: Inconsistency in input xml document"
                          << std::endl
                          << "tag " << S(thisS) << " in input document contains"
                          << " multiple instances of tag " << S(nameS)
                          << std::endl
                          << "but it does not have a maxOccurs=\"*\" attribute "
                          << "in the template." << std::endl;
                     exit(1);
                  }
                  else if (repCount > rep) {
                     std::cerr
                          << "xml-hddm: Inconsistency in input xml document"
                          << std::endl
                          << "tag " << S(nameS) << " in the template has "
                          << "maxOccurs=" << rep << std::endl
                          << "but the input document contains more than "
                          << rep << " instances." << std::endl;
                     exit(1);
                  }
               }
            }
         }

         int size = (int)ofsbuf.tellp();
         if (explicit_repeat_count && rep > 1)
         {
            ofx << (int32_t)((size > 0)? size+sizeof(int) : sizeof(int))
                << (int32_t)repCount;
         }
         else
         {
            ofx << (int32_t)((size > 0)? size : 0);
         }
         if (size > 0)
         {
            ofs << ofsbuf.str();
         }
      }
   }
}
