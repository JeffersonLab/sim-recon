/*
 *  xml-hddm :	tool that reads in a plain-text xml data document
 *              and writes it out as a hddm (Hall D Data Model)
 *              following the template provided.
 *
 *  Version 1.2 - Richard Jones, December 2005.
 *  - Updated code to use STL strings and vectors instead of old c-style
 *    pre-allocated arrays and strXXX functions.
 *  - Moved functions into classes grouped by function for better clarity.
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
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <rpc/rpc.h>
#include <unistd.h>

#include <fstream>
#include <string>
#include <sstream>

#include "particleType.h"

#define X(XString) XString.unicode_str()
#define S(XString) XString.c_str()

class HDDMmaker
{
 public:
   FILE* ofd; 
   HDDMmaker() {};
   ~HDDMmaker() {};
   
   void constructDocument(DOMElement* el);
   void outputStream(DOMElement* thisEl, DOMElement* modelEl, XDR* xdrs);
   inline void writeOut(const XString& string)
   {
      fwrite(string.c_str(),1,string.size(),ofd);
   }
};

void usage()
{
   cerr << "\nUsage:\n"
        << "    xml-hddm [-o <filename>] -t <template> [input file]\n\n"
        << "Options:\n"
        <<  "    -t <template>    read template from <template>\n"
        <<  "    -o <filename>    write to hddm file <filename>"
        << endl;
}


int main(int argC, char* argV[])
{
   ifstream* ifs;
   XString templFilename;
   XString outFilename;

   try
   {
      XMLPlatformUtils::Initialize();
   }
   catch (const XMLException* toCatch)
   {
      XString msg(toCatch->getMessage());
      cerr << "hddm-xml: Error during initialization! :\n"
           << S(msg) << endl;
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
      builder.ofd = fopen(fname.c_str(),"w");
   }
   else
   {
      builder.ofd = stdout;
   }

#if defined OLD_STYLE_XERCES_PARSER
   DOMDocument* document = parseInputDocument(templFilename.c_str(),true);
#else
   DOMDocument* document = buildDOMDocument(templFilename.c_str(),true);
#endif
   if (document == 0)
   {
      cerr << "xml-hddm : Error parsing template HDDM document, "
           << "cannot continue" << endl;
      return 1;
   }

   DOMElement* rootEl = document->getDocumentElement();
   XString rootS(rootEl->getTagName());
   if (rootS != "HDDM")
   {
      cerr << "xml-hddm error: root element of input document is "
           << "\"" << S(rootS) << "\", expected \"HDDM\""
           << endl;
      exit(1);
   }

   builder.constructDocument(rootEl);

   XString xmlFile;
   if (argInd == argC)
   {
      ifs = new ifstream(0);
   }
   else if (argInd == argC - 1)
   {
      xmlFile = XString(argV[argInd]);
      ifs = new ifstream(xmlFile.c_str());
   }
   else
   {
      usage();
      return 1;
   }

   if (! ifs->is_open())
   {
      cerr << "xml-hddm: Error opening input stream " << xmlFile << endl;
      exit(1);
   }

   XString xmlHeader;
   XString docHeader;
   XString line;
   if (std::getline(*ifs,line))
   {
      if (line.substr(0,5) != "<?xml")
      {
         cerr << "xml-hddm: Error reading input stream " << xmlFile << endl;
         cerr << "Input file does not appear to be an xml document!" << endl;
         exit(1);
      }
   }
   else
   {
      cerr << "xml-hddm: Error reading from input stream " << xmlFile << endl;
      exit(1);
   }
   xmlHeader = line;
   if (std::getline(*ifs,line) && line.substr(0,5) != "<HDDM")
   {
      cerr << "xml-hddm: Input document tag is not HDDM!"
           << endl;
      exit(1);
   }
   docHeader = line;

   XDR* xdrs = new XDR;
   xdrstdio_create(xdrs,builder.ofd,XDR_ENCODE);
   while (getline(*ifs,line))
   {
      if (line.size() > 500000)
      {
         cerr << "xml-hddm: line too long in input document" << endl;
         exit(1);
      }

      XString text(line);

      std::stringstream tmpFileStr;
      tmpFileStr << "tmp" << getpid();
      ofstream ofs(tmpFileStr.str().c_str());
      if (! ofs.is_open())
      {
         cerr << "xml-hddm: Error opening temp file " << tmpFileStr << endl;
         exit(2);
      }
      ofs << xmlHeader << endl;
      ofs << docHeader << endl;

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
               cerr << "xml-hddm: tag too long in input document" << endl;
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
            ofs << text.substr(0,end+1) << endl;
            text.erase(0,end+1);
            end = 0;
         }
         else
         {
            XString endTag;
            endTag = "</" + text.substr(start+1,
                            text.find_first_of(" \t",start)-start-1);
            while (text.find(endTag) == XString::npos)
            {
               ofs << text << endl;
               std::getline(*ifs,text);
            }
            ofs << text << endl;
         }
         ofs << "</HDDM>" << endl;
         ofs.close();

#if defined OLD_STYLE_XERCES_PARSER
         document = parseInputDocument(tmpFile.str().c_str(),false);
#else
         document = buildDOMDocument(tmpFileStr.str().c_str(),false);
#endif
         if (document == 0)
         {
            cerr << "xml-hddm : Error parsing HDDM document, "
                 << "cannot continue" << endl;
            return 1;
         }
         unlink(tmpFileStr.str().c_str());

         DOMElement* thisEl = document->getDocumentElement();

	 int base = xdr_getpos(xdrs);
	 int first=0;
	 xdr_int(xdrs,&first);
	 first = xdr_getpos(xdrs);
         builder.outputStream(thisEl,rootEl,xdrs);
	 int last = xdr_getpos(xdrs);
	 int size = last - first;
	 xdr_setpos(xdrs,base);
	 xdr_int(xdrs,&size);
	 xdr_setpos(xdrs,last);
      }
   }
   xdr_destroy(xdrs);

   XMLPlatformUtils::Terminate();
   return 0;
}

/* Generate the xml document template in normal form */

void HDDMmaker::constructDocument(DOMElement* el)
{
   static int indent = 0;
   for (int n = 0; n < indent; n++)
   {
      writeOut("  ");
   }
   
   XString tagS(el->getTagName());
   writeOut("<"), writeOut(S(tagS));
   DOMNamedNodeMap* attrList = el->getAttributes();
   int attrListLength = attrList->getLength();
   for (int a = 0; a < attrListLength; a++)
   {
      DOMNode* node = attrList->item(a);
      XString nameS(node->getNodeName());
      XString valueS(node->getNodeValue());
      writeOut(" "), writeOut(S(nameS)),
      writeOut("=\""), writeOut(S(valueS)), writeOut("\"");
   }

   DOMNodeList* contList = el->getChildNodes();
   int contListLength = contList->getLength();
   if (contListLength > 0)
   {
      writeOut(">"), writeOut("\n");
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
         writeOut("  ");
      }
      writeOut("</"), writeOut(S(tagS)), writeOut(">\n");
   }
   else
   {
      writeOut(" />\n");
   }
}

/* Generate the output binary stream according the HDDM template */

void HDDMmaker::outputStream(DOMElement* thisEl, DOMElement* modelEl, XDR* xdrs)
{
   XString modelS(modelEl->getTagName());
   XString thisS(thisEl->getTagName());

   DOMNamedNodeMap* modelAttrList = modelEl->getAttributes();
   int modelAttrListLength = modelAttrList->getLength();
   DOMNamedNodeMap* thisAttrList = thisEl->getAttributes();
   int thisAttrListLength = thisAttrList->getLength();
   XString minAttS("minOccurs");
   XString maxAttS("maxOccurs");
   XString minS(modelEl->getAttribute(X(minAttS)));
   XString maxS(modelEl->getAttribute(X(maxAttS)));
   int expectAttrCount = modelAttrList->getLength()
                         - (minS == ""? 0 : 1)
                         - (maxS == ""? 0 : 1);
   if (thisAttrListLength != expectAttrCount)
   {
      cerr << "xml-hddm: Inconsistency in input xml document" << endl
           << "tag " << S(thisS) << " in input document with "
           << thisAttrListLength << " attributes " << endl
           << "matched to tag " << S(modelS) << " in hddm template "
           << "with " << expectAttrCount << " attributes." << endl;
      exit(1);
   }
   for (int a = 0; a < modelAttrListLength; a++)
   {
      XString attrS(modelAttrList->item(a)->getNodeName());
      XString typeS(modelAttrList->item(a)->getNodeValue());
      XString valueS(thisEl->getAttribute(X(attrS)));
      if (attrS == "maxOccurs" || attrS == "minOccurs")
      {
         continue;
      }
      else if (valueS == "")
      {
         cerr << "xml-hddm: Inconsistency in input xml document" << endl
              << "tag " << S(thisS) << " in input document is missing "
              << "attribute " << S(attrS) << endl;
         exit(1);
      }
      if (typeS == "int")
      {
         int val;
         sscanf(S(valueS),"%d",&val);
         xdr_int(xdrs,&val);
      }
      if (typeS == "long")
      {
         long long val;
         sscanf(S(valueS),"%lld",&val);
#ifndef XDR_LONGLONG_MISSING
         xdr_longlong_t(xdrs,&val);
#else
         int* ival = (int*)&val;
# if __BIG_ENDIAN__
         xdr_int(xdrs,&ival[0]);
         xdr_int(xdrs,&ival[1]);
# else
         xdr_int(xdrs,&ival[1]);
         xdr_int(xdrs,&ival[0]);
# endif
#endif
      }
      else if (typeS == "float")
      {
         float val;
         sscanf(S(valueS),"%g",&val);
         xdr_float(xdrs,&val);
      }
      else if (typeS == "double")
      {
         double val;
         sscanf(S(valueS),"%Lg",&val);
         xdr_double(xdrs,&val);
      }
      else if (typeS == "boolean")
      {
         int val;
         sscanf(S(valueS),"%d",&val);
         xdr_int(xdrs,&val);
      }
      else if (typeS == "Particle_t")
      {
         int val;
         for (val = 1; val < 99; val++)
         {
            if (valueS == ParticleType((Particle_t)val))
            {
               break;
            }
         }
         xdr_int(xdrs,&val);
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
	 int start=0;
         int base = xdr_getpos(xdrs);
	 xdr_int(xdrs,&start);
         start = xdr_getpos(xdrs);
         DOMElement* model = (DOMElement*) mode;
         XString modelS(model->getTagName());
	 XString repAttS("maxOccurs");
         XString repS(model->getAttribute(X(repAttS)));
	 int rep = (repS == "unbounded")? 9999 : atoi(S(repS));
         int repCount=0;
         if (rep > 1)
         {
            xdr_int(xdrs,&repCount);
         }
         repCount = 0;
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
                  outputStream(instanceEl,model,xdrs);
                  if (repCount++ && (rep == 0))
                  {
                     cerr << "xml-hddm: Inconsistency in input xml document"
                          << endl
                          << "tag " << S(thisS) << " in input document contains"
                          << " multiple instances of tag " << S(nameS) << endl
                          << "but it does not have a maxOccurs=\"*\" attribute "
                          << "in the template." << endl;
                     exit(1);
                  }
               }
            }
         }
	 int end = xdr_getpos(xdrs);
         int size = end - start;
	 xdr_setpos(xdrs,base);
	 xdr_int(xdrs,&size);
	 if (rep > 1)
	 {
	    xdr_int(xdrs,&repCount);
	 }
	 xdr_setpos(xdrs,end);
      }
   }
}
