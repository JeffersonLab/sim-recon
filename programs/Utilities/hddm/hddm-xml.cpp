/*
 *  hddm-xml :	tool that reads in a HDDM document (Hall D Data Model)
 *		and translates it into plain-text xml.
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
 *  Original version - Richard Jones, June 4 2001.
 *
 *
 *  Programmer's Notes:
 *  -------------------
 * 1. The output from hddm-xml is a well-formed xml document.
 *
 * 2. The hddm stream contains a xml header that functions as a prototype
 *    of the xml output.
 *
 * 3. This tool can read any hddm stream.  No recompilation is required.
 *
 * 4. The code has been tested with the xerces-c DOM implementation from
 *    Apache, and is intended to be used with the xerces-c library.
 *
 * 5. Output is sent by default to stdout and can be changed with the
 *    -o option.
 */

//#define VERBOSE_HDDM_LOGGING 1

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>

#include "XParsers.hpp"
#include "XString.hpp"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <rpc/rpc.h>
#include <unistd.h>
#include <xstream/xdr.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include "particleType.h"


#define X(str) XString(str).unicode_str()
#define S(str) str.c_str()

using namespace xercesc;

class XMLmaker
{
 public:
   std::ofstream xout; 

   XMLmaker() {};
   ~XMLmaker() {};

   void writeXML(const XString& s);
   void constructXML(xstream::xdr::istream& ifx, DOMElement* el, int depth);
};

void usage()
{
   std::cerr
        << "\nUsage:\n"
        << "    hddm-xml [-o <filename>] [HDDM file]\n\n"
        << "Options:\n"
        <<  "    -o <filename>	write to <filename>.xml"
        << std::endl;
}


int main(int argC, char* argV[])
{
   XString xFilename;

   try
   {
      XMLPlatformUtils::Initialize();
   }
   catch (const XMLException* toCatch)
   {
      XString msg(toCatch->getMessage());
      std::cerr
           << "hddm-xml: Error during initialization! :\n"
           << S(msg) << std::endl;
      return 1;
   }

   int reqcount=0;
   int argInd;
   for (argInd = 1; argInd < argC; argInd++)
   {
      if (argV[argInd][0] != '-')
      {
         break;
      }
      else if (strcmp(argV[argInd],"-o") == 0)
      {
         xFilename = argV[++argInd];
      }
      else if (sscanf(argV[argInd],"%d",&reqcount))
      {
	 reqcount = 1-reqcount;
      }
      else
      {
         usage();
         return 1;
      }
   }

   XString hddmFile;
   std::istream* ifs;
   if (argInd == argC)
   {
      ifs = &std::cin;
   }
   else if (argInd == argC - 1)
   {
      hddmFile = XString(argV[argInd]);
      ifs = new std::ifstream(hddmFile.c_str());
   }
   else
   {
      usage();
      return 1;
   }
   if (!ifs->good())
   {
      std::cerr
           << "hddm-xml: Error opening input stream " << hddmFile << std::endl;
      exit(1);
   }
   std::ostringstream tmpFileStr;
   tmpFileStr << "tmp" << getpid();
   std::ofstream ofs(tmpFileStr.str().c_str());
   if (! ofs.is_open())
   {
      std::cerr
           << "hddm-xml: Error opening temp file " << tmpFileStr << std::endl;
      exit(2);
   }

   ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
   XString xmlHeader;
   XString line;
   if (getline(*ifs,line))
   {
      if (line.substr(0,5) == "<?xml")
      {
         std::cerr
              << "hddm-xml: Error reading input stream " << hddmFile
              << std::endl;
         std::cerr
              << "Input file appears to be an xml document!" << std::endl;
         exit(1);
      }
      else if (line.substr(0,5) == "<HDDM")
      {
         xmlHeader = line + "\n";
         ofs << line;
      }
      else
      {
         std::cerr
              << "hddm-xml: Input stream does not contain valid hddm header"
              << std::endl;
         exit(1);
      }
   }
   else
   {
      std::cerr
           << "hddm-xml: Error reading from input stream " << hddmFile 
           << std::endl;
      exit(1);
   }
   while (getline(*ifs,line))
   {
      ofs << line;
      if (line == "</HDDM>")
      {
         break;
      }
   }
   ofs.close();

#if defined OLD_STYLE_XERCES_PARSER
   DOMDocument* document = parseInputDocument(tmpFileStr.str().c_str(),false);
#else
   DOMDocument* document = buildDOMDocument(tmpFileStr.str().c_str(),false);
#endif
   if (document == 0)
   {
      std::cerr
           << "hddm-xml : Error parsing HDDM document, "
           << "cannot continue" << std::endl;
      return 1;
   }
   unlink(tmpFileStr.str().c_str());

   DOMElement* rootEl = document->getDocumentElement();
   XString rootS(rootEl->getTagName());
   if (rootS != "HDDM")
   {
      std::cerr
           << "hddm-xml error: root element of input document is "
           << "\"" << S(rootS) << "\", expected \"HDDM\""
           << std::endl;
      return 1;
   }

   XMLmaker builder;
   if (xFilename.size())
   {
      XString fname(xFilename + ".xml");
      builder.xout.open(fname.c_str());
   }
  
   builder.writeXML("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
   builder.writeXML(xmlHeader);

   xstream::xdr::istream ifx(*ifs);
   while (--reqcount && ifs->good())
   {
      DOMNodeList* contList = rootEl->getChildNodes();
      int contLength = contList->getLength();
      int tsize;
      ifx >> tsize;
#ifdef VERBOSE_HDDM_LOGGING
      XString tnameS(rootEl->getTagName());
      std::cerr << "hddm-xml : tag " << S(tnameS)
                << " found with size " << tsize
                << std::endl;
#endif
      if (tsize <= 0)
      {
         break;
      }
      for (int c = 0; c < contLength; c++)
      {
         DOMNode* cont = contList->item(c);
         short type = cont->getNodeType();
         if (type == DOMNode::ELEMENT_NODE)
         {
            DOMElement* contEl = (DOMElement*) cont;
            int size;
            ifx >> size;
#ifdef VERBOSE_HDDM_LOGGING
            XString cnameS(contEl->getTagName());
            std::cerr << "hddm-xml : tag " << S(cnameS)
                      << " found with size " << size
                      << std::endl;
#endif
            if (size > 0)
            {
               builder.constructXML(ifx,contEl,1);
            }
            else {
               XString repS(contEl->getAttribute(X("minOccurs")));
               int rep = (repS == "")? 1 : atoi(S(repS));
               if (size != 0) {
                  std::cerr << "hddm-xml : weird size, continue? [y/n] ";
                  std::string ans;
                  std::cin >> ans;
                  if (ans[0] != 'y' && ans[0] != 'Y') {
                     exit(5);
                  }
               }
            }
         }
      }
   }

   builder.writeXML("</HDDM>\n");

   if (ifs != &std::cin)
   {
      ((std::ifstream*)ifs)->close();
   }
   XMLPlatformUtils::Terminate();
   return 0;
}

/* write a string to xml output stream, either stdout or a file */

void XMLmaker::writeXML(const XString& s)
{
      if (xout.is_open())
      {
         xout << s;
      }
      else
      {
         std::cout << s;
      }
}

/* Generate the output xml document according the DOM;
 * at entry the buffer pointer bp points the the word after the word count
 */

void XMLmaker::constructXML(xstream::xdr::istream& ifx,
                            DOMElement* el, int depth)
{
   XString repS(el->getAttribute(X("maxOccurs")));
   int rep = (repS == "unbounded")? 9999 :
             (repS == "")? 1 :
             atoi(S(repS));
   if (rep > 1)
   {
      ifx >> rep;
   }

   for (int r = 0; r < rep; r++)
   {
      XString tagS(el->getTagName());
      for (int d = 0; d < depth; d++)
      {
         writeXML("  ");
      }
      writeXML("<");
      writeXML(S(tagS));
      DOMNamedNodeMap* attrList = el->getAttributes();
      int listLength = attrList->getLength();
      for (int a = 0; a < listLength; a++)
      {
         XString nameS(attrList->item(a)->getNodeName());
         XString typeS(attrList->item(a)->getNodeValue());
         std::ostringstream attrStr;
         if (typeS == "int")
         {
            int32_t value;
	    ifx >> value;
            attrStr << " " << nameS << "=\"" << value << "\"";
         }
	 else if (typeS == "long")
         {
            int64_t value;
            ifx >> value;
            attrStr << " " << nameS << "=\"" << value << "\"";
         }
         else if (typeS == "float")
         {
            float value;
            ifx >> value;
            attrStr << " " << nameS << "=\"" << value << "\"";
         }
         else if (typeS == "double")
         {
            double value;
            ifx >> value;
            attrStr << " " << nameS << "=\"" << value << "\"";
         }
         else if (typeS == "boolean")
         {
            bool_t value;
            ifx >> value;
            attrStr << " " << nameS << "=\"" << value << "\"";
         }
         else if (typeS == "Particle_t")
         {
            int32_t value;
            ifx >> value;
            attrStr << " " << nameS << "=\"" << ParticleType((Particle_t)value) << "\"";
         }
         else if (typeS == "string" || typeS == "anyURI")
         {
            std::string value;
            ifx >> value;
            attrStr << " " << nameS << "=\"" << value << "\"";
         }
         else if (nameS == "minOccurs" || nameS == "maxOccurs")
         {
            ;
         }
         else
         {
            attrStr << " " << nameS << "=\"" << typeS << "\"";
         }
         writeXML(attrStr.str());
      }

      DOMNodeList* contList = el->getChildNodes();
      int contLength = contList->getLength();
      if (contLength > 1)
      {
         writeXML(">\n");
      }
      else
      {
         writeXML(" />\n");
      }

      for (int c = 0; c < contLength; c++)
      {
         DOMNode* cont = contList->item(c);
         short type = cont->getNodeType();
         if (type == DOMNode::ELEMENT_NODE)
         {
            DOMElement* contEl = (DOMElement*) cont;
            int size;
            ifx >> size;
#ifdef VERBOSE_HDDM_LOGGING
            XString cnameS(contEl->getTagName());
            std::cerr << "hddm-xml : tag " << S(cnameS)
                      << " found with size " << size
                      << std::endl;
#endif
            if (size > 0) {
               constructXML(ifx,contEl,depth +1);
            }
            else {
               XString irepS(contEl->getAttribute(X("minOccurs")));
               int irep = (irepS == "")? 1 : atoi(S(irepS));
               if (size != 0) {
                  std::cerr << "hddm-xml : weird size, continue? [y/n] ";
                  std::string ans;
                  std::cin >> ans;
                  if (ans[0] != 'y' && ans[0] != 'Y') {
                     exit(5);
                  }
               }
            }
         }
      }

      if (contLength > 1)
      {
         for (int d = 0; d < depth; d++)
         {
            writeXML("  ");
         }
         XString endTag("</"+tagS+">\n");
         writeXML(endTag);
      }
   }
}
