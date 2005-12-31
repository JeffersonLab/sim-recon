/*
 *  hddm-xml :	tool that reads in a HDDM document (Hall D Data Model)
 *		and translates it into plain-text xml.
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

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>

#include "XParsers.hpp"
#include "XString.hpp"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <rpc/rpc.h>
#include <string.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include "particleType.h"

#define X(XString) XString.unicode_str()
#define S(XString) XString.c_str()

class XMLmaker
{
 public:
   ofstream xout; 

   XMLmaker() {};
   ~XMLmaker() {};

   void writeXML(const XString& s);
   void constructXML(XDR* xdrs, DOMElement* el, int depth);
};

void usage()
{
   cerr << "\nUsage:\n"
        << "    hddm-xml [-o <filename>] [HDDM file]\n\n"
        << "Options:\n"
        <<  "    -o <filename>	write to <filename>.xml"
        << endl;
}

int getline(std::stringstream& buf, FILE* stream)
{
   int count;
   for (count=0; 1; ++count)
   {
      char c = fgetc(stream);
      if (c == EOF) break;
      buf << c;
      if (c == '\n') break;
   }
   return count;
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
      cerr << "hddm-xml: Error during initialization! :\n"
           << S(msg) << endl;
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
   FILE* ifp;
   if (argInd == argC)
   {
      ifp = stdin;
   }
   else if (argInd == argC - 1)
   {
      hddmFile = XString(argV[argInd]);
      ifp = fopen(hddmFile.c_str(),"r");
   }
   else
   {
      usage();
      return 1;
   }
   if (!ifp)
   {
      cerr << "hddm-xml: Error opening input stream " << hddmFile << endl;
      exit(1);
   }
   std::ostringstream tmpFileStr;
   tmpFileStr << "tmp" << getpid();
   ofstream ofs(tmpFileStr.str().c_str());
   if (! ofs.is_open())
   {
      cerr << "hddm-xml: Error opening temp file " << tmpFileStr << endl;
      exit(2);
   }

   ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
   XString xmlHeader;
   std::stringstream buf;
   if (getline(buf,ifp))
   {
      XString line;
      std::getline(buf,line);
      if (line.substr(0,5) == "<?xml")
      {
         cerr << "hddm-xml: Error reading input stream " << hddmFile << endl;
         cerr << "Input file appears to be an xml document!" << endl;
         exit(1);
      }
      else if (line.substr(0,5) == "<HDDM")
      {
         xmlHeader = line + "\n";
         ofs << line;
      }
      else
      {
         cerr << "hddm-xml: Input stream does not contain valid hddm header"
              << endl;
         exit(1);
      }
   }
   else
   {
      cerr << "hddm-xml: Error reading from input stream " << hddmFile << endl;
      exit(1);
   }
   while (getline(buf,ifp))
   {
      XString line;
      std::getline(buf,line);
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
      cerr << "hddm-xml : Error parsing HDDM document, "
           << "cannot continue" << endl;
      return 1;
   }
   unlink(tmpFileStr.str().c_str());

   DOMElement* rootEl = document->getDocumentElement();
   XString rootS(rootEl->getTagName());
   if (rootS != "HDDM")
   {
      cerr << "hddm-xml error: root element of input document is "
           << "\"" << S(rootS) << "\", expected \"HDDM\""
           << endl;
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

   XDR* xdrs = new XDR;
   xdrstdio_create(xdrs,ifp,XDR_DECODE);
   int icount;
   while (--reqcount && xdr_int(xdrs,&icount))
   {
      DOMNodeList* contList = rootEl->getChildNodes();
      int contLength = contList->getLength();
      for (int c = 0; c < contLength; c++)
      {
         DOMNode* cont = contList->item(c);
         short type = cont->getNodeType();
         if (type == DOMNode::ELEMENT_NODE)
         {
            DOMElement* contEl = (DOMElement*) cont;
            int size;
	    xdr_int(xdrs,&size);
            if (size > 0)
            {
               builder.constructXML(xdrs,contEl,1);
            }
         }
      }
   }

   builder.writeXML("</HDDM>\n");

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
         cout << s;
      }
}

/* Generate the output xml document according the DOM;
 * at entry the buffer pointer bp points the the word after the word count
 */

void XMLmaker::constructXML(XDR* xdrs, DOMElement* el, int depth)
{
   XString repAttS("maxOccurs");
   XString repS(el->getAttribute(X(repAttS)));
   int rep = (repS == "unbounded")? 9999 : atoi(S(repS));
   if (rep > 1)
   {
      xdr_int(xdrs,&rep);
   }
   else
   {
      rep = 1;
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
            int value;
	    xdr_int(xdrs,&value);
            attrStr << " " << nameS << "=\"" << value << "\"";
         }
	 else if (typeS == "long")
         {
            long long value;
#ifndef XDR_LONGLONG_MISSING
	    xdr_longlong_t(xdrs,&value);
#else
            int* ival = (int*)&value;
# if __BIG_ENDIAN__
	    xdr_int(xdrs,&ival[0]);
	    xdr_int(xdrs,&ival[1]);
# else
	    xdr_int(xdrs,&ival[0]);
	    xdr_int(xdrs,&ival[1]);
# endif
#endif
            attrStr << " " << nameS << "=\"" << value << "\"";
         }
         else if (typeS == "float")
         {
            float value;
	    xdr_float(xdrs,&value);
            attrStr << " " << nameS << "=\"" << value << "\"";
         }
         else if (typeS == "double")
         {
            double value;
	    xdr_double(xdrs,&value);
            attrStr << " " << nameS << "=\"" << value << "\"";
         }
         else if (typeS == "boolean")
         {
            bool_t value;
	    xdr_bool(xdrs,&value);
            attrStr << " " << nameS << "=\"" << value << "\"";
         }
         else if (typeS == "Particle_t")
         {
            Particle_t value;
            xdr_int(xdrs,(int*)&value);
            attrStr << " " << nameS << "=\"" << ParticleType(value) << "\"";
         }
         else if (typeS == "string" || typeS == "anyURI")
         {
            char* value = new char [999];
            xdr_string(xdrs,&value,999);
            attrStr << " " << nameS << "=\"" << value << "\"";
	    delete [] value;
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
	    xdr_int(xdrs,&size);
            if (size > 0)
            {
               constructXML(xdrs,contEl,depth +1);
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
