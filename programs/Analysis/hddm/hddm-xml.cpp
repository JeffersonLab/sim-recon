/*
 *  hddm-xml :	tool that reads in a HDDM document (Hall D Data Model)
 *		and translates it into plain-text xml.
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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE true
#endif

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>

#include "XParsers.hpp"
#include "XString.hpp"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <rpc/xdr.h>
#include <string.h>
#include <unistd.h>

#include <fstream>

#include "hddm-xml.hpp"
#include "particleType.h"

#define X(XString) XString.unicodeForm()
#define S(XString) XString.localForm()

ofstream xout; 

void usage()
{
   cerr << "\nUsage:\n"
        << "    hddm-xml [-o <filename>] [HDDM file]\n\n"
        << "Options:\n"
        <<  "    -o <filename>	write to <filename>.xml"
        << endl;
}

/* write a string to xml output stream, either stdout or a file */

void writeXML(const char* s)
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

void constructXML(XDR* xdrs, DOMElement* el, int depth)
{
   XString repAttS("maxOccurs");
   XString repS(el->getAttribute(X(repAttS)));
   int rep = (repS.equals("unbounded"))? 9999 : atoi(S(repS));
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
         char attrStr[500];
         if (typeS.equals("int"))
         {
            int value;
	    xdr_int(xdrs,&value);
            sprintf(attrStr," %s=\"%d\"",S(nameS),value);
         }
	 else if (typeS.equals("long"))
         {
            long long value;
	    xdr_longlong_t(xdrs,&value);
            sprintf(attrStr," %s=\"%lld\"",S(nameS),value);
         }
         else if (typeS.equals("float"))
         {
            float value;
	    xdr_float(xdrs,&value);
            sprintf(attrStr," %s=\"%g\"",S(nameS),value);
         }
         else if (typeS.equals("double"))
         {
            double value;
	    xdr_double(xdrs,&value);
            sprintf(attrStr," %s=\"%g\"",S(nameS),value);
         }
         else if (typeS.equals("bool"))
         {
            bool_t value;
	    xdr_bool(xdrs,&value);
            sprintf(attrStr," %s=\"%d\"",S(nameS),value);
         }
         else if (typeS.equals("Particle_t"))
         {
            Particle_t value;
            xdr_int(xdrs,(int*)&value);
            sprintf(attrStr," %s=\"%s\"",S(nameS),ParticleType(value));
         }
         else if (typeS.equals("string") || typeS.equals("anyURI"))
         {
            char* value = new char [999];
            xdr_string(xdrs,&value,999);
            sprintf(attrStr," %s=\"%s\"",S(nameS),value);
	    delete [] value;
         }
         else if (nameS.equals("minOccurs") || nameS.equals("maxOccurs"))
         {
            attrStr[0] = 0;
         }
         else
         {
            sprintf(attrStr," %s=\"%s\"",S(nameS),S(typeS));
         }
         writeXML(attrStr);
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
         char endTag[500];
         sprintf(endTag,"</%s>\n",S(tagS));
         writeXML(endTag);
      }
   }
}

int main(int argC, char* argV[])
{
   char* xFilename = 0;

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
      else if (strcmp(argV[argInd],"-o") == 0)
      {
         xFilename = argV[++argInd];
      }
      else
      {
         usage();
         return 1;
      }
   }

   char* hddmFile;
   FILE* ifd;
   if (argInd == argC)
   {
      ifd = stdin;
      hddmFile = new char[1];
      *hddmFile = 0;
   }
   else if (argInd == argC - 1)
   {
      hddmFile = argV[argInd];
      ifd = fopen(hddmFile,"r");
   }
   else
   {
      usage();
      return 1;
   }
   if (!ifd)
   {
      cerr << "hddm-xml: Error opening input stream " << hddmFile << endl;
      exit(1);
   }
   int pid;
   sscanf(getenv("$$"),"%d",&pid);
   char tmpFile[30];
   sprintf(tmpFile,"tmp%d",pid);
   ofstream ofs(tmpFile);
   if (! ofs.is_open())
   {
      cerr << "hddm-xml: Error opening temp file " << tmpFile << endl;
      exit(2);
   }

   ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
   char xmlHeader[500];
   size_t lineSize = 500;
   char* line = new char [500];
   if (getline(&line,&lineSize,ifd))
   {
      if (strstr(line,"<?xml") != 0)
      {
         cerr << "hddm-xml: Error reading input stream " << hddmFile << endl;
         cerr << "Input file appears to be an xml document!" << endl;
         exit(1);
      }
      else if (strstr(line,"<HDDM") == line)
      {
         strncpy(xmlHeader,line,500);
         ofs << line << endl;
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
   while (getline(&line,&lineSize,ifd))
   {
      ofs << line << endl;
      if (strstr(line,"</HDDM>") != 0)
      {
         break;
      }
   }
   ofs.close();
   delete [] line;

#if defined OLD_STYLE_XERCES_PARSER
   DOMDocument* document = parseInputDocument(tmpFile,false);
#else
   DOMDocument* document = buildDOMDocument(tmpFile,false);
#endif
   if (document == 0)
   {
      cerr << "hddm-xml : Error parsing HDDM document, "
           << "cannot continue" << endl;
      return 1;
   }
   unlink(tmpFile);

   DOMElement* rootEl = document->getDocumentElement();
   XString rootS(rootEl->getTagName());
   if (!rootS.equals("HDDM"))
   {
      cerr << "hddm-xml error: root element of input document is "
           << "\"" << S(rootS) << "\", expected \"HDDM\""
           << endl;
      return 1;
   }

   if (xFilename)
   {
      char fname[500];
      sprintf(fname,"%s.xml",xFilename);
      xout.open(fname);
   }
  
   writeXML("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
   writeXML(xmlHeader);

   XDR* xdrs = new XDR;
   xdrstdio_create(xdrs,ifd,XDR_DECODE);
   int icount;
   while (xdr_int(xdrs,&icount))
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
               constructXML(xdrs,contEl,1);
            }
         }
      }
   }

   writeXML("</HDDM>\n");

   XMLPlatformUtils::Terminate();
   return 0;
}
