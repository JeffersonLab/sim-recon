/*
 *  hddm-xml :	tool that reads in a HDDM document (Hall D Data Model)
 *		and translates it into plain-text xml.
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

#include <util/PlatformUtils.hpp>
#include <sax/SAXException.hpp>
#include <sax/SAXParseException.hpp>
#include <parsers/DOMParser.hpp>
#include <dom/DOM_DOMException.hpp>
#include <dom/DOM_NamedNodeMap.hpp>

#include <assert.h>
#include <fstream.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "hddm-xml.hpp"
#include "particleType.h"

ofstream xout; 

/* inline functions to unpack values from input stream */

inline int unpack_int(int* &bp)
{
   return *(bp++);
}

inline float unpack_float(int* &bp)
{
   return *(float*)(bp++);
}

inline double unpack_double(int* &bp)
{
   double val;
   int* pval = (int*) &val;
   pval[0] = *(bp++);
   pval[1] = *(bp++);
   return val;
}

/* write a string to xml output stream, either stdout or a file */

void writeXML(char* s)
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
   at entry the buffer pointer bp points the the word after the word count */

void constructXML(int* bp, DOM_Element el, int depth)
{
   DOMString rep = el.getAttribute("repeat");
   int repeats;
   if (rep != 0)
   {
      repeats = unpack_int(bp);
   }
   else
   {
      repeats = 1;
   }

   for (int r = 0; r < repeats; r++)
   {
      char* tagStr = el.getTagName().transcode();
      for (int d = 0; d < depth; d++)
      {
         writeXML("  ");
      }
      writeXML("<");
      writeXML(tagStr);
      DOM_NamedNodeMap attrList = el.getAttributes();
      int listLength = attrList.getLength();
      for (int a = 0; a < listLength; a++)
      {
         DOMString nameS = attrList.item(a).getNodeName();
         char* name = nameS.transcode();
         DOMString typeS = attrList.item(a).getNodeValue();
         char attrStr[500];
         if (typeS.equals("int"))
         {
            int value = unpack_int(bp);
            sprintf(attrStr," %s=\"%d\"",name,value);
         }
         else if (typeS.equals("float"))
         {
            float value = unpack_float(bp);
            sprintf(attrStr," %s=\"%g\"",name,value);
         }
         else if (typeS.equals("double"))
         {
            double value = unpack_double(bp);
            sprintf(attrStr," %s=\"%g\"",name,value);
         }
         else if (typeS.equals("bool"))
         {
            int value = unpack_int(bp);
            sprintf(attrStr," %s=\"%d\"",name,value);
         }
         else if (typeS.equals("Particle_t"))
         {
            Particle_t value = (Particle_t)unpack_int(bp);
            sprintf(attrStr," %s=\"%s\"",name,ParticleType(value));
         }
         else if (nameS.equals("repeat"))
         {
            attrStr[0] = 0;
         }
         else
         {
            char* value = typeS.transcode();
            sprintf(attrStr," %s=\"%s\"",name,value);
            delete [] value;
         }
         writeXML(attrStr);
         delete [] name;
      }

      DOM_NodeList contList = el.getChildNodes();
      int contLength = contList.getLength();
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
         DOM_Node cont = contList.item(c);
         short type = cont.getNodeType();
         if (type == ELEMENT_NODE)
         {
            DOM_Element contEl = (DOM_Element&) cont;
            int size = unpack_int(bp);
            if (size > 0)
            {
               constructXML(bp,contEl,depth +1);
               bp += size;
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
         sprintf(endTag,"</%s>\n",tagStr);
         writeXML(endTag);
      }
      delete [] tagStr;
   }
}

void usage()
{
   cerr << "\nUsage:\n"
        << "    hddm-xml [-o <filename>] [HDDM file]\n\n"
        << "Options:\n"
        <<  "    -o <filename>	write to <filename>.xml"
        << endl;
}


int main(int argC, char* argV[])
{
   char* xFilename = 0;

   try
   {
      XMLPlatformUtils::Initialize();
   }
   catch (const XMLException& toCatch)
   {
      cerr << "hddm-xml: Error during initialization! :\n"
           << StrX(toCatch.getMessage()) << endl;
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
   ifstream* ifs;
   if (argInd == argC)
   {
      ifs = new ifstream(0);
      hddmFile = new char[1];
      *hddmFile = 0;
   }
   else if (argInd == argC - 1)
   {
      hddmFile = argV[argInd];
      ifs = new ifstream(hddmFile);
   }
   else
   {
      usage();
      return 1;
   }
   if (! *ifs)
   {
      cerr << "hddm-xml: Error opening input stream " << hddmFile << endl;
      exit(1);
   }
   char* tmpFile = tmpnam(0);
   ofstream ofs(tmpFile);
   if (! ofs)
   {
      cerr << "hddm-xml: Error opening temp file " << tmpFile << endl;
      exit(2);
   }

   ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
   char line[500];
   char xmlHeader[500];
   if (ifs->getline(line,500))
   {
      if (strstr(line,"<?xml") != 0)
      {
         cerr << "hddm-xml: Error reading input stream " << hddmFile << endl;
         cerr << "Input file appears to be an xml document!" << endl;
         exit(1);
      }
      ofs << line << endl;
   }
   else
   {
      cerr << "hddm-xml: Error reading from input stream " << hddmFile << endl;
      exit(1);
   }
   if (ifs->getline(line,500) && (strstr(line,"<HDDM") == line))
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
   while (ifs->getline(line,500))
   {
      ofs << line << endl;
      if (strstr(line,"</HDDM>") != 0)
      {
         break;
      }
   }
   ofs.close();

   DOMParser parser;
   parser.setValidationScheme(DOMParser::Val_Never);
   parser.setCreateEntityReferenceNodes(false);
   parser.setDoNamespaces(false);

   MyOwnErrorHandler errorHandler;
   parser.setErrorHandler(&errorHandler);

   try
   {
      parser.parse(tmpFile);
      unlink(tmpFile);
   }
   catch (const XMLException& toCatch)
   {
      cerr << "\nhddm-xml: Error during parsing: '" << tmpFile << "'\n"
           << "Exception message is:  "
           << StrX(toCatch.getMessage()) << "\n" << endl;
      return -1;
   }
   catch (const DOM_DOMException& toCatch)
   {
      cerr << "\nhddm-xml: Error during parsing: '" << tmpFile << "'\n"
           << "Exception message is:  "
           << toCatch.msg.transcode() << "\n" << endl;
      XMLPlatformUtils::Terminate();
      return 4;
   }
   catch (...)
   {
      cerr << "\nhddm-xml: Unexpected exception during parsing: '"
           << tmpFile << "'\n";
      XMLPlatformUtils::Terminate();
      return 4;
   }

   if (errorHandler.getSawErrors())
   {
      cerr << "\nErrors occured, no output available\n" << endl;
   }

   DOM_Document doc = parser.getDocument();
   DOM_Element rootEl = doc.getDocumentElement();
   char* rootStr = rootEl.getTagName().transcode();
   if (strcmp(rootStr,"HDDM") != 0)
   {
      cerr << "hddm-xml error: root element of input document is "
           << "\"" << rootStr << "\", expected \"HDDM\""
           << endl;
      delete [] rootStr;
      return 1;
   }
   delete [] rootStr;

   if (xFilename)
   {
      char fname[500];
      sprintf(fname,"%s.xml",xFilename);
      xout.open(fname);
   }
  
   writeXML("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
   writeXML(xmlHeader);
   writeXML("\n");

   int buff[100000];
   int icount;
   while (ifs->read((char*)&icount,sizeof(int)))
   {
      ifs->read((char*)buff,icount*sizeof(int));
      DOM_NodeList contList = rootEl.getChildNodes();
      int contLength = contList.getLength();
      int* bp = buff;
      for (int c = 0; c < contLength; c++)
      {
         DOM_Node cont = contList.item(c);
         short type = cont.getNodeType();
         if (type == ELEMENT_NODE)
         {
            DOM_Element contEl = (DOM_Element&) cont;
            int size = unpack_int(bp);
            if (size > 0)
            {
               constructXML(bp,contEl,1);
            }
         }
      }
   }

   writeXML("</HDDM>\n");

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
   cerr << "\nhddm-xml: Error at file " << StrX(e.getSystemId())
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << StrX(e.getMessage()) << endl;
}

void MyOwnErrorHandler::fatalError(const SAXParseException& e)
{
   fSawErrors = true;
   cerr << "\nhddm-xml: Fatal Error at file " << StrX(e.getSystemId())
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << StrX(e.getMessage()) << endl;
}

void MyOwnErrorHandler::warning(const SAXParseException& e)
{
   cerr << "\nhddm-xml: Warning at file " << StrX(e.getSystemId())
        << ", line " << e.getLineNumber()
        << ", char " << e.getColumnNumber()
        << "\n  Message: " << StrX(e.getMessage()) << endl;
}

void MyOwnErrorHandler::resetErrors()
{
}
