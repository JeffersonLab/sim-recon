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
 * 2. The hddm file contains a xml header that functions as a prototype
 *    of the xml output.
 *
 * 3. This tool can read any hddm file.  No recompilation is required.
 *
 * 3. The code has been tested with the xerces-c DOM implementation from
 *    Apache, and is intended to be used with the xerces-c library.
 *
 * 4. Output is sent by default to stdout and can be changed with the
 *    -o option.
 */

#include <util/PlatformUtils.hpp>
#include <sax/SAXException.hpp>
#include <sax/SAXParseException.hpp>
#include <parsers/DOMParser.hpp>
#include <dom/DOM_DOMException.hpp>
#include <dom/DOM_NamedNodeMap.hpp>

DOM_Element* tagList[100000];
char* tagListName[100000];
int tagListStackPtr[100000];
int tagListLength = 0;
int tagStack[100000];
int tagStackLength = 0;

#include <assert.h>
#include <fstream.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "hddm-xml.hpp"
#include "particleType.h"

ofstream xout; 


int constructGroup(DOM_Element& el)
{
   char* tagStr = el.getTagName().transcode();
   int t;
   for (t=0; t < tagListLength; t++)
   {
      if (strcmp(tagStr,tagListName[t]) == 0) break;
   }

   if (t < tagListLength)
   {
      DOM_NamedNodeMap oldAttr = tagList[t]->getAttributes();
      DOM_NamedNodeMap newAttr = el.getAttributes();
      int listLength = oldAttr.getLength();
      if (newAttr.getLength() != listLength)
      {
         cerr << "hddm-c error: inconsistent usage of tag "
              << "\"" << tagStr << "\" in xml document." << endl;
         exit(1);
      }
      for (int n = 0; n < listLength; n++)
      {
         char* name = oldAttr.item(n).getNodeName().transcode();
         char* value = oldAttr.item(n).getNodeValue().transcode();
         short type = oldAttr.item(n).getNodeType();
         char* was = tagList[t]->getAttribute(name).transcode();
         char* is = el.getAttribute(name).transcode();
         if (strcmp(was,is) != 0)
         {
            cerr << "hddm-c error: inconsistent usage of attribute "
                 << "\"" << name << "\" in tag "
                 << "\"" << tagStr << "\" in xml document." << endl;
            exit(1);
         }
         delete [] name;
         delete [] value;
         delete [] was;
         delete [] is;
      }
      DOM_NodeList oldList = tagList[t]->getChildNodes();
      DOM_NodeList newList = el.getChildNodes();
      listLength = oldList.getLength();
      if (newList.getLength() != listLength)
      {
         cerr << "hddm-c error: inconsistent usage of tag "
              << "\"" << tagStr << "\" in xml document." << endl;
         exit(1);
      }
      for (int n = 0; n < listLength; n++)
      {
         DOM_Node cont = oldList.item(n);
         char* name = cont.getNodeName().transcode();
         short type = cont.getNodeType();
         if (type == ELEMENT_NODE)
         {
            DOM_NodeList contList = el.getElementsByTagName(name);
            if (contList.getLength() != 1)
            {
                cerr << "hddm-c error: inconsistent usage of tag "
                     << "\"" << tagStr << "\" in xml document." << endl;
                exit(1);
            }
         }
         delete [] name;
      }
      delete [] tagStr;
      return t;
   }

   tagList[t] = new DOM_Element(el);
   tagListName[t] = new char [strlen(tagStr)+1];
   strcpy(tagListName[t],tagStr);
   tagListLength++;

   DOM_NodeList contList = el.getChildNodes();
   int contLength = contList.getLength();
   for (int c = 0; c < contLength; c++)
   {
      DOM_Node cont = contList.item(c);
      short type = cont.getNodeType();
      if (type == ELEMENT_NODE)
      {
         DOM_Element contEl = (DOM_Element&) cont;
         constructGroup(contEl);
      }
   }

   tagListStackPtr[t] = tagStackLength;
   tagStack[tagStackLength++] = t;		// this tag
   int& repFlag = tagStack[tagStackLength++];	// repeat flag
   int& wcount = tagStack[tagStackLength++];	// words of data
   int& pcount = tagStack[tagStackLength++];	// words of pointers

   wcount = 0;
   DOM_NamedNodeMap varList = el.getAttributes();
   int varCount = varList.getLength();
   for (int v = 0; v < varCount; v++)
   {
      DOM_Node var = varList.item(v);
      char* typeStr = var.getNodeValue().transcode();
      if (strcmp(typeStr,"int") == 0)
      {
         wcount += 1;
      }
      else if (strcmp(typeStr,"float") == 0)
      {
         wcount += 1;
      }
      else if (strcmp(typeStr,"double") == 0)
      {
         wcount += 2;
      }
      else if (strcmp(typeStr,"bool") == 0)
      {
         wcount += 1;
      }
      else if (strcmp(typeStr,"Particle_t") == 0)
      {
         wcount += 1;
      }
      else
      {
         /* ignore attributes with unrecognized values */
      }
      delete [] typeStr;
   }

   pcount = 0;
   for (int c = 0; c < contLength; c++)
   {
      DOM_Node cont = contList.item(c);
      short type = cont.getNodeType();
      if (type == ELEMENT_NODE)
      {
         DOM_Element contEl = (DOM_Element&) cont;
         tagStack[tagStackLength++] = constructGroup(contEl);
         ++pcount;
      }
   }

   repFlag = 0;
   DOMString rep = el.getAttribute("repeat");
   if (rep != 0)
   {
      repFlag = 1;
   }

   delete [] tagStr;
   return t;
}


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


void constructXML(int* sp, int t)
{
   int s = tagListStackPtr[t];
   assert (tagStack[s++] == t);
   int repFlag = tagStack[s++];
   int wcount = tagStack[s++];
   int pcount = tagStack[s++];

   int repeats;
   if (repFlag)
   {
      repeats = *(sp++);
   }
   else
   {
      repeats = 1;
   }

   static int depth = 0;

   for (int i = 0; i < repeats; i++)
   {
      if (depth == 0 || t > 0)
      {
         for (int d = 0; d < depth; d++)
         {
            writeXML("  ");
         }
         writeXML("<");
         writeXML(tagListName[t]);
         DOM_NamedNodeMap attrList = tagList[t]->getAttributes();
         int listLength = attrList.getLength();
         for (int n = 0; n < listLength; n++)
         {
            char* name = attrList.item(n).getNodeName().transcode();
            char* typeStr = attrList.item(n).getNodeValue().transcode();
            char attrStr[500];
            if (strcmp(typeStr,"int") == 0)
            {
               int value = *(sp++);
               sprintf(attrStr," %s=\"%d\"",name,value);
            }
            else if (strcmp(typeStr,"float") == 0)
            {
               float value = *(float*)(sp++);
               sprintf(attrStr," %s=\"%f\"",name,value);
            }
            else if (strcmp(typeStr,"double") == 0)
            {
               double value = *(double*)(sp++);
               ++sp;
               sprintf(attrStr," %s=\"%f\"",name,value);
            }
            else if (strcmp(typeStr,"bool") == 0)
            {
               int value = *(sp++);
               sprintf(attrStr," %s=\"%d\"",name,value);
            }
            else if (strcmp(typeStr,"Particle_t") == 0)
            {
               Particle_t value = *(Particle_t*)(sp++);
               sprintf(attrStr," %s=\"%s\"",name,ParticleType(value));
            }
            else if (strcmp(name,"repeat") == 0)
            {
               attrStr[0] = 0;
            }
            else
            {
               sprintf(attrStr," %s=\"%s\"",name,typeStr);
            }
            writeXML(attrStr);
            delete [] name;
            delete [] typeStr;
         }
         if (pcount > 0)
         {
            writeXML(">\n");
         }
         else
         {
            writeXML(" />\n");
         }
         ++depth;
      }

      if (pcount > 0)
      {
         int ss = s;
         for (int p = 0; p < pcount; p++)
         {
            int tt = tagStack[ss++];      
            int wc = *(sp++);
            if (wc > 0)
            {
               constructXML(sp,tt);
               sp += wc;
            }
         }
         if (depth > 1)
         {
            --depth;
            for (int d = 0; d < depth; d++)
            {
               writeXML("  ");
            }
            char endTag[500];
            sprintf(endTag,"</%s>\n",tagListName[t]);
            writeXML(endTag);
         }
      }
      else
      {
         --depth;
      }
   }
}

void usage()
{
   cerr << "\nUsage:\n"
        << "    hddm-xml [-o <filename>] {HDDM file}\n\n"
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
         cerr << "Unknown option \'" << argV[argInd]
              << "\', ignoring it\n" << endl;
      }
   }

   if (argInd != argC - 1)
   {
      usage();
      return 1;
   }
   char* hddmFile = argV[argInd];
   ifstream ifs(hddmFile);
   if (! ifs)
   {
      cerr << "hddm-xml: Error opening input file " << hddmFile << endl;
      exit(1);
   }
   char* tmpFile = tmpnam(0);
   ofstream ofs(tmpFile);
   if (! ofs)
   {
      cerr << "hddm-xml: Error opening temp file " << tmpFile << endl;
      exit(2);
   }

   ofs << "<?xml version=\"1.0\"?>";
   char line[500];
   if (ifs.getline(line,500))
   {
      ofs << line << endl;
   }
   else
   {
      cerr << "hddm-xml: Error reading from input file " << hddmFile << endl;
      exit(1);
   }
   if (ifs.getline(line,500) && (strstr(line,"<HDDM") != 0))
   {
      ofs << line << endl;
   }
   else
   {
      cerr << "hddm-xml: Input file does not contain valid hddm header" << endl;
      exit(1);
   }
   while (ifs.getline(line,500))
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
           << "Exception message is:  \n"
           << StrX(toCatch.getMessage()) << "\n" << endl;
      return -1;
   }
   catch (const DOM_DOMException& toCatch)
   {
      cerr << "\nhddm-xml: Error during parsing: '" << tmpFile << "'\n"
           << "Exception message is:  \n"
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

   constructGroup(rootEl);

   if (xFilename)
   {
      char fname[500];
      sprintf(fname,"%s.xml",xFilename);
      xout.open(fname);
   }
  
   writeXML("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");

   int buff[100000];
   int icount;
   while (ifs.read((char*)&icount,sizeof(int)))
   {
      ifs.read((char*)buff,icount*sizeof(int));
      constructXML(buff,0);
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
