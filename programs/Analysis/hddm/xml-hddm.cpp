/*
 *  xml-hddm :	tool that reads in a plain-text xml data document
 *              and writes it out as a hddm (Hall D Data Model)
 *              following the data model provided.
 *
 *  Original version - Richard Jones, October 1 2001.
 *
 *
 *  Programmer's Notes:
 *  -------------------
 * 1. The output from xml-hddm is a valid hddm data stream.
 *
 * 2. Two inputs are required: the input xml data document, and a hddm
 *    data model describing the structure of the data in the document.
 *    Both documents must be well-formed.  In addition, the data document
 *    must conform to the hddm specification document.  Only if both of
 *    these requirements are met is it possible to insure that the data
 *    document can be expressed in a hddm stream.
 *
 * 3. There is an ambiguity in the interpretation of a xml data document
 *    that conforms to a hddm data model with more than one tag directly
 *    descending from the <HDDM> document tag.  Since the data model does
 *    not treat the order of content tags at a given level as significant
 *    in distinguishing one data model from another, it is not possible
 *    just looking at the sequence of level-1 tags in the document to tell
 *    where one "event" (in-memory structure) stops and the next one begins.
 *    Thus it is decided that each level-1 tag constititutes a new event.
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

#include "xml-hddm.hpp"
#include "particleType.h"

ofstream out; 

/* inline functions to pack values into the output stream */

inline void pack_int(int* &bp, int val)
{
   *(bp++) = val;
}

inline void pack_float(int* &bp, float val)
{
   *(float*)(bp++) = val;
}

inline void pack_double(int* &bp, double val)
{
   int* pval = (int*) bp;
   *(bp++) = pval[0];
   *(bp++) = pval[1];
}

inline void writeOut(char* string)
{
   if (out.is_open())
   {
      out.write(string,strlen(string));
   }
   else
   {
      cout << string;
   }
}

/* Generate the xml document model in normal form */

void constructDocument(DOM_Element& el)
{
   static int indent = 0;
   for (int n = 0; n < indent; n++)
   {
      writeOut("  ");
   }
   
   char* tagStr = el.getTagName().transcode();
   writeOut("<"), writeOut(tagStr);
   DOM_NamedNodeMap attrList = el.getAttributes();
   int attrListLength = attrList.getLength();
   for (int a = 0; a < attrListLength; a++)
   {
      DOM_Node node = attrList.item(a);
      char* name = node.getNodeName().transcode();
      char* value = node.getNodeValue().transcode();
      writeOut(" "), writeOut(name),
      writeOut("=\""), writeOut(value), writeOut("\"");
      delete [] name;
      delete [] value;
   }

   DOM_NodeList contList = el.getChildNodes();
   int contListLength = contList.getLength();
   if (contListLength > 0)
   {
      writeOut(">"), writeOut("\n");
      indent++;
      for (int c = 0; c < contListLength; c++)
      {
         DOM_Node node = contList.item(c);
         if (node.getNodeType() == ELEMENT_NODE)
         {
            DOM_Element contEl = (DOM_Element&) node;
            constructDocument(contEl);
         }
      }
      indent--;
      for (int n = 0; n < indent; n++)
      {
         writeOut("  ");
      }
      writeOut("</"), writeOut(tagStr), writeOut(">\n");
   }
   else
   {
      writeOut(" />\n");
   }
   delete [] tagStr;
}

/* Generate the output binary stream according the HDDM model */

void outputStream(DOM_Element thisEl, DOM_Element modelEl, int* &bp)
{
   char* modelStr = modelEl.getTagName().transcode();
   char* thisStr = thisEl.getTagName().transcode();

   DOM_NamedNodeMap modelAttrList = modelEl.getAttributes();
   int modelAttrListLength = modelAttrList.getLength();
   DOM_NamedNodeMap thisAttrList = thisEl.getAttributes();
   int thisAttrListLength = thisAttrList.getLength();
   DOMString rep = modelEl.getAttribute("repeat");
   int expectAttrCount = modelAttrListLength - ((rep == 0)? 0 : 1);
   if (thisAttrListLength != expectAttrCount)
   {
      cerr << "xml-hddm: Inconsistency in input xml document" << endl
           << "tag " << thisStr << " in input document with "
           << thisAttrListLength << " attributes " << endl
           << "matched to tag" << modelStr << " in hddm model "
           << "with " << expectAttrCount << " attributes." << endl;
      exit(1);
   }
   for (int a = 0; a < modelAttrListLength; a++)
   {
      DOMString attrS = modelAttrList.item(a).getNodeName();
      DOMString typeS = modelAttrList.item(a).getNodeValue();
      DOMString valueS = thisEl.getAttribute(attrS);
      if (attrS.equals("repeat"))
      {
         continue;
      }
      else if (valueS == 0)
      {
         char* attrStr = attrS.transcode();
         cerr << "xml-hddm: Inconsistency in input xml document" << endl
              << "tag " << thisStr << " in input document is missing "
              << "attribute " << attrStr << endl;
         exit(1);
      }
      char* value = valueS.transcode();
      if (typeS.equals("int"))
      {
         int val = atoi(value);
         pack_int(bp,val);
      }
      else if (typeS.equals("float"))
      {
         float val = atof(value);
         pack_float(bp,val);
      }
      else if (typeS.equals("double"))
      {
         double val = atof(value);
         pack_double(bp,val);
      }
      else if (typeS.equals("bool"))
      {
         int val = atoi(value);
         pack_int(bp,val);
      }
      else if (typeS.equals("Particle_t"))
      {
         int val;
         for (val = 1; val < 99; val++)
         {
            if (strcmp(value,ParticleType((Particle_t)val)) == 0)
            {
               break;
            }
         }
         pack_int(bp,val);
      }
      delete [] value;
   }

   DOM_NodeList thisList = thisEl.getChildNodes();
   int thisListLength = thisList.getLength();
   DOM_NodeList modelList = modelEl.getChildNodes();
   int modelListLength = modelList.getLength();
   for (int m = 0; m < modelListLength; m++)
   {
      DOM_Node mode = modelList.item(m);
      short type = mode.getNodeType();
      if (type == ELEMENT_NODE)
      {
         int* bpStart = bp++;
         DOM_Element model = (DOM_Element&) mode;
         DOMString modelS = model.getTagName();
         DOMString re = model.getAttribute("repeat");
         int repCounter;
         int* repCount;
         if (re != 0)
         {
            repCount = bp++;
         }
         else
         {
             repCount = &repCounter;
         }
         *repCount = 0;
         for (int i = 0; i < thisListLength; i++)
         {
            DOM_Node instance = thisList.item(i);
            short type = instance.getNodeType();
            if (type == ELEMENT_NODE)
            {
               DOM_Element instancel = (DOM_Element&) instance;
               DOMString nameS = instancel.getTagName();
               if (nameS.equals(modelS))
               {
                  outputStream(instancel,model,bp);
                  if ((*repCount)++ && (re == 0))
                  {
                     char* name = nameS.transcode();
                     cerr << "xml-hddm: Inconsistency in input xml document"
                          << endl
                          << "tag " << thisStr << " in input document contains"
                          << " multiple instances of tag " << name << endl
                          << "but it does not have a repeat=\"*\" attribute "
                          << "in the data model." << endl;
                     exit(1);
                  }
               }
            }
         }
         int size = bp - bpStart;
         *bpStart = size -1;
      }
   }
   delete [] modelStr;
   delete [] thisStr;
}

void usage()
{
   cerr << "\nUsage:\n"
        << "    xml-hddm [-o <filename>] -m <model> [input file]\n\n"
        << "Options:\n"
        <<  "    -m <model file>  read hddm model from <model>.xml\n"
        <<  "    -o <filename>    write to <filename>.hddm"
        << endl;
}


int main(int argC, char* argV[])
{
   ifstream* ifs;
   char* modelFilename = 0;
   char* outFilename = 0;

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
      else if (strcmp(argV[argInd],"-m") == 0)
      {
         modelFilename = argV[++argInd];
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

   if (modelFilename == 0)
   {
      usage();
      exit(1);
   }

   if (outFilename)
   {
      char fname[500];
      sprintf(fname,"%s.xml",outFilename);
      out.open(fname);
   }

   DOMParser modelParser;
   modelParser.setValidationScheme(DOMParser::Val_Never);
   modelParser.setCreateEntityReferenceNodes(false);
   modelParser.setDoNamespaces(false);

   MyOwnErrorHandler errorHandler;
   modelParser.setErrorHandler(&errorHandler);

   try
   {
      char fname[500];
      sprintf(fname,"%s.xml",modelFilename);
      modelParser.parse(fname);
      ifs = new ifstream(fname);
   }

   catch (const XMLException& toCatch)
   {
      cerr << "\nxml-hddm: Error during parsing: '" << modelFilename
           << ".xml'\n"
           << "Exception message is: "
           << StrX(toCatch.getMessage()) << "\n" << endl;
      return -1;
   }
   catch (const DOM_DOMException& toCatch)
   {
      cerr << "\nxml-hddm: Error during parsing: '" << modelFilename
           << ".xml'\n"
           << "Exception message is: "
           << toCatch.msg.transcode() << "\n" << endl;
      XMLPlatformUtils::Terminate();
      return 4;
   }
   catch (...)
   {
      cerr << "\nxml-hddm: Unexpected exception during parsing: '"
           << modelFilename << ".xml'\n";
      XMLPlatformUtils::Terminate();
      return 4;
   }

   DOM_Document modelDoc = modelParser.getDocument();
   DOM_Element rootEl = modelDoc.getDocumentElement();
   char* rootStr = rootEl.getTagName().transcode();
   if (strcmp(rootStr,"HDDM") != 0)
   {
      cerr << "xml-hddm error: root element of input document is "
           << "\"" << rootStr << "\", expected \"HDDM\""
           << endl;
      exit(1);
   }
   delete [] rootStr;

   writeOut("\n");
   constructDocument(rootEl);

   char* hddmFile;
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
      cerr << "xml-hddm: Error opening input stream " << hddmFile << endl;
      exit(1);
   }

   char xmlHeader[500];
   char docHeader[500];
   char* line = new char[500000];
   if (ifs->getline(line,500))
   {
      if (strstr(line,"<?xml") == 0)
      {
         cerr << "xml-hddm: Error reading input stream " << hddmFile << endl;
         cerr << "Input file does not appear to be an xml document!" << endl;
         exit(1);
      }
   }
   else
   {
      cerr << "xml-hddm: Error reading from input stream " << hddmFile << endl;
      exit(1);
   }
   strncpy(xmlHeader,line,500);
   if (ifs->getline(line,500) && (strstr(line,"<HDDM") != line))
   {
      cerr << "xml-hddm: Input document tag is not HDDM!"
           << endl;
      exit(1);
   }
   strncpy(docHeader,line,500);

   int* buffer = new int[1000000];
   char* tmpFile = tmpnam(0);
   ofstream ofs;

   while (ifs->getline(line,500000))
   {
      if (strlen(line) > 499990)
      {
         cerr << "xml-hddm: line too long in input document" << endl;
         exit(1);
      }

      char* text = line;
      DOMParser parser;
      parser.setValidationScheme(DOMParser::Val_Never);
      parser.setDoNamespaces(false);
      MyOwnErrorHandler errorHandler;
      parser.setErrorHandler(&errorHandler);

      ofs.open(tmpFile);
      if (! ofs.is_open())
      {
         cerr << "xml-hddm: Error opening temp file " << tmpFile << endl;
         exit(2);
      }
      ofs << xmlHeader << endl;
      ofs << docHeader << endl;

      while (text)
      {
         if ((text = strstr(line,"<")) == 0)
         {
            break;
         }
         else if (strstr(text,"</") == text)
         {
            break;
         }
         char* end;
         while ((end = strstr(text,">")) == 0)
         {
            int len = strlen(line);
            if (len > 400000)
            {
               cerr << "xml-hddm: tag too long in input document" << endl;
               exit(1);
            }
            else
            {
               ifs->getline(&line[len],100000);
            }
         }
         if (strstr(end-1,"/>") == end-1)
         {
            *end = 0;
            ofs << text << ">" << endl;
            text = end + 1;
         }
         else
         {
            char endTag[500];
            sprintf(endTag,"</%s>",strtok(text+1," \t"));
            text[strlen(text)] = ' ';
            while ((end = strstr(text,endTag)) == 0)
            {
               ofs << text << endl;
               ifs->getline(line,500000);
               text = line;
            }
            *end = 0;
            ofs << text << endTag << endl;
            text = end + strlen(endTag);
         }
         ofs << "</HDDM>" << endl;
         ofs.close();

         try
         {
            parser.parse(tmpFile);
            unlink(tmpFile);
         }
         catch (const XMLException& toCatch)
         {
            cerr << "\nxml-hddm: Error during parsing: '" << tmpFile << "'\n"
                 << "Exception message is:  "
                 << StrX(toCatch.getMessage()) << "\n" << endl;
            return -1;
         }
         catch (const DOM_DOMException& toCatch)
         {
            cerr << "\nxml-hddm: Error during parsing: '" << tmpFile << "'\n"
                 << "Exception message is:  "
                 << toCatch.msg.transcode() << "\n" << endl;
            XMLPlatformUtils::Terminate();
            return 4;
         }
         catch (...)
         {
            cerr << "\nxml-hddm: Unexpected exception during parsing: '"
                 << tmpFile << "'\n";
            XMLPlatformUtils::Terminate();
            return 4;
         }

         if (errorHandler.getSawErrors())
         {
            cerr << "\nErrors occured, no output available\n" << endl;
         }

         DOM_Document doc = parser.getDocument();
         DOM_Element thisEl = doc.getDocumentElement();

         int* bp = buffer;
         outputStream(thisEl,rootEl,++bp);
         int size = bp - buffer;
         *buffer = size -1;
         if (out.is_open())
         {
            out.write(buffer,size * sizeof(int));
         }
         else
         {
            cout.write(buffer,size * sizeof(int));
         }
      }
   }

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
