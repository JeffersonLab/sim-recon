/*
 *  hddm-xml :	tool that reads in a HDDM document (Hall D Data Model)
 *		and translates it into plain-text xml.
 *
 *  Version 1.3 - Richard Jones, July 2014.
 *  - Added support for input hddm streams with additional features
 *    provided through the c++ API, including on-the-fly compression with
 *    zlib and bzlib2, and per-record crc32 integrity checks.
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

// #define VERBOSE_HDDM_LOGGING 1
#define BAD_CRC_IS_ONLY_WARNING 1

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>

#include "XParsers.hpp"
#include "XString.hpp"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <rpc/rpc.h>
#include <unistd.h>
#include <xstream/z.h>
#include <xstream/bz.h>
#include <xstream/xdr.h>
#include <xstream/digest.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include "particleType.h"


#define X(str) XString(str).unicode_str()
#define S(str) str.c_str()

using namespace xercesc;

int explicit_repeat_count = 1;

class XMLmaker
{
 public:
   std::ofstream xout; 

   XMLmaker() {};
   ~XMLmaker() {};

   void writeXML(const XString& s);
   void constructXML(xstream::xdr::istream *ifx, DOMElement* el,
                     int size, int depth);
};

class istreambuffer : public std::streambuf {
 public:
   istreambuffer(char* buffer, std::streamsize bufferLength) {
      setg(buffer, buffer, buffer + bufferLength);
   }

   std::streampos tellg() {
      return gptr() - eback();
   }

   void seekg(std::streampos pos) {
      reset();
      gbump(pos);
   }

   int size() {
      return egptr() - gptr();
   }

   void reset() {
      char *gbegin = eback();
      char *gend = egptr();
      setg(gbegin, gbegin, gend);
   }

   char *getbuf() {
      return eback();
   }
};

class ostreambuffer : public std::streambuf {
 public:
   ostreambuffer(char* buffer, std::streamsize bufferLength) {
      setp(buffer, buffer + bufferLength);
   }

   std::streampos tellp() {
      return pptr() - pbase();
   }

   void seekp(std::streampos pos) {
      reset();
      pbump(pos);
   }

   int size() {
      return pptr() - pbase();
   }

   void reset() {
      char *pbegin = pbase();
      char *pend = epptr();
      setp(pbegin, pend);
   }

   char *getbuf() {
      return pbase();
   }
};

void usage()
{
   std::cerr
        << "\nUsage:\n"
        << "    hddm-xml [-n <count>] [-o <filename>] [HDDM file]\n\n"
        << "Options:\n"
        <<  "    -o <filename>  write to <filename>.xml"
        <<  "    -n <count>	    limit output to <count> records"
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

   int reqcount=-1;
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
      else if (strcmp(argV[argInd],"-n") == 0)
      {
         if (!sscanf(argV[++argInd],"%d",&reqcount))
         {
            usage();
            return 1;
         }
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
           << "hddm-xml: Error opening temp file " << tmpFileStr.str() << std::endl;
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

   int event_buffer_size;
   char *event_buffer = new char[event_buffer_size = 1000000];
   istreambuffer *isbuf = new istreambuffer(event_buffer,event_buffer_size);
   xstream::xdr::istream *ifx = new xstream::xdr::istream(isbuf);
   int integrity_check_mode = 0;
   int compression_mode = 0;
   while (reqcount && ifs->good())
   {
      DOMNodeList* contList = rootEl->getChildNodes();
      int contLength = contList->getLength();
      int tsize;
      ifs->read(event_buffer,4);
      if (ifs->eof()) {
         break;
      }
      isbuf->reset();
      *ifx >> tsize;
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
      else if (tsize == 1) {
         int size, format, flags;
         ifs->read(event_buffer+4,4);
         *ifx >> size;
         ifs->read(event_buffer+8,size);
         *ifx >> format >> flags;
         int compression_flags = flags & 0xf0;
         int integrity_flags = flags & 0x0f;
         std::streambuf *fin_sb = 0;
         xstream::z::istreambuf *zin_sb = 0;
         xstream::bz::istreambuf *bzin_sb = 0;
         int *leftovers = new int[100];
         int sizeof_leftovers = sizeof(int[100]);
         leftovers[0] = 0;
         if (compression_flags == compression_mode) {
            fin_sb = ifs->rdbuf();
         }
         else if (size == 8 && format == 0 && compression_flags == 0x10) {
            if (compression_mode == 0x20) {
               bzin_sb = (xstream::bz::istreambuf*)ifs->rdbuf();
            }
            compression_mode = compression_flags;
            zin_sb = new xstream::z::istreambuf(ifs->rdbuf(),
                                                leftovers, sizeof_leftovers);
            ifs->rdbuf(zin_sb);
            if (bzin_sb != 0)
               delete bzin_sb;
         }
         else if (size == 8 && format == 0 && compression_flags == 0x20) {
            if (compression_mode == 0x10) {
               zin_sb = (xstream::z::istreambuf*)ifs->rdbuf();
            }
            compression_mode = compression_flags;
            bzin_sb = new xstream::bz::istreambuf(ifs->rdbuf(),
                                                  leftovers, sizeof_leftovers);
            ifs->rdbuf(bzin_sb);
            if (zin_sb != 0)
               delete zin_sb;
         }
         else {
            if (compression_mode == 0x20) {
               bzin_sb = (xstream::bz::istreambuf*)ifs->rdbuf();
               fin_sb = bzin_sb->get_streambuf();
            }
            else if (compression_mode == 0x10) {
               zin_sb = (xstream::z::istreambuf*)ifs->rdbuf();
               fin_sb = zin_sb->get_streambuf();
            }
            compression_mode = compression_flags;
            ifs->rdbuf(fin_sb);
            if (zin_sb != 0)
               delete zin_sb;
            if (bzin_sb != 0)
               delete bzin_sb;
         }
         if (size == 8 && format == 0 && integrity_flags == 0x0) {
            integrity_check_mode = 0;
         }
         else if (size == 8 && format == 0 && integrity_flags == 0x1) {
            integrity_check_mode = 1;
         }
         else {
            std::cerr << "hddm-xml error: unrecognized stream modifier"
                         " encountered, this stream is no longer readable."
                      << std::endl;
            break;
         }
         continue;
      }
      else if (tsize+4 > event_buffer_size) {
         delete ifx;
         delete isbuf;
         char *new_buffer = new char[event_buffer_size = tsize+1000];
         isbuf = new istreambuffer(new_buffer,event_buffer_size);
         ifx = new xstream::xdr::istream(isbuf);
         memcpy(new_buffer,event_buffer,4);
         *ifx >> tsize;
         delete[] event_buffer;
         event_buffer = new_buffer;
      }
      ifs->read(event_buffer+4,tsize);
      --reqcount;

      if (integrity_check_mode == 1) {
         char crcbuf[10];
         istreambuffer sbuf(crcbuf,10);
         xstream::xdr::istream xstr(&sbuf);
         unsigned int recorded_crc;
         ifs->read(crcbuf,4);
         xstr >> recorded_crc;
         xstream::digest::crc32 crc;
         std::ostream out(&crc);
         out.write(event_buffer,tsize+4);
         out.flush();
         if (crc.digest() != recorded_crc) {
#if BAD_CRC_IS_ONLY_WARNING
            static int bad_crc_warning_needed = true;
            char errmsg[] =
                 "WARNING: data integrity crc check failed on input.\n"
                 "This may be the result of a bug in the xstream library\n"
                 "if you are analyzing a data file that was generated by\n"
                 "code prior to svn rev 18530. If this concerns you, \n"
                 "regenerate the file using a newer build of the sim-recon\n"
                 "tools and it should go away.\n";
            if (bad_crc_warning_needed) {
               std::cerr << errmsg << std::endl;
               bad_crc_warning_needed = false;
            }
#else
            std::cerr << "hddm-xml error: crc32 check error on input stream"
                         " encountered, this stream is no longer readable."
                      << std::endl;
            break;
#endif
         }
      }

      for (int c = 0; c < contLength; c++)
      {
         DOMNode* cont = contList->item(c);
         short type = cont->getNodeType();
         if (type == DOMNode::ELEMENT_NODE)
         {
            DOMElement* contEl = (DOMElement*) cont;
            int size;
            *ifx >> size;
#ifdef VERBOSE_HDDM_LOGGING
            XString cnameS(contEl->getTagName());
            std::cerr << "hddm-xml : top-level tag " << S(cnameS)
                      << " found with size " << size
                      << std::endl;
#endif
            if (size > 0)
            {
               builder.constructXML(ifx,contEl,size,1);
            }
            else {
               XString repS(contEl->getAttribute(X("minOccurs")));
               int rep = (repS == "")? 1 : atoi(S(repS));
               if (rep != 0) {
                  XString conameS(contEl->getTagName());
                  std::cerr << "hddm-xml warning: top-level tag " << S(conameS)
                            << " found with zero size "
                            << "inside an event with size " << tsize
                            << " continue? [y/n] ";
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

void XMLmaker::constructXML(xstream::xdr::istream *ifx,
                            DOMElement* el, int size, int depth)
{
   XString tagS(el->getTagName());
   XString repS(el->getAttribute(X("maxOccurs")));
   int rep = (repS == "unbounded")? INT_MAX :
             (repS == "")? 1 :
             atoi(S(repS));
   if (explicit_repeat_count && rep > 1)
   {
      *ifx >> rep;
      size -= 4;
   }

   int r;
   for (r = 0; r < rep && size > 0; r++)
   {
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
	    *ifx >> value;
            size -= 4;
            attrStr << " " << nameS << "=\"" << value << "\"";
         }
	 else if (typeS == "long")
         {
            int64_t value;
            *ifx >> value;
            size -= 8;
            attrStr << " " << nameS << "=\"" << value << "\"";
         }
         else if (typeS == "float")
         {
            float value;
            *ifx >> value;
            size -= 4;
            attrStr << " " << nameS << "=\"" << value << "\"";
         }
         else if (typeS == "double")
         {
            double value;
            *ifx >> value;
            size -= 8;
            attrStr << " " << nameS << "=\"" << value << "\"";
         }
         else if (typeS == "boolean")
         {
            bool_t value;
            *ifx >> value;
            size -= 4;
            attrStr << " " << nameS << "=\"" << value << "\"";
         }
         else if (typeS == "Particle_t")
         {
            int32_t value;
            *ifx >> value;
            size -= 4;
            attrStr << " " << nameS << "=\"" << ParticleType((Particle_t)value) << "\"";
         }
         else if (typeS == "string" || typeS == "anyURI")
         {
            std::string value;
            *ifx >> value;
            int strsize = value.size();
            size -= strsize + 4 + ((strsize % 4)? 4-(strsize % 4) : 0);
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
            int csize;
            *ifx >> csize;
            size -= 4;
#ifdef VERBOSE_HDDM_LOGGING
            XString cnameS(contEl->getTagName());
            std::cerr << "hddm-xml : tag " << S(cnameS)
                      << " found with size " << csize
                      << std::endl;
#endif
            if (csize > 0) {
               constructXML(ifx,contEl,csize,depth +1);
               size -= csize;
            }
#ifdef VERBOSE_HDDM_LOGGING
            else {
               XString irepS(contEl->getAttribute(X("minOccurs")));
               int irep = (irepS == "")? 1 : atoi(S(irepS));
               if (irep != 0) {
                  XString conameS(contEl->getTagName());
                  std::cerr << "hddm-xml warning: tag " << S(conameS)
                            << " found with zero size, "
                            << "continue? [y/n] ";
                  std::string ans;
                  std::cin >> ans;
                  if (ans[0] != 'y' && ans[0] != 'Y') {
                     exit(5);
                  }
               }
            }
#endif
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
   if (size != 0) {
      std::cerr << "hddm-xml : size mismatch in tag " << S(tagS)
                << ", remainder is " << size
                << ", cannot continue." << std::endl;
      exit(5);
   }
   else if (explicit_repeat_count && r != rep) {
      std::cerr << "hddm-xml : repeat count mismatch in tag " << S(tagS)
                << ", expected " << rep << " but saw " << r
                << ", cannot continue." << std::endl;
      exit(5);
   }
}
