/*
 * hddm-root :  tool that reads in a HDDM document (Hall D Data Model)
 *              and copies the contents into root trees for easy browsing.
 *
 * author: richard.t.jones at uconn.edu
 * version: january 2, 2017
 *
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
#include <xstream/z.h>
#include <xstream/bz.h>
#include <xstream/xdr.h>
#include <xstream/digest.h>

#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include <TTree.h>
#include <TFile.h>

#include "particleType.h"


#define X(str) XString(str).unicode_str()
#define S(str) str.c_str()

using namespace xercesc;

int explicit_repeat_count = 1;

void usage()
{
   std::cerr
        << "\nUsage:\n"
        << "    hddm-root [-o <filename>] [HDDM file]\n\n"
        << "Options:\n"
        <<  "    -o <filename>    write to <filename>.root"
        << std::endl;
}

struct attribute_t {
   XString type;
   union payload {
      int i;
      float f;
      double d;
      long int l;
      Particle_t p;
      char s[32];
   };
   payload datum;

   attribute_t() {}
   attribute_t(XString varname)
    : type("") { datum.l = 0; }
   attribute_t(XString varname, XString vartype) 
    : type(vartype) { datum.l = 0; }
   attribute_t(const attribute_t &src)
    : type(src.type)
   {
      if (src.type == "string") {
         strncpy(datum.s, src.datum.s, 30);
      }
      else {
         datum.l = src.datum.l;
      }
   }
   attribute_t operator=(const attribute_t &src) {
      attribute_t copy(src);
      return copy;
   }
   ~attribute_t() {}
};

typedef std::map<XString,attribute_t> attribute_table;

class TreeMaker
{
 public:
   TreeMaker(XString filename) {
      fRootFile = new TFile(S(filename), "recreate");
   }
   ~TreeMaker() {
      std::map<XString,TTree*>::iterator iter;
      for (iter = fTrees.begin(); iter != fTrees.end(); ++iter)
         iter->second->Write();
      delete fRootFile;
   }

   void build(const DOMElement* elem);
   void filltrees(xstream::xdr::istream *ifx, DOMElement* el, int size);

 private:
   TreeMaker(const TreeMaker &src) {}
   TreeMaker operator=(const TreeMaker &src) {
      TreeMaker copy(*this);
      return copy;
   }

 protected:
   TFile *fRootFile;
   std::map<XString,TTree*> fTrees;
   attribute_table fColumns;
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

int main(int argC, char* argV[])
{
   XString rootFilename;

   try
   {
      XMLPlatformUtils::Initialize();
   }
   catch (const XMLException* toCatch)
   {
      XString msg(toCatch->getMessage());
      std::cerr
           << "hddm-root: Error during initialization! :\n"
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
         rootFilename = argV[++argInd];
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
           << "hddm-root: Error opening input stream " << hddmFile << std::endl;
      exit(1);
   }
   std::ostringstream tmpFileStr;
   tmpFileStr << "tmp" << getpid();
   std::ofstream ofs(tmpFileStr.str().c_str());
   if (! ofs.is_open())
   {
      std::cerr
           << "hddm-root: Error opening temp file " << tmpFileStr.str() << std::endl;
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
              << "hddm-root: Error reading input stream " << hddmFile
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
              << "hddm-root: Input stream does not contain valid hddm header"
              << std::endl;
         exit(1);
      }
   }
   else
   {
      std::cerr
           << "hddm-root: Error reading from input stream " << hddmFile 
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
           << "hddm-root : Error parsing HDDM document, "
           << "cannot continue" << std::endl;
      return 1;
   }
   unlink(tmpFileStr.str().c_str());

   DOMElement* rootEl = document->getDocumentElement();
   XString rootS(rootEl->getTagName());
   if (rootS != "HDDM")
   {
      std::cerr
           << "hddm-root error: root element of input document is "
           << "\"" << S(rootS) << "\", expected \"HDDM\""
           << std::endl;
      return 1;
   }

   // open root file for output and initialize the trees for writing
   TreeMaker builder(rootFilename);
   builder.build(rootEl);
  
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
      std::cerr << "hddm-root : tag " << S(tnameS)
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
            std::cerr << "hddm-root error: unrecognized stream modifier"
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
            std::cerr << "hddm-root error: crc32 check error on input stream"
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
            std::cerr << "hddm-root : top-level tag " << S(cnameS)
                      << " found with size " << size
                      << std::endl;
#endif
            if (size > 0)
            {
               builder.filltrees(ifx,contEl,size);
            }
            else {
               XString repS(contEl->getAttribute(X("minOccurs")));
               int rep = (repS == "")? 1 : atoi(S(repS));
               if (rep != 0) {
                  XString conameS(contEl->getTagName());
                  std::cerr << "hddm-root warning: top-level tag " << S(conameS)
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

   if (ifs != &std::cin)
   {
      ((std::ifstream*)ifs)->close();
   }

   XMLPlatformUtils::Terminate();
   return 0;
}

void TreeMaker::build(const DOMElement* elem)
{
   // Recursively create TTree objects to hold the contents of the hddm model
   // in the form of a row/column table, like a relational database model.

   XString elemS(elem->getTagName());
   if (elemS != "HDDM") {
      TTree *tree = new TTree(S(elemS), S(XString(elemS + " tree")));
      attribute_table::iterator iter;
      for (iter = fColumns.begin(); iter != fColumns.end(); ++iter) {
         XString nameS = iter->first;
         XString typeS = iter->second.type;
         if (typeS == "int" || typeS == "boolean" || typeS == "Particle_t")
            tree->Branch(S(nameS), (void*)&fColumns[nameS].datum.i,
                                                 S(XString(nameS + "/I")));
         else if (typeS == "long")
            tree->Branch(S(nameS), (void*)&fColumns[nameS].datum.l,
                                                 S(XString(nameS + "/L")));
         else if (typeS == "float")
            tree->Branch(S(nameS), (void*)&fColumns[nameS].datum.f,
                                                 S(XString(nameS + "/F")));
         else if (typeS == "double")
            tree->Branch(S(nameS), (void*)&fColumns[nameS].datum.d,
                                                 S(XString(nameS + "/D")));
         else if (typeS == "string" || typeS == "anyURI")
            tree->Branch(S(nameS), (void*)&fColumns[nameS].datum.s,
                                                 S(XString(nameS + "/C")));
      }
      DOMNamedNodeMap* attrList = elem->getAttributes();
      int attrListLength = attrList->getLength();
      for (int a = 0; a < attrListLength; a++) {
         DOMNode* node = attrList->item(a);
         XString nameS(node->getNodeName());
         XString typeS(node->getNodeValue());
         XString colnameS(elemS + "_" + nameS);
         fColumns[colnameS] = attribute_t(typeS);
         if (typeS == "int" || typeS == "boolean" || typeS == "Particle_t")
            tree->Branch(S(nameS), (void*)&fColumns[colnameS].datum.i,
                                                   S(XString(nameS + "/I")));
         else if (typeS == "long")
            tree->Branch(S(nameS), (void*)&fColumns[colnameS].datum.l,
                                                   S(XString(nameS + "/L")));
         else if (typeS == "float")
            tree->Branch(S(nameS), (void*)&fColumns[colnameS].datum.f,
                                                   S(XString(nameS + "/F")));
         else if (typeS == "double")
            tree->Branch(S(nameS), (void*)&fColumns[colnameS].datum.d,
                                                   S(XString(nameS + "/D")));
         else if (typeS == "string" || typeS == "anyURI")
            tree->Branch(S(nameS), (void*)&fColumns[colnameS].datum.s,
                                                   S(XString(nameS + "/C")));
      }
      fTrees[elemS] = tree;
   }

   DOMNodeList* contList = elem->getChildNodes();
   int contLength = contList->getLength();
   for (int c = 0; c < contLength; c++) {
      DOMNode* cont = contList->item(c);
      short type = cont->getNodeType();
      if (type == DOMNode::ELEMENT_NODE) {
         DOMElement* contEl = (DOMElement*) cont;
         build(contEl);
      }
   }
}

void TreeMaker::filltrees(xstream::xdr::istream *ifx, DOMElement* el, int size)
{
   XString tagS(el->getTagName());
   XString repS(el->getAttribute(X("maxOccurs")));
   int rep = (repS == "unbounded")? INT_MAX :
             (repS == "")? 1 :
             atoi(S(repS));
   if (explicit_repeat_count && rep > 1) {
      *ifx >> rep;
      size -= 4;
   }

   int r;
   for (r = 0; r < rep && size > 0; r++) {
      DOMNamedNodeMap* attrList = el->getAttributes();
      int listLength = attrList->getLength();
      for (int a = 0; a < listLength; a++)
      {
         XString nameS(attrList->item(a)->getNodeName());
         XString typeS(attrList->item(a)->getNodeValue());
         XString colnameS(tagS + "_" + nameS);
         if (typeS == "int" || typeS == "boolean" || typeS == "Particle_t")
         {
            *ifx >> fColumns[colnameS].datum.i;
            size -= 4;
         }
         else if (typeS == "long")
         {
            *ifx >> fColumns[colnameS].datum.l;
            size -= 8;
         }
         else if (typeS == "float")
         {
            *ifx >> fColumns[colnameS].datum.f;
            size -= 4;
         }
         else if (typeS == "double")
         {
            *ifx >> fColumns[colnameS].datum.d;
            size -= 8;
         }
         else if (typeS == "string" || typeS == "anyURI")
         {
            std::string value;
            *ifx >> value;
            strncpy(fColumns[colnameS].datum.s, value.c_str(), 30);
            int strsize = value.size();
            size -= strsize + 4 + ((strsize % 4)? 4-(strsize % 4) : 0);
         }
      }
      fTrees[tagS]->Fill();

      DOMNodeList* contList = el->getChildNodes();
      int contLength = contList->getLength();
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
            std::cerr << "hddm-root : tag " << S(cnameS)
                      << " found with size " << csize
                      << std::endl;
#endif
            if (csize > 0) {
               filltrees(ifx,contEl,csize);
               size -= csize;
            }
#ifdef VERBOSE_HDDM_LOGGING
            else {
               XString irepS(contEl->getAttribute(X("minOccurs")));
               int irep = (irepS == "")? 1 : atoi(S(irepS));
               if (irep != 0) {
                  XString conameS(contEl->getTagName());
                  std::cerr << "hddm-root warning: tag " << S(conameS)
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
   }
   if (size != 0) {
      std::cerr << "hddm-root : size mismatch in tag " << S(tagS)
                << ", remainder is " << size
                << ", cannot continue." << std::endl;
      exit(5);
   }
   else if (explicit_repeat_count && r != rep) {
      std::cerr << "hddm-root : repeat count mismatch in tag " << S(tagS)
                << ", expected " << rep << " but saw " << r
                << ", cannot continue." << std::endl;
      exit(5);
   }
}
