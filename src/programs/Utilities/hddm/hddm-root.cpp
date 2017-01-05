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
int write_xml_output_to_stdout = 0;

void usage()
{
   std::cerr
        << "\nUsage:\n"
        << "    hddm-root [-x] [-n <count>] [-o <filename>] [HDDM file]\n\n"
        << "Options:\n"
        <<  "    -o <filename>  write to output root file <filename>\n"
        <<  "    -n <count>     limit output to <count> rows\n"
        <<  "    -x             write xml to stdout"
        << " (in addition to root file output)\n"
        << std::endl;
}

typedef xstream::xdr::istream ixstream;

class attribute_t {
 protected:
   attribute_t() : fName(""), fType("") {}
   attribute_t(XString name) : fName(name), fType("") {}
   attribute_t(XString name, XString type) : fName(name), fType(type) {}
   virtual ~attribute_t() {}

 public:
   virtual void reset() = 0;
   virtual void *address() = 0;
   virtual std::string toString() = 0;
   virtual int read(ixstream *ifx) = 0;
   virtual XString get_name() { return fName; }
   virtual XString get_type() { return fType; }

 private:
   attribute_t(const attribute_t &src);
   attribute_t &operator=(const attribute_t &src);

 protected:
   XString fName;
   XString fType;
};

class int_attribute_t : public attribute_t {
 public:
   int_attribute_t() : attribute_t("", "int"), value(0) {}
   int_attribute_t(XString name) : attribute_t(name, "int"), value(0) {}
   virtual ~int_attribute_t() {}

   int_attribute_t &operator=(const int_attribute_t &src) {
      fName = src.fName;
      value = src.value;
      return *this;
   }

   virtual void reset() {
      value = 0;
   }
   virtual void set(int val) {
      value = val;
   }
   virtual void *address() {
      return &value;
   }
   virtual std::string toString() {
      std::stringstream str;
      str << value;
      return str.str();
   }
   virtual int read(ixstream *ifx) {
      *ifx >> value;
      return 4;
   }

   int value;
};

class boolean_attribute_t : public attribute_t {
 public:
   boolean_attribute_t() : attribute_t("", "boolean"), value(0) {}
   boolean_attribute_t(XString name) : attribute_t(name, "boolean"), value(0) {}
   virtual ~boolean_attribute_t() {}

   boolean_attribute_t &operator=(const boolean_attribute_t &src) {
      fName = src.fName;
      value = src.value;
      return *this;
   }

   virtual void reset() {
      value = 0;
   }
   virtual void set(int val) {
      value = val;
   }
   virtual void *address() {
      return &value;
   }
   virtual std::string toString() {
      std::stringstream str;
      str << value;
      return str.str();
   }
   virtual int read(ixstream *ifx) {
      *ifx >> value;
      return 4;
   }

   int value;
};

class Particle_attribute_t : public attribute_t {
 public:
   Particle_attribute_t() : attribute_t("", "Particle_t"), value(Unknown) {}
   Particle_attribute_t(XString name) : attribute_t(name, "Particle_t"), 
                                        value(Unknown) {}
   virtual ~Particle_attribute_t() {}

   Particle_attribute_t &operator=(const Particle_attribute_t &src) {
      fName = src.fName;
      value = src.value;
      return *this;
   }

   virtual void reset() {
      value = Unknown;
   }
   virtual void set(Particle_t val) {
      value = val;
   }
   virtual void *address() {
      return &value;
   }
   virtual std::string toString() {
      std::stringstream str;
      str << ParticleType(value);
      return str.str();
   }
   virtual int read(ixstream *ifx) {
      int val;
      *ifx >> val;
      value = (Particle_t)val;
      return 4;
   }

   Particle_t value;
};

class long_attribute_t : public attribute_t {
 public:
   long_attribute_t() : attribute_t("", "long"), value(0) {}
   long_attribute_t(XString name) : attribute_t(name, "long"), value(0) {}
   virtual ~long_attribute_t() {}

   long_attribute_t &operator=(const long_attribute_t &src) {
      fName = src.fName;
      value = src.value;
      return *this;
   }

   virtual void reset() {
      value = 0;
   }
   virtual void set(long int val) {
      value = val;
   }
   virtual void *address() {
      return &value;
   }
   virtual std::string toString() {
      std::stringstream str;
      str << value;
      return str.str();
   }
   virtual int read(ixstream *ifx) {
      *ifx >> value;
      return 8;
   }

   long int value;
};

class float_attribute_t : public attribute_t {
 public:
   float_attribute_t() : attribute_t("", "float"), value(0) {}
   float_attribute_t(XString name) : attribute_t(name, "float"), value(0) {}
   virtual ~float_attribute_t() {}

   float_attribute_t &operator=(const float_attribute_t &src) {
      fName = src.fName;
      value = src.value;
      return *this;
   }

   virtual void reset() {
      value = 0;
   }
   virtual void set(float val) {
      value = val;
   }
   virtual void *address() {
      return &value;
   }
   virtual std::string toString() {
      std::stringstream str;
      str << value;
      return str.str();
   }
   virtual int read(ixstream *ifx) {
      *ifx >> value;
      return 4;
   }

   float value;
};

class double_attribute_t : public attribute_t {
 public:
   double_attribute_t() : attribute_t("", "double"), value(0) {}
   double_attribute_t(XString name) : attribute_t(name, "double"), value(0) {}
   ~double_attribute_t() {}

   double_attribute_t &operator=(const double_attribute_t &src) {
      fName = src.fName;
      value = src.value;
      return *this;
   }

   virtual void reset() {
      value = 0;
   }
   virtual void set(double val) {
      value = val;
   }
   virtual void *address() {
      return &value;
   }
   virtual std::string toString() {
      std::stringstream str;
      str << value;
      return str.str();
   }
   virtual int read(ixstream *ifx) {
      *ifx >> value;
      return 8;
   }

   double value;
};

class string_attribute_t : public attribute_t {
 public:
   string_attribute_t() : attribute_t("", "string") { reset(); }
   string_attribute_t(XString name) : attribute_t(name, "string") { reset(); }
   ~string_attribute_t() {}

   string_attribute_t &operator=(const string_attribute_t &src) {
      fName = src.fName;
      strncpy(value, src.value, 80);
      return *this;
   }

   virtual void reset() {
      strncpy(value, "", 80);
   }
   virtual void set(char *val) {
      strncpy(value, val, 80);
   }
   virtual void *address() {
      return &value;
   }
   virtual std::string toString() {
      return std::string(value);
   }
   virtual int read(ixstream *ifx) {
      std::string val;
      *ifx >> val;
      strncpy(value, val.c_str(), 80);
      return (val.size() + 7) / 4 * 4;
   }

   char value[80];
};

class anyURI_attribute_t : public attribute_t {
 public:
   anyURI_attribute_t() : attribute_t("", "anyURI") { reset(); }
   anyURI_attribute_t(XString name) : attribute_t(name, "anyURI") { reset(); }
   ~anyURI_attribute_t() {}

   anyURI_attribute_t &operator=(const anyURI_attribute_t &src) {
      fName = src.fName;
      strncpy(value, src.value, 80);
      return *this;
   }

   virtual void reset() {
      strncpy(value, "", 80);
   }
   virtual void set(char *val) {
      strncpy(value, val, 80);
   }
   virtual void *address() {
      return &value;
   }
   virtual std::string toString() {
      return std::string(value);
   }
   virtual int read(ixstream *ifx) {
      std::string val;
      *ifx >> val;
      strncpy(value, val.c_str(), 80);
      return (val.size() + 7) / 4 * 4;
   }

   char value[80];
};

class constant_attribute_t : public attribute_t {
 public:
   constant_attribute_t() : attribute_t("", "constant"), value(0)
   { reset(); }
   constant_attribute_t(XString name) : attribute_t(name, "constant"), value(0)
   { reset(); }
   ~constant_attribute_t() {}

   void reset() {
      if (value)
         delete value;
      value = new char[4];
      strncpy(value, "", 4);
   }
   void set(const char *str) {
      if (value)
         delete value;
      int size = (strlen(str) + 7) / 4 * 4;
      value = new char[size];
      strncpy(value, str, size);
   }

   constant_attribute_t &operator=(const constant_attribute_t &src) {
      fName = src.fName;
      set(src.value);
      return *this;
   }

   virtual void *address() {
      return value;
   }
   virtual std::string toString() {
      return std::string(value);
   }
   virtual int read(ixstream *ifx) {
      return 0;
   }

   char *value;
};

class element_t {
 public:
   element_t(TTree *tree)
    : fKey(0), fTree(tree), fRepeats(0) {}
   element_t(TTree *tree, XString name)
    : fKey(0), fTree(tree), fName(name), fRepeats(0) {}
   ~element_t() {}

   void add_attribute(attribute_t *attr) {
      fAttributes.push_back(attr);
   }
   void add_element(element_t *elem) {
      fElements.push_back(elem);
   }
   void add_key(int_attribute_t *attr) {
      fKey = attr;
   }
   void set_repeating() {
      fRepeats = 1;
   }

   int read(ixstream *ifx) {
      int size;
      *ifx >> size;
      if (size == 0)
         return 4;
      int seen;
      int reps;
      if (fRepeats) {
         *ifx >> reps;
         seen = 4;
      }
      else {
         seen = 0;
      }

      static int indent = 0;
      if (write_xml_output_to_stdout) {
         if (indent == 0) {
            std::cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
                      << std::endl
                      << "<HDDM class=\"s\" version=\"1.0\" "
                      << "xmlns=\"http://www.gluex.org/hddm\">"
                      << std::endl;
            ++indent;
         }
      }
      while (seen < size) {
         if (write_xml_output_to_stdout) {
            for (int i=0; i < indent; ++i)
               std::cout << "  ";
            std::cout << "<" << fName;
         }
         std::list<attribute_t*>::iterator ater;
         for (ater = fAttributes.begin(); ater != fAttributes.end(); ++ater) {
            seen += (*ater)->read(ifx);
            if (write_xml_output_to_stdout) {
               std::cout << " " << (*ater)->get_name() << "=\"" 
                         << (*ater)->toString() << "\"";
            }
         }
         if (fElements.size() == 0) {
            if (write_xml_output_to_stdout) {
               std::cout << " />" << std::endl;
            }
            if (fTree)
               fTree->Fill();
            if (fKey)
               ++fKey->value;
         }
         else {
            if (write_xml_output_to_stdout) {
               std::cout << ">" << std::endl;
               ++indent;
            }
            std::list<element_t*>::iterator eter;
            for (eter = fElements.begin(); eter != fElements.end(); ++eter) {
               seen += (*eter)->read(ifx);
            }
            if (write_xml_output_to_stdout) {
               --indent;
               for (int i=0; i < indent; ++i) {
                  std::cout << "  ";
               }
               std::cout << "</" << fName << ">" << std::endl;
            }
            if (fTree)
               fTree->Fill();
         }
         --reps;
      }
      assert (seen == size);
      if (fRepeats)
         assert (reps == 0);
      return size + 4;
   }

   std::list<attribute_t*> fAttributes;
   std::list<element_t*> fElements;
   int_attribute_t *fKey;
   TTree *fTree;
   XString fName;
   int fRepeats;

 private:
   element_t(const element_t &src);
   element_t &operator=(const element_t &src);
};

typedef std::map<XString,XString> attribute_list;
typedef std::map<XString,attribute_t*> attribute_table;

class TreeMaker
{
 public:
   TreeMaker(XString filename) {
      fRootFile = new TFile(S(filename), "recreate");
   }
   ~TreeMaker() {
      delete fRootFile;
   }

   void build(const DOMElement* elem, element_t *parent_element,
                                      attribute_list columns);
   int filltrees(ixstream *ifx, element_t *parent_element);
   int savetrees(element_t *parent_element);

 private:
   TreeMaker(const TreeMaker &src) {}
   TreeMaker operator=(const TreeMaker &src) {
      TreeMaker copy(*this);
      return copy;
   }

 protected:
   TFile *fRootFile;
   element_t *fTopElement;
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
      else if (strcmp(argV[argInd],"-x") == 0)
      {
         write_xml_output_to_stdout = 1;
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
   std::ostringstream doc;
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
   doc << std::endl;
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
         ofs << line << std::endl;
         doc << line << std::endl;
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
      ofs << line << std::endl;
      doc << line << std::endl;
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
   attribute_list columns;
   element_t root_element(0);
   builder.build(rootEl, &root_element, columns);
  
   int event_buffer_size;
   char *event_buffer = new char[event_buffer_size = 1000000];
   istreambuffer *isbuf = new istreambuffer(event_buffer,event_buffer_size);
   ixstream *ifx = new ixstream(isbuf);
   int integrity_check_mode = 0;
   int compression_mode = 0;
   while (reqcount && ifs->good())
   {
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
         ifx = new ixstream(isbuf);
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
         ixstream xstr(&sbuf);
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
      builder.filltrees(ifx, &root_element);
   }
   if (write_xml_output_to_stdout) {
      std::cout << "</HDDM>" << std::endl;
   }
   builder.savetrees(&root_element);
   TNamed docString("document_metadata_XML", doc.str().c_str());
   docString.Write();

   if (ifs != &std::cin)
   {
      ((std::ifstream*)ifs)->close();
   }

   XMLPlatformUtils::Terminate();
   return 0;
}

void TreeMaker::build(const DOMElement* elem, element_t *parent_element,
                                              attribute_list columns)
{
   // Recursively create TTree objects to hold the contents of the hddm model
   // in the form of a row/column table, like a relational database model.

   XString elemS(elem->getTagName());
   TTree *tree = new TTree(S(elemS), S(XString(elemS + " tree")));
   element_t *this_element = new element_t(tree, elemS);
   XString keyS("HDDM_MASTER_ORDERING_KEY");
   if (fColumns.find(keyS) == fColumns.end()) {
      int_attribute_t *key = new int_attribute_t("key");
      fColumns[keyS] = key;
   }
   tree->Branch("key", fColumns[keyS]->address(), "key/I");
   this_element->add_key((int_attribute_t*)fColumns[keyS]);
   XString repS(elem->getAttribute(X("maxOccurs")));
   int rep = (repS == "unbounded")? INT_MAX :
             (repS == "")? 1 :
             atoi(S(repS));
   if (explicit_repeat_count && rep > 1)
      this_element->set_repeating();
   attribute_list::iterator iter;
   for (iter = columns.begin(); iter != columns.end(); ++iter) {
      XString colS = iter->first;
      XString nameS = iter->second;
      XString typeS = fColumns[colS]->get_type();
      if (typeS == "int" || typeS == "boolean" || typeS == "Particle_t")
         tree->Branch(S(nameS), fColumns[colS]->address(),
                                              S(XString(nameS + "/I")));
      else if (typeS == "long")
         tree->Branch(S(nameS), fColumns[colS]->address(),
                                              S(XString(nameS + "/L")));
      else if (typeS == "float")
         tree->Branch(S(nameS), fColumns[colS]->address(),
                                              S(XString(nameS + "/F")));
      else if (typeS == "double")
         tree->Branch(S(nameS), fColumns[colS]->address(),
                                              S(XString(nameS + "/D")));
      else if (typeS == "string" || typeS == "anyURI")
         tree->Branch(S(nameS), fColumns[colS]->address(),
                                              S(XString(nameS + "/C")));
      else {
         tree->Branch(S(nameS), fColumns[colS]->address(),
                                              S(XString(nameS + "/C")));
      }
   }
    
   DOMNamedNodeMap* attrList = elem->getAttributes();
   int attrListLength = attrList->getLength();
   for (int a = 0; a < attrListLength; a++) {
      DOMNode* node = attrList->item(a);
      XString nameS(node->getNodeName());
      XString typeS(node->getNodeValue());
      XString colS(elemS + "_" + nameS);
      if (columns.find(nameS) != columns.end())
         nameS = colS;
      if (typeS == "int") {
         fColumns[colS] = new int_attribute_t(nameS);
         tree->Branch(S(nameS), fColumns[colS]->address(),
                                                S(XString(nameS + "/I")));
      }
      else if (typeS == "boolean") {
         fColumns[colS] = new boolean_attribute_t(nameS);
         tree->Branch(S(nameS), fColumns[colS]->address(),
                                                S(XString(nameS + "/I")));
      }
      else if (typeS == "Particle_t") {
         fColumns[colS] = new Particle_attribute_t(nameS);
         tree->Branch(S(nameS), fColumns[colS]->address(),
                                                S(XString(nameS + "/I")));
      }
      else if (typeS == "long") {
         fColumns[colS] = new long_attribute_t(nameS);
         tree->Branch(S(nameS), fColumns[colS]->address(),
                                                S(XString(nameS + "/L")));
      }
      else if (typeS == "float") {
         fColumns[colS] = new float_attribute_t(nameS);
         tree->Branch(S(nameS), fColumns[colS]->address(),
                                                S(XString(nameS + "/F")));
      }
      else if (typeS == "double") {
         fColumns[colS] = new double_attribute_t(nameS);
         tree->Branch(S(nameS), fColumns[colS]->address(),
                                                S(XString(nameS + "/D")));
      }
      else if (typeS == "string") {
         fColumns[colS] = new string_attribute_t(nameS);
         tree->Branch(S(nameS), fColumns[colS]->address(),
                                                S(XString(nameS + "/C")));
      }
      else if (typeS == "anyURI") {
         fColumns[colS] = new anyURI_attribute_t(nameS);
         tree->Branch(S(nameS), fColumns[colS]->address(),
                                                S(XString(nameS + "/C")));
      }
      else if (nameS == "minOccurs" || nameS == "maxOccurs") {
         continue;
      }
      else {
         fColumns[colS] = new constant_attribute_t(nameS);
         tree->Branch(S(nameS), fColumns[colS]->address(),
                                                S(XString(nameS + "/C")));
         ((constant_attribute_t*)fColumns[colS])->set(S(typeS));
      }
      columns[colS] = nameS;
      this_element->add_attribute(fColumns[colS]);
   }
   parent_element->add_element(this_element);

   DOMNodeList* contList = elem->getChildNodes();
   int contLength = contList->getLength();
   for (int c = 0; c < contLength; c++) {
      DOMNode* cont = contList->item(c);
      short type = cont->getNodeType();
      if (type == DOMNode::ELEMENT_NODE) {
         DOMElement* contEl = (DOMElement*) cont;
         build(contEl, this_element, columns);
      }
   }
}

int TreeMaker::filltrees(ixstream *ifx, element_t *parent_element)
{
   element_t *HDDMelement = *parent_element->fElements.begin();
   std::list<element_t*>::iterator iter;
   int size = 0;
   for (iter = HDDMelement->fElements.begin();
        iter != HDDMelement->fElements.end();
        ++iter)
   {
      size += (*iter)->read(ifx);
#ifdef VERBOSE_HDDM_LOGGING
      XString cnameS((*iter)->fName);
      std::cerr << "hddm-root : top-level tag " << S(cnameS)
                << " found with size " << size
                << std::endl;
#endif
   }
   return size;
}

int TreeMaker::savetrees(element_t *parent_element)
{
   int count = 0;
   std::list<element_t*>::iterator iter;
   for (iter = parent_element->fElements.begin();
        iter != parent_element->fElements.end();
        ++iter)
   {
      if ((*iter)->fTree) {
         (*iter)->fTree->Write();
         ++count;
      }
      count += savetrees(*iter);
   }
   return count;
}
