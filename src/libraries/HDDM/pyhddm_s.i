// author: Sean Dobbs (s-dobbs@northwestern.edu), 2015
//
// This file defines the interface for accessing MC HDDM files through python
// It is a known problem that python 2 does not play well with C++ iostreams
// HDDM IO is streamed, so a wrapper class is needed.  I tried various methods
// that were more idiomatically python, but they didn't work right for whatever reason.
// The method used below is just a standard hand-rolled wrapper class.

%module pyhddm_s
%{
#include "hddm_s.h"
#include "hddm_s.hpp"
%}

// SWIG does not like this variable for some reason
%ignore HDDM_s_DocumentString;

//%include "hddm_s.h"
%include "hddm_s.hpp"
%include "std_string.i"

//// proxy classes for accessing streamed IO in python
%inline %{
  #include <iostream>
  #include <fstream>

  class hddm_istream_proxy {
  public:
    hddm_istream_proxy( const std::string& fname ) {
      open( fname );
    }

    virtual ~hddm_istream_proxy() {
      delete ifs;
      delete hddm_ifs;
    }

    void open( const std::string& fname ) {
      filename = fname;

      // error check?
      ifs = new std::ifstream(fname.c_str());
      if (ifs && ifs->is_open()) 
	hddm_ifs = new hddm_s::istream(*ifs);
      else  {
	ifs = NULL;
	// maybe should throw an exception?
	return;
      }
    }
    
    void close() {
	delete hddm_ifs;
	delete ifs;
	ifs = NULL;
	hddm_ifs = NULL;
    }

    void reset() {
      close();
      open(filename);
    }

    bool eof() {
      return ifs->eof();
    }
    

    bool read() {
      // don't read in anything if the file doesn't exist
      if(ifs == NULL)
	return false;

      // stop reading once we hit the end of file
      // and close everything out
      if(ifs->eof()) {
	delete hddm_ifs;
	delete ifs;
	ifs = NULL;
	hddm_ifs = NULL;
	return false;
      }

      (*hddm_ifs) >> record;

      // skip over comments
      while (true) {
	hddm_s::PhysicsEvent &re
	  = record.getPhysicsEvent();
	int runnum = re.getRunNo();
	int eventnum = re.getEventNo();

	if (runnum == 0 && eventnum == 0) {
	  (*hddm_ifs) >> record;
	} else {
	  break;
	}
      }
      
      return true;
    }

    void skip(int nskip) {
      hddm_ifs->skip(nskip);
    }

    //////////////// get methods ////////////////////////////////

    int getRunNumber() {
      hddm_s::PhysicsEvent &re
	= record.getPhysicsEvent();
      return re.getRunNo();
    }

    int getEventNumber() {
      hddm_s::PhysicsEvent &re
	= record.getPhysicsEvent();
      return re.getEventNo();
    }

    // placeholder function
    // fix this to 1 right now
    long getUid() {
      return 1L;
    }

    hddm_s::PhysicsEvent &getEvent() {
      return record.getPhysicsEvent();
    }

    std::string getFilename() {
      return filename;
    }

  private:
    hddm_s::HDDM record;
    std::ifstream *ifs;
    hddm_s::istream *hddm_ifs;

    std::string filename;
  };
  

%}
