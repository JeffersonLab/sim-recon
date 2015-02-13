// author: Sean Dobbs (s-dobbs@northwestern.edu), 15 April 2014
//
// This file defines the interface for accessing REST files through python
// It is a known problem that python 2 does not play well with C++ iostreams
// HDDM IO is streamed, so a wrapper class is needed.  I tried various methods
// that were more idiomatically python, but they didn't work right for whatever reason.
// The method used below is just a standard hand-rolled wrapper class.

%module pyhddm_r
%{
#include "hddm_r.h"
#include "hddm_r.hpp"
%}

// SWIG does not like this variable for some reason
%ignore HDDM_r_DocumentString;

//%include "hddm_r.h"
%include "hddm_r.hpp"
%include "std_string.i"

//// proxy classes for accessing streamed IO in python
%inline %{
  #include <iostream>
  #include <fstream>

  class hddm_istream_proxy {
  public:
    hddm_istream_proxy( const std::string& fname ) {
      ifs = NULL;
      hddm_ifs = NULL;
      open( fname );
    }

    virtual ~hddm_istream_proxy() {
      if(hddm_ifs != NULL)
	delete hddm_ifs;
      //if(ifs != NULL)
      //delete ifs;
    }

    void open( const std::string& fname ) {
      filename = fname;

      // error check?
      ifs = new std::ifstream(fname.c_str());
      if (ifs && ifs->is_open())  {
	hddm_ifs = new hddm_r::istream(*ifs);
      } else {
	ifs = NULL;
	// maybe should throw an exception?
	return;
      }
    }
    
    void close() {
      if(hddm_ifs != NULL) {
	delete hddm_ifs;
	hddm_ifs = NULL;
      }
      if(ifs != NULL) {
	//delete ifs;
	ifs = NULL;
      }
    }

    void reset() {
      //ifs->seekg(0);
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
	//delete hddm_ifs;
	//delete ifs;
	//ifs = NULL;
	//hddm_ifs = NULL;
	return false;
      }

      (*hddm_ifs) >> record;

      // skip over comments
      /*
      while (true) {
	hddm_r::ReconstructedPhysicsEvent &re
	  = record.getReconstructedPhysicsEvent();
	int runnum = re.getRunNo();
	int eventnum = re.getEventNo();

	if (runnum == 0 && eventnum == 0) {
	  (*hddm_ifs) >> record;
	} else {
	  break;
	}
      }
      */

      return true;
    }

    void skip(int nskip) {
      hddm_ifs->skip(nskip);
    }

    //////////////// get methods ////////////////////////////////

    int getRunNumber() {
      hddm_r::ReconstructedPhysicsEvent &re
	= record.getReconstructedPhysicsEvent();
      return re.getRunNo();
    }

    int getEventNumber() {
      hddm_r::ReconstructedPhysicsEvent &re
	= record.getReconstructedPhysicsEvent();
      return re.getEventNo();
    }

    // placeholder function
    // fix this to 1 right now
    long getUid() {
      return 1L;
    }

    hddm_r::ReconstructedPhysicsEvent &getEvent() {
      return record.getReconstructedPhysicsEvent();
    }

    std::string getFilename() {
      return filename;
    }

  private:
    hddm_r::HDDM record;
    std::ifstream *ifs;
    hddm_r::istream *hddm_ifs;

    std::string filename;
  };
  

%}
