
			Build Note for GXtwist
	    		    Richard Jones
    		          January 25, 2006

GXtwist is the GEANT3-based simulation program for the GlueX tagging
spectrometer and associated electron beam line.  The geometry is defined
in a xml document (HDDS-1.0 schema) and translated into a Geant3 volume
description using the tool hdds-geant.  This document is intended as a
quick-start guide for building and using the hdds tools.  For more
information and for discussion of features and bugs, please go to
http://portal.gluex.org and look for "forums".

1) Since you are reading this, you have already done the first step.
   From your svn working source directory you typed:

   halld> svn checkout https://halldsvn.jlab.org/repos/trunk/sim-recon/src

   To build the gxtwist executable you also need to check out the hdds
   module. The following line places this at a level parallel with src.

   halld> svn checkout https://halldsvn.jlab.org/repos/trunk/hdds

2) Download the xerces-c xml library from xml.apache.org and unpack
   it somewhere on your system or GlueX working area.  Getting the sources
   and doing the build yourself (next step) makes sure that you have a
   working installation for your configuration.

3) Build the xerces xml library for your system.

   This is pretty simple.  The instructions are found on the xml.apache.org
   web site.  There are just three steps: define XERCESCROOT to point to the
   base directory where you unpacked xerces-c, then runConfigure and gmake.
   The result is a shared library in the directory xerces-c/lib.

4) Somewhere, preferably at the top of your working source directory you
   should make a script called setup that sets up some environment
   variables that are needed to locate the CERN libraries on your system.
   What that looks like depends on your shell, but for tcsh it looks like:

   halld> cat setup
   ...
   setenv CERN <your local cernlib path, contains pro, old, new, 2000...>
   setenv CERN_LEVEL <pro or whatever you use>
   setenv CERN_ROOT ${CERN}/${CERN_LEVEL}
   ...

   For the bash or ksh shells you should use the export command instead of
   setenv.  You will need to source this file before every session, (or
   invoke it with the . operator for ksh or bash).

5) Now go to gxtwist and build the interactive version.

   halld> cd ../gxtwist
   halld> make gxtwist++

   Now find the start of the long string of errors that was just produced by
   make, and find out what you did wrong in steps 1..5.  Iterate until the
   package builds without errors.

6) Start up interactive Geant and plot the detector.

   halld> ./gxtwist++
   ... lots of output
   GEANT> next; dcut hill y 0 1 10 .01 .01

7) Install the field map in your working directory and generate events

   The field map is not stored in svn as a part of the code repository
   because it is too large (over 200MB).  Before you can start tracking
   particles through the tagger you need to fetch a copy of the file
   taggerBfield-XXX.map into your working directory.  There are two
   varieties of the field map, one with the quadrupole magnet on and
   one with it off.  They are available at the following web address.

   http://zeus.phys.uconn.edu/halld/tagger/simulation/taggerBfield-quad.map
   http://zeus.phys.uconn.edu/halld/tagger/simulation/taggerBfield-noquad.map

Now you are on your way.  If you ever want to fiddle with the geometry,
you need to know where to find it.  The hdds description of the tagger
is in the hdds directory inside this folder, not in the upper-level hdds
project in svn.  It is checked out automatically when you check out gxtwist.

Richard Jones
Storrs, Connecticut
richard.t.jones@uconn.edu
