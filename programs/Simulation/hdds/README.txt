
			Build Notes for HDDS Tools
	    		    Richard Jones
	    		  September 10, 2003

This document is intended as a quick-start guide for building and using
the hdds tools.  For more information and for discussion of features and
bugs, please go to http://zeus.phys.cuonn.edu/forums/gluex.org

1) Since you are reading this, you have already done the first step.
   From your cvs source directory you typed:

   halld> cvs checkout hdds
   
2) Download the xerces-c xml library from xml.apache.org and unpack
   it somewhere on your system or GlueX working area.  Getting the sources
   and doing the build yourself (next step) makes sure that you have a
   working installation for your configuration.

3) Build the xerces xml library for your system.

   This is pretty simple.  The instructions are found on the xml.apache.org
   web site.  There are just three steps: define XERCESCROOT to point to the
   base directory where you unpacked xerces-c, then runConfigure and gmake.
   The result is a shared library in the directory xerces-c/lib.

4) Somewhere, perferrably at the top of your cvs source directory you
   should make a script called setup that sets up some environment
   variables that are needed to locate the XML libraries on your system.
   What that looks like depends on your shell, but for tcsh it looks like:

   halld> cat setup
   setenv CVS_RSH ssh
   setenv CVSROOT <userid>@login1.jlab.org:/halld/cvsroot
   setenv CERN <your local cernlib path, contains pro, old, new, 2000...>
   setenv CERN_ROOT ${CERN}/<pro or whatever you use>
   setenv HALLD_ROOT <your local source directory path>
   setenv PATH ${PATH}:${HALLD_ROOT}/bin.Linux:${CERN_ROOT}/bin
   setenv BUILDS ${HALLD_ROOT}
   setenv HALLDLIB ${HALLD_ROOT}/lib.Linux
   setenv QQ_DIR ${HALLD_ROOT}/libmcfast/qq_v9_2b_Linux+2.2
   setenv STDHEP_DIR ${HALLD_ROOT}/libmcfast/stdhep_v4_08_Linux+2.2
   setenv XERCESCROOT <your Xerces-C installation path>
   setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/${XERCESCROOT/lib
   setenv JAVAROOT <your java j2re or j2sdk installation path>

   For the bash or ksh shells you should use the export command instead of
   setenv.  You will need to source this file before every session, (or
   invoke it with the . operator for ksh or bash).  Note that the JAVAROOT
   definition is not necessary if you do not intend to use java tools.

6) Now if everything was installed correctly, you can build the hdds tools
   by going to the hdds directory you created in step 1 and doing make.

   halld> cd hdds
   halld> make hdds-geant hdds-mcfast hdds-root

7) The tool hdds-geant will run immediately.  If you invoke it by the
   command "make hddsGeant3.f" then it will use hdds-geant to create the
   above file and then copy it into ../HDGeant where it is needed for the
   compilation of hdgeant.  The tool hdds-mcfast requires that the cvs
   package libmcfast has been checked out in your working area, and that
   the environment variable MCFAST_DIR points to libmcfast/mcfast* where
   mcfast* is the long name of the production mcfast library.  If you
   invoke hdds-mcfast by the command "make hddsMCfast.db" then it will
   use hdds-mcfast to create the above file and then copy it into
   ../HDFast where it is needed by HDFast at run time.

8) For information about the design and use of the hddm tools, see the
   documentation in index.html.
