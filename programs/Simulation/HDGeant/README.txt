
			Build Notes for HDGeant
	    		    Richard Jones
	    		    July 10, 2001
                      (updated October 30, 2003)

1) Since you are reading this, you have already done the first step.
   From your cvs source directory you typed:

   halld> cvs checkout HDGeant

2) Download the xerces-c xml library from xml.apache.org and unpack
   it somewhere on your system or GlueX working area.  Getting the sources
   and doing the build yourself (next step) makes sure that you have a
   working installation for your configuration.

3) Build the xerces xml library for your system.

   This is pretty simple.  The instructions are found on the xml.apache.org
   web site.  There are just three steps: define XERCESCROOT to point to the
   base directory where you unpacked xerces-c, then runConfigure and gmake.
   The result is a shared library in the directory xerces-c/lib.

4) Somewhere, perferably at the top of your cvs source directory you
   should make a script called setup that sets up some environment
   variables that are needed to locate the CERN libraries on your system.
   What that looks like depends on your shell, but for tcsh it looks like:

   halld> cat setup
   setenv CVS_RSH ssh
   setenv CVSROOT <userid>@login1.jlab.org:/halld/cvsroot
   setenv CERN <your local cernlib path, contains pro, old, new, 2000...>
   setenv CERN_LEVEL <pro or whatever you use>
   setenv CERN_ROOT ${CERN}/${CERN_LEVEL}
   setenv HALLD_ROOT <your local source directory path>
   setenv PATH ${HALLD_ROOT}/bin.Linux:${CERN_ROOT}/bin:${PATH}
   setenv BUILDS ${HALLD_ROOT}
   setenv HALLDLIB ${HALLD_ROOT}/lib.Linux
   setenv QQ_DIR ${HALLD_ROOT}/libmcfast/qq_v9_2b_Linux+2.2
   setenv STDHEP_DIR ${HALLD_ROOT}/libmcfast/stdhep_v4_08_Linux+2.2
   setenv XERCESCROOT <your installation path to xerces-c>
   setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/${XERCESCROOT}/lib

   For the bash or ksh shells you should use the export command instead of
   setenv.  You will need to source this file before every session, (or
   invoke it with the . operator for ksh or bash).

5) Now you need to check out the geometry package and build it.

   halld> cvs checkout hdds
   halld> cd hdds
   halld> make hddsGeant3.f

6) Now go to HDGeant and build the interactive version.

   halld> cd ../HDGeant
   halld> make hdgeant++.x

   Now find the start of the long string of errors that was just produced by
   make, and find out what you did wrong in steps 1..5.  Iterate until the
   package builds without errors.

7) Start up interactive Geant and plot the detector.

   halld> ./hdgeant++.x
   ... lots of output
   GEANT> exec ray#init
halld
10000
10000
10000
   GEANT> exec ray#draw

Now you are on your way...  Right now why don't you stop and send an email
to richard.t.jones@uconn.edu letting me know how things went.

Richard Jones
Storrs, Connecticut


Special note for users of Redhat 9
----------------------------------

During the link step you may get an error that the library gcc_s is not
found on your system.  This seems to be a fault of g77 under gcc 3.2.
The following workaround was verified to work by Ed Brash.

Building the non-interactive hdgeant requires a library called libgcc_s.a
(appears to be a g77 requirement).  This library is no longer present in
gcc3.2. The fix is to make a softlink from libgcc.a to libgcc_s.a in the
/usr/lib/gcc-lib/i386-redhat-linux/3.2 directory. 


Mailing lists
-------------

Please post problem reports, suggestions and fixes to the appropriate
forum on the GlueX message board.  The moderator of the list will respond
to any requests or forward them to the correct party.

http://zeus.phys.uconn.edu/forums/gluex.org
