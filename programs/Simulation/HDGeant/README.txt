
			Build Notes for HDGeant
	    		    Richard Jones
	    		    July 10, 2001

1) Since you are reading this, you have already done the first step.
   From your cvs source directory you typed:

   halld> cvs checkout HDGeant
   
2) Somewhere, perferrably at the top of your cvs source directory you
   should make a script called setup that sets up some environment
   variables that are needed to locate the CERN libraries on your system.
   What that looks like depends on your shell, but for tcsh it looks like:

   halld> cat setup
   setenv CVSROOT :pserver:<userid>@jlabs1.jlab.org:/halld/cvsroot
   setenv CERN <your local cernlib path, contains pro, old, new, 2000...>
   setenv CERN_ROOT ${CERN}/<pro or whatever you use>
   setenv HALLD_ROOT <your local source directory path>
   setenv PATH ${PATH}:${HALLD_ROOT}/bin.Linux:${CERN_ROOT}/bin
   setenv BUILDS ${HALLD_ROOT}
   setenv HALLDLIB ${HALLD_ROOT}/lib.Linux
   setenv QQ_DIR ${HALLD_ROOT}/libmcfast/qq_v9_2b_Linux+2.2
   setenv STDHEP_DIR ${HALLD_ROOT}/libmcfast/stdhep_v4_08_Linux+2.2
   setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/${HALLD_ROOT}/xerces-c/lib

   For the bash or ksh shells you should use the export command instead of
   setenv.  You will need to source this file before every session, (or
   invoke it with the . operator for ksh or bash).

3) Download the xerces xml library from cvs.

   halld> cd ${HALLD_ROOT}
   halld> cvs checkout xerces-c

4) Now you need to check out the geometry package and build it.

   halld> cvs checkout hdds
   halld> cd hdds
   halld> make hddsGeant3.f

5) Now go to HDGeant and build the interactive version.

   halld> cd ../HDGeant
   halld> make hdgeant++.x

   Now find the start of the long string of errors that was just produced by
   make, and find out what you did wrong in steps 1..5.  Iterate until the
   package builds without errors.

6) Start up interactive Geant and plot the detector.

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


