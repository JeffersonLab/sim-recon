
			Build Notes for HDDM Tools
	    		    Richard Jones
	    		  September 10, 2003

This document is intended as a quick-start guide for building and using
the hdds tools.  For more information and for discussion of features and
bugs, please go to http://portal.gluex.org and look for "forums".

1) Since you are reading this, you have already done the first step.
   From your cvs source directory you typed:

   halld> cvs checkout hddm
   
2) Download the xerces-c xml library from xml.apache.org and unpack
   it somewhere on your system or GlueX working area.  Getting the sources
   and doing the build yourself (next step) makes sure that you have a
   working installation for your configuration.

3) Build the xerces xml library for your system.

   This is pretty simple.  The instructions are found on the xml.apache.org
   web site.  There are just three steps: define XERCESCROOT to point to the
   base directory where you unpacked xerces-c, then runConfigure and gmake.
   The result is a shared library in the directory xerces-c/lib.

4) Download and build the Xalan xml tool from xml.apache.org and unpack
   it somewhere on your system or GlueX working area.  If you have java
   installed on your system, downloading Xalan-j is simple and does not
   require a build step.  Otherwise download Xalan-c and build it in a
   similar way to the one described under step 3.

5) Somewhere, perferrably at the top of your cvs source directory you
   should make a script called setup that sets up some environment
   variables that are needed to locate the XML libraries on your system.
   What that looks like depends on your shell, but for tcsh it looks like:

   halld> cat setup
   setenv CVS_RSH ssh
   setenv CVSROOT <userid>@login1.jlab.org:/halld/cvsroot
   setenv CERN <your local cernlib root, contains pro, old, new, 2000...>
   setenv CERN_VERSION <your current cernlib release: pro, old, new, 2000...>
   setenv CERN_ROOT ${CERN}/${CERN_VERSION}
   setenv HALLD_ROOT <your local source directory path>
   setenv PATH ${PATH}:${HALLD_ROOT}/bin.Linux:${CERN_ROOT}/bin
   setenv XERCESCROOT <your Xerces-C installation path>
   setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/${XERCESCROOT/lib
   setenv JAVAROOT <your java j2re or j2sdk installation path>
   setenv XALANJROOT <your xalan-j installation path>
   setenv CLASSPATH=${XALANJROOT}/bin/xalan.jar

   For the bash or ksh shells you should use the export command instead of
   setenv.  You will need to source this file before every session, (or
   invoke it with the . operator for ksh or bash).  Note that the JAVAROOT
   definition is not necessary if you do not intend to use java tools.
   If you do use java, you should also make sure that the java command is
   in your path, so that the command java is found.

6) Check out the xstream package from the halld cvs repository and build
   and install it in your build area in the usual way.

   halld> cd <directory-containing-hddm>
   halld> cvs checkout xstream
   halld> make install -C xstream

7) Now if everything was installed correctly, you can build the hddm tools
   by going to the hddm directory you created in step 1 and doing make.

   halld> cd hddm
   halld> make hddm-c hddm-xml xml-hddm

8) The tools hddm-schema, schema-hddm and xml-xml are written in perl and
   require Xerces-Perl from xml.apache.org to be installed if you want to
   use them.  Installation instructions are found at xml.apache.org.  By
   default these tools use the Xalan-j package.  If you installed the c++
   version of Xalan instead, go into hddm-schema and schema-hddm and
   comment out the line where xalan is started under java and remove the
   comment from the line where xalan-c is invoked.

9) For information about the design and use of the hddm tools, see the
   documentation in index.html.  For more information about the xstream
   binary streams interface, see http://xstream.sourceforge.net .
