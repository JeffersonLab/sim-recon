This is JANA integration plugin that allows to use CCDB from JANA

Steps to run JANA with CCDB:



PREPARE CCDB
============

1. Install CCDB. Follow the instructions of file install.linux.txt located in 
   the $CCDB_HOME directory

2. Don't forget to setup CCDB environment in a terminal you going to run JANA with CCDB in



ACCESS TO DATABASE
==================

The simpliest thing to try CCDB would be to download sqlite file tith database.

SqLite file could be downloaded from:
https://halldweb1.jlab.org/ccdb/gluex_calib.sqlite



PREPARE JANA
============

1) copy $CCDB_HOME/janaccdb to Jana plugins directory

     cp -r $CCDB_HOME/janaccdb $JANA_HOME/src/plugins


2) modify $JANA_HOME/src/plugins/Makefile, add 'janaccdb' to DIRS string

     nano $JANA_HOME/src/plugins/Makefile
     DIRS := janadot janactl ... janaccdb


3) Rebuild JANA 



RUN ANALYSIS
============

1) Set JANA_CALIB_URL to mysql or sqlite connection

      export JANA_CALIB_URL=sqlite:///path/to/gluex.sqlite


2) Now it is possible to run analysis specifying --plugin=janaccdb, like:

      hd_root --plugin=phys_tree --plugin=janaccdb hdgeant_smeared.hddm



CCDB JANA PLUGIN DEBUG OUTPUT
=============================
If one builds JANA with CCDB_DEBUG_OUTPUT preprocessor definition:

     make CPPFLAGS=-DCCDB_DEBUG_OUTPUT

This makes janaccdb plugin to print explicit information to jout



MORE DOCUMENTATION
==================

More documentation may be found in $CCDB_HOME/doc directory... with time...
