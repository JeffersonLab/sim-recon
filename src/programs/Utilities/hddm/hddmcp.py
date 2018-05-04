#!/bin/env python
#
# hddmcp.py - python alternative to the old hddmcp, multithreaded.
#
# author: richard.t.jones at uconn.edu
# version: march 23, 2018

import threading
import hddm_s
import sys

def usage():
   print "Usage: hddmcp.py [-n #threads] <file1.hddm> <file2.hddm> ... <outfile.hddm>"
   sys.exit(1)

def copy(fin):
   for r in fin:
      fout.write(r)

nthreads = 1
infiles = []
iarg = 1
while iarg < len(sys.argv):
   arg = sys.argv[iarg]
   iarg += 1
   if arg == '-n':
      nthreads = int(sys.argv[iarg])
      iarg += 1
      continue
   elif arg[0:2] == '-n':
      nthreads = int(arg[2:])
      continue
   infiles.append(arg)
outfile = infiles.pop()
if len(infiles) < 1:
   usage()
fout = hddm_s.ostream(outfile)
fout.compression = hddm_s.k_bz2_compression

for infile in infiles:
   print "starting", infile
   fin = hddm_s.istream(infile)
   threads = []
   for i in range(0, nthreads):
      t = threading.Thread(target=copy, args=(fin,))
      t.start()
      threads.append(t)
   for t in threads:
      t.join()
   print "finished", infile
