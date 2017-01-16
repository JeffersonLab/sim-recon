#!/bin/env python
#
# resort.py - reads in a hddm_s file and writes out the events in
#             strict order of increasing event number.
#
# author: richard.t.jones at uconn.edu
# version: december 31, 2016

import hddm_s
import sys

if len(sys.argv) != 3:
   print "Usage: resort.py <infile.hddm> <outfile.hddm>"
   sys.exit(1)
else:
   out = hddm_s.ostream(sys.argv[2])

cache = {}
nextno = 1
geom_pending = 0

for r in hddm_s.istream(sys.argv[1]):
   for geom in r.getGeometrys():
      geom_pending = 1
      md5simul = geom.md5simulation
      md5smear = geom.md5smear
      md5recons = geom.md5reconstruction
      r.deleteGeometrys()
   thisno = r.getPhysicsEvent().eventNo
   if thisno == nextno:
      if geom_pending:
         geom = r.addGeometrys()
         geom[0].md5simulation = md5simul
         geom[0].md5smear = md5smear
         geom[0].md5reconstruction = md5recons
         geom_pending = 0
      out.write(r)
      nextno += 1
      while nextno in cache:
         out.write(cache[nextno])
         del cache[nextno]
         nextno += 1
   else:
      cache[thisno] = r
