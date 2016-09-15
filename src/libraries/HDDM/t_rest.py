#!/usr/bin/env python
#
# t_rest.py - standard suite of tests for hddm i/o performance,
#             with the aim of validating the correct behavior
#             and benchmarking the performance of the library
#             in single/multi-threaded mode, with and without
#             compression, in sequential and random access.
#
# author: richard.t.jones at uconn.edu
# version: august 12, 2016
#
# usage: t_rest.py [-n <events>] [-p <threads>] <input_restfile.hddm>
#
# notes:
# 1) The input file is read only once, and copies are made in the local
#    directory in different compression modes. There should be enough
#    space in the local directory to hold six copies of the entire
#    input file, or its first N events if option -n is given: one
#    single-threaded output and one multi-threaded output file for
#    each compression mode: bz2, zlib, and none.

from __future__ import absolute_import, division, print_function, unicode_literals

maxevents = -1
multithreads = 4

import os
import sys
import time
import argparse

import threading
import hddm_r

lock = threading.Lock()

def event_reader(tid, fin):
   global postable
   for rec in fin:
      evno = rec.getReconstructedPhysicsEvent().eventNo
      if evno > 0 and evno % 10000 == 0:
         print("   event", evno, tid)
      postable[evno] = fin.position

def event_writer(tid, fin, fout):
   for rec in fin:
      evno = rec.getReconstructedPhysicsEvent().eventNo
      if evno > 0 and evno % 10000 == 0:
         print("   event", evno, tid)
      fout.write(rec)

def random_reader(tid, fin, post, verbose):
   for evno in sorted(post):
      fin.position = post[evno]
      rec = fin.read()
      if evno != rec.getReconstructedPhysicsEvent().eventNo:
         lock.acquire()
         print("bad event lookup for event", evno, "in file", fin)
         print("requested", evno, "but got", rec.getReconstructedPhysicsEvent().eventNo)
         lock.release()
         break
      elif fin.position != post[evno]:
         lock.acquire()
         print("bad record position for event", evno, "in file", fin)
         print("requested", post[evno], "but got", fin.position)
         lock.release()
         break
      elif evno > 0 and evno % 10000 == 0:
         print("   event", evno)

# main processing thread begins here

if __name__ == '__main__' :
   parser = argparse.ArgumentParser()
   parser.add_argument("input_restfile", help="input hddm file in REST format")
   parser.add_argument("-n", "--Nevents", 
                       help="number of events to read from input hddm file")
   parser.add_argument("-p", "--Nthreads", 
                       help="number of i/o threads in multi-threading tests")
   args = parser.parse_args()
   if args.Nevents:
      maxevents = int(args.Nevents)
   if args.Nthreads:
      multithreads = int(args.Nthreads)

   rest_z = "t_rest_z.hddm"
   rest_bz2 = "t_rest_bz2.hddm"
   rest_none = "t_rest_none.hddm"

   print("test 0: copying records from input file")

   fout0 = hddm_r.ostream(rest_none)
   fout1 = hddm_r.ostream(rest_z)
   fout1.compression = hddm_r.k_z_compression
   fout2 = hddm_r.ostream(rest_bz2)
   fout2.compression = hddm_r.k_bz2_compression
   fin = hddm_r.istream(args.input_restfile)
   for rec in fin:
      evno = rec.getReconstructedPhysicsEvent().eventNo
      if evno > 0 and evno % 10000 == 0:
         print("   event", evno)
      fout0.write(rec)
      fout1.write(rec)
      fout2.write(rec)
      if fin.recordsRead == maxevents:
         break
   del fout0
   del fout1
   del fout2

   print("test 1: reading single-threaded from input file, bz2 compression")

   fin = hddm_r.istream(rest_bz2)
   postable_bz2 = {}
   t0 = time.time()
   for rec in fin:
      evno = rec.getReconstructedPhysicsEvent().eventNo
      if evno > 0 and evno % 10000 == 0:
         print("   event", evno)
      postable_bz2[evno] = fin.position
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")

   print("test 2: reading single-threaded from input file, zlib compression")

   fin = hddm_r.istream(rest_z)
   postable_z = {}
   t0 = time.time()
   for rec in fin:
      evno = rec.getReconstructedPhysicsEvent().eventNo
      if evno > 0 and evno % 10000 == 0:
         print("   event", evno)
      postable_z[evno] = fin.position
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")

   print("test 3: reading single-threaded from input file, no compression")

   fin = hddm_r.istream(rest_none)
   postable_none = {}
   t0 = time.time()
   for rec in fin:
      evno = rec.getReconstructedPhysicsEvent().eventNo
      if evno > 0 and evno % 10000 == 0:
         print("   event", evno)
      postable_none[evno] = fin.position
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")
   tread = tlapse

   print("test 4: writing single-threaded to output file, bz2 compression")

   fout = hddm_r.ostream("." + rest_bz2)
   fout.compression = hddm_r.k_bz2_compression
   t0 = time.time()
   for rec in hddm_r.istream(rest_none):
      evno = rec.getReconstructedPhysicsEvent().eventNo
      if evno > 0 and evno % 10000 == 0:
         print("   event", evno)
      fout.write(rec)
   tlapse = time.time() - t0 - tread + 1e-29
   print("   ", fout.recordsWritten, "events written in", tlapse, "seconds,",
         fout.recordsWritten / tlapse, "ev/s,", fout.bytesWritten / tlapse / 1e6, "MB/s")
   del fout

   print("test 5: writing single-threaded to output file, zlib compression")

   fout = hddm_r.ostream("." + rest_z)
   fout.compression = hddm_r.k_z_compression
   t0 = time.time()
   for rec in hddm_r.istream(rest_none):
      evno = rec.getReconstructedPhysicsEvent().eventNo
      if evno > 0 and evno % 10000 == 0:
         print("   event", evno)
      fout.write(rec)
   tlapse = time.time() - t0 - tread + 1e-29
   print("   ", fout.recordsWritten, "events written in", tlapse, "seconds,",
         fout.recordsWritten / tlapse, "ev/s,", fout.bytesWritten / tlapse / 1e6, "MB/s")
   del fout

   print("test 6: writing single-threaded to output file, no compression")

   fout = hddm_r.ostream("." + rest_none)
   t0 = time.time()
   for rec in hddm_r.istream(rest_none):
      evno = rec.getReconstructedPhysicsEvent().eventNo
      if evno > 0 and evno % 10000 == 0:
         print("   event", evno)
      fout.write(rec)
   tlapse = time.time() - t0 - tread + 1e-29
   print("   ", fout.recordsWritten, "events written in", tlapse, "seconds,",
         fout.recordsWritten / tlapse, "ev/s,", fout.bytesWritten / tlapse / 1e6, "MB/s")
   del fout

   print("test 7: reading multi-threaded from input file, bz2 compression")

   fin = hddm_r.istream("." + rest_bz2)
   postable = {}
   threads = []
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = event_reader, args=(tid, fin))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")
   postable_dot_bz2 = postable

   print("test 8: reading multi-threaded from input file, zlib compression")

   fin = hddm_r.istream("." + rest_z)
   postable = {}
   threads = []
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = event_reader, args=(tid, fin))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")
   postable_dot_z = postable

   print("test 9: reading multi-threaded from input file, no compression")

   fin = hddm_r.istream("." + rest_none)
   postable = {}
   threads = []
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = event_reader, args=(tid, fin))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")
   tread = tlapse
   postable_dot_none = postable

   print("test 10: writing multi-threaded to output file, bz2 compression")

   fin = hddm_r.istream(rest_none)
   fout = hddm_r.ostream(".." + rest_bz2)
   fout.compression = hddm_r.k_bz2_compression
   threads = []
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = event_writer, args=(tid, fin, fout))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 - tread + 1e-30
   print("   ", fout.recordsWritten, "events written out in", tlapse, "seconds,",
         fout.recordsWritten / tlapse, "ev/s,", fout.bytesWritten / tlapse / 1e6, "MB/s")
   del fout

   print("test 11: writing multi-threaded to output file, zlib compression")

   fin = hddm_r.istream(rest_none)
   fout = hddm_r.ostream(".." + rest_z)
   fout.compression = hddm_r.k_z_compression
   threads = []
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = event_writer, args=(tid, fin, fout))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 - tread + 1e-30
   print("   ", fout.recordsWritten, "events written out in", tlapse, "seconds,",
         fout.recordsWritten / tlapse, "ev/s,", fout.bytesWritten / tlapse / 1e6, "MB/s")
   del fout

   print("test 12: writing multi-threaded to output file, no compression")

   fin = hddm_r.istream(rest_none)
   fout = hddm_r.ostream(".." + rest_none)
   threads = []
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = event_writer, args=(tid, fin, fout))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 - tread + 1e-30
   print("   ", fout.recordsWritten, "events written out in", tlapse, "seconds,",
         fout.recordsWritten / tlapse, "ev/s,", fout.bytesWritten / tlapse / 1e6, "MB/s")
   del fout

   print("test 13: reading from multi-threaded output file, bz2 compression")

   fin = hddm_r.istream(".." + rest_bz2)
   postable = {}
   threads = []
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = event_reader, args=(tid, fin))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")
   postable_dotdot_bz2 = postable

   print("test 14: reading from multi-threaded output file, zlib compression")

   fin = hddm_r.istream(".." + rest_z)
   postable = {}
   threads = []
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = event_reader, args=(tid, fin))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")
   postable_dotdot_z = postable

   print("test 15: reading from multi-threaded output file, no compression")

   fin = hddm_r.istream(".." + rest_none)
   postable = {}
   threads = []
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = event_reader, args=(tid, fin))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")
   postable_dotdot_none = postable

   print("test 16: single-threaded random access reads, bz2 compression")

   fin = hddm_r.istream(rest_bz2)
   t0 = time.time()
   for evno in sorted(postable_bz2):
      fin.position = postable_bz2[evno]
      rec = fin.read()
      if evno != rec.getReconstructedPhysicsEvent().eventNo:
         print("bad event lookup for event", evno, "in file", rest_bz2)
         print("requested", evno, "but got", rec.getReconstructedPhysicsEvent().eventNo)
         sys.exit(1)
      elif fin.position != postable_bz2[evno]:
         print("bad record position for event", evno, "in file", rest_bz2)
         print("requested", postable_bz2[evno], "but got", fin.position)
         sys.exit(1)
      elif evno > 0 and evno % 1000 == 0:
         print("   event", evno)
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")

   print("test 17: single-threaded random access reads, z compression")

   fin = hddm_r.istream(rest_z)
   t0 = time.time()
   for evno in sorted(postable_z):
      fin.position = postable_z[evno]
      rec = fin.read()
      if evno != rec.getReconstructedPhysicsEvent().eventNo:
         print("bad event lookup for event", evno, "in file", rest_z)
         print("requested", evno, "but got", rec.getReconstructedPhysicsEvent().eventNo)
         sys.exit(1)
      elif fin.position != postable_z[evno]:
         print("bad record position for event", evno, "in file", rest_z)
         print("requested", postable_z[evno], "but got", fin.position)
         sys.exit(1)
      elif evno > 0 and evno % 1000 == 0:
         print("   event", evno)
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")

   print("test 18: single-threaded random access reads, no compression")

   fin = hddm_r.istream(rest_none)
   t0 = time.time()
   for evno in sorted(postable_none):
      fin.position = postable_none[evno]
      rec = fin.read()
      if evno != rec.getReconstructedPhysicsEvent().eventNo:
         print("bad event lookup for event", evno, "in file", rest_none)
         print("requested", evno, "but got", rec.getReconstructedPhysicsEvent().eventNo)
         sys.exit(1)
      elif fin.position != postable_none[evno]:
         print("bad record position for event", evno, "in file", rest_none)
         print("requested", postable_none[evno], "but got", fin.position)
         sys.exit(1)
      elif evno > 0 and evno % 1000 == 0:
         print("   event", evno)
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")

   print("test 19: multi-threaded random access reads, bz2 compression")

   fin = hddm_r.istream(rest_bz2)
   post = postable_bz2
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = random_reader, args=(tid, fin, post, 0))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")

   print("test 20: multi-threaded random access reads, z compression")

   fin = hddm_r.istream(rest_z)
   post = postable_z
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = random_reader, args=(tid, fin, post, 0))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")

   print("test 21: multi-threaded random access reads, no compression")

   fin = hddm_r.istream(rest_none)
   post = postable_none
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = random_reader, args=(tid, fin, post, 0))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")

   print("test 22: multi-threaded random access reads, mt-bz2 compression")

   fin = hddm_r.istream(".." + rest_bz2)
   post = postable_dotdot_bz2
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = random_reader, args=(tid + 10, fin, post, 0))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")

   print("test 23: multi-threaded random access reads, mt-z compression")

   fin = hddm_r.istream(".." + rest_z)
   post = postable_dotdot_z
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = random_reader, args=(tid + 10, fin, post, 0))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")

   print("test 24: multi-threaded random access reads, mt-no compression")

   fin = hddm_r.istream(".." + rest_none)
   post = postable_dotdot_none
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = random_reader, args=(tid + 10, fin, post, 0))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")

   print("test 13 repeated: reading from multi-threaded output file, bz2 compression")

   fin = hddm_r.istream(".." + rest_bz2)
   postable = {}
   threads = []
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = event_reader, args=(tid, fin))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")
   postable_dotdot_bz2 = postable

   print("test 14 repeated: reading from multi-threaded output file, zlib compression")

   fin = hddm_r.istream(".." + rest_z)
   postable = {}
   threads = []
   t0 = time.time()
   for tid in range(0, multithreads):
      thread = threading.Thread(target = event_reader, args=(tid, fin))
      threads.append(thread)
      thread.start()
   for thread in threads:
      thread.join()
   tlapse = time.time() - t0 + 1e-30
   print("   ", fin.recordsRead, "events read in", tlapse, "seconds,",
         fin.recordsRead / tlapse, "ev/s,", fin.bytesRead / tlapse / 1e6, "MB/s")
   postable_dotdot_z = postable

