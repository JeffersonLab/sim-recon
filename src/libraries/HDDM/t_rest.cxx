//
// t_rest - standard suite of tests for hddm i/o performance,
//          with the aim of validating the correct behavior
//          and benchmarking the performance of the library
//          in single/multi-threaded mode, with and without
//          compression, in sequential and random access.
//
// author: richard.t.jones at uconn.edu
// version: august 12, 2016
//
// usage: t_rest [-n <events>] [-p <threads>] <input_restfile.hddm>
//
// notes:
// 1) The input file is read only once, and copies are made in the local
//    directory in different compression modes. There should be enough
//    space in the local directory to hold six copies of the entire
//    input file, or its first N events if option -n is given: one
//    single-threaded output and one multi-threaded output file for
//    each compression mode: bz2, zlib, and none.
//

#include <iostream>
#include <HDDM/hddm_r.hpp>
#include <fstream>
#include <map>

#include <pthread.h>
#include <unistd.h>
#include <time.h>

int max_test = 99;
int maxevents = -1;
int multithreads = 4;
float tread = 0;

pthread_mutex_t lock;

std::map<int, hddm_r::streamposition> postable;

void *event_reader(void *arg) {
   int tid = *((int**)arg)[0];
   hddm_r::istream &fin = *((hddm_r::istream**)arg)[1];
   hddm_r::HDDM rec;
   while (fin >> rec) {
      int evno = rec.getReconstructedPhysicsEvent().getEventNo();
      if (evno > 0 && evno % 10000 == 0)
         printf("   event %d %d\n", evno, tid);
      pthread_mutex_lock(&lock);
      postable[evno] = fin.getPosition();
      pthread_mutex_unlock(&lock);
   }
   return 0;
}

void *event_writer(void *arg) {
   int tid = *((int**)arg)[0];
   hddm_r::istream &fin = *((hddm_r::istream**)arg)[1];
   hddm_r::ostream &fout = *((hddm_r::ostream**)arg)[2];
   hddm_r::HDDM rec;
   while (fin >> rec) {
      int evno = rec.getReconstructedPhysicsEvent().getEventNo();
      if (evno > 0 && evno % 10000 == 0)
         printf("   event %d %d\n", evno, tid);
      fout << rec;
   }
   return 0;
}

void *random_reader(void *arg) {
   int tid = *((int**)arg)[0];
   hddm_r::istream &fin = *((hddm_r::istream**)arg)[1];
   std::map<int, hddm_r::streamposition> &post =
    *((std::map<int, hddm_r::streamposition>**)arg)[2];
   std::map<int, hddm_r::streamposition>::iterator iter;
   for (iter = post.begin(); iter != post.end(); ++iter) {
      fin.setPosition(iter->second);
      hddm_r::HDDM rec;
      fin >> rec;
      if (iter->first != rec.getReconstructedPhysicsEvent().getEventNo()) {
         pthread_mutex_lock(&lock);
         printf("bad event lookup for event %d in thread %d: ", iter->first, tid);
         printf("requested %d but got %ld\n", iter->first,
                rec.getReconstructedPhysicsEvent().getEventNo());
         printf("requested (%ld,%d,%d) and got (%ld,%d,%d)\n", 
                post[iter->first].block_start,
                post[iter->first].block_offset,
                post[iter->first].block_status,
                fin.getPosition().block_start,
                fin.getPosition().block_offset,
                fin.getPosition().block_status);
         pthread_mutex_unlock(&lock);
         break;
      }
      else if (fin.getPosition() != post[iter->first]) {
         pthread_mutex_lock(&lock);
         printf("bad record position for event %d: ", iter->first);
         printf("requested (%ld,%d,%d) but got (%ld,%d,%d)\n", 
                post[iter->first].block_start,
                post[iter->first].block_offset,
                post[iter->first].block_status,
                fin.getPosition().block_start,
                fin.getPosition().block_offset,
                fin.getPosition().block_status);
         pthread_mutex_unlock(&lock);
         break;
      }
      else if (iter->first > 0 && iter->first % 10000 == 0) {
         printf("   event %d\n", iter->first);
      }
   }
   return 0;
}

void usage()
{
   printf("usage: t_rest [options] <input_hddm_file>\n"
          " where options may include any of the following:\n"
          "   -n <count> : stop after reading <count> events from the\n"
          "                input hddm file, default is the entire file\n"
          "   -p <count> : use <count> i/o threads in multi-thread tests,\n"
          "                default is 4\n"
          "   -h, --help : display this help message\n"
          "\n");
   exit(0);
}

int main(int argc, const char *argv[])
{
   int narg;
   for (narg=1; narg < argc; ++narg) {
      if (strncmp(argv[narg], "-n", 2) == 0) {
         maxevents = atoi(argv[++narg]);
      }
      else if (strncmp(argv[narg], "--Nevents", 9) == 0) {
         maxevents = atoi(argv[++narg]);
      }
      else if (strncmp(argv[narg], "-p", 2) == 0) {
         multithreads = atoi(argv[++narg]);
      }
      else if (strncmp(argv[narg], "-Nthreads", 9) == 0) {
         multithreads = atoi(argv[++narg]);
      }
      else if (strncmp(argv[narg], "-h", 2) == 0) {
         usage();
      }
      else if (strncmp(argv[narg], "-?", 2) == 0) {
         usage();
      }
      else if (strncmp(argv[narg], "--help", 6) == 0) {
         usage();
      }
      else {
         break;
      }
   }
   if (narg == argc) {
      usage();
   }

   std::string rest_z = "t_rest_z.hddm";
   std::string rest_bz2 = "t_rest_bz2.hddm";
   std::string rest_none = "t_rest_none.hddm";

   std::map<int, hddm_r::streamposition> postable_bz2;
   std::map<int, hddm_r::streamposition> postable_z;
   std::map<int, hddm_r::streamposition> postable_none;

   pthread_mutex_init(&lock, 0);

   printf("test 0: copying records from input file\n");

   if (max_test > 0) {
      std::ifstream ifs(argv[narg]);
      hddm_r::istream fin(ifs);
      std::ofstream ofs0(rest_none.c_str());
      std::ofstream ofs1(rest_z.c_str());
      std::ofstream ofs2(rest_bz2.c_str());
      hddm_r::ostream fout0(ofs0);
      hddm_r::ostream fout1(ofs1);
      fout1.setCompression(hddm_r::k_z_compression);
      hddm_r::ostream fout2(ofs2);
      fout2.setCompression(hddm_r::k_bz2_compression);
      hddm_r::HDDM rec;
      while (fin >> rec) {
         int evno = rec.getReconstructedPhysicsEvent().getEventNo();
         if (evno > 0 && evno % 10000 == 0)
            printf("   event %d\n", evno);
         fout0 << rec;
         fout1 << rec;
         fout2 << rec;
         if (fin.getRecordsRead() == maxevents)
            break;
      }
   }

   printf("test 1: reading single-threaded from input file, bz2 compression\n");

   if (max_test > 1) {
      std::ifstream ifs(rest_bz2);
      hddm_r::istream fin(ifs);
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      hddm_r::HDDM rec;
      while (fin >> rec) {
         int evno = rec.getReconstructedPhysicsEvent().getEventNo();
         if (evno > 0 && evno % 10000 == 0)
            printf("   event %d\n", evno);
         postable_bz2[evno] = fin.getPosition();
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }

   printf("test 2: reading single-threaded from input file, zlib compression\n");

   if (max_test > 2) {
      std::ifstream ifs(rest_z);
      hddm_r::istream fin(ifs);
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      hddm_r::HDDM rec;
      while (fin >> rec) {
         int evno = rec.getReconstructedPhysicsEvent().getEventNo();
         if (evno > 0 && evno % 10000 == 0)
            printf("   event %d\n", evno);
         postable_z[evno] = fin.getPosition();
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }

   printf("test 3: reading single-threaded from input file, no compression\n");

   if (max_test > 3) {
      std::ifstream ifs(rest_none);
      hddm_r::istream fin(ifs);
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      hddm_r::HDDM rec;
      while (fin >> rec) {
         int evno = rec.getReconstructedPhysicsEvent().getEventNo();
         if (evno > 0 && evno % 10000 == 0)
            printf("   event %d\n", evno);
         postable_none[evno] = fin.getPosition();
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
      tread = tlapse;
   }

   printf("test 4: writing single-threaded to output file, bz2 compression\n");

   if (max_test > 4) {
      std::ifstream ifs(rest_none);
      hddm_r::istream fin(ifs);
      std::ofstream ofs("." + rest_bz2);
      hddm_r::ostream fout(ofs);
      fout.setCompression(hddm_r::k_bz2_compression);
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      hddm_r::HDDM rec;
      while (fin >> rec) {
         int evno = rec.getReconstructedPhysicsEvent().getEventNo();
         if (evno > 0 && evno % 10000 == 0)
            printf("   event %d\n", evno);
         fout << rec;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      tlapse -= tread;
      printf("   %d events written in %f seconds,", fout.getRecordsWritten(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fout.getRecordsWritten() / tlapse, 
             fout.getBytesWritten() / tlapse / 1e6);
   }

   printf("test 5: writing single-threaded to output file, zlib compression\n");

   if (max_test > 5) {
      std::ifstream ifs(rest_none);
      hddm_r::istream fin(ifs);
      std::ofstream ofs("." + rest_z);
      hddm_r::ostream fout(ofs);
      fout.setCompression(hddm_r::k_z_compression);
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      hddm_r::HDDM rec;
      while (fin >> rec) {
         int evno = rec.getReconstructedPhysicsEvent().getEventNo();
         if (evno > 0 && evno % 10000 == 0)
            printf("   event %d\n", evno);
         fout << rec;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      tlapse -= tread;
      printf("   %d events written in %f seconds,", fout.getRecordsWritten(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fout.getRecordsWritten() / tlapse, 
             fout.getBytesWritten() / tlapse / 1e6);
   }

   printf("test 6: writing single-threaded to output file, no compression\n");

   if (max_test > 6) {
      std::ifstream ifs(rest_none);
      hddm_r::istream fin(ifs);
      std::ofstream ofs("." + rest_none);
      hddm_r::ostream fout(ofs);
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      hddm_r::HDDM rec;
      while (fin >> rec) {
         int evno = rec.getReconstructedPhysicsEvent().getEventNo();
         if (evno > 0 && evno % 10000 == 0)
            printf("   event %d\n", evno);
         fout << rec;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      tlapse -= tread;
      printf("   %d events written in %f seconds,", fout.getRecordsWritten(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fout.getRecordsWritten() / tlapse, 
             fout.getBytesWritten() / tlapse / 1e6);
   }

   printf("test 7: reading multi-threaded from input file, bz2 compression\n");

   if (max_test > 7) {
      std::ifstream ifs("." + rest_bz2);
      hddm_r::istream fin(ifs);
      postable.clear();
      std::list<pthread_t*> threads;
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      for (int tid=0; tid < multithreads; ++tid) {
         pthread_t *thread = new pthread_t;
         void *args[2] = {&tid, &fin};
         pthread_create(thread, 0, event_reader, args);
         threads.push_back(thread);
      }
      std::list<pthread_t*>::iterator iter;
      for (iter = threads.begin(); iter != threads.end(); ++iter) {
         void *retval;
         pthread_join(**iter, &retval);
         delete *iter;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }
   std::map<int, hddm_r::streamposition> postable_dot_bz2(postable);

   printf("test 8: reading multi-threaded from input file, zlib compression\n");

   if (max_test > 8) {
      std::ifstream ifs("." + rest_z);
      hddm_r::istream fin(ifs);
      postable.clear();
      std::list<pthread_t*> threads;
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      for (int tid=0; tid < multithreads; ++tid) {
         pthread_t *thread = new pthread_t;
         void *args[2] = {&tid, &fin};
         pthread_create(thread, 0, event_reader, args);
         threads.push_back(thread);
      }
      std::list<pthread_t*>::iterator iter;
      for (iter = threads.begin(); iter != threads.end(); ++iter) {
         void *retval;
         pthread_join(**iter, &retval);
         delete *iter;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }
   std::map<int, hddm_r::streamposition> postable_dot_z(postable);

   printf("test 9: reading multi-threaded from input file, no compression\n");

   if (max_test > 9) {
      std::ifstream ifs("." + rest_none);
      hddm_r::istream fin(ifs);
      postable.clear();
      std::list<pthread_t*> threads;
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      for (int tid=0; tid < multithreads; ++tid) {
         pthread_t *thread = new pthread_t;
         void *args[2] = {&tid, &fin};
         pthread_create(thread, 0, event_reader, args);
         threads.push_back(thread);
      }
      std::list<pthread_t*>::iterator iter;
      for (iter = threads.begin(); iter != threads.end(); ++iter) {
         void *retval;
         pthread_join(**iter, &retval);
         delete *iter;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }
   std::map<int, hddm_r::streamposition> postable_dot_none(postable);

   printf("test 10: writing multi-threaded to output file, bz2 compression\n");

   if (max_test > 10) {
      std::ifstream ifs(rest_none);
      hddm_r::istream fin(ifs);
      std::ofstream ofs(".." + rest_bz2);
      hddm_r::ostream fout(ofs);
      fout.setCompression(hddm_r::k_bz2_compression);
      std::list<pthread_t*> threads;
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      for (int tid=0; tid < multithreads; ++tid) {
         pthread_t *thread = new pthread_t;
         void *args[3] = {&tid, &fin, &fout};
         pthread_create(thread, 0, event_writer, args);
         threads.push_back(thread);
      }
      std::list<pthread_t*>::iterator iter;
      for (iter = threads.begin(); iter != threads.end(); ++iter) {
         void *retval;
         pthread_join(**iter, &retval);
         delete *iter;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      tlapse -= tread;
      printf("   %d events written in %f seconds,", fout.getRecordsWritten(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fout.getRecordsWritten() / tlapse, 
             fout.getBytesWritten() / tlapse / 1e6);
   }

   printf("test 11: writing multi-threaded to output file, zlib compression\n");

   if (max_test > 11) {
      std::ifstream ifs(rest_none);
      hddm_r::istream fin(ifs);
      std::ofstream ofs(".." + rest_z);
      hddm_r::ostream fout(ofs);
      fout.setCompression(hddm_r::k_z_compression);
      std::list<pthread_t*> threads;
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      for (int tid=0; tid < multithreads; ++tid) {
         pthread_t *thread = new pthread_t;
         void *args[3] = {&tid, &fin, &fout};
         pthread_create(thread, 0, event_writer, args);
         threads.push_back(thread);
      }
      std::list<pthread_t*>::iterator iter;
      for (iter = threads.begin(); iter != threads.end(); ++iter) {
         void *retval;
         pthread_join(**iter, &retval);
         delete *iter;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      tlapse -= tread;
      printf("   %d events written in %f seconds,", fout.getRecordsWritten(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fout.getRecordsWritten() / tlapse, 
             fout.getBytesWritten() / tlapse / 1e6);
   }

   printf("test 12: writing multi-threaded to output file, no compression\n");

   if (max_test > 12) {
      std::ifstream ifs(rest_none);
      hddm_r::istream fin(ifs);
      std::ofstream ofs(".." + rest_none);
      hddm_r::ostream fout(ofs);
      std::list<pthread_t*> threads;
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      for (int tid=0; tid < multithreads; ++tid) {
         pthread_t *thread = new pthread_t;
         void *args[3] = {&tid, &fin, &fout};
         pthread_create(thread, 0, event_writer, args);
         threads.push_back(thread);
      }
      std::list<pthread_t*>::iterator iter;
      for (iter = threads.begin(); iter != threads.end(); ++iter) {
         void *retval;
         pthread_join(**iter, &retval);
         delete *iter;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      tlapse -= tread;
      printf("   %d events written in %f seconds,", fout.getRecordsWritten(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fout.getRecordsWritten() / tlapse, 
             fout.getBytesWritten() / tlapse / 1e6);
   }

   printf("test 13: reading from multi-threaded output file, bz2 compression\n");

   if (max_test > 13) {
      std::ifstream ifs(".." + rest_bz2);
      hddm_r::istream fin(ifs);
      postable.clear();
      std::list<pthread_t*> threads;
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      for (int tid=0; tid < multithreads; ++tid) {
         pthread_t *thread = new pthread_t;
         void *args[2] = {&tid, &fin};
         pthread_create(thread, 0, event_reader, args);
         threads.push_back(thread);
      }
      std::list<pthread_t*>::iterator iter;
      for (iter = threads.begin(); iter != threads.end(); ++iter) {
         void *retval;
         pthread_join(**iter, &retval);
         delete *iter;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }
   std::map<int, hddm_r::streamposition> postable_dotdot_bz2(postable);

   printf("test 14: reading from multi-threaded output file, zlib compression\n");

   if (max_test > 14) {
      std::ifstream ifs(".." + rest_z);
      hddm_r::istream fin(ifs);
      postable.clear();
      std::list<pthread_t*> threads;
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      for (int tid=0; tid < multithreads; ++tid) {
         pthread_t *thread = new pthread_t;
         void *args[2] = {&tid, &fin};
         pthread_create(thread, 0, event_reader, args);
         threads.push_back(thread);
      }
      std::list<pthread_t*>::iterator iter;
      for (iter = threads.begin(); iter != threads.end(); ++iter) {
         void *retval;
         pthread_join(**iter, &retval);
         delete *iter;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }
   std::map<int, hddm_r::streamposition> postable_dotdot_z(postable);

   printf("test 15: reading from multi-threaded output file, no compression\n");

   if (max_test > 15) {
      std::ifstream ifs(".." + rest_none);
      hddm_r::istream fin(ifs);
      postable.clear();
      std::list<pthread_t*> threads;
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      for (int tid=0; tid < multithreads; ++tid) {
         pthread_t *thread = new pthread_t;
         void *args[2] = {&tid, &fin};
         pthread_create(thread, 0, event_reader, args);
         threads.push_back(thread);
      }
      std::list<pthread_t*>::iterator iter;
      for (iter = threads.begin(); iter != threads.end(); ++iter) {
         void *retval;
         pthread_join(**iter, &retval);
         delete *iter;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }
   std::map<int, hddm_r::streamposition> postable_dotdot_none(postable);

   printf("test 16: single-threaded random access reads, bz2 compression\n");

   if (max_test > 16) {
      std::ifstream ifs(rest_bz2);
      hddm_r::istream fin(ifs);
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      std::map<int, hddm_r::streamposition>::iterator iter;
      for (iter = postable_bz2.begin(); iter != postable_bz2.end(); ++iter) {
         fin.setPosition(iter->second);
         hddm_r::HDDM rec;
         fin >> rec;
         if (iter->first != rec.getReconstructedPhysicsEvent().getEventNo()) {
            printf("bad event lookup for event %d in file %s: ", 
                   iter->first, rest_bz2.c_str());
            printf("requested %d but got %ld\n", iter->first, 
                   rec.getReconstructedPhysicsEvent().getEventNo());
            printf("requested (%ld,%d,%d) and got (%ld,%d,%d)\n",
                   iter->second.block_start,
                   iter->second.block_offset,
                   iter->second.block_status,
                   fin.getPosition().block_start,
                   fin.getPosition().block_offset,
                   fin.getPosition().block_status);
            exit(1);
         }
         else if (fin.getPosition() != iter->second) {
            printf("bad record position for event %d in file %s: ",
                   iter->first, rest_bz2.c_str());
            printf("requested (%ld,%d,%d) but got (%ld,%d,%d)\n",
                   iter->second.block_start,
                   iter->second.block_offset,
                   iter->second.block_status,
                   fin.getPosition().block_start,
                   fin.getPosition().block_offset,
                   fin.getPosition().block_status);
            exit(1);
         }
         else if (iter->first > 0 && iter->first % 1000 == 0 ) {
            printf("   event %d\n", iter->first);
         }
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }

   printf("test 17: single-threaded random access reads, z compression\n");

   if (max_test > 17) {
      std::ifstream ifs(rest_z);
      hddm_r::istream fin(ifs);
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      std::map<int, hddm_r::streamposition>::iterator iter;
      for (iter = postable_z.begin(); iter != postable_z.end(); ++iter) {
         fin.setPosition(iter->second);
         hddm_r::HDDM rec;
         fin >> rec;
         if (iter->first != rec.getReconstructedPhysicsEvent().getEventNo()) {
            printf("bad event lookup for event %d in file %s: ", 
                   iter->first, rest_z.c_str());
            printf("requested %d but got %ld\n", iter->first, 
                   rec.getReconstructedPhysicsEvent().getEventNo());
            printf("requested (%ld,%d,%d) and got (%ld,%d,%d)\n",
                   iter->second.block_start,
                   iter->second.block_offset,
                   iter->second.block_status,
                   fin.getPosition().block_start,
                   fin.getPosition().block_offset,
                   fin.getPosition().block_status);
            exit(1);
         }
         else if (fin.getPosition() != iter->second) {
            printf("bad record position for event %d in file %s: ",
                   iter->first, rest_z.c_str());
            printf("requested (%ld,%d,%d) but got (%ld,%d,%d)\n",
                   iter->second.block_start,
                   iter->second.block_offset,
                   iter->second.block_status,
                   fin.getPosition().block_start,
                   fin.getPosition().block_offset,
                   fin.getPosition().block_status);
            exit(1);
         }
         else if (iter->first > 0 && iter->first % 1000 == 0 ) {
            printf("   event %d\n", iter->first);
         }
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }

   printf("test 18: single-threaded random access reads, no compression\n");

   if (max_test > 18) {
      std::ifstream ifs(rest_none);
      hddm_r::istream fin(ifs);
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      std::map<int, hddm_r::streamposition>::iterator iter;
      for (iter = postable_none.begin(); iter != postable_none.end(); ++iter) {
         fin.setPosition(iter->second);
         hddm_r::HDDM rec;
         fin >> rec;
         if (iter->first != rec.getReconstructedPhysicsEvent().getEventNo()) {
            printf("bad event lookup for event %d in file %s: ", 
                   iter->first, rest_none.c_str());
            printf("requested %d but got %ld\n", iter->first, 
                   rec.getReconstructedPhysicsEvent().getEventNo());
            printf("requested (%ld,%d,%d) and got (%ld,%d,%d)\n",
                   iter->second.block_start,
                   iter->second.block_offset,
                   iter->second.block_status,
                   fin.getPosition().block_start,
                   fin.getPosition().block_offset,
                   fin.getPosition().block_status);
            exit(1);
         }
         else if (fin.getPosition() != iter->second) {
            printf("bad record position for event %d in file %s: ",
                   iter->first, rest_none.c_str());
            printf("requested (%ld,%d,%d) but got (%ld,%d,%d)\n",
                   iter->second.block_start,
                   iter->second.block_offset,
                   iter->second.block_status,
                   fin.getPosition().block_start,
                   fin.getPosition().block_offset,
                   fin.getPosition().block_status);
            exit(1);
         }
         else if (iter->first > 0 && iter->first % 1000 == 0 ) {
            printf("   event %d\n", iter->first);
         }
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }

   printf("test 19: multi-threaded random access reads, bz2 compression\n");

   if (max_test > 19) {
      std::ifstream ifs(rest_bz2);
      hddm_r::istream fin(ifs);
      std::list<pthread_t*> threads;
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      for (int tid=0; tid < multithreads; ++tid) {
         pthread_t *thread = new pthread_t;
         void *args[3] = {&tid, &fin, &postable_bz2};
         pthread_create(thread, 0, random_reader, args);
         threads.push_back(thread);
      }
      std::list<pthread_t*>::iterator iter;
      for (iter = threads.begin(); iter != threads.end(); ++iter) {
         void *retval;
         pthread_join(**iter, &retval);
         delete *iter;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }

   printf("test 20: multi-threaded random access reads, z compression\n");

   if (max_test > 20) {
      std::ifstream ifs(rest_z);
      hddm_r::istream fin(ifs);
      std::list<pthread_t*> threads;
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      for (int tid=0; tid < multithreads; ++tid) {
         pthread_t *thread = new pthread_t;
         void *args[3] = {&tid, &fin, &postable_z};
         pthread_create(thread, 0, random_reader, args);
         threads.push_back(thread);
      }
      std::list<pthread_t*>::iterator iter;
      for (iter = threads.begin(); iter != threads.end(); ++iter) {
         void *retval;
         pthread_join(**iter, &retval);
         delete *iter;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }

   printf("test 21: multi-threaded random access reads, no compression\n");

   if (max_test > 21) {
      std::ifstream ifs(rest_none);
      hddm_r::istream fin(ifs);
      std::list<pthread_t*> threads;
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      for (int tid=0; tid < multithreads; ++tid) {
         pthread_t *thread = new pthread_t;
         void *args[3] = {&tid, &fin, &postable_none};
         pthread_create(thread, 0, random_reader, args);
         threads.push_back(thread);
      }
      std::list<pthread_t*>::iterator iter;
      for (iter = threads.begin(); iter != threads.end(); ++iter) {
         void *retval;
         pthread_join(**iter, &retval);
         delete *iter;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }

   printf("test 22: multi-threaded random access reads, mt-bz2 compression\n");

   if (max_test > 22) {
      std::ifstream ifs(".." + rest_bz2);
      hddm_r::istream fin(ifs);
      std::list<pthread_t*> threads;
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      for (int tid=0; tid < multithreads; ++tid) {
         pthread_t *thread = new pthread_t;
         void *args[3] = {&tid, &fin, &postable_dotdot_bz2};
         pthread_create(thread, 0, random_reader, args);
         threads.push_back(thread);
      }
      std::list<pthread_t*>::iterator iter;
      for (iter = threads.begin(); iter != threads.end(); ++iter) {
         void *retval;
         pthread_join(**iter, &retval);
         delete *iter;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }

   printf("test 23: multi-threaded random access reads, mt-z compression\n");

   if (max_test > 23) {
      std::ifstream ifs(".." + rest_z);
      hddm_r::istream fin(ifs);
      std::list<pthread_t*> threads;
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      for (int tid=0; tid < multithreads; ++tid) {
         pthread_t *thread = new pthread_t;
         void *args[3] = {&tid, &fin, &postable_dotdot_z};
         pthread_create(thread, 0, random_reader, args);
         threads.push_back(thread);
      }
      std::list<pthread_t*>::iterator iter;
      for (iter = threads.begin(); iter != threads.end(); ++iter) {
         void *retval;
         pthread_join(**iter, &retval);
         delete *iter;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }

   printf("test 24: multi-threaded random access reads, mt-no compression\n");

   if (max_test > 24) {
      std::ifstream ifs(".." + rest_none);
      hddm_r::istream fin(ifs);
      std::list<pthread_t*> threads;
      struct timespec t0;
      clock_gettime(CLOCK_REALTIME, &t0);
      for (int tid=0; tid < multithreads; ++tid) {
         pthread_t *thread = new pthread_t;
         void *args[3] = {&tid, &fin, &postable_dotdot_none};
         pthread_create(thread, 0, random_reader, args);
         threads.push_back(thread);
      }
      std::list<pthread_t*>::iterator iter;
      for (iter = threads.begin(); iter != threads.end(); ++iter) {
         void *retval;
         pthread_join(**iter, &retval);
         delete *iter;
      }
      struct timespec t1;
      clock_gettime(CLOCK_REALTIME, &t1);
      float tlapse = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)*1e-9;
      printf("   %d events read in %f seconds,", fin.getRecordsRead(), tlapse);
      printf(" %f ev/s, %f MB/s\n",  fin.getRecordsRead() / tlapse, 
             fin.getBytesRead() / tlapse / 1e6);
   }
}
