# 15 "timel.c"
# 1 "/cern/pro/include/kerngen/pilot.h" 1
# 16 "timel.c" 2
# 27 "timel.c"
# 1 "/usr/include/sys/types.h" 1 3
# 26 "/usr/include/sys/types.h" 3
# 1 "/usr/include/features.h" 1 3
# 283 "/usr/include/features.h" 3
# 1 "/usr/include/sys/cdefs.h" 1 3
# 284 "/usr/include/features.h" 2 3
# 311 "/usr/include/features.h" 3
# 1 "/usr/include/gnu/stubs.h" 1 3
# 312 "/usr/include/features.h" 2 3
# 27 "/usr/include/sys/types.h" 2 3



# 1 "/usr/include/bits/types.h" 1 3
# 26 "/usr/include/bits/types.h" 3
# 1 "/usr/include/features.h" 1 3
# 27 "/usr/include/bits/types.h" 2 3


# 1 "/usr/lib/gcc-lib/i386-redhat-linux/2.96/include/stddef.h" 1 3
# 199 "/usr/lib/gcc-lib/i386-redhat-linux/2.96/include/stddef.h" 3
typedef unsigned int size_t;
# 30 "/usr/include/bits/types.h" 2 3


typedef unsigned char __u_char;
typedef unsigned short __u_short;
typedef unsigned int __u_int;
typedef unsigned long __u_long;

__extension__ typedef unsigned long long int __u_quad_t;
__extension__ typedef long long int __quad_t;
# 49 "/usr/include/bits/types.h" 3
typedef signed char __int8_t;
typedef unsigned char __uint8_t;
typedef signed short int __int16_t;
typedef unsigned short int __uint16_t;
typedef signed int __int32_t;
typedef unsigned int __uint32_t;

__extension__ typedef signed long long int __int64_t;
__extension__ typedef unsigned long long int __uint64_t;

typedef __quad_t *__qaddr_t;

typedef __u_quad_t __dev_t;
typedef __u_int __uid_t;
typedef __u_int __gid_t;
typedef __u_long __ino_t;
typedef __u_int __mode_t;
typedef __u_int __nlink_t;
typedef long int __off_t;
typedef __quad_t __loff_t;
typedef int __pid_t;
typedef int __ssize_t;
typedef __u_long __rlim_t;
typedef __u_quad_t __rlim64_t;
typedef __u_int __id_t;

typedef struct
  {
    int __val[2];
  } __fsid_t;


typedef int __daddr_t;
typedef char *__caddr_t;
typedef long int __time_t;
typedef unsigned int __useconds_t;
typedef long int __suseconds_t;
typedef long int __swblk_t;

typedef long int __clock_t;


typedef int __clockid_t;


typedef int __timer_t;






typedef int __key_t;


typedef unsigned short int __ipc_pid_t;



typedef long int __blksize_t;




typedef long int __blkcnt_t;
typedef __quad_t __blkcnt64_t;


typedef __u_long __fsblkcnt_t;
typedef __u_quad_t __fsblkcnt64_t;


typedef __u_long __fsfilcnt_t;
typedef __u_quad_t __fsfilcnt64_t;


typedef __u_quad_t __ino64_t;


typedef __loff_t __off64_t;


typedef long int __t_scalar_t;
typedef unsigned long int __t_uscalar_t;


typedef int __intptr_t;


typedef unsigned int __socklen_t;




# 1 "/usr/include/bits/pthreadtypes.h" 1 3
# 23 "/usr/include/bits/pthreadtypes.h" 3
# 1 "/usr/include/bits/sched.h" 1 3
# 68 "/usr/include/bits/sched.h" 3
struct __sched_param
  {
    int __sched_priority;
  };
# 24 "/usr/include/bits/pthreadtypes.h" 2 3


struct _pthread_fastlock
{
  long int __status;
  int __spinlock;

};



typedef struct _pthread_descr_struct *_pthread_descr;





typedef struct __pthread_attr_s
{
  int __detachstate;
  int __schedpolicy;
  struct __sched_param __schedparam;
  int __inheritsched;
  int __scope;
  size_t __guardsize;
  int __stackaddr_set;
  void *__stackaddr;
  size_t __stacksize;
} pthread_attr_t;



typedef struct
{
  struct _pthread_fastlock __c_lock;
  _pthread_descr __c_waiting;
} pthread_cond_t;



typedef struct
{
  int __dummy;
} pthread_condattr_t;


typedef unsigned int pthread_key_t;





typedef struct
{
  int __m_reserved;
  int __m_count;
  _pthread_descr __m_owner;
  int __m_kind;
  struct _pthread_fastlock __m_lock;
} pthread_mutex_t;



typedef struct
{
  int __mutexkind;
} pthread_mutexattr_t;



typedef int pthread_once_t;
# 140 "/usr/include/bits/pthreadtypes.h" 3
typedef unsigned long int pthread_t;
# 144 "/usr/include/bits/types.h" 2 3
# 31 "/usr/include/sys/types.h" 2 3


typedef __u_char u_char;
typedef __u_short u_short;
typedef __u_int u_int;
typedef __u_long u_long;
typedef __quad_t quad_t;
typedef __u_quad_t u_quad_t;
typedef __fsid_t fsid_t;


typedef __loff_t loff_t;



typedef __ino_t ino_t;
# 58 "/usr/include/sys/types.h" 3
typedef __dev_t dev_t;




typedef __gid_t gid_t;




typedef __mode_t mode_t;




typedef __nlink_t nlink_t;




typedef __uid_t uid_t;





typedef __off_t off_t;
# 96 "/usr/include/sys/types.h" 3
typedef __pid_t pid_t;




typedef __id_t id_t;




typedef __ssize_t ssize_t;




typedef __daddr_t daddr_t;
typedef __caddr_t caddr_t;



typedef __key_t key_t;
# 126 "/usr/include/sys/types.h" 3
# 1 "/usr/include/time.h" 1 3
# 67 "/usr/include/time.h" 3
# 1 "/usr/include/bits/types.h" 1 3
# 68 "/usr/include/time.h" 2 3


typedef __time_t time_t;
# 79 "/usr/include/time.h" 3
# 1 "/usr/include/bits/types.h" 1 3
# 80 "/usr/include/time.h" 2 3


typedef __clockid_t clockid_t;
# 91 "/usr/include/time.h" 3
# 1 "/usr/include/bits/types.h" 1 3
# 92 "/usr/include/time.h" 2 3


typedef __timer_t timer_t;
# 127 "/usr/include/sys/types.h" 2 3
# 140 "/usr/include/sys/types.h" 3
# 1 "/usr/lib/gcc-lib/i386-redhat-linux/2.96/include/stddef.h" 1 3
# 141 "/usr/include/sys/types.h" 2 3



typedef unsigned long int ulong;
typedef unsigned short int ushort;
typedef unsigned int uint;
# 184 "/usr/include/sys/types.h" 3
typedef int int8_t __attribute__ ((__mode__ (__QI__)));
typedef int int16_t __attribute__ ((__mode__ (__HI__)));
typedef int int32_t __attribute__ ((__mode__ (__SI__)));
typedef int int64_t __attribute__ ((__mode__ (__DI__)));


typedef unsigned int u_int8_t __attribute__ ((__mode__ (__QI__)));
typedef unsigned int u_int16_t __attribute__ ((__mode__ (__HI__)));
typedef unsigned int u_int32_t __attribute__ ((__mode__ (__SI__)));
typedef unsigned int u_int64_t __attribute__ ((__mode__ (__DI__)));

typedef int register_t __attribute__ ((__mode__ (__word__)));
# 206 "/usr/include/sys/types.h" 3
# 1 "/usr/include/endian.h" 1 3
# 22 "/usr/include/endian.h" 3
# 1 "/usr/include/features.h" 1 3
# 23 "/usr/include/endian.h" 2 3
# 37 "/usr/include/endian.h" 3
# 1 "/usr/include/bits/endian.h" 1 3
# 38 "/usr/include/endian.h" 2 3
# 207 "/usr/include/sys/types.h" 2 3


# 1 "/usr/include/sys/select.h" 1 3
# 25 "/usr/include/sys/select.h" 3
# 1 "/usr/include/features.h" 1 3
# 26 "/usr/include/sys/select.h" 2 3


# 1 "/usr/include/bits/types.h" 1 3
# 29 "/usr/include/sys/select.h" 2 3


# 1 "/usr/include/bits/select.h" 1 3
# 32 "/usr/include/sys/select.h" 2 3


# 1 "/usr/include/bits/sigset.h" 1 3
# 23 "/usr/include/bits/sigset.h" 3
typedef int __sig_atomic_t;




typedef struct
  {
    unsigned long int __val[(1024 / (8 * sizeof (unsigned long int)))];
  } __sigset_t;
# 35 "/usr/include/sys/select.h" 2 3



typedef __sigset_t sigset_t;




# 1 "/usr/include/time.h" 1 3
# 106 "/usr/include/time.h" 3
struct timespec
  {
    long int tv_sec;
    long int tv_nsec;
  };
# 44 "/usr/include/sys/select.h" 2 3

# 1 "/usr/include/bits/time.h" 1 3
# 63 "/usr/include/bits/time.h" 3
# 1 "/usr/include/bits/types.h" 1 3
# 64 "/usr/include/bits/time.h" 2 3



struct timeval
  {
    __time_t tv_sec;
    __suseconds_t tv_usec;
  };
# 46 "/usr/include/sys/select.h" 2 3



typedef long int __fd_mask;







typedef struct
  {






    __fd_mask __fds_bits[1024 / (8 * sizeof (__fd_mask))];


  } fd_set;






typedef __fd_mask fd_mask;
# 96 "/usr/include/sys/select.h" 3
extern int select (int __nfds, fd_set *__restrict __readfds,
                   fd_set *__restrict __writefds,
                   fd_set *__restrict __exceptfds,
                   struct timeval *__restrict __timeout) ;
# 210 "/usr/include/sys/types.h" 2 3


# 1 "/usr/include/sys/sysmacros.h" 1 3
# 213 "/usr/include/sys/types.h" 2 3
# 224 "/usr/include/sys/types.h" 3
typedef __blkcnt_t blkcnt_t;



typedef __fsblkcnt_t fsblkcnt_t;



typedef __fsfilcnt_t fsfilcnt_t;
# 28 "timel.c" 2
# 1 "/usr/include/time.h" 1 3
# 28 "/usr/include/time.h" 3
# 1 "/usr/include/features.h" 1 3
# 29 "/usr/include/time.h" 2 3
# 38 "/usr/include/time.h" 3
# 1 "/usr/lib/gcc-lib/i386-redhat-linux/2.96/include/stddef.h" 1 3
# 39 "/usr/include/time.h" 2 3



# 1 "/usr/include/bits/time.h" 1 3
# 39 "/usr/include/bits/time.h" 3
# 1 "/usr/include/bits/types.h" 1 3
# 40 "/usr/include/bits/time.h" 2 3
extern long int __sysconf (int);
# 43 "/usr/include/time.h" 2 3
# 56 "/usr/include/time.h" 3
# 1 "/usr/include/bits/types.h" 1 3
# 57 "/usr/include/time.h" 2 3


typedef __clock_t clock_t;
# 118 "/usr/include/time.h" 3
struct tm
{
  int tm_sec;
  int tm_min;
  int tm_hour;
  int tm_mday;
  int tm_mon;
  int tm_year;
  int tm_wday;
  int tm_yday;
  int tm_isdst;


  long int tm_gmtoff;
  __const char *tm_zone;




};




struct itimerspec
  {
    struct timespec it_interval;
    struct timespec it_value;
  };


struct sigevent;
# 163 "/usr/include/time.h" 3
extern clock_t clock (void) ;


extern time_t time (time_t *__timer) ;


extern double difftime (time_t __time1, time_t __time0)
             __attribute__ ((__const__));


extern time_t mktime (struct tm *__tp) ;





extern size_t strftime (char *__restrict __s, size_t __maxsize,
                        __const char *__restrict __format,
                        __const struct tm *__restrict __tp) ;
# 194 "/usr/include/time.h" 3
extern struct tm *gmtime (__const time_t *__timer) ;



extern struct tm *localtime (__const time_t *__timer) ;




extern struct tm *gmtime_r (__const time_t *__restrict __timer,
                            struct tm *__restrict __tp) ;



extern struct tm *localtime_r (__const time_t *__restrict __timer,
                               struct tm *__restrict __tp) ;




extern char *asctime (__const struct tm *__tp) ;


extern char *ctime (__const time_t *__timer) ;






extern char *asctime_r (__const struct tm *__restrict __tp,
                        char *__restrict __buf) ;


extern char *ctime_r (__const time_t *__restrict __timer,
                      char *__restrict __buf) ;




extern char *__tzname[2];
extern int __daylight;
extern long int __timezone;




extern char *tzname[2];



extern void tzset (void) ;



extern int daylight;
extern long int timezone;





extern int stime (__const time_t *__when) ;
# 271 "/usr/include/time.h" 3
extern time_t timegm (struct tm *__tp) ;


extern time_t timelocal (struct tm *__tp) ;


extern int dysize (int __year) __attribute__ ((__const__));





extern int nanosleep (__const struct timespec *__requested_time,
                      struct timespec *__remaining) ;



extern int clock_getres (clockid_t __clock_id, struct timespec *__res) ;


extern int clock_gettime (clockid_t __clock_id, struct timespec *__tp) ;


extern int clock_settime (clockid_t __clock_id, __const struct timespec *__tp)
            ;
# 309 "/usr/include/time.h" 3
extern int timer_create (clockid_t __clock_id,
                         struct sigevent *__restrict __evp,
                         timer_t *__restrict __timerid) ;


extern int timer_delete (timer_t __timerid) ;


extern int timer_settime (timer_t __timerid, int __flags,
                          __const struct itimerspec *__restrict __value,
                          struct itimerspec *__restrict __ovalue) ;


extern int timer_gettime (timer_t __timerid, struct itimerspec *__value)
            ;


extern int timer_getoverrun (timer_t __timerid) ;
# 29 "timel.c" 2
# 1 "/usr/include/sys/times.h" 1 3
# 26 "/usr/include/sys/times.h" 3
# 1 "/usr/include/features.h" 1 3
# 27 "/usr/include/sys/times.h" 2 3


# 1 "/usr/include/time.h" 1 3
# 30 "/usr/include/sys/times.h" 2 3





struct tms
  {
    clock_t tms_utime;
    clock_t tms_stime;

    clock_t tms_cutime;
    clock_t tms_cstime;
  };






extern clock_t times (struct tms *__buffer) ;
# 30 "timel.c" 2
# 1 "/usr/include/sys/param.h" 1 3
# 22 "/usr/include/sys/param.h" 3
# 1 "/usr/lib/gcc-lib/i386-redhat-linux/2.96/include/limits.h" 1 3
# 11 "/usr/lib/gcc-lib/i386-redhat-linux/2.96/include/limits.h" 3
# 1 "/usr/lib/gcc-lib/i386-redhat-linux/2.96/include/syslimits.h" 1 3






# 1 "/usr/lib/gcc-lib/i386-redhat-linux/2.96/include/limits.h" 1 3
# 130 "/usr/lib/gcc-lib/i386-redhat-linux/2.96/include/limits.h" 3
# 1 "/usr/include/limits.h" 1 3
# 26 "/usr/include/limits.h" 3
# 1 "/usr/include/features.h" 1 3
# 27 "/usr/include/limits.h" 2 3
# 144 "/usr/include/limits.h" 3
# 1 "/usr/include/bits/posix1_lim.h" 1 3
# 126 "/usr/include/bits/posix1_lim.h" 3
# 1 "/usr/include/bits/local_lim.h" 1 3
# 36 "/usr/include/bits/local_lim.h" 3
# 1 "/usr/include/linux/limits.h" 1 3
# 37 "/usr/include/bits/local_lim.h" 2 3
# 127 "/usr/include/bits/posix1_lim.h" 2 3
# 145 "/usr/include/limits.h" 2 3



# 1 "/usr/include/bits/posix2_lim.h" 1 3
# 149 "/usr/include/limits.h" 2 3
# 131 "/usr/lib/gcc-lib/i386-redhat-linux/2.96/include/limits.h" 2 3
# 8 "/usr/lib/gcc-lib/i386-redhat-linux/2.96/include/syslimits.h" 2 3
# 12 "/usr/lib/gcc-lib/i386-redhat-linux/2.96/include/limits.h" 2 3
# 23 "/usr/include/sys/param.h" 2 3
# 1 "/usr/include/linux/limits.h" 1 3
# 24 "/usr/include/sys/param.h" 2 3
# 1 "/usr/include/linux/param.h" 1 3



# 1 "/usr/include/asm/param.h" 1 3
# 5 "/usr/include/linux/param.h" 2 3
# 25 "/usr/include/sys/param.h" 2 3
# 42 "/usr/include/sys/param.h" 3
# 1 "/usr/include/sys/types.h" 1 3
# 43 "/usr/include/sys/param.h" 2 3
# 31 "timel.c" 2
# 1 "/usr/include/sys/time.h" 1 3
# 22 "/usr/include/sys/time.h" 3
# 1 "/usr/include/features.h" 1 3
# 23 "/usr/include/sys/time.h" 2 3

# 1 "/usr/include/bits/types.h" 1 3
# 25 "/usr/include/sys/time.h" 2 3

# 1 "/usr/include/time.h" 1 3
# 27 "/usr/include/sys/time.h" 2 3

# 1 "/usr/include/bits/time.h" 1 3
# 29 "/usr/include/sys/time.h" 2 3

# 1 "/usr/include/sys/select.h" 1 3
# 31 "/usr/include/sys/time.h" 2 3


typedef __suseconds_t suseconds_t;
# 56 "/usr/include/sys/time.h" 3
struct timezone
  {
    int tz_minuteswest;
    int tz_dsttime;
  };

typedef struct timezone *__restrict __timezone_ptr_t;
# 72 "/usr/include/sys/time.h" 3
extern int gettimeofday (struct timeval *__restrict __tv,
                         __timezone_ptr_t __tz) ;




extern int settimeofday (__const struct timeval *__tv,
                         __const struct timezone *__tz) ;





extern int adjtime (__const struct timeval *__delta,
                    struct timeval *__olddelta) ;




enum __itimer_which
  {

    ITIMER_REAL = 0,


    ITIMER_VIRTUAL = 1,



    ITIMER_PROF = 2

  };



struct itimerval
  {

    struct timeval it_interval;

    struct timeval it_value;
  };




typedef int __itimer_which_t;




extern int getitimer (__itimer_which_t __which,
                      struct itimerval *__value) ;




extern int setitimer (__itimer_which_t __which,
                      __const struct itimerval *__restrict __new,
                      struct itimerval *__restrict __old) ;



extern int utimes (__const char *__file, __const struct timeval __tvp[2])
            ;
# 32 "timel.c" 2
# 1 "/usr/include/sys/resource.h" 1 3
# 22 "/usr/include/sys/resource.h" 3
# 1 "/usr/include/features.h" 1 3
# 23 "/usr/include/sys/resource.h" 2 3


# 1 "/usr/include/bits/resource.h" 1 3
# 24 "/usr/include/bits/resource.h" 3
# 1 "/usr/include/bits/types.h" 1 3
# 25 "/usr/include/bits/resource.h" 2 3







enum __rlimit_resource
{

  RLIMIT_CPU = 0,



  RLIMIT_FSIZE = 1,



  RLIMIT_DATA = 2,



  RLIMIT_STACK = 3,



  RLIMIT_CORE = 4,






  RLIMIT_RSS = 5,



  RLIMIT_NOFILE = 7,
  RLIMIT_OFILE = RLIMIT_NOFILE,




  RLIMIT_AS = 9,



  RLIMIT_NPROC = 6,



  RLIMIT_MEMLOCK = 8,



  RLIMIT_LOCKS = 10,


  RLIMIT_NLIMITS = 11,
  RLIM_NLIMITS = RLIMIT_NLIMITS


};
# 107 "/usr/include/bits/resource.h" 3
typedef __rlim_t rlim_t;







struct rlimit
  {

    rlim_t rlim_cur;

    rlim_t rlim_max;
  };
# 134 "/usr/include/bits/resource.h" 3
enum __rusage_who
{

  RUSAGE_SELF = 0,



  RUSAGE_CHILDREN = -1,



  RUSAGE_BOTH = -2

};


# 1 "/usr/include/bits/time.h" 1 3
# 151 "/usr/include/bits/resource.h" 2 3


struct rusage
  {

    struct timeval ru_utime;

    struct timeval ru_stime;

    long int ru_maxrss;


    long int ru_ixrss;

    long int ru_idrss;

    long int ru_isrss;


    long int ru_minflt;

    long int ru_majflt;

    long int ru_nswap;


    long int ru_inblock;

    long int ru_oublock;

    long int ru_msgsnd;

    long int ru_msgrcv;

    long int ru_nsignals;



    long int ru_nvcsw;


    long int ru_nivcsw;
  };







enum __priority_which
{
  PRIO_PROCESS = 0,

  PRIO_PGRP = 1,

  PRIO_USER = 2

};
# 26 "/usr/include/sys/resource.h" 2 3
# 43 "/usr/include/sys/resource.h" 3
typedef int __rlimit_resource_t;
typedef int __rusage_who_t;
typedef int __priority_which_t;





extern int getrlimit (__rlimit_resource_t __resource,
                      struct rlimit *__rlimits) ;
# 71 "/usr/include/sys/resource.h" 3
extern int setrlimit (__rlimit_resource_t __resource,
                      __const struct rlimit *__rlimits) ;
# 89 "/usr/include/sys/resource.h" 3
extern int getrusage (__rusage_who_t __who, struct rusage *__usage) ;





extern int getpriority (__priority_which_t __which, id_t __who) ;



extern int setpriority (__priority_which_t __which, id_t __who, int __prio)
            ;
# 33 "timel.c" 2
# 58 "timel.c"
struct tms tps;
static float timlim;
static time_t timstart, timlast;
static int tml_init = 1;
float deftim = 999.;
# 78 "timel.c"
static void time_st(timl)
float timl;
{
    times(&tps);
    timlim = timl;
    timstart = tps.tms_utime+tps.tms_cutime+tps.tms_stime+tps.tms_cstime;
    timlast = timstart;
    tml_init = 0;
    return;
}

static void time_init()
{
        struct rlimit rlimit;
        float maxtime;

        maxtime=deftim;

        if (getrlimit(RLIMIT_CPU, &rlimit)==0) {
                if ( rlimit.rlim_cur != ((unsigned long int)(~0UL)) )
                   maxtime = (float) rlimit.rlim_cur;
        }

        time_st(maxtime);
        return;
}

void timest_(timl)
float *timl;
{
 struct rlimit rlimit;
 float maxtime;

 if (tml_init != 0) {


    maxtime = *timl;
    if (getrlimit(RLIMIT_CPU, &rlimit)==0) {
           maxtime = (float) rlimit.rlim_cur;
           maxtime = ( maxtime > *timl ) ? *timl : maxtime;
    }
    time_st(maxtime);
 }
 return;
}
void timex_(tx)





float *tx;
{
   time_t timnow;
   if (tml_init) {
       time_init();
       *tx = 0.;
   }
   else {
       times(&tps);
       timnow = tps.tms_utime+tps.tms_cutime+tps.tms_stime+tps.tms_cstime;
       *tx = (float) (timnow - timstart) / 100;
   }
   return;
}

void timed_(td)





float *td;
{
   time_t timnow;
   if (tml_init) {
       time_init();
       *td = timlim;
   }
   else {
       times(&tps);
       timnow = tps.tms_utime+tps.tms_cutime+tps.tms_stime+tps.tms_cstime;
       *td = (float) (timnow - timlast) / 100;
       timlast = timnow;
   }
   return;
}

void timel_(tl)





float *tl;
{
   time_t timnow;
   if (tml_init) {
       time_init();
       *tl = timlim;
   }
   else {
       times(&tps);
       timnow = tps.tms_utime+tps.tms_cutime+tps.tms_stime+tps.tms_cstime;
       *tl = timlim - (float) (timnow - timstart) / 100;
   }
   return;
}
