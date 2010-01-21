/*
 * ntypes.h
 *
*/

#ifndef ntypesH_INCLUDED
#define ntypesH_INCLUDED

static const char sccsid_ntypesH[] = "@(#)ntypes.h\t5.2\tCreated 7/27/97 19:02:56, \tcompiled "__DATE__;

/*
 * This include file defines the basic types
 * used to define the data structures written on E852 data tapes
 * and elsewhere.
 *
 * int8, uint8    -  8 bit wide MSB integers (more commonly known as 'char')
 * int16, uint16  - 16 bit wide MSB integers (usually 'short')
 * int32, uint32  - 32 bit wide MSB integers (usually 'int')
 * int64, uint64  - 64 bit wide MSB integers (on some systems 'long', on others 'long long')
 *
 * float32   - 32 bit wide IEEE-754 floating point numbers (usually C 'float',  FORTRAN 'REAL')
 * float64   - 64 bit wide IEEE-754 floating point numbers (usually C 'double', FORTRAN 'DOUBLE PRECISION')
 * 
 * MSB means "big-endian", like SGI-MIPS, RS6000 and others.
 * LSB means "little-endian", like Intel x86, WindowsNT and others.
 *
 * Read 'man math' for more details on IEEE-754 floating point numbers.
 *
*/

#ifdef sgi

/*
 * SGI-only definitions 
*/

typedef int                   int32;
typedef unsigned int         uint32;
typedef long long             int64;
typedef unsigned long long   uint64;

#else

/*
 * these definitions work for rs6000 under AIX 3.2.5, untested for others
 * The IBM xlc compiler seem to understand 'long long', but it is not documented, so...
*/

typedef int                   int32;
typedef unsigned int         uint32;
typedef long long             int64;
typedef unsigned long long   uint64;

#endif

typedef short             int16;
typedef unsigned short   uint16;
typedef char               int8;
typedef unsigned char     uint8;

typedef float    float32;
typedef double   float64;

/* These are the 'float' structures ... */

typedef struct { float32 x,y,z; }               vector3_t;
typedef struct { float32 t; vector3_t space; }  vector4_t;
typedef struct { float32 rho,phi,z;}            vector3cyl_t;
typedef struct { vector3cyl_t vert[8];}         cylcoord_t;
typedef struct { vector3_t    vert[8];}         cartcoord_t;

/* ... and their 'double' siblings */

typedef struct { float64 x,y,z; } dvector3_t;
typedef struct { float64 t; dvector3_t space; } dvector4_t;
typedef struct { float64 rho,phi,z;     } dvector3cyl_t;
typedef struct { dvector3cyl_t vert[8]; } dcylcoord_t;
typedef struct { dvector3_t    vert[8]; } dcartcoord_t;


/*
 * Kludge to get rid of the Particle_t without breaking everybody's code
*/

#include <particleType.h>


#endif
/* end file */
