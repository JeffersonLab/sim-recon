/*
 This file was generated from the file grkuta.F from the GEANT3.2.1 package.
 It was converted into C using the f2c program and the f2c.h header file
 was pasted directly into it in place of the #include.
 
 This was done because it makes the DANA swimmer agree very closely with
 the GEANT3 swimmer. It is not intended that this routine will survive in
 this form, but will be converted to a method of DMagneticFieldStepper
 in a more readable form.
 
 12/26/2006  -David Lawrence
*/




/* grkuta.F -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

/* f2c.h  --  Standard Fortran to C header file */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."

	- From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */

#ifndef F2C_INCLUDE
#define F2C_INCLUDE

typedef long int integer;
typedef unsigned long int uinteger;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;
#ifdef INTEGER_STAR_8	/* Adjust for integer*8. */
typedef long long longint;		/* system-dependent */
typedef unsigned long long ulongint;	/* system-dependent */
#define qbit_clear(a,b)	((a) & ~((ulongint)1 << (b)))
#define qbit_set(a,b)	((a) |  ((ulongint)1 << (b)))
#endif

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

/* I/O stuff */

#ifdef f2c_i2
/* for -i2 */
typedef short flag;
typedef short ftnlen;
typedef short ftnint;
#else
typedef long int flag;
typedef long int ftnlen;
typedef long int ftnint;
#endif

/*external read, write*/
typedef struct
{	flag cierr;
	ftnint ciunit;
	flag ciend;
	char *cifmt;
	ftnint cirec;
} cilist;

/*internal read, write*/
typedef struct
{	flag icierr;
	char *iciunit;
	flag iciend;
	char *icifmt;
	ftnint icirlen;
	ftnint icirnum;
} icilist;

/*open*/
typedef struct
{	flag oerr;
	ftnint ounit;
	char *ofnm;
	ftnlen ofnmlen;
	char *osta;
	char *oacc;
	char *ofm;
	ftnint orl;
	char *oblnk;
} olist;

/*close*/
typedef struct
{	flag cerr;
	ftnint cunit;
	char *csta;
} cllist;

/*rewind, backspace, endfile*/
typedef struct
{	flag aerr;
	ftnint aunit;
} alist;

/* inquire */
typedef struct
{	flag inerr;
	ftnint inunit;
	char *infile;
	ftnlen infilen;
	ftnint	*inex;	/*parameters in standard's order*/
	ftnint	*inopen;
	ftnint	*innum;
	ftnint	*innamed;
	char	*inname;
	ftnlen	innamlen;
	char	*inacc;
	ftnlen	inacclen;
	char	*inseq;
	ftnlen	inseqlen;
	char 	*indir;
	ftnlen	indirlen;
	char	*infmt;
	ftnlen	infmtlen;
	char	*inform;
	ftnint	informlen;
	char	*inunf;
	ftnlen	inunflen;
	ftnint	*inrecl;
	ftnint	*innrec;
	char	*inblank;
	ftnlen	inblanklen;
} inlist;

#define VOID void

union Multitype {	/* for multiple entry points */
	integer1 g;
	shortint h;
	integer i;
	/* longint j; */
	real r;
	doublereal d;
	complex c;
	doublecomplex z;
	};

typedef union Multitype Multitype;

/*typedef long int Long;*/	/* No longer used; formerly in Namelist */

struct Vardesc {	/* for Namelist */
	char *name;
	char *addr;
	ftnlen *dims;
	int  type;
	};
typedef struct Vardesc Vardesc;

struct Namelist {
	char *name;
	Vardesc **vars;
	int nvars;
	};
typedef struct Namelist Namelist;

#define loc_abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)loc_abs(x)
//#define min(a,b) ((a) <= (b) ? (a) : (b))
//#define max(a,b) ((a) >= (b) ? (a) : (b))
//#define dmin(a,b) (doublereal)min(a,b)
//#define dmax(a,b) (doublereal)max(a,b)
#define bit_test(a,b)	((a) >> (b) & 1)
#define bit_clear(a,b)	((a) & ~((uinteger)1 << (b)))
#define bit_set(a,b)	((a) |  ((uinteger)1 << (b)))

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef int /* Unknown procedure type */ (*U_fp)(...);
typedef shortint (*J_fp)(...);
typedef integer (*I_fp)(...);
typedef real (*R_fp)(...);
typedef doublereal (*D_fp)(...), (*E_fp)(...);
typedef /* Complex */ VOID (*C_fp)(...);
typedef /* Double Complex */ VOID (*Z_fp)(...);
typedef logical (*L_fp)(...);
typedef shortlogical (*K_fp)(...);
typedef /* Character */ VOID (*H_fp)(...);
typedef /* Subroutine */ int (*S_fp)(...);
#else
typedef int /* Unknown procedure type */ (*U_fp)();
typedef shortint (*J_fp)();
typedef integer (*I_fp)();
typedef real (*R_fp)();
typedef doublereal (*D_fp)(), (*E_fp)();
typedef /* Complex */ VOID (*C_fp)();
typedef /* Double Complex */ VOID (*Z_fp)();
typedef logical (*L_fp)();
typedef shortlogical (*K_fp)();
typedef /* Character */ VOID (*H_fp)();
typedef /* Subroutine */ int (*S_fp)();
#endif
/* E_fp is for real functions when -R is not specified */
typedef VOID C_f;	/* complex function */
typedef VOID H_f;	/* character function */
typedef VOID Z_f;	/* double complex function */
typedef doublereal E_f;	/* real function with -R not specified */

/* undef any lower-case symbols that your C compiler predefines, e.g.: */

#ifndef Skip_f2c_Undefs
#undef cray
#undef gcos
#undef mc68010
#undef mc68020
#undef mips
#undef pdp11
#undef sgi
#undef sparc
#undef sun
#undef sun2
#undef sun3
#undef sun4
#undef u370
#undef u3b
#undef u3b2
#undef u3b5
#undef unix
#undef vax
#endif
#endif


/* $Id: grkuta.F,v 1.1.1.1 1995/10/24 10:21:42 cernlib Exp $ */

/* $Log: grkuta.F,v $ */
/* Revision 1.1.1.1  1995/10/24 10:21:42  cernlib */
/* Geant */

#include "HDGEOMETRY/DMagneticFieldMap.h"
//#include "cmath"
//using namespace std;

/* Builtin functions */
extern "C" {double sqrt(double), sin(double);}

/* #include "geant321/pilot.h" */
/* CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani */
/* -- Author : */
/* Subroutine */ int grkuta_(double *charge, double *step, double *vect, double *vout, const DMagneticFieldMap *bfield)
{
    /* System generated locals */
    double d__1, d__2, d__3;
    /*static*/ double equiv_2[3], equiv_5[3];


    /* Local variables */
    /*static*/ double a, b, c__;
    /*static*/ /*double f[4];*/
    /*static*/ double *f=&vout[7]; /* not that vout is decremented below so indexes start at 1 like FORTRAN */
    /*static*/ double h__;
    /*static*/ integer j;
#define x (equiv_2)
#define y (equiv_2 + 1)
#define z__ (equiv_2 + 2)
    /*static*/ double f1, f2, f3, h2, f4, h4, g1, g2, g3, g4, g5, g6, at, bt, 
	    ct, ph, hp, tl;
#define xt (equiv_5)
#define yt (equiv_5 + 1)
#define zt (equiv_5 + 2)
    /*static*/ double ph2, cba, rho, est, tet, hxp[3], dxt, dyt, dzt, ang2;
#define xyz (equiv_2)
    /*static*/ double rho1;
    /*static*/ integer iter;
    /*static*/ double cost;
    /*static*/ integer ncut;
    /*static*/ double pinv, rest, sint;
#define xyzt (equiv_5)
    extern /* Subroutine */ int gufld_(double *, double *);
    /*static*/ double hnorm, secxs[4], secys[4], seczs[4];

/* . */
/* .    ****************************************************************** */
/* .    *                                                                * */
/* .    *  Runge-Kutta method for tracking a particle through a magnetic * */
/* .    *  field. Uses Nystroem algorithm (See Handbook Nat. Bur. of     * */
/* .    *  Standards, procedure 25.5.20)                                 * */
/* .    *                                                                * */
/* .    *  Input parameters                                              * */
/* .    *       CHARGE    Particle charge                                * */
/* .    *       STEP      Step size                                      * */
/* .    *       VECT      Initial co-ords,direction cosines,momentum     * */
/* .    *  Output parameters                                             * */
/* .    *       VOUT      Output co-ords,direction cosines,momentum      * */
/* .    *  User routine called                                           * */
/* .    *       CALL GUFLD(X,F)                                          * */
/* .    *                                                                * */
/* .    *    ==>Called by : <USER>, GUSWIM                               * */
/* .    *       Authors    R.Brun, M.Hansroul  *********                 * */
/* .    *                  V.Perevoztchikov (CUT STEP implementation)    * */
/* .    *                                                                * */
/* .    *                                                                * */
/* .    ****************************************************************** */
/* . */

/* . */
/* .    ------------------------------------------------------------------ */
/* . */
/*             This constant is for units CM,GEV/C and KGAUSS */

    /* Parameter adjustments */
    --vout;
    --vect;

    /* Function Body */
    iter = 0;
    ncut = 0;
    for (j = 1; j <= 7; ++j) {
	vout[j] = vect[j];
/* L10: */
    }
    pinv = *charge * 2.9979251e-3 / vect[7];
    tl = 0.f;
    h__ = *step;


L20:
    rest = *step - tl;
    if (loc_abs(h__) > loc_abs(rest)) {
	h__ = rest;
    }
    //gufld_(&vout[1], f);
    bfield->GetField(vout[1], vout[2], vout[3], f[0], f[1], f[2]);


/*             Start of integration */

    *x = vout[1];
    *y = vout[2];
    *z__ = vout[3];
    a = vout[4];
    b = vout[5];
    c__ = vout[6];

    h2 = h__ * .5;
    h4 = h2 * .5;
    ph = pinv * h__;
    ph2 = ph * .5;
    secxs[0] = (b * f[2] - c__ * f[1]) * ph2;
    secys[0] = (c__ * f[0] - a * f[2]) * ph2;
    seczs[0] = (a * f[1] - b * f[0]) * ph2;
/* Computing 2nd power */
    d__1 = secxs[0];
/* Computing 2nd power */
    d__2 = secys[0];
/* Computing 2nd power */
    d__3 = seczs[0];
    ang2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    if (ang2 > 9.86960440109) {
	goto L40;
    }
    dxt = h2 * a + h4 * secxs[0];
    dyt = h2 * b + h4 * secys[0];
    dzt = h2 * c__ + h4 * seczs[0];
    *xt = *x + dxt;
    *yt = *y + dyt;
    *zt = *z__ + dzt;

/*              Second intermediate point */

    est = loc_abs(dxt) + loc_abs(dyt) + loc_abs(dzt);
    if (est > h__) {
	goto L30;
    }
    //gufld_(xyzt, f);
	 bfield->GetField(xyzt[0], xyzt[1], xyzt[2], f[0], f[1], f[2]);
    at = a + secxs[0];
    bt = b + secys[0];
    ct = c__ + seczs[0];

    secxs[1] = (bt * f[2] - ct * f[1]) * ph2;
    secys[1] = (ct * f[0] - at * f[2]) * ph2;
    seczs[1] = (at * f[1] - bt * f[0]) * ph2;
    at = a + secxs[1];
    bt = b + secys[1];
    ct = c__ + seczs[1];
    secxs[2] = (bt * f[2] - ct * f[1]) * ph2;
    secys[2] = (ct * f[0] - at * f[2]) * ph2;
    seczs[2] = (at * f[1] - bt * f[0]) * ph2;
    dxt = h__ * (a + secxs[2]);
    dyt = h__ * (b + secys[2]);
    dzt = h__ * (c__ + seczs[2]);
    *xt = *x + dxt;
    *yt = *y + dyt;
    *zt = *z__ + dzt;
    at = a + secxs[2] * 2.;
    bt = b + secys[2] * 2.;
    ct = c__ + seczs[2] * 2.;

    est = loc_abs(dxt) + loc_abs(dyt) + loc_abs(dzt);
    if (est > loc_abs(h__) * 2.f) {
	goto L30;
    }
    //gufld_(xyzt, f);
	 bfield->GetField(xyzt[0], xyzt[1], xyzt[2], f[0], f[1], f[2]);

    *z__ += (c__ + (seczs[0] + seczs[1] + seczs[2]) * .33333333333333331) * 
	    h__;
    *y += (b + (secys[0] + secys[1] + secys[2]) * .33333333333333331) * h__;
    *x += (a + (secxs[0] + secxs[1] + secxs[2]) * .33333333333333331) * h__;

    secxs[3] = (bt * f[2] - ct * f[1]) * ph2;
    secys[3] = (ct * f[0] - at * f[2]) * ph2;
    seczs[3] = (at * f[1] - bt * f[0]) * ph2;
    a += (secxs[0] + secxs[3] + (secxs[1] + secxs[2]) * 2.) * 
	    .33333333333333331;
    b += (secys[0] + secys[3] + (secys[1] + secys[2]) * 2.) * 
	    .33333333333333331;
    c__ += (seczs[0] + seczs[3] + (seczs[1] + seczs[2]) * 2.) * 
	    .33333333333333331;

    est = (d__1 = secxs[0] + secxs[3] - (secxs[1] + secxs[2]), loc_abs(d__1)) + (
	    d__2 = secys[0] + secys[3] - (secys[1] + secys[2]), loc_abs(d__2)) + (
	    d__3 = seczs[0] + seczs[3] - (seczs[1] + seczs[2]), loc_abs(d__3));

    if (est > 1e-4 && loc_abs(h__) > 1e-4f) {
	goto L30;
    }
    ++iter;
    ncut = 0;
/*               If too many iterations, go to HELIX */
    if (iter > 1992) {
	goto L40;
    }

    tl += h__;
    if (est < 3.1250000000000001e-6) {
	h__ *= 2.;
    }
    cba = 1. / sqrt(a * a + b * b + c__ * c__);
    vout[1] = *x;
    vout[2] = *y;
    vout[3] = *z__;
    vout[4] = cba * a;
    vout[5] = cba * b;
    vout[6] = cba * c__;
    rest = *step - tl;
    if (*step < 0.f) {
	rest = -rest;
    }
    if (rest > dabs(*step) * 1e-5f) {
	goto L20;
    }

    goto L999;

/* *              CUT STEP */
L30:
    ++ncut;
/*               If too many cuts , go to HELIX */
    if (ncut > 11) {
	goto L40;
    }
    h__ *= .5;
    goto L20;

/* *              ANGLE TOO BIG, USE HELIX */
L40:
    f1 = f[0];
    f2 = f[1];
    f3 = f[2];
/* Computing 2nd power */
    d__1 = f1;
/* Computing 2nd power */
    d__2 = f2;
/* Computing 2nd power */
    d__3 = f3;
    f4 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    rho = -f4 * pinv;
    tet = rho * *step;
    if (tet != 0.f) {
	hnorm = 1. / f4;
	f1 *= hnorm;
	f2 *= hnorm;
	f3 *= hnorm;

	hxp[0] = f2 * vect[6] - f3 * vect[5];
	hxp[1] = f3 * vect[4] - f1 * vect[6];
	hxp[2] = f1 * vect[5] - f2 * vect[4];
	hp = f1 * vect[4] + f2 * vect[5] + f3 * vect[6];

	rho1 = 1. / rho;
	sint = sin(tet);
/* Computing 2nd power */
	d__1 = sin(tet * .5);
	cost = d__1 * d__1 * 2.;

	g1 = sint * rho1;
	g2 = cost * rho1;
	g3 = (tet - sint) * hp * rho1;
	g4 = -cost;
	g5 = sint;
	g6 = cost * hp;
	vout[1] = vect[1] + (g1 * vect[4] + g2 * hxp[0] + g3 * f1);
	vout[2] = vect[2] + (g1 * vect[5] + g2 * hxp[1] + g3 * f2);
	vout[3] = vect[3] + (g1 * vect[6] + g2 * hxp[2] + g3 * f3);
	vout[4] = vect[4] + (g4 * vect[4] + g5 * hxp[0] + g6 * f1);
	vout[5] = vect[5] + (g4 * vect[5] + g5 * hxp[1] + g6 * f2);
	vout[6] = vect[6] + (g4 * vect[6] + g5 * hxp[2] + g6 * f3);

    } else {
	vout[1] = vect[1] + *step * vect[4];
	vout[2] = vect[2] + *step * vect[5];
	vout[3] = vect[3] + *step * vect[6];

    }

L999:
    return 0;
} /* grkuta_ */

#undef xyzt
#undef xyz
#undef zt
#undef yt
#undef xt
#undef z__
#undef y
#undef x


