#ifndef _DMatrixSIMD_
#define _DMatrixSIMD_

#include <math.h>
#include <emmintrin.h> // Header file for SSE2 SIMD instructions
#include <iostream>
#include <iomanip>
using namespace std;

#include "DMatrix2x1.h"
#include "DMatrix2x2.h"
#include "DMatrix4x2.h"
#include "DMatrix4x4.h"
#include "DMatrix2x4.h"

// Matrix multiplication:  (4x2) x (2x4)
inline DMatrix4x4 operator*(const DMatrix4x2 &m1,
			    const DMatrix2x4 &m2){
  __m128d a11=_mm_set1_pd(m2(0,0));
  __m128d a12=_mm_set1_pd(m2(0,1)); 
  __m128d a13=_mm_set1_pd(m2(0,2));
  __m128d a14=_mm_set1_pd(m2(0,3));
  __m128d a21=_mm_set1_pd(m2(1,0));
  __m128d a22=_mm_set1_pd(m2(1,1)); 
  __m128d a23=_mm_set1_pd(m2(1,2));
  __m128d a24=_mm_set1_pd(m2(1,3));
  return DMatrix4x4(_mm_add_pd(_mm_mul_pd(m1.GetV(0,0),a11),
			       _mm_mul_pd(m1.GetV(0,1),a21)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(0,0),a12),
			       _mm_mul_pd(m1.GetV(0,1),a22)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(0,0),a13),
			       _mm_mul_pd(m1.GetV(0,1),a23)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(0,0),a14),
			       _mm_mul_pd(m1.GetV(0,1),a24)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(1,0),a11),
			       _mm_mul_pd(m1.GetV(1,1),a21)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(1,0),a12),
			       _mm_mul_pd(m1.GetV(1,1),a22)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(1,0),a13),
			       _mm_mul_pd(m1.GetV(1,1),a23)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(1,0),a14),
			       _mm_mul_pd(m1.GetV(1,1),a24)));
}

// Matrix multiplication:  (2x2) x (2x4)
inline DMatrix2x4 operator*(const DMatrix2x2 &m1,const DMatrix2x4 &m2){
  return DMatrix2x4(_mm_add_pd(_mm_mul_pd(m1.GetV(0),_mm_set1_pd(m2(0,0))),
			       _mm_mul_pd(m1.GetV(1),_mm_set1_pd(m2(1,0)))),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(0),_mm_set1_pd(m2(0,1))),
			       _mm_mul_pd(m1.GetV(1),_mm_set1_pd(m2(1,1)))),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(0),_mm_set1_pd(m2(0,2))),
			       _mm_mul_pd(m1.GetV(1),_mm_set1_pd(m2(1,2)))), 
		    _mm_add_pd(_mm_mul_pd(m1.GetV(0),_mm_set1_pd(m2(0,3))),
			       _mm_mul_pd(m1.GetV(1),_mm_set1_pd(m2(1,3)))));
}


#include "DMatrix5x1.h"
#include "DMatrix5x2.h"
#include "DMatrix5x5.h"
#include "DMatrix2x5.h"
#include "DMatrix1x5.h"

// Multiply a 5x1 matrix by a 1x5 matrix
inline DMatrix5x5 operator*(const DMatrix5x1 &m1,const DMatrix1x5 &m2){
  __m128d b1=_mm_set1_pd(m2(0));
  __m128d b2=_mm_set1_pd(m2(1));
  __m128d b3=_mm_set1_pd(m2(2));
  __m128d b4=_mm_set1_pd(m2(3));
  __m128d b5=_mm_set1_pd(m2(4));
  return DMatrix5x5(_mm_mul_pd(m1.GetV(0),b1),_mm_mul_pd(m1.GetV(0),b2),
		    _mm_mul_pd(m1.GetV(0),b3),_mm_mul_pd(m1.GetV(0),b4),
		    _mm_mul_pd(m1.GetV(0),b5),
		    _mm_mul_pd(m1.GetV(1),b1),_mm_mul_pd(m1.GetV(1),b2),
		    _mm_mul_pd(m1.GetV(1),b3),_mm_mul_pd(m1.GetV(1),b4),
		    _mm_mul_pd(m1.GetV(1),b5),
		    _mm_mul_pd(m1.GetV(2),b1),_mm_mul_pd(m1.GetV(2),b2),
		    _mm_mul_pd(m1.GetV(2),b3),_mm_mul_pd(m1.GetV(2),b4),
		    _mm_mul_pd(m1.GetV(2),b5));
}

// Multiply a 5x2 matrix by a 2x5 matrix
inline DMatrix5x5 operator*(const DMatrix5x2 &m1,const DMatrix2x5 &m2){
  __m128d a11=_mm_set1_pd(m2(0,0));
  __m128d a12=_mm_set1_pd(m2(0,1)); 
  __m128d a13=_mm_set1_pd(m2(0,2));
  __m128d a14=_mm_set1_pd(m2(0,3));
  __m128d a15=_mm_set1_pd(m2(0,4));
  __m128d a21=_mm_set1_pd(m2(1,0));
  __m128d a22=_mm_set1_pd(m2(1,1)); 
  __m128d a23=_mm_set1_pd(m2(1,2));
  __m128d a24=_mm_set1_pd(m2(1,3)); 
  __m128d a25=_mm_set1_pd(m2(1,4));

  return DMatrix5x5(_mm_add_pd(_mm_mul_pd(m1.GetV(0,0),a11),
			       _mm_mul_pd(m1.GetV(0,1),a21)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(0,0),a12),
			       _mm_mul_pd(m1.GetV(0,1),a22)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(0,0),a13),
			       _mm_mul_pd(m1.GetV(0,1),a23)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(0,0),a14),
			       _mm_mul_pd(m1.GetV(0,1),a24)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(0,0),a15),
			       _mm_mul_pd(m1.GetV(0,1),a25)),
		    
		    _mm_add_pd(_mm_mul_pd(m1.GetV(1,0),a11),
			       _mm_mul_pd(m1.GetV(1,1),a21)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(1,0),a12),
			       _mm_mul_pd(m1.GetV(1,1),a22)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(1,0),a13),
			       _mm_mul_pd(m1.GetV(1,1),a23)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(1,0),a14),
			       _mm_mul_pd(m1.GetV(1,1),a24)), 
		    _mm_add_pd(_mm_mul_pd(m1.GetV(1,0),a15),
			       _mm_mul_pd(m1.GetV(1,1),a25)),

		    _mm_add_pd(_mm_mul_pd(m1.GetV(2,0),a11),
			       _mm_mul_pd(m1.GetV(2,1),a21)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(2,0),a12),
			       _mm_mul_pd(m1.GetV(2,1),a22)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(2,0),a13),
			       _mm_mul_pd(m1.GetV(2,1),a23)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(2,0),a14),
			       _mm_mul_pd(m1.GetV(2,1),a24)), 
		    _mm_add_pd(_mm_mul_pd(m1.GetV(2,0),a15),
			       _mm_mul_pd(m1.GetV(2,1),a25)));		    
}

#endif // _DMatrixSIMD_
