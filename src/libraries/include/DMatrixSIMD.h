#ifndef _DMatrixSIMD_
#define _DMatrixSIMD_

#include <math.h>
#ifdef USE_SIMD
#include <emmintrin.h> // Header file for SSE2 SIMD instructions
#endif
#ifdef USE_SSE3
#include <pmmintrin.h> // Header file for SSE3 SIMD instructions
#endif
#include <iostream>
#include <iomanip>
using namespace std;

#include "DMatrix2x1.h"
#include "DMatrix2x2.h"
#include "DMatrix3x2.h"
#include "DMatrix3x3.h"
#include "DMatrix2x3.h"
#include "DMatrix4x2.h"
#include "DMatrix4x4.h"
#include "DMatrix2x4.h"


#ifndef USE_SIMD

// Matrix multiplication:  (3x2) x (2x3)
inline DMatrix3x3 operator*(const DMatrix3x2 &m1,const DMatrix2x3 &m2){
  return DMatrix3x3(m1(0,0)*m2(0,0)+m1(0,1)*m2(1,0),
		    m1(0,0)*m2(0,1)+m1(0,1)*m2(1,1),
		    m1(0,0)*m2(0,2)+m1(0,1)*m2(1,2),
		    m1(1,0)*m2(0,0)+m1(1,1)*m2(1,0),
		    m1(1,0)*m2(0,1)+m1(1,1)*m2(1,1),
		    m1(1,0)*m2(0,2)+m1(1,1)*m2(1,2),
		    m1(2,0)*m2(0,0)+m1(2,1)*m2(1,0),
		    m1(2,0)*m2(0,1)+m1(2,1)*m2(1,1),
		    m1(2,0)*m2(0,2)+m1(2,1)*m2(1,2)
		    );
}

// Matrix multiplication:  (2x2) x (2x3)
inline DMatrix2x3 operator*(const DMatrix2x2 &m1,const DMatrix2x3 &m2){
  return DMatrix2x3(m1(0,0)*m2(0,0)+m1(0,1)*m2(1,0),
		    m1(0,0)*m2(0,1)+m1(0,1)*m2(1,1),
		    m1(0,0)*m2(0,2)+m1(0,1)*m2(1,2),
		    m1(1,0)*m2(0,0)+m1(1,1)*m2(1,0),
		    m1(1,0)*m2(0,1)+m1(1,1)*m2(1,1),
		    m1(1,0)*m2(0,2)+m1(1,1)*m2(1,2)
		    );
  

}

#else

// Matrix multiplication:  (3x2) x (2x3)
inline DMatrix3x3 operator*(const DMatrix3x2 &m1,
			    const DMatrix2x3 &m2){
  __m128d a11=_mm_set1_pd(m2(0,0));
  __m128d a12=_mm_set1_pd(m2(0,1)); 
  __m128d a13=_mm_set1_pd(m2(0,2));
  __m128d a21=_mm_set1_pd(m2(1,0));
  __m128d a22=_mm_set1_pd(m2(1,1)); 
  __m128d a23=_mm_set1_pd(m2(1,2));
  return DMatrix3x3(_mm_add_pd(_mm_mul_pd(m1.GetV(0,0),a11),
			       _mm_mul_pd(m1.GetV(0,1),a21)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(0,0),a12),
			       _mm_mul_pd(m1.GetV(0,1),a22)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(0,0),a13),
			       _mm_mul_pd(m1.GetV(0,1),a23)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(1,0),a11),
			       _mm_mul_pd(m1.GetV(1,1),a21)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(1,0),a12),
			       _mm_mul_pd(m1.GetV(1,1),a22)),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(1,0),a13),
			       _mm_mul_pd(m1.GetV(1,1),a23)));
}

// Matrix multiplication:  (2x2) x (2x3)
inline DMatrix2x3 operator*(const DMatrix2x2 &m1,const DMatrix2x3 &m2){
  return DMatrix2x3(_mm_add_pd(_mm_mul_pd(m1.GetV(0),_mm_set1_pd(m2(0,0))),
			       _mm_mul_pd(m1.GetV(1),_mm_set1_pd(m2(1,0)))),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(0),_mm_set1_pd(m2(0,1))),
			       _mm_mul_pd(m1.GetV(1),_mm_set1_pd(m2(1,1)))),
		    _mm_add_pd(_mm_mul_pd(m1.GetV(0),_mm_set1_pd(m2(0,2))),
			       _mm_mul_pd(m1.GetV(1),_mm_set1_pd(m2(1,2)))));
}



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
#endif


#include "DMatrix5x1.h"
#include "DMatrix5x2.h"
#include "DMatrix5x5.h"
#include "DMatrix2x5.h"
#include "DMatrix1x5.h"

#ifndef USE_SIMD

// Multiply a 5x1 matrix by its transpose 
inline DMatrix5x5 MultiplyTranspose(const DMatrix5x1 &m1){
  DMatrix5x5 temp;
  for (unsigned int i=0;i<5;i++){
    for (unsigned int j=0;j<5;j++){
      temp(i,j)=m1(i)*m1(j);      
    }
  }
  return temp;
}

// Multiply a 5x1 matrix by a 1x5 matrix
inline DMatrix5x5 operator*(const DMatrix5x1 &m1,const DMatrix1x5 &m2){
  DMatrix5x5 temp;
  for (unsigned int i=0;i<5;i++){
    for (unsigned int j=0;j<5;j++){
      temp(i,j)=m1(i)*m2(j);      
    }
  }
  return temp;
}

// Multiply a 5x2 matrix by a 2x5 matrix
inline DMatrix5x5 operator*(const DMatrix5x2 &m1,const DMatrix2x5 &m2){
  DMatrix5x5 temp;
  for (unsigned int i=0;i<5;i++){
    for (unsigned int j=0;j<5;j++){
      temp(i,j)=m1(i,0)*m2(0,j)+m1(i,1)*m2(1,j);      
    }
  }
  return temp;

}

#else

// Multiply a 5x1 matrix by its transpose 
inline DMatrix5x5 MultiplyTranspose(const DMatrix5x1 &m1){
  __m128d b1=_mm_set1_pd(m1(0));
  __m128d b2=_mm_set1_pd(m1(1));
  __m128d b3=_mm_set1_pd(m1(2));
  __m128d b4=_mm_set1_pd(m1(3));
  __m128d b5=_mm_set1_pd(m1(4));
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

#endif
#endif // _DMatrixSIMD_
