#ifndef _DMatrixSIMD_
#define _DMatrixSIMD_

#include <math.h>
#ifdef USE_SSE2
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
#include "DMatrix3x1.h"
#include "DMatrix3x2.h"
#include "DMatrix3x3.h"
#include "DMatrix2x3.h"
#include "DMatrix1x3.h"
#include "DMatrix1x2.h"
// We are not currently using any of the 4x2,4x4,etc matrices in the main code
#include "DMatrix4x1.h"
#include "DMatrix4x2.h"
#include "DMatrix4x4.h"
#include "DMatrix2x4.h"
#include "DMatrix1x4.h"
#include "align_16.h"


#ifndef USE_SSE2

// Matrix multiplication:  (2x1) x (1x2)
inline DMatrix2x2 operator*(const DMatrix2x1 &m1,const DMatrix1x2 &m2){
  return DMatrix2x2(m1(0)*m2(0),m1(0)*m2(1),
		    m1(1)*m2(0),m1(1)*m2(1)
		    );
}


// Matrix multiplication:  (3x1) x (1x3)
inline DMatrix3x3 operator*(const DMatrix3x1 &m1,const DMatrix1x3 &m2){
  return DMatrix3x3(m1(0)*m2(0),m1(0)*m2(1),m1(0)*m2(2),
		    m1(1)*m2(0),m1(1)*m2(1),m1(1)*m2(2),
		    m1(2)*m2(0),m1(2)*m2(1),m1(2)*m2(2)
		    );
}


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

// Matrix multiplication:  (4x2) x (2x4)
inline DMatrix4x4 operator*(const DMatrix4x2 &m1,const DMatrix2x4 &m2){
  return DMatrix4x4(m1(0,0)*m2(0,0)+m1(0,1)*m2(1,0),
		    m1(0,0)*m2(0,1)+m1(0,1)*m2(1,1),
		    m1(0,0)*m2(0,2)+m1(0,1)*m2(1,2),
		    m1(0,0)*m2(0,3)+m1(0,1)*m2(1,3),

		    m1(1,0)*m2(0,0)+m1(1,1)*m2(1,0),
		    m1(1,0)*m2(0,1)+m1(1,1)*m2(1,1),
		    m1(1,0)*m2(0,2)+m1(1,1)*m2(1,2),
		    m1(1,0)*m2(0,3)+m1(1,1)*m2(1,3),

		    m1(2,0)*m2(0,0)+m1(2,1)*m2(1,0),
		    m1(2,0)*m2(0,1)+m1(2,1)*m2(1,1),
		    m1(2,0)*m2(0,2)+m1(2,1)*m2(1,2),
		    m1(2,0)*m2(0,3)+m1(2,1)*m2(1,3),    

		    m1(3,0)*m2(0,0)+m1(3,1)*m2(1,0),
		    m1(3,0)*m2(0,1)+m1(3,1)*m2(1,1),
		    m1(3,0)*m2(0,2)+m1(3,1)*m2(1,2),
		    m1(3,0)*m2(0,3)+m1(3,1)*m2(1,3)
		    );
}


#else

// Matrix multiplication:  (3x2) x (2x3)
inline DMatrix3x3 operator*(const DMatrix3x2 &m1,
			    const DMatrix2x3 &m2){
  ALIGNED_16_BLOCK_WITH_PTR(__m128d, 6, p)
  __m128d &a11=p[0];
  __m128d &a12=p[1];
  __m128d &a13=p[2];
  __m128d &a21=p[3];
  __m128d &a22=p[4];
  __m128d &a23=p[5];
  a11=_mm_set1_pd(m2(0,0));
  a12=_mm_set1_pd(m2(0,1)); 
  a13=_mm_set1_pd(m2(0,2));
  a21=_mm_set1_pd(m2(1,0));
  a22=_mm_set1_pd(m2(1,1)); 
  a23=_mm_set1_pd(m2(1,2));
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

//----------------------------------------------------------------------------
// NOTE:  We are not currently using any of the 2x4,4x2, etc matrices in the 
// code.  As a consequence I am commenting out the following routines until 
// I get around to writing the non-SIMD versions of this code... whenever that 
// may be...
// Matrix multiplication:  (4x2) x (2x4)

inline DMatrix4x4 operator*(const DMatrix4x2 &m1,
			    const DMatrix2x4 &m2){
  ALIGNED_16_BLOCK_WITH_PTR(__m128d, 8, p)
    __m128d &a11=p[0];
  __m128d &a12=p[1];
  __m128d &a13=p[2];
  __m128d &a14=p[3];
  __m128d &a21=p[4];
  __m128d &a22=p[5];
  __m128d &a23=p[6];
  __m128d &a24=p[7];
  a11=_mm_set1_pd(m2(0,0));
  a12=_mm_set1_pd(m2(0,1)); 
  a13=_mm_set1_pd(m2(0,2));
  a14=_mm_set1_pd(m2(0,3));
  a21=_mm_set1_pd(m2(1,0));
  a22=_mm_set1_pd(m2(1,1)); 
  a23=_mm_set1_pd(m2(1,2));
  a24=_mm_set1_pd(m2(1,3));
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

// end of un-used matrix block...
#endif


#include "DMatrix5x1.h"
#include "DMatrix5x2.h"
#include "DMatrix5x5.h"
#include "DMatrix2x5.h"
#include "DMatrix1x5.h"

#ifndef USE_SSE2

// Multiply a 4x1 matrix by a 1x4 matrix
inline DMatrix4x4 operator*(const DMatrix4x1 &m1,const DMatrix1x4 &m2){
  DMatrix4x4 temp;
  for (unsigned int i=0;i<4;i++){
    for (unsigned int j=0;j<4;j++){
      temp(i,j)=m1(i)*m2(j);      
    }
  }
  return temp;
}



// Find the tranpose of a 5x2 matrix
inline DMatrix2x5 Transpose(const DMatrix5x2 &M){
  return DMatrix2x5(M(0,0),M(1,0),M(2,0),M(3,0),M(4,0),
		    M(0,1),M(1,1),M(2,1),M(3,1),M(4,1));
}

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

// Multiply a 4x1 matrix by a 1x4 matrix
inline DMatrix4x4 operator*(const DMatrix4x1 &m1,const DMatrix1x4 &m2){
  ALIGNED_16_BLOCK_WITH_PTR(__m128d, 4, p)
  __m128d &b1=p[0];
  __m128d &b2=p[1];
  __m128d &b3=p[2];
  __m128d &b4=p[3];
  b1=_mm_set1_pd(m2(0));
  b2=_mm_set1_pd(m2(1));
  b3=_mm_set1_pd(m2(2));
  b4=_mm_set1_pd(m2(3));
  return DMatrix4x4(_mm_mul_pd(m1.GetV(0),b1),_mm_mul_pd(m1.GetV(0),b2),
		    _mm_mul_pd(m1.GetV(0),b3),_mm_mul_pd(m1.GetV(0),b4),
		    _mm_mul_pd(m1.GetV(1),b1),_mm_mul_pd(m1.GetV(1),b2),
		    _mm_mul_pd(m1.GetV(1),b3),_mm_mul_pd(m1.GetV(1),b4)
		    );
}




// Find the tranpose of a 5x2 matrix
inline DMatrix2x5 Transpose(const DMatrix5x2 &M){
  return DMatrix2x5(_mm_setr_pd(M(0,0),M(0,1)),_mm_setr_pd(M(1,0),M(1,1)),
		    _mm_setr_pd(M(2,0),M(2,1)),_mm_setr_pd(M(3,0),M(3,1)),
		    _mm_setr_pd(M(4,0),M(4,1)));
}


// Multiply a 5x1 matrix by its transpose 
inline DMatrix5x5 MultiplyTranspose(const DMatrix5x1 &m1){
  ALIGNED_16_BLOCK_WITH_PTR(__m128d, 5, p)
  __m128d &b1=p[0];
  __m128d &b2=p[1];
  __m128d &b3=p[2];
  __m128d &b4=p[3];
  __m128d &b5=p[4];
  b1=_mm_set1_pd(m1(0));
  b2=_mm_set1_pd(m1(1));
  b3=_mm_set1_pd(m1(2));
  b4=_mm_set1_pd(m1(3));
  b5=_mm_set1_pd(m1(4));
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
  ALIGNED_16_BLOCK_WITH_PTR(__m128d, 5, p)
  __m128d &b1=p[0];
  __m128d &b2=p[1];
  __m128d &b3=p[2];
  __m128d &b4=p[3];
  __m128d &b5=p[4];
  b1=_mm_set1_pd(m2(0));
  b2=_mm_set1_pd(m2(1));
  b3=_mm_set1_pd(m2(2));
  b4=_mm_set1_pd(m2(3));
  b5=_mm_set1_pd(m2(4));
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
  ALIGNED_16_BLOCK_WITH_PTR(__m128d, 10, p)
  __m128d &a11=p[0];
  __m128d &a12=p[1];
  __m128d &a13=p[2];
  __m128d &a14=p[3];
  __m128d &a15=p[4];
  __m128d &a21=p[5];
  __m128d &a22=p[6];
  __m128d &a23=p[7];
  __m128d &a24=p[8];
  __m128d &a25=p[9];
  a11=_mm_set1_pd(m2(0,0));
  a12=_mm_set1_pd(m2(0,1)); 
  a13=_mm_set1_pd(m2(0,2));
  a14=_mm_set1_pd(m2(0,3));
  a15=_mm_set1_pd(m2(0,4));
  a21=_mm_set1_pd(m2(1,0));
  a22=_mm_set1_pd(m2(1,1)); 
  a23=_mm_set1_pd(m2(1,2));
  a24=_mm_set1_pd(m2(1,3)); 
  a25=_mm_set1_pd(m2(1,4));
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
