#include <align_16.h>

#ifndef USE_SSE2
#include "DMatrixDSym.h"
// Matrix class without SIMD instructions

class DMatrix4x4{
 public:
  DMatrix4x4(){
    for (unsigned int i=0;i<4;i++){
      for (unsigned int j=0;j<4;j++){
	mA[i][j]=0.;
      }
    }
  } 
  DMatrix4x4(double c11, double c12, double c13, double c14,
	     double c21, double c22, double c23, double c24,
	     double c31, double c32, double c33, double c34,
	     double c41, double c42, double c43, double c44
	     ){
    mA[0][0]=c11;
    mA[0][1]=c12;
    mA[0][2]=c13;
    mA[0][3]=c14;
    mA[1][0]=c21;
    mA[1][1]=c22;
    mA[1][2]=c23; 
    mA[1][3]=c24;
    mA[2][0]=c31;
    mA[2][1]=c32;
    mA[2][2]=c33;    
    mA[2][3]=c34;
    mA[3][0]=c41;
    mA[3][1]=c42;
    mA[3][2]=c43;
    mA[3][3]=c44;

  }
  DMatrix4x4(const DMatrix2x2 &m1,const DMatrix2x2 &m2,
	     const DMatrix2x2 &m3,const DMatrix2x2 &m4){
    mA[0][0]=m1(0,0);
    mA[0][1]=m1(0,1);
    mA[1][0]=m1(1,0);
    mA[1][1]=m1(1,1);
    mA[0][2]=m2(0,0);
    mA[0][3]=m2(0,1);
    mA[1][2]=m2(1,0);
    mA[1][3]=m2(1,1);
    mA[2][0]=m3(0,0);
    mA[2][1]=m3(0,1);
    mA[3][0]=m3(1,0);
    mA[3][1]=m3(1,1);
    mA[2][2]=m4(0,0);
    mA[2][3]=m4(0,1);
    mA[3][2]=m4(1,0);
    mA[3][3]=m4(1,1);
  }

  ~DMatrix4x4(){};

  double &operator() (int row, int col){
    return mA[row][col];
  } 
  double operator() (int row, int col) const{
    return mA[row][col];
  }
  // Assignment
  DMatrix4x4 &operator=(const DMatrix4x4 &m1){
    for (unsigned int i=0;i<4;i++){
      for (unsigned int j=0;j<4;j++){
	mA[i][j]=m1(i,j);
      }
    }    
    return *this;
  }




  // unary minus
  DMatrix4x4 operator-(){
    return DMatrix4x4(-mA[0][0],-mA[0][1],-mA[0][2],-mA[0][3],
		      -mA[1][0],-mA[1][1],-mA[1][2],-mA[1][3],
		      -mA[2][0],-mA[2][1],-mA[2][2],-mA[2][3],
		      -mA[3][0],-mA[3][1],-mA[3][2],-mA[3][3]
		      );
  }
  // Matrix Addition
  DMatrix4x4 operator+(const DMatrix4x4 &m2){
    return DMatrix4x4(mA[0][0]+m2(0,0),mA[0][1]+m2(0,1),mA[0][2]+m2(0,2),mA[0][3]+m2(0,3),
		      mA[1][0]+m2(1,0),mA[1][1]+m2(1,1),mA[1][2]+m2(1,2),mA[1][3]+m2(1,3),
		      mA[2][0]+m2(2,0),mA[2][1]+m2(2,1),mA[2][2]+m2(2,2),mA[2][3]+m2(2,3),
		      mA[3][0]+m2(3,0),mA[3][1]+m2(3,1),mA[3][2]+m2(3,2),mA[3][3]+m2(3,3)
		      );
  }


  // Matrix subtraction
  DMatrix4x4 operator-(const DMatrix4x4 &m2){
    return DMatrix4x4(mA[0][0]-m2(0,0),mA[0][1]-m2(0,1),mA[0][2]-m2(0,2),mA[0][3]-m2(0,3),
		      mA[1][0]-m2(1,0),mA[1][1]-m2(1,1),mA[1][2]-m2(1,2),mA[1][3]-m2(1,3),
		      mA[2][0]-m2(2,0),mA[2][1]-m2(2,1),mA[2][2]-m2(2,2),mA[2][3]-m2(2,3),
		      mA[3][0]-m2(3,0),mA[3][1]-m2(3,1),mA[3][2]-m2(3,2),mA[3][3]-m2(3,3)
		      );
  }

  // Matrix multiplication:  (4x4) x (4x4)
  DMatrix4x4 operator*(const DMatrix4x4 &m2){
    return DMatrix4x4(mA[0][0]*m2(0,0)+mA[0][1]*m2(1,0)+mA[0][2]*m2(2,0)+mA[0][3]*m2(3,0),
		      mA[0][0]*m2(0,1)+mA[0][1]*m2(1,1)+mA[0][2]*m2(2,1)+mA[0][3]*m2(3,1),
		      mA[0][0]*m2(0,2)+mA[0][1]*m2(1,2)+mA[0][2]*m2(2,2)+mA[0][3]*m2(3,2),
		      mA[0][0]*m2(0,3)+mA[0][1]*m2(1,3)+mA[0][2]*m2(2,3)+mA[0][3]*m2(3,3),

		      mA[1][0]*m2(0,0)+mA[1][1]*m2(1,0)+mA[1][2]*m2(2,0)+mA[1][3]*m2(3,0),
		      mA[1][0]*m2(0,1)+mA[1][1]*m2(1,1)+mA[1][2]*m2(2,1)+mA[1][3]*m2(3,1),
		      mA[1][0]*m2(0,2)+mA[1][1]*m2(1,2)+mA[1][2]*m2(2,2)+mA[1][3]*m2(3,2),
		      mA[1][0]*m2(0,3)+mA[1][1]*m2(1,3)+mA[1][2]*m2(2,3)+mA[1][3]*m2(3,3),

		      mA[2][0]*m2(0,0)+mA[2][1]*m2(1,0)+mA[2][2]*m2(2,0)+mA[2][3]*m2(3,0),
		      mA[2][0]*m2(0,1)+mA[2][1]*m2(1,1)+mA[2][2]*m2(2,1)+mA[2][3]*m2(3,1),
		      mA[2][0]*m2(0,2)+mA[2][1]*m2(1,2)+mA[2][2]*m2(2,2)+mA[2][3]*m2(3,2),
		      mA[2][0]*m2(0,3)+mA[2][1]*m2(1,3)+mA[2][2]*m2(2,3)+mA[2][3]*m2(3,3),

		      mA[3][0]*m2(0,0)+mA[3][1]*m2(1,0)+mA[3][2]*m2(2,0)+mA[3][3]*m2(3,0),
		      mA[3][0]*m2(0,1)+mA[3][1]*m2(1,1)+mA[3][2]*m2(2,1)+mA[3][3]*m2(3,1),
		      mA[3][0]*m2(0,2)+mA[3][1]*m2(1,2)+mA[3][2]*m2(2,2)+mA[3][3]*m2(3,2),
		      mA[3][0]*m2(0,3)+mA[3][1]*m2(1,3)+mA[3][2]*m2(2,3)+mA[3][3]*m2(3,3)
		      );
  }

  // Matrix multiplication:  (4x4) x (4x2)
  DMatrix4x2 operator*(const DMatrix4x2 &m2){
    return DMatrix4x2(mA[0][0]*m2(0,0)+mA[0][1]*m2(1,0)+mA[0][2]*m2(2,0)+mA[0][3]*m2(3,0),
		      mA[0][0]*m2(0,1)+mA[0][1]*m2(1,1)+mA[0][2]*m2(2,1)+mA[0][3]*m2(3,1),

		      mA[1][0]*m2(0,0)+mA[1][1]*m2(1,0)+mA[1][2]*m2(2,0)+mA[1][3]*m2(3,0),
		      mA[1][0]*m2(0,1)+mA[1][1]*m2(1,1)+mA[1][2]*m2(2,1)+mA[1][3]*m2(3,1),

		      mA[2][0]*m2(0,0)+mA[2][1]*m2(1,0)+mA[2][2]*m2(2,0)+mA[2][3]*m2(3,0),
		      mA[2][0]*m2(0,1)+mA[2][1]*m2(1,1)+mA[2][2]*m2(2,1)+mA[2][3]*m2(3,1),  

		      mA[3][0]*m2(0,0)+mA[3][1]*m2(1,0)+mA[3][2]*m2(2,0)+mA[3][3]*m2(3,0),
		      mA[3][0]*m2(0,1)+mA[3][1]*m2(1,1)+mA[3][2]*m2(2,1)+mA[3][3]*m2(3,1)
		      
		      );
  } 
  
  
  // Matrix multiplication:  (4x4) x (4x1)
  DMatrix4x1 operator*(const DMatrix4x1 &m2){
    return DMatrix4x1(mA[0][0]*m2(0)+mA[0][1]*m2(1)+mA[0][2]*m2(2)+mA[0][3]*m2(3),
		      mA[1][0]*m2(0)+mA[1][1]*m2(1)+mA[1][2]*m2(2)+mA[1][3]*m2(3),
		      mA[2][0]*m2(0)+mA[2][1]*m2(1)+mA[2][2]*m2(2)+mA[2][3]*m2(3),
		      mA[3][0]*m2(0)+mA[3][1]*m2(1)+mA[3][2]*m2(2)+mA[3][3]*m2(3)
		      );
  } 
  



  DMatrix4x4 Invert(){
    DMatrix2x2 F(mA[0][0],mA[0][1],mA[1][0],mA[1][1]);
    DMatrix2x2 Finv=F.Invert();
    DMatrix2x2 G(mA[0][2],mA[0][3],mA[1][2],mA[1][3]); 
    DMatrix2x2 H(mA[2][0],mA[2][1],mA[3][0],mA[3][1]);
    DMatrix2x2 J(mA[2][2],mA[2][3],mA[3][2],mA[3][3]);
    DMatrix2x2 Jinv=J.Invert();
    DMatrix2x2 FF=(F-G*Jinv*H).Invert();
    DMatrix2x2 JJ=(J-H*Finv*G).Invert();
    return DMatrix4x4(FF,-FF*G*Jinv,-JJ*H*Finv,JJ);     
  }
  
  // Find the transpose of this matrix
  DMatrix4x4 Transpose(){
    DMatrix4x4 temp;
    for (unsigned int i=0;i<4;i++){
      for (unsigned int j=0;j<4;j++){
	temp(i,j)=mA[j][i];
      }
    }
    return temp;
  }

  DMatrixDSym GetSub(unsigned int lowerBound, unsigned int upperBound){
     if (upperBound >= lowerBound) return DMatrixDSym();
     DMatrixDSym subMatrix(upperBound - lowerBound);
     for (unsigned int i=lowerBound; i <= upperBound; i++){
        for (unsigned int j=lowerBound; j <= upperBound; j++){
           subMatrix(i,j) = mA[i][j];
        }
     }
     return subMatrix;
  }

  bool IsPosDef(){
     if(mA[0][0] > 0. &&
           GetSub(0,1).Determinant() > 0. && GetSub(0,2).Determinant() > 0. &&
           GetSub(0,3).Determinant() > 0.)
        return true;
     else return false;
  }

  void Print(){
     cout << "DMatrix4x4:" <<endl;
     cout << "     |      0    |      1    |      2    |      3    |" <<endl;
     cout << "------------------------------------------------------" <<endl;

     for (unsigned int i=0;i<4;i++){
        cout << "   " << i <<" |";
        for (unsigned int j=0;j<4;j++){
           cout <<setw(11)<<setprecision(4)<< mA[i][j] <<" "; 
        } 
        cout << endl;
     }      
  }


 private:
  double mA[4][4];

};

#else

class DMatrix4x4{
   public:
      DMatrix4x4()
         : mA( ALIGNED_16_BLOCK_PTR(union dvec, 4, mA) )
      {
         for (unsigned int i=0;i<4;i++){
            mA[i].v[0]=_mm_setzero_pd();
            mA[i].v[1]=_mm_setzero_pd();
         }
      }
      DMatrix4x4(double c11, double c12, double c13, double c14,
            double c21, double c22, double c23, double c24,
            double c31, double c32, double c33, double c34,
            double c41, double c42, double c43, double c44)
         : mA( ALIGNED_16_BLOCK_PTR(union dvec, 4, mA) )
      {
         mA[0].v[0]=_mm_setr_pd(c11,c21);
         mA[0].v[1]=_mm_setr_pd(c31,c41);
         mA[1].v[0]=_mm_setr_pd(c12,c22);
         mA[1].v[1]=_mm_setr_pd(c32,c42);
         mA[2].v[0]=_mm_setr_pd(c13,c23);
         mA[2].v[1]=_mm_setr_pd(c33,c43);
         mA[3].v[0]=_mm_setr_pd(c14,c24);
         mA[3].v[1]=_mm_setr_pd(c34,c44);	
      }

      DMatrix4x4(__m128d m11, __m128d m12, __m128d m13, __m128d m14,
            __m128d m21, __m128d m22, __m128d m23, __m128d m24)
         : mA( ALIGNED_16_BLOCK_PTR(union dvec, 4, mA) )
      {
         mA[0].v[0]=m11;
         mA[1].v[0]=m12;
         mA[2].v[0]=m13;
         mA[3].v[0]=m14; 
         mA[0].v[1]=m21;
         mA[1].v[1]=m22;
         mA[2].v[1]=m23;
         mA[3].v[1]=m24;
      }
      DMatrix4x4(const DMatrix2x2 &m1,const DMatrix2x2 &m2,
            const DMatrix2x2 &m3,const DMatrix2x2 &m4)
         : mA( ALIGNED_16_BLOCK_PTR(union dvec, 4, mA) )
      {
         mA[0].v[0]=m1.GetV(0);
         mA[1].v[0]=m1.GetV(1);
         mA[0].v[1]=m3.GetV(0);
         mA[1].v[1]=m3.GetV(1);
         mA[2].v[0]=m2.GetV(0);
         mA[3].v[0]=m2.GetV(1);
         mA[2].v[1]=m4.GetV(0);
         mA[3].v[1]=m4.GetV(1);
      }
      DMatrix4x4(const DMatrix4x4& dm)
         : mA( ALIGNED_16_BLOCK_PTR(union dvec, 4, mA) )
      {
         mA[0].v[0]=dm.mA[0].v[0];
         mA[1].v[0]=dm.mA[1].v[0];
         mA[2].v[0]=dm.mA[2].v[0];
         mA[3].v[0]=dm.mA[3].v[0];
         mA[0].v[1]=dm.mA[0].v[1];
         mA[1].v[1]=dm.mA[1].v[1];
         mA[2].v[1]=dm.mA[2].v[1];
         mA[3].v[1]=dm.mA[3].v[1];
      }
      ~DMatrix4x4(){}; 

      __m128d GetV(int pair, int col) const{
         return mA[col].v[pair];
      }
      void SetV(int pair, int col, __m128d v){
         mA[col].v[pair]=v;
      }

      double &operator() (int row, int col){
         return mA[col].d[row];
      } 
      double operator() (int row, int col) const{
         return mA[col].d[row];
      }

      DMatrix4x4 operator-(){
         ALIGNED_16_BLOCK_WITH_PTR(__m128d, 1, p)
            __m128d &zero=p[0];
         zero=_mm_setzero_pd();
         return DMatrix4x4(_mm_sub_pd(zero,GetV(0,0)),
               _mm_sub_pd(zero,GetV(0,1)),
               _mm_sub_pd(zero,GetV(0,2)),
               _mm_sub_pd(zero,GetV(0,3)),
               _mm_sub_pd(zero,GetV(1,0)),
               _mm_sub_pd(zero,GetV(1,1)),
               _mm_sub_pd(zero,GetV(1,2)),
               _mm_sub_pd(zero,GetV(1,3)));
      }

      DMatrix4x4 operator-(const DMatrix4x4 &m2){
         return DMatrix4x4(_mm_sub_pd(GetV(0,0),m2.GetV(0,0)),
               _mm_sub_pd(GetV(0,1),m2.GetV(0,1)),
               _mm_sub_pd(GetV(0,2),m2.GetV(0,2)),
               _mm_sub_pd(GetV(0,3),m2.GetV(0,3)),
               _mm_sub_pd(GetV(1,0),m2.GetV(1,0)),
               _mm_sub_pd(GetV(1,1),m2.GetV(1,1)),
               _mm_sub_pd(GetV(1,2),m2.GetV(1,2)),
               _mm_sub_pd(GetV(1,3),m2.GetV(1,3)));
      }

      DMatrix4x4 operator+(const DMatrix4x4 &m2){
         return DMatrix4x4(_mm_add_pd(GetV(0,0),m2.GetV(0,0)),
               _mm_add_pd(GetV(0,1),m2.GetV(0,1)),
               _mm_add_pd(GetV(0,2),m2.GetV(0,2)),
               _mm_add_pd(GetV(0,3),m2.GetV(0,3)),
               _mm_add_pd(GetV(1,0),m2.GetV(1,0)),
               _mm_add_pd(GetV(1,1),m2.GetV(1,1)),
               _mm_add_pd(GetV(1,2),m2.GetV(1,2)),
               _mm_add_pd(GetV(1,3),m2.GetV(1,3)));
      }


      DMatrix4x2 operator*(const DMatrix4x2 &m2){
         ALIGNED_16_BLOCK_WITH_PTR(__m128d, 8, p);
         __m128d &m11=p[0];
         __m128d &m12=p[1];
         __m128d &m21=p[2];
         __m128d &m22=p[3];
         __m128d &m31=p[4];
         __m128d &m32=p[5];
         __m128d &m41=p[6];
         __m128d &m42=p[7];
         m11=_mm_set1_pd(m2(0,0));
         m12=_mm_set1_pd(m2(0,1)); 
         m21=_mm_set1_pd(m2(1,0));
         m22=_mm_set1_pd(m2(1,1));  
         m31=_mm_set1_pd(m2(2,0));
         m32=_mm_set1_pd(m2(2,1)); 
         m41=_mm_set1_pd(m2(3,0));
         m42=_mm_set1_pd(m2(3,1));
         return 
            DMatrix4x2(_mm_add_pd(_mm_mul_pd(GetV(0,0),m11),
                     _mm_add_pd(_mm_mul_pd(GetV(0,1),m21),
                        _mm_add_pd(_mm_mul_pd(GetV(0,2),m31),
                           _mm_mul_pd(GetV(0,3),m41)))),
                  _mm_add_pd(_mm_mul_pd(GetV(0,0),m12),
                     _mm_add_pd(_mm_mul_pd(GetV(0,1),m22),
                        _mm_add_pd(_mm_mul_pd(GetV(0,2),m32),
                           _mm_mul_pd(GetV(0,3),m42)))),
                  _mm_add_pd(_mm_mul_pd(GetV(1,0),m11),
                     _mm_add_pd(_mm_mul_pd(GetV(1,1),m21),
                        _mm_add_pd(_mm_mul_pd(GetV(1,2),m31),
                           _mm_mul_pd(GetV(1,3),m41)))),
                  _mm_add_pd(_mm_mul_pd(GetV(1,0),m12),
                     _mm_add_pd(_mm_mul_pd(GetV(1,1),m22),
                        _mm_add_pd(_mm_mul_pd(GetV(1,2),m32),
                           _mm_mul_pd(GetV(1,3),m42)))));

      }


      DMatrix4x4 &operator=(const DMatrix4x4 &m1){
         for (unsigned int i=0;i<4;i++){
            mA[i].v[0]=m1.GetV(0,i);
            mA[i].v[1]=m1.GetV(1,i);
         }    
         return *this;
      }

      DMatrix4x4 Invert(){
         DMatrix2x2 F(GetV(0,0),GetV(0,1));
         DMatrix2x2 Finv=F.Invert();
         DMatrix2x2 G(GetV(0,2),GetV(0,3)); 
         DMatrix2x2 H(GetV(1,0),GetV(1,1));
         DMatrix2x2 J(GetV(1,2),GetV(1,3));
         DMatrix2x2 Jinv=J.Invert();
         DMatrix2x2 FF=(F-G*Jinv*H).Invert();
         DMatrix2x2 JJ=(J-H*Finv*G).Invert();
         return DMatrix4x4(FF,-FF*G*Jinv,-JJ*H*Finv,JJ);     
      }

      // Find the transpose of this matrix
      DMatrix4x4 Transpose(){
#define SWAP(i,j,k,m) _mm_setr_pd(mA[(j)].d[(i)],mA[(m)].d[(k)])

         return DMatrix4x4(SWAP(0,0,0,1),SWAP(1,0,1,1),SWAP(2,0,2,1),SWAP(3,0,3,1),
               SWAP(0,2,0,3),SWAP(1,2,1,3),SWAP(2,2,2,3),SWAP(3,2,3,3));

      }

      // Matrix multiplication:  (4x4) x (4x1)
      DMatrix4x1 operator*(const DMatrix4x1 &m2){
         ALIGNED_16_BLOCK_WITH_PTR(__m128d, 4, p)
            __m128d &a1=p[0];
         __m128d &a2=p[1];
         __m128d &a3=p[2];
         __m128d &a4=p[3];
         a1=_mm_set1_pd(m2(0));
         a2=_mm_set1_pd(m2(1));
         a3=_mm_set1_pd(m2(2));
         a4=_mm_set1_pd(m2(3));
         return 
            DMatrix4x1(_mm_add_pd(_mm_mul_pd(GetV(0,0),a1),
                     _mm_add_pd(_mm_mul_pd(GetV(0,1),a2),
                        _mm_add_pd(_mm_mul_pd(GetV(0,2),a3),
                           _mm_mul_pd(GetV(0,3),a4)))),
                  _mm_add_pd(_mm_mul_pd(GetV(1,0),a1),
                     _mm_add_pd(_mm_mul_pd(GetV(1,1),a2),
                        _mm_add_pd(_mm_mul_pd(GetV(1,2),a3),
                           _mm_mul_pd(GetV(1,3),a4)))));
      }

      // Matrix multiplication: (4x4) x (4x4)
      DMatrix4x4 operator*(const DMatrix4x4 &m2){
         struct dvec1{
            __m128d v[4][4];
         };
         ALIGNED_16_BLOCK_WITH_PTR(struct dvec1, 1, p)
            struct dvec1 &temp=p[0];
         for (unsigned int i=0;i<4;i++){
            for (unsigned int j=0;j<4;j++){
               temp.v[i][j]=_mm_set1_pd(m2(i,j));
            }
         }
         // Preprocessor macro for multiplying two __m128d elements together
#define MUL(i,j,k) _mm_mul_pd(GetV((i),(j)),temp.v[(j)][(k)])

         return DMatrix4x4(_mm_add_pd(MUL(0,0,0),
                  _mm_add_pd(MUL(0,1,0),
                     _mm_add_pd(MUL(0,2,0),
                        MUL(0,3,0)))),
               _mm_add_pd(MUL(0,0,1),
                  _mm_add_pd(MUL(0,1,1),
                     _mm_add_pd(MUL(0,2,1),
                        MUL(0,3,1)))),
               _mm_add_pd(MUL(0,0,2),
                  _mm_add_pd(MUL(0,1,2),
                     _mm_add_pd(MUL(0,2,2),
                        MUL(0,3,2)))),
               _mm_add_pd(MUL(0,0,3),
                  _mm_add_pd(MUL(0,1,3),
                     _mm_add_pd(MUL(0,2,3),
                        MUL(0,3,3)))),  
               _mm_add_pd(MUL(1,0,0),
                  _mm_add_pd(MUL(1,1,0),
                     _mm_add_pd(MUL(1,2,0),
                        MUL(1,3,0)))),
               _mm_add_pd(MUL(1,0,1),
                  _mm_add_pd(MUL(1,1,1),
                     _mm_add_pd(MUL(1,2,1),
                        MUL(1,3,1)))),
               _mm_add_pd(MUL(1,0,2),
                     _mm_add_pd(MUL(1,1,2),
                        _mm_add_pd(MUL(1,2,2),
                           MUL(1,3,2)))),
               _mm_add_pd(MUL(1,0,3),
                     _mm_add_pd(MUL(1,1,3),
                        _mm_add_pd(MUL(1,2,3),
                           MUL(1,3,3)))));  
      }


      void Print(){
         cout << "DMatrix4x4:" <<endl;
         cout << "         |      0    |      1    |      2    |      3    |" <<endl;
         cout << "----------------------------------------------------------" <<endl;

         for (unsigned int i=0;i<4;i++){
            for (unsigned int j=0;j<4;j++){
               cout << mA[j].d[i] <<" "; 
            } 
            cout << endl;
         }      
      }

   private:
      union dvec{
         __m128d v[2];
         double d[4];
      };
      ALIGNED_16_BLOCK(union dvec, 4, mA)
};
#endif
