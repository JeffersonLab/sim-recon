#include <align_16.h>

#ifndef USE_SSE2

#include "DMatrixDSym.h"
// Matrix class without SIMD instructions

class DMatrix5x5{
   public:
      DMatrix5x5(){
         for (unsigned int i=0;i<5;i++)
            for (unsigned int j=0;j<5;j++)
               mA[i][j]=0.;
      }
      // Copy constructor
      DMatrix5x5(const DMatrix5x5 &m2){
         for (unsigned int i=0;i<5;i++){
            for (unsigned int j=0;j<5;j++){
               mA[i][j]=m2(i,j);
            }
         }
      }
      ~DMatrix5x5(){};

      // Constructor using block matrices from matrix inversion 
      DMatrix5x5(const DMatrix2x2 &A,const DMatrix2x3 &B,
            const DMatrix3x2 &C,const DMatrix3x3 &D){
         mA[0][0]=A(0,0);
         mA[0][1]=A(0,1);
         mA[1][0]=A(1,0);
         mA[1][1]=A(1,1);
         for (unsigned int i=0;i<2;i++){
            for (unsigned int j=0;j<3;j++){
               mA[i][j+2]=B(i,j);
            }
         }
         for (unsigned int i=0;i<3;i++){
            for (unsigned int j=0;j<2;j++){
               mA[i+2][j]=C(i,j);
               mA[i+2][j+2]=D(i,j);
            }
            mA[i+2][4]=D(i,2);
         }
      }


      // Access by indices
      double &operator() (int row, int col){
         return mA[row][col];
      } 
      double operator() (int row, int col) const{
         return mA[row][col];
      }

      // Assignment operator
      DMatrix5x5 &operator=(const DMatrix5x5 &m1){
         for (unsigned int i=0;i<5;i++){
            for (unsigned int j=0;j<5;j++){
               mA[i][j]=m1(i,j);
            }
         }    
         return *this;
      }

      // Find the transpose of this matrix
      DMatrix5x5 Transpose(){
         DMatrix5x5 temp;
         for (unsigned int i=0;i<5;i++){
            for (unsigned int j=0;j<5;j++){
               temp(i,j)=mA[j][i];
            }
         }
         return temp;
      }
      // Matrix addition
      DMatrix5x5 operator+(const DMatrix5x5 &m2) const{
         DMatrix5x5 temp;

         for (unsigned int row=0;row<5;row++){
            for (unsigned int col=0;col<5;col++){
               temp(row,col)=mA[row][col]+m2(row,col);
            }
         }
         return temp;
      }
      DMatrix5x5 &operator+=(const DMatrix5x5 &m2){
         for (unsigned int i=0;i<5;i++){
            for (unsigned int j=0;j<5;j++){
               mA[i][j]+=m2(i,j);
            }
         }
         return *this;
      }

      // Method for adding symmetric matrices
      DMatrix5x5 AddSym(const DMatrix5x5 &m2) const{
         DMatrix5x5 temp;
         for (unsigned int i=0;i<5;i++){
            for (unsigned int j=i;j<5;j++){
               temp(i,j)=mA[i][j]+m2(i,j);
               temp(j,i)=temp(i,j);
            }
         }
         return temp;
      } 

      // Method for subtracting  symmetric matrices
      DMatrix5x5 SubSym(const DMatrix5x5 &m2) const{
         DMatrix5x5 temp;
         for (unsigned int i=0;i<5;i++){
            for (unsigned int j=i;j<5;j++){
               temp(i,j)=mA[i][j]-m2(i,j);
               temp(j,i)=temp(i,j);
            }
         }
         return temp;
      }

      // Matrix subtraction
      DMatrix5x5 operator-(const DMatrix5x5 &m2) const{
         DMatrix5x5 temp;

         for (unsigned int row=0;row<5;row++){
            for (unsigned int col=0;col<5;col++){
               temp(row,col)=mA[row][col]-m2(row,col);
            }
         }
         return temp;
      }
      DMatrix5x5 &operator-=(const DMatrix5x5 &m2){
         for (unsigned int i=0;i<5;i++){
            for (unsigned int j=0;j<5;j++){
               mA[i][j]-=m2(i,j);
            }
         }
         return *this;
      }
      // Scale a matrix by a double
      DMatrix5x5 &operator*=(const double c){
         for (unsigned int i=0;i<5;i++){
            for (unsigned int j=0;j<5;j++){
               mA[i][j]*=c;
            }
         }
         return *this;
      }

      // Matrix multiplication:  (5x5) x (5x1)
      DMatrix5x1 operator*(const DMatrix5x1 &m2){
         return DMatrix5x1(mA[0][0]*m2(0)+mA[0][1]*m2(1)+mA[0][2]*m2(2)
               +mA[0][3]*m2(3)+mA[0][4]*m2(4),
               mA[1][0]*m2(0)+mA[1][1]*m2(1)+mA[1][2]*m2(2)
               +mA[1][3]*m2(3)+mA[1][4]*m2(4),
               mA[2][0]*m2(0)+mA[2][1]*m2(1)+mA[2][2]*m2(2)
               +mA[2][3]*m2(3)+mA[2][4]*m2(4),
               mA[3][0]*m2(0)+mA[3][1]*m2(1)+mA[3][2]*m2(2)
               +mA[3][3]*m2(3)+mA[3][4]*m2(4),  
               mA[4][0]*m2(0)+mA[4][1]*m2(1)+mA[4][2]*m2(2)
               +mA[4][3]*m2(3)+mA[4][4]*m2(4));

      }

      // Matrix multiplication:  (5x5) x (5x2)
      DMatrix5x2 operator*(const DMatrix5x2 &m2){
         return DMatrix5x2(mA[0][0]*m2(0,0)+mA[0][1]*m2(1,0)+mA[0][2]*m2(2,0)
               +mA[0][3]*m2(3,0)+mA[0][4]*m2(4,0),
               mA[0][0]*m2(0,1)+mA[0][1]*m2(1,1)+mA[0][2]*m2(2,1)
               +mA[0][3]*m2(3,1)+mA[0][4]*m2(4,1),

               mA[1][0]*m2(0,0)+mA[1][1]*m2(1,0)+mA[1][2]*m2(2,0)
               +mA[1][3]*m2(3,0)+mA[1][4]*m2(4,0),
               mA[1][0]*m2(0,1)+mA[1][1]*m2(1,1)+mA[1][2]*m2(2,1)
               +mA[1][3]*m2(3,1)+mA[1][4]*m2(4,1),

               mA[2][0]*m2(0,0)+mA[2][1]*m2(1,0)+mA[2][2]*m2(2,0)
               +mA[2][3]*m2(3,0)+mA[2][4]*m2(4,0),
               mA[2][0]*m2(0,1)+mA[2][1]*m2(1,1)+mA[2][2]*m2(2,1)
               +mA[2][3]*m2(3,1)+mA[2][4]*m2(4,1),

               mA[3][0]*m2(0,0)+mA[3][1]*m2(1,0)+mA[3][2]*m2(2,0)
               +mA[3][3]*m2(3,0)+mA[3][4]*m2(4,0),  
               mA[3][0]*m2(0,1)+mA[3][1]*m2(1,1)+mA[3][2]*m2(2,1)
               +mA[3][3]*m2(3,1)+mA[3][4]*m2(4,1),  

               mA[4][0]*m2(0,0)+mA[4][1]*m2(1,0)+mA[4][2]*m2(2,0)
                  +mA[4][3]*m2(3,0)+mA[4][4]*m2(4,0),
               mA[4][0]*m2(0,1)+mA[4][1]*m2(1,1)+mA[4][2]*m2(2,1)
                  +mA[4][3]*m2(3,1)+mA[4][4]*m2(4,1)
                  );


      }
      // Matrix multiplication: (5x5) x (5x5)
      DMatrix5x5 operator*(const DMatrix5x5 &m2){
         DMatrix5x5 temp;

         for (unsigned int j=0;j<5;j++){
            for (unsigned int i=0;i<5;i++){
               double temp2=0.;
               for (unsigned int k=0;k<5;k++){
                  temp2+=mA[j][k]*m2(k,i);
               }
               temp(j,i)=temp2;
            }
         }
         return temp;
      }

      // Zero out the matrix
      DMatrix5x5 Zero(){
         for (unsigned int i=0;i<5;i++){
            for (unsigned int j=0;j<5;j++){
               mA[i][j]=0.;
            }
         }
         return *this;
      }


      // Matrix inversion by blocks for a symmetric matrix
      DMatrix5x5 InvertSym(){
         DMatrix2x2 A(mA[0][0],mA[0][1],mA[1][0],mA[1][1]);
         DMatrix3x2 C(mA[2][0],mA[2][1],mA[3][0],mA[3][1],mA[4][0],mA[4][1]);
         DMatrix3x2 CAinv=C*A.Invert();
         DMatrix3x3 D(mA[2][2],mA[2][3],mA[2][4],mA[3][2],mA[3][3],mA[3][4],
               mA[4][2],mA[4][3],mA[4][4]);
         DMatrix2x3 B(mA[0][2],mA[0][3],mA[0][4],mA[1][2],mA[1][3],mA[1][4]);
         DMatrix2x3 BDinv=B*D.InvertSym();
         DMatrix2x2 AA=(A-BDinv*C).Invert();
         DMatrix3x3 DD=(D-CAinv*B).InvertSym();
         return DMatrix5x5(AA,-AA*BDinv,-DD*CAinv,DD);
      }


      // The following code performs the matrix operation ABA^T, where B is a symmetric matrix.   The end 
      // result is also a symmetric matrix.
      DMatrix5x5 SandwichMultiply(const DMatrix5x5 &A){
         DMatrix5x5 temp;
         double sums[5];
         // First row/column of ACA^T
         sums[0]= mA[0][0]*A(0,0);
         sums[0]+=mA[0][1]*A(0,1); 
         sums[0]+=mA[0][2]*A(0,2);
         sums[0]+=mA[0][3]*A(0,3);
         sums[0]+=mA[0][4]*A(0,4);

         sums[1]= mA[1][0]*A(0,0);
         sums[1]+=mA[1][1]*A(0,1); 
         sums[1]+=mA[1][2]*A(0,2);
         sums[1]+=mA[1][3]*A(0,3);
         sums[1]+=mA[1][4]*A(0,4);

         sums[2]= mA[2][0]*A(0,0);
         sums[2]+=mA[2][1]*A(0,1); 
         sums[2]+=mA[2][2]*A(0,2);
         sums[2]+=mA[2][3]*A(0,3);
         sums[2]+=mA[2][4]*A(0,4);

         sums[3]= mA[3][0]*A(0,0);
         sums[3]+=mA[3][1]*A(0,1); 
         sums[3]+=mA[3][2]*A(0,2);
         sums[3]+=mA[3][3]*A(0,3);
         sums[3]+=mA[3][4]*A(0,4);

         sums[4]= mA[4][0]*A(0,0);
         sums[4]+=mA[4][1]*A(0,1); 
         sums[4]+=mA[4][2]*A(0,2);
         sums[4]+=mA[4][3]*A(0,3);
         sums[4]+=mA[4][4]*A(0,4);

         temp(0,0)+=A(0,0)*sums[0];
         temp(0,0)+=A(0,1)*sums[1];
         temp(0,0)+=A(0,2)*sums[2];
         temp(0,0)+=A(0,3)*sums[3];
         temp(0,0)+=A(0,4)*sums[4];

         temp(1,0)+=A(1,0)*sums[0];
         temp(1,0)+=A(1,1)*sums[1];
         temp(1,0)+=A(1,2)*sums[2];
         temp(1,0)+=A(1,3)*sums[3];
         temp(1,0)+=A(1,4)*sums[4];

         temp(2,0)+=A(2,0)*sums[0];
         temp(2,0)+=A(2,1)*sums[1];
         temp(2,0)+=A(2,2)*sums[2];
         temp(2,0)+=A(2,3)*sums[3];
         temp(2,0)+=A(2,4)*sums[4];

         temp(3,0)+=A(3,0)*sums[0];
         temp(3,0)+=A(3,1)*sums[1];
         temp(3,0)+=A(3,2)*sums[2];
         temp(3,0)+=A(3,3)*sums[3];
         temp(3,0)+=A(3,4)*sums[4];

         temp(4,0)+=A(4,0)*sums[0];
         temp(4,0)+=A(4,1)*sums[1];
         temp(4,0)+=A(4,2)*sums[2];
         temp(4,0)+=A(4,3)*sums[3];
         temp(4,0)+=A(4,4)*sums[4];

         temp(0,1)=temp(1,0);
         temp(0,2)=temp(2,0);
         temp(0,3)=temp(3,0);
         temp(0,4)=temp(4,0);

         // Second row/column of ACA^T
         sums[0] =mA[0][0]*A(1,0);
         sums[0]+=mA[0][1]*A(1,1);
         sums[0]+=mA[0][2]*A(1,2);
         sums[0]+=mA[0][3]*A(1,3);
         sums[0]+=mA[0][4]*A(1,4);

         sums[1]=mA[1][0]*A(1,0);
         sums[1]+=mA[1][1]*A(1,1);
         sums[1]+=mA[1][2]*A(1,2);
         sums[1]+=mA[1][3]*A(1,3);
         sums[1]+=mA[1][4]*A(1,4);

         sums[2]=mA[2][0]*A(1,0);
         sums[2]+=mA[2][1]*A(1,1);
         sums[2]+=mA[2][2]*A(1,2);
         sums[2]+=mA[2][3]*A(1,3);
         sums[2]+=mA[2][4]*A(1,4);

         sums[3]=mA[3][0]*A(1,0);
         sums[3]+=mA[3][1]*A(1,1);
         sums[3]+=mA[3][2]*A(1,2);
         sums[3]+=mA[3][3]*A(1,3);
         sums[3]+=mA[3][4]*A(1,4);

         sums[4]=mA[4][0]*A(1,0);
         sums[4]+=mA[4][1]*A(1,1);
         sums[4]+=mA[4][2]*A(1,2);
         sums[4]+=mA[4][3]*A(1,3);
         sums[4]+=mA[4][4]*A(1,4);

         temp(1,1)+=A(1,0)*sums[0];
         temp(1,1)+=A(1,1)*sums[1];
         temp(1,1)+=A(1,2)*sums[2];
         temp(1,1)+=A(1,3)*sums[3];
         temp(1,1)+=A(1,4)*sums[4];

         temp(2,1)+=A(2,0)*sums[0];
         temp(2,1)+=A(2,1)*sums[1];
         temp(2,1)+=A(2,2)*sums[2];
         temp(2,1)+=A(2,3)*sums[3];
         temp(2,1)+=A(2,4)*sums[4];

         temp(3,1)+=A(3,0)*sums[0];
         temp(3,1)+=A(3,1)*sums[1];
         temp(3,1)+=A(3,2)*sums[2];
         temp(3,1)+=A(3,3)*sums[3];
         temp(3,1)+=A(3,4)*sums[4];

         temp(4,1)+=A(4,0)*sums[0];
         temp(4,1)+=A(4,1)*sums[1];
         temp(4,1)+=A(4,2)*sums[2];
         temp(4,1)+=A(4,3)*sums[3];
         temp(4,1)+=A(4,4)*sums[4]; 

         temp(1,2)=temp(2,1);
         temp(1,3)=temp(3,1);
         temp(1,4)=temp(4,1);

         // Third row/column of ACA^T
         sums[0] =mA[0][0]*A(2,0);
         sums[0]+=mA[0][1]*A(2,1);
         sums[0]+=mA[0][2]*A(2,2);
         sums[0]+=mA[0][3]*A(2,3);
         sums[0]+=mA[0][4]*A(2,4);             

         sums[1] =mA[1][0]*A(2,0);
         sums[1]+=mA[1][1]*A(2,1);
         sums[1]+=mA[1][2]*A(2,2);
         sums[1]+=mA[1][3]*A(2,3);
         sums[1]+=mA[1][4]*A(2,4);

         sums[2] =mA[2][0]*A(2,0);
         sums[2]+=mA[2][1]*A(2,1);
         sums[2]+=mA[2][2]*A(2,2);
         sums[2]+=mA[2][3]*A(2,3);
         sums[2]+=mA[2][4]*A(2,4);

         sums[3] =mA[3][0]*A(2,0);
         sums[3]+=mA[3][1]*A(2,1);
         sums[3]+=mA[3][2]*A(2,2);
         sums[3]+=mA[3][3]*A(2,3);
         sums[3]+=mA[3][4]*A(2,4);

         sums[4] =mA[4][0]*A(2,0);
         sums[4]+=mA[4][1]*A(2,1);
         sums[4]+=mA[4][2]*A(2,2);
         sums[4]+=mA[4][3]*A(2,3);
         sums[4]+=mA[4][4]*A(2,4);

         temp(2,2)+=A(2,0)*sums[0]; 
         temp(2,2)+=A(2,1)*sums[1];
         temp(2,2)+=A(2,2)*sums[2]; 
         temp(2,2)+=A(2,3)*sums[3]; 
         temp(2,2)+=A(2,4)*sums[4];

         temp(3,2)+=A(3,0)*sums[0]; 
         temp(3,2)+=A(3,1)*sums[1];
         temp(3,2)+=A(3,2)*sums[2]; 
         temp(3,2)+=A(3,3)*sums[3]; 
         temp(3,2)+=A(3,4)*sums[4];

         temp(4,2)+=A(4,0)*sums[0]; 
         temp(4,2)+=A(4,1)*sums[1];
         temp(4,2)+=A(4,2)*sums[2]; 
         temp(4,2)+=A(4,3)*sums[3]; 
         temp(4,2)+=A(4,4)*sums[4];

         temp(2,3)=temp(3,2);
         temp(2,4)=temp(4,2);

         // Fourth row/column of ACA^T
         sums[0] =mA[0][0]*A(3,0);
         sums[0]+=mA[0][1]*A(3,1);
         sums[0]+=mA[0][2]*A(3,2);
         sums[0]+=mA[0][3]*A(3,3);
         sums[0]+=mA[0][4]*A(3,4);         

         sums[1] =mA[1][0]*A(3,0);
         sums[1]+=mA[1][1]*A(3,1);
         sums[1]+=mA[1][2]*A(3,2);
         sums[1]+=mA[1][3]*A(3,3);
         sums[1]+=mA[1][4]*A(3,4);

         sums[2] =mA[2][0]*A(3,0);
         sums[2]+=mA[2][1]*A(3,1);
         sums[2]+=mA[2][2]*A(3,2);
         sums[2]+=mA[2][3]*A(3,3);
         sums[2]+=mA[2][4]*A(3,4);

         sums[3] =mA[3][0]*A(3,0);
         sums[3]+=mA[3][1]*A(3,1);
         sums[3]+=mA[3][2]*A(3,2);
         sums[3]+=mA[3][3]*A(3,3);
         sums[3]+=mA[3][4]*A(3,4);

         sums[4] =mA[4][0]*A(3,0);
         sums[4]+=mA[4][1]*A(3,1);
         sums[4]+=mA[4][2]*A(3,2);
         sums[4]+=mA[4][3]*A(3,3);
         sums[4]+=mA[4][4]*A(3,4);

         temp(3,3)+=A(3,0)*sums[0];
         temp(3,3)+=A(3,1)*sums[1];
         temp(3,3)+=A(3,2)*sums[2];
         temp(3,3)+=A(3,3)*sums[3];
         temp(3,3)+=A(3,4)*sums[4];

         temp(4,3)+=A(4,0)*sums[0];
         temp(4,3)+=A(4,1)*sums[1];
         temp(4,3)+=A(4,2)*sums[2];
         temp(4,3)+=A(4,3)*sums[3];
         temp(4,3)+=A(4,4)*sums[4];

         temp(3,4)=temp(4,3);

         // Last entry ACA^T[4][4]  
         sums[0] =mA[0][0]*A(4,0);
         sums[0]+=mA[0][1]*A(4,1);
         sums[0]+=mA[0][2]*A(4,2);
         sums[0]+=mA[0][3]*A(4,3);
         sums[0]+=mA[0][4]*A(4,4);

         sums[1] =mA[1][0]*A(4,0);
         sums[1]+=mA[1][1]*A(4,1);
         sums[1]+=mA[1][2]*A(4,2);
         sums[1]+=mA[1][3]*A(4,3);
         sums[1]+=mA[1][4]*A(4,4);

         sums[2] =mA[2][0]*A(4,0);
         sums[2]+=mA[2][1]*A(4,1);
         sums[2]+=mA[2][2]*A(4,2);
         sums[2]+=mA[2][3]*A(4,3);
         sums[2]+=mA[2][4]*A(4,4);

         sums[3] =mA[3][0]*A(4,0);
         sums[3]+=mA[3][1]*A(4,1);
         sums[3]+=mA[3][2]*A(4,2);
         sums[3]+=mA[3][3]*A(4,3);
         sums[3]+=mA[3][4]*A(4,4);

         sums[4] =mA[4][0]*A(4,0);
         sums[4]+=mA[4][1]*A(4,1);
         sums[4]+=mA[4][2]*A(4,2);
         sums[4]+=mA[4][3]*A(4,3);
         sums[4]+=mA[4][4]*A(4,4);

         temp(4,4)+=A(4,0)*sums[0];
         temp(4,4)+=A(4,1)*sums[1];
         temp(4,4)+=A(4,2)*sums[2];
         temp(4,4)+=A(4,3)*sums[3];
         temp(4,4)+=A(4,4)*sums[4];

         return temp;
      }



      // The following code performs the matrix operation ABA^T, where B is a symmetric matrix.   The end 
      // result is also a symmetric matrix.
      DMatrix5x5 SandwichMultiply2(const DMatrix5x5 &A){
         DMatrix5x5 temp;
         double sums[5]={0};
         // First row/column of ACA^T
         for (unsigned int i=0;i<5;i++){
            for (unsigned int k=0;k<5;k++){
               sums[i]+=mA[i][k]*A(0,k);
            }
         }
         for (unsigned int i=0;i<5;i++){
            for (unsigned int k=0;k<5;k++){
               temp(i,0)+=A(i,k)*sums[k];
            }
            temp(0,i)=temp(i,0);
         }
         // Second row/column of ACA^T
         for (unsigned int i=0;i<5;i++){
            sums[i]=0.;
            for (unsigned int k=0;k<5;k++){
               sums[i]+=mA[i][k]*A(1,k);
            }
         }    
         for (unsigned int i=1;i<5;i++){
            for (unsigned int k=0;k<5;k++){
               temp(i,1)+=A(i,k)*sums[k];
            }
            temp(1,i)=temp(i,1);
         } 
         // Third row/column of ACA^T
         for (unsigned int i=0;i<5;i++){
            sums[i]=0.;
            for (unsigned int k=0;k<5;k++){
               sums[i]+=mA[i][k]*A(2,k);
            }
         }      
         for (unsigned int i=2;i<5;i++){
            for (unsigned int k=0;k<5;k++){
               temp(i,2)+=A(i,k)*sums[k];
            }
            temp(2,i)=temp(i,2);
         } 
         // Fourth row/column of ACA^T
         for (unsigned int i=0;i<5;i++){
            sums[i]=0.;
            for (unsigned int k=0;k<5;k++){
               sums[i]+=mA[i][k]*A(3,k);
            }
         }      
         for (unsigned int i=3;i<5;i++){
            for (unsigned int k=0;k<5;k++){
               temp(i,3)+=A(i,k)*sums[k];
            }
            temp(3,i)=temp(i,3);
         }
         // Last entry ACA^T[4][4]   
         for (unsigned int i=0;i<5;i++){
            sums[i]=0.;
            for (unsigned int k=0;k<5;k++){
               sums[i]+=mA[i][k]*A(4,k);
            }
         }
         for (unsigned int k=0;k<5;k++){
            temp(4,4)+=A(4,k)*sums[k];
         }
         return temp;
      }

      double SandwichMultiply(const DMatrix5x1 &A){
         double a1=A(0),a2=A(1),a3=A(2),a4=A(3),a5=A(4);
         return a1*(a1*mA[0][0]+2.*a2*mA[0][1]+2.*a3*mA[0][2]+2.*a4*mA[0][3]
               +2.*a5*mA[0][4])
            +a2*(a2*mA[1][1]+2.*a3*mA[1][2]+2.*a4*mA[1][3]+2.*a5*mA[1][4])
            +a3*(a3*mA[2][2]+2.*a4*mA[2][3]+2.*a5*mA[2][4])
            +a4*(a4*mA[3][3]+2.*a5*mA[3][4])
            +a5*a5*mA[4][4];

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
               GetSub(0,3).Determinant() > 0. && GetSub(0,4).Determinant() > 0.)
                  return true;
         else return false;
      }

      void Print(){
         cout << "DMatrix5x5:" <<endl;
         cout << "     |      0    |      1    |      2    |      3    |      4    |" <<endl;
         cout << "----------------------------------------------------------------------" <<endl;

         for (unsigned int i=0;i<5;i++){
            cout <<"   "<< i << " |";
            for (unsigned int j=0;j<5;j++){
               cout << setw(11)<<setprecision(6)<<mA[i][j] <<" "; 
            } 
            cout << endl;
         }      
      }

   private:
      double mA[5][5];

};


// Scale 5x5 matrix by a floating point number
inline DMatrix5x5 operator*(const double c,const DMatrix5x5 &M){
   DMatrix5x5 temp;
   for (unsigned int i=0;i<5;i++){
      for (unsigned int j=0;j<5;j++){
         temp(i,j)=c*M(i,j);
      }
   }
   return temp;
}


#else

// Matrix class with SIMD instructions

class DMatrix5x5{
   public:
      DMatrix5x5()
         : mA( ALIGNED_16_BLOCK_PTR(union dvec, 5, mA) )
      {
         for (unsigned int i=0;i<5;i++){
            for (unsigned int j=0;j<3;j++){
               mA[i].v[j]=_mm_setzero_pd();
            }
         }
      }
      DMatrix5x5( __m128d m11, __m128d m12, __m128d m13, __m128d m14, __m128d m15,
            __m128d m21, __m128d m22, __m128d m23, __m128d m24, __m128d m25,
            __m128d m31, __m128d m32, __m128d m33, __m128d m34, __m128d m35)
         : mA( ALIGNED_16_BLOCK_PTR(union dvec, 5, mA) )
      {
         mA[0].v[0]=m11;
         mA[1].v[0]=m12;  
         mA[2].v[0]=m13;
         mA[3].v[0]=m14;  
         mA[4].v[0]=m15;
         mA[0].v[1]=m21;
         mA[1].v[1]=m22;  
         mA[2].v[1]=m23;
         mA[3].v[1]=m24;  
         mA[4].v[1]=m25;
         mA[0].v[2]=m31;
         mA[1].v[2]=m32;  
         mA[2].v[2]=m33;
         mA[3].v[2]=m34;  
         mA[4].v[2]=m35;
      }
      // Constructor for a symmetric matrix
      DMatrix5x5( __m128d m11, __m128d m12, __m128d m13, __m128d m14, __m128d m15,
            __m128d m23, __m128d m24, __m128d m25,
            __m128d m35)
         : mA( ALIGNED_16_BLOCK_PTR(union dvec, 5, mA) )
      {
         mA[0].v[0]=m11;
         mA[1].v[0]=m12;  
         mA[2].v[0]=m13;
         mA[3].v[0]=m14;  
         mA[4].v[0]=m15;
         mA[0].v[1]=_mm_setr_pd(mA[2].d[0],mA[3].d[0]);
         mA[1].v[1]=_mm_setr_pd(mA[2].d[1],mA[3].d[1]);  
         mA[2].v[1]=m23;
         mA[3].v[1]=m24;  
         mA[4].v[1]=m25;
         mA[0].v[2]=_mm_setr_pd(mA[4].d[0],0.);
         mA[1].v[2]=_mm_setr_pd(mA[4].d[1],0.);
         mA[2].v[2]=_mm_setr_pd(mA[4].d[2],0.);
         mA[3].v[2]=_mm_setr_pd(mA[4].d[3],0.);
         mA[4].v[2]=m35;
      }
      // Constructor for symmetric matrix by elements
      DMatrix5x5(double C11, double C12, const double C13, double C14, double C15,
            double C22, const double C23, double C24, double C25,
            double C33, double C34, double C35,
            double C44, double C45,
            double C55)
         : mA( ALIGNED_16_BLOCK_PTR(union dvec, 5, mA) )
      {
         mA[0].v[0]=_mm_setr_pd(C11,C12);
         mA[1].v[0]=_mm_setr_pd(C12,C22);
         mA[2].v[0]=_mm_setr_pd(C13,C23);
         mA[3].v[0]=_mm_setr_pd(C14,C24);
         mA[4].v[0]=_mm_setr_pd(C15,C25);
         mA[0].v[1]=_mm_setr_pd(C13,C14);
         mA[1].v[1]=_mm_setr_pd(C23,C24);
         mA[2].v[1]=_mm_setr_pd(C33,C34);
         mA[3].v[1]=_mm_setr_pd(C34,C44);
         mA[4].v[1]=_mm_setr_pd(C35,C45);
         mA[0].v[2]=_mm_setr_pd(C15,0);
         mA[1].v[2]=_mm_setr_pd(C25,0);
         mA[2].v[2]=_mm_setr_pd(C35,0);
         mA[3].v[2]=_mm_setr_pd(C45,0);
         mA[4].v[2]=_mm_setr_pd(C55,0);

      }

      // Constructor using block matrices from matrix inversion 
      DMatrix5x5(const DMatrix2x2 &A,const DMatrix2x3 &B,
            const DMatrix3x2 &C,const DMatrix3x3 &D)
         : mA( ALIGNED_16_BLOCK_PTR(union dvec, 5, mA) )
      {
         mA[0].v[0]=A.GetV(0);
         mA[1].v[0]=A.GetV(1);
         for (unsigned int i=0;i<3;i++){
            mA[2+i].v[0]=B.GetV(i);
            for (unsigned int j=0;j<2;j++){
               mA[2+i].v[j+1]=D.GetV(j,i);
            }
         }
         mA[0].v[1]=C.GetV(0,0);
         mA[0].v[2]=C.GetV(1,0);
         mA[1].v[1]=C.GetV(0,1);
         mA[1].v[2]=C.GetV(1,1);
      }
      ~DMatrix5x5(){};

      __m128d GetV(int pair, int col) const{
         return mA[col].v[pair];
      }
      void SetV(int pair, int col, __m128d v){
         mA[col].v[pair]=v;
      }

      // Copy constructor
      DMatrix5x5(const DMatrix5x5 &m2)
         : mA( ALIGNED_16_BLOCK_PTR(union dvec, 5, mA) )
      {
         for (unsigned int i=0;i<5;i++){
            for (unsigned int j=0;j<3;j++){
               mA[i].v[j]=m2.GetV(j,i);
            }
         }
      }

      // return a column of the 5x5 matrix as a DMatrix5x1 object 
      DMatrix5x1 GetColumn(int col){
         return DMatrix5x1(mA[col].v[0],mA[col].v[1],mA[col].v[2]);

      }
      // return the trace of the matrix
      double Trace(){
         return(mA[0].d[0]+mA[1].d[1]+mA[2].d[2]+mA[3].d[3]+mA[4].d[4]);
      }


      // Access by indices
      double &operator() (int row, int col){
         return mA[col].d[row];
      } 
      double operator() (int row, int col) const{
         return mA[col].d[row];
      }

      // Assignment operator
      DMatrix5x5 &operator=(const DMatrix5x5 &m1){
         for (unsigned int i=0;i<5;i++){
            mA[i].v[0]=m1.GetV(0,i);
            mA[i].v[1]=m1.GetV(1,i);
            mA[i].v[2]=m1.GetV(2,i);
         }    
         return *this;
      }

      // Matrix addition
      DMatrix5x5 operator+(const DMatrix5x5 &m2) const{
#define ADD(i,j) _mm_add_pd(GetV((i),(j)),m2.GetV((i),(j)))

         return DMatrix5x5(ADD(0,0),ADD(0,1),ADD(0,2),ADD(0,3),ADD(0,4),ADD(1,0),ADD(1,1),ADD(1,2),
               ADD(1,3),ADD(1,4),ADD(2,0),ADD(2,1),ADD(2,2),ADD(2,3),ADD(2,4));
      }
      DMatrix5x5 &operator+=(const DMatrix5x5 &m2){
         for (unsigned int i=0;i<5;i++){
            for (unsigned int j=0;j<3;j++){
               mA[i].v[j]=ADD(j,i);
            }
         }
         return *this;
      }
      // Method for adding symmetric matrices
      DMatrix5x5 AddSym(const DMatrix5x5 &m2) const{
         return DMatrix5x5(ADD(0,0),ADD(0,1),ADD(0,2),ADD(0,3),ADD(0,4),ADD(1,2),
               ADD(1,3),ADD(1,4),ADD(2,4));
      }

      // Matrix subtraction
      DMatrix5x5 operator-(const DMatrix5x5 &m2) const{
#define SUB(i,j) _mm_sub_pd(GetV((i),(j)),m2.GetV((i),(j)))

         return DMatrix5x5(SUB(0,0),SUB(0,1),SUB(0,2),SUB(0,3),SUB(0,4),SUB(1,0),SUB(1,1),SUB(1,2),
               SUB(1,3),SUB(1,4),SUB(2,0),SUB(2,1),SUB(2,2),SUB(2,3),SUB(2,4));
      }
      // method for subtracting a symmetric matrix from another symmetric matrix
      DMatrix5x5 SubSym(const DMatrix5x5 &m2) const{
         return DMatrix5x5(SUB(0,0),SUB(0,1),SUB(0,2),SUB(0,3),SUB(0,4),SUB(1,2),
               SUB(1,3),SUB(1,4),SUB(2,4));
      }

      // Matrix multiplication:  (5x5) x (5x1)
      DMatrix5x1 operator*(const DMatrix5x1 &m2){
         ALIGNED_16_BLOCK_WITH_PTR(__m128d, 5, p)
            __m128d &a1=p[0];
         __m128d &a2=p[1];
         __m128d &a3=p[2];
         __m128d &a4=p[3];
         __m128d &a5=p[4];
         a1=_mm_set1_pd(m2(0));
         a2=_mm_set1_pd(m2(1));
         a3=_mm_set1_pd(m2(2));
         a4=_mm_set1_pd(m2(3));
         a5=_mm_set1_pd(m2(4));
         return 
            DMatrix5x1(_mm_add_pd(_mm_mul_pd(GetV(0,0),a1),
                     _mm_add_pd(_mm_mul_pd(GetV(0,1),a2),
                        _mm_add_pd(_mm_mul_pd(GetV(0,2),a3),
                           _mm_add_pd(_mm_mul_pd(GetV(0,3),a4),
                              _mm_mul_pd(GetV(0,4),a5))))),
                  _mm_add_pd(_mm_mul_pd(GetV(1,0),a1),
                     _mm_add_pd(_mm_mul_pd(GetV(1,1),a2),
                        _mm_add_pd(_mm_mul_pd(GetV(1,2),a3),
                           _mm_add_pd(_mm_mul_pd(GetV(1,3),a4),
                              _mm_mul_pd(GetV(1,4),a5))))),
                  _mm_add_pd(_mm_mul_pd(GetV(2,0),a1),
                     _mm_add_pd(_mm_mul_pd(GetV(2,1),a2),
                        _mm_add_pd(_mm_mul_pd(GetV(2,2),a3),
                           _mm_add_pd(_mm_mul_pd(GetV(2,3),a4),
                              _mm_mul_pd(GetV(2,4),a5))))));


      }

      // Matrix multiplication:  (5x5) x (5x2)
      DMatrix5x2 operator*(const DMatrix5x2 &m2){
         ALIGNED_16_BLOCK_WITH_PTR(__m128d, 10, p)
            __m128d &m11=p[0];
         __m128d &m12=p[1];
         __m128d &m21=p[2];
         __m128d &m22=p[3];
         __m128d &m31=p[4];
         __m128d &m32=p[5];
         __m128d &m41=p[6];
         __m128d &m42=p[7];
         __m128d &m51=p[8];
         __m128d &m52=p[9];
         m11=_mm_set1_pd(m2(0,0));
         m12=_mm_set1_pd(m2(0,1)); 
         m21=_mm_set1_pd(m2(1,0));
         m22=_mm_set1_pd(m2(1,1));  
         m31=_mm_set1_pd(m2(2,0));
         m32=_mm_set1_pd(m2(2,1)); 
         m41=_mm_set1_pd(m2(3,0));
         m42=_mm_set1_pd(m2(3,1)); 
         m51=_mm_set1_pd(m2(4,0));
         m52=_mm_set1_pd(m2(4,1));
         return 
            DMatrix5x2(_mm_add_pd(_mm_mul_pd(GetV(0,0),m11),
                     _mm_add_pd(_mm_mul_pd(GetV(0,1),m21),
                        _mm_add_pd(_mm_mul_pd(GetV(0,2),m31),
                           _mm_add_pd(_mm_mul_pd(GetV(0,3),m41),
                              _mm_mul_pd(GetV(0,4),m51))))),
                  _mm_add_pd(_mm_mul_pd(GetV(0,0),m12),
                     _mm_add_pd(_mm_mul_pd(GetV(0,1),m22),
                        _mm_add_pd(_mm_mul_pd(GetV(0,2),m32),
                           _mm_add_pd(_mm_mul_pd(GetV(0,3),m42),
                              _mm_mul_pd(GetV(0,4),m52))))),
                  _mm_add_pd(_mm_mul_pd(GetV(1,0),m11),
                     _mm_add_pd(_mm_mul_pd(GetV(1,1),m21),
                        _mm_add_pd(_mm_mul_pd(GetV(1,2),m31),
                           _mm_add_pd(_mm_mul_pd(GetV(1,3),m41),
                              _mm_mul_pd(GetV(1,4),m51))))),
                  _mm_add_pd(_mm_mul_pd(GetV(1,0),m12),
                     _mm_add_pd(_mm_mul_pd(GetV(1,1),m22),
                        _mm_add_pd(_mm_mul_pd(GetV(1,2),m32),
                           _mm_add_pd(_mm_mul_pd(GetV(1,3),m42),
                              _mm_mul_pd(GetV(1,4),m52))))), 
                  _mm_add_pd(_mm_mul_pd(GetV(2,0),m11),
                     _mm_add_pd(_mm_mul_pd(GetV(2,1),m21),
                        _mm_add_pd(_mm_mul_pd(GetV(2,2),m31),
                           _mm_add_pd(_mm_mul_pd(GetV(2,3),m41),
                              _mm_mul_pd(GetV(2,4),m51))))),
                  _mm_add_pd(_mm_mul_pd(GetV(2,0),m12),
                        _mm_add_pd(_mm_mul_pd(GetV(2,1),m22),
                           _mm_add_pd(_mm_mul_pd(GetV(2,2),m32),
                              _mm_add_pd(_mm_mul_pd(GetV(2,3),m42),
                                 _mm_mul_pd(GetV(2,4),m52))))));
      }

      // The following code performs the matrix operation ABA^T,where B is a 5x5 matrix, and A is 5x1
      double SandwichMultiply(const DMatrix5x1 &A){
         double a1=A(0),a2=A(1),a3=A(2),a4=A(3),a5=A(4);
         return a1*(a1*mA[0].d[0]+2.*a2*mA[0].d[1]+2.*a3*mA[0].d[2]+2.*a4*mA[0].d[3]
               +2.*a5*mA[0].d[4])
            +a2*(a2*mA[1].d[1]+2.*a3*mA[1].d[2]+2.*a4*mA[1].d[3]+2.*a5*mA[1].d[4])
            +a3*(a3*mA[2].d[2]+2.*a4*mA[2].d[3]+2.*a5*mA[2].d[4])
            +a4*(a4*mA[3].d[3]+2.*a5*mA[3].d[4])
            +a5*a5*mA[4].d[4];

      }



#ifdef USE_SSE3

      // The following code performs the matrix operation ABA^T, where B is a symmetric matrix
      DMatrix5x5 SandwichMultiply(const DMatrix5x5 &A){
         //__m128d A1112=_mm_setr_pd(A(0,0),A(0,1));
         //__m128d A1314=_mm_setr_pd(A(0,2),A(0,3));
#define A1112 _mm_setr_pd(A(0,0),A(0,1))
#define A1314 _mm_setr_pd(A(0,2),A(0,3))
         //__m128d A15=_mm_set1_pd(A(0,4));
         union dvec1{
            __m128d v[3];
            double d[6];
         };
         ALIGNED_16_BLOCK_WITH_PTR(union dvec1, 1, p1)
            union dvec1 &temp=p1[0];

         union dvec2{
            __m128d v;
            double d[2];
         };
         ALIGNED_16_BLOCK_WITH_PTR(union dvec2, 1, p2)
            union dvec2 &temp2=p2[0];
         temp2.v=_mm_set1_pd(A(0,4));

         // BA^T column 1
         temp.v[0]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,0),A1112),
                     _mm_mul_pd(GetV(0,1),A1112)),
                  _mm_hadd_pd(_mm_mul_pd(GetV(1,0),A1314),
                     _mm_mul_pd(GetV(1,1),A1314))),
               _mm_mul_pd(GetV(0,4),temp2.v));
         temp.v[1]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,2),A1112),
                     _mm_mul_pd(GetV(0,3),A1112)),
                  _mm_hadd_pd(_mm_mul_pd(GetV(1,2),A1314),
                     _mm_mul_pd(GetV(1,3),A1314))),
               _mm_mul_pd(GetV(1,4),temp2.v));
         temp.v[2] =_mm_set1_pd(mA[4].d[0]*A(0,0)+mA[4].d[1]*A(0,1)+mA[4].d[2]*A(0,2)+mA[4].d[3]*A(0,3)+mA[4].d[4]*A(0,4));



         // C=ABA^T elements C11,C12
#define A2122 _mm_setr_pd(A(1,0),A(1,1))
#define A2324 _mm_setr_pd(A(1,2),A(1,3))
         //    __m128d A2122=_mm_setr_pd(A(1,0),A(1,1));
         //__m128d A2324=_mm_setr_pd(A(1,2),A(1,3));
         temp2.v=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A1112,temp.v[0]),_mm_mul_pd(A2122,temp.v[0])),
                  _mm_hadd_pd(_mm_mul_pd(A1314,temp.v[1]),_mm_mul_pd(A2324,temp.v[1]))),
               _mm_mul_pd(A.GetV(0,4),temp.v[2]));
         double C11=temp2.d[0];
         double C12=temp2.d[1];

         // C=ABA^T elements C13,C14
#define A3132 _mm_setr_pd(A(2,0),A(2,1))
#define A3334 _mm_setr_pd(A(2,2),A(2,3))
#define A4142 _mm_setr_pd(A(3,0),A(3,1))
#define A4344 _mm_setr_pd(A(3,2),A(3,3))
         //__m128d A3132=_mm_setr_pd(A(2,0),A(2,1));
         //__m128d A3334=_mm_setr_pd(A(2,2),A(2,3));
         //__m128d A4142=_mm_setr_pd(A(3,0),A(3,1));
         //__m128d A4344=_mm_setr_pd(A(3,2),A(3,3));
         temp2.v=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A3132,temp.v[0]),_mm_mul_pd(A4142,temp.v[0])),
                  _mm_hadd_pd(_mm_mul_pd(A3334,temp.v[1]),_mm_mul_pd(A4344,temp.v[1]))),
               _mm_mul_pd(A.GetV(1,4),temp.v[2]));
         double C13=temp2.d[0];
         double C14=temp2.d[1];

         // BA^T column 2
         //#define A25 _mm_set1_pd(A(1,4))
         temp2.v=_mm_set1_pd(A(1,4));
         temp.v[0]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,0),A2122),
                     _mm_mul_pd(GetV(0,1),A2122)),
                  _mm_hadd_pd(_mm_mul_pd(GetV(1,0),A2324),
                     _mm_mul_pd(GetV(1,1),A2324))),
               _mm_mul_pd(GetV(0,4),temp2.v));
         temp.v[1]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,2),A2122),
                     _mm_mul_pd(GetV(0,3),A2122)),
                  _mm_hadd_pd(_mm_mul_pd(GetV(1,2),A2324),
                     _mm_mul_pd(GetV(1,3),A2324))),
               _mm_mul_pd(GetV(1,4),temp2.v));
         temp.v[2]=_mm_set1_pd(mA[4].d[0]*A(1,0)+mA[4].d[1]*A(1,1)+mA[4].d[2]*A(1,2)+mA[4].d[3]*A(1,3)+mA[4].d[4]*A(1,4));

         // C=ABA^T elements C22,C23
#define A2535 _mm_setr_pd(A(1,4),A(2,4))
         temp2.v=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A2122,temp.v[0]),_mm_mul_pd(A3132,temp.v[0])),
                  _mm_hadd_pd(_mm_mul_pd(A2324,temp.v[1]),_mm_mul_pd(A3334,temp.v[1]))),
               _mm_mul_pd(A2535,temp.v[2]));
         double C22=temp2.d[0];
         double C23=temp2.d[1];

         // C=ABA^T elements C24,C25
#define A5152 _mm_setr_pd(A(4,0),A(4,1))
#define A5354 _mm_setr_pd(A(4,2),A(4,3)) 
#define A4555 _mm_setr_pd(A(3,4),A(4,4))
         //__m128d A5152=_mm_setr_pd(A(4,0),A(4,1));
         //__m128d A5354=_mm_setr_pd(A(4,2),A(4,3)); 
         //__m128d A4555=_mm_setr_pd(A(3,4),A(4,4));
         temp2.v=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A4142,temp.v[0]),_mm_mul_pd(A5152,temp.v[0])),
                  _mm_hadd_pd(_mm_mul_pd(A4344,temp.v[1]),_mm_mul_pd(A5354,temp.v[1]))),
               _mm_mul_pd(A4555,temp.v[2]));  
         double C24=temp2.d[0];
         double C25=temp2.d[1];

         // BA^T column 3
         //#define A35 _mm_set1_pd(A(2,4))
         temp2.v=_mm_set1_pd(A(2,4));
         temp.v[0]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,0),A3132),
                     _mm_mul_pd(GetV(0,1),A3132)),
                  _mm_hadd_pd(_mm_mul_pd(GetV(1,0),A3334),
                     _mm_mul_pd(GetV(1,1),A3334))),
               _mm_mul_pd(GetV(0,4),temp2.v));
         temp.v[1]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,2),A3132),
                     _mm_mul_pd(GetV(0,3),A3132)),
                  _mm_hadd_pd(_mm_mul_pd(GetV(1,2),A3334),
                     _mm_mul_pd(GetV(1,3),A3334))),
               _mm_mul_pd(GetV(1,4),temp2.v));
         temp.v[2]=_mm_set1_pd(mA[4].d[0]*A(2,0)+mA[4].d[1]*A(2,1)+mA[4].d[2]*A(2,2)+mA[4].d[3]*A(2,3)+mA[4].d[4]*A(2,4));


         // C=ABA^T elements C33,C34,C35
#define A3545 _mm_setr_pd(A(2,4),A(3,4))
         temp2.v=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A3132,temp.v[0]),_mm_mul_pd(A4142,temp.v[0])),
                  _mm_hadd_pd(_mm_mul_pd(A3334,temp.v[1]),_mm_mul_pd(A4344,temp.v[1]))),
               _mm_mul_pd(A3545,temp.v[2]));
         double C33=temp2.d[0];
         double C34=temp2.d[1];
         double C35=A(4,0)*temp.d[0]+A(4,1)*temp.d[1]+A(4,2)*temp.d[2]+A(4,3)*temp.d[3]+A(4,4)*temp.d[4];

         // BA^T column 4
         //#define A45 _mm_set1_pd(A(3,4))
         temp2.v=_mm_set1_pd(A(3,4));
         temp.v[0]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,0),A4142),
                     _mm_mul_pd(GetV(0,1),A4142)),
                  _mm_hadd_pd(_mm_mul_pd(GetV(1,0),A4344),
                     _mm_mul_pd(GetV(1,1),A4344))),
               _mm_mul_pd(GetV(0,4),temp2.v));
         temp.v[1]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,2),A4142),
                     _mm_mul_pd(GetV(0,3),A4142)),
                  _mm_hadd_pd(_mm_mul_pd(GetV(1,2),A4344),
                     _mm_mul_pd(GetV(1,3),A4344))),
               _mm_mul_pd(GetV(1,4),temp2.v));
         temp.v[2]=_mm_set1_pd(mA[4].d[0]*A(3,0)+mA[4].d[1]*A(3,1)+mA[4].d[2]*A(3,2)+mA[4].d[3]*A(3,3)+mA[4].d[4]*A(3,4));

         // C=ABA^T elements C44,C45
         temp2.v=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A4142,temp.v[0]),_mm_mul_pd(A5152,temp.v[0])),
                  _mm_hadd_pd(_mm_mul_pd(A4344,temp.v[1]),_mm_mul_pd(A5354,temp.v[1]))),
               _mm_mul_pd(A4555,temp.v[2]));
         double C44=temp2.d[0];
         double C45=temp2.d[1];

         // BA^T column 5
         //#define A55 _mm_set1_pd(A(4,4))
         temp2.v=_mm_set1_pd(A(4,4));
         temp.v[0]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,0),A5152),
                     _mm_mul_pd(GetV(0,1),A5152)),
                  _mm_hadd_pd(_mm_mul_pd(GetV(1,0),A5354),
                     _mm_mul_pd(GetV(1,1),A5354))),
               _mm_mul_pd(GetV(0,4),temp2.v));
         temp.v[1]=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(GetV(0,2),A5152),
                     _mm_mul_pd(GetV(0,3),A5152)),
                  _mm_hadd_pd(_mm_mul_pd(GetV(1,2),A5354),
                     _mm_mul_pd(GetV(1,3),A5354))),
               _mm_mul_pd(GetV(1,4),temp2.v));
         temp.v[2]=_mm_set1_pd(mA[4].d[0]*A(4,0)+mA[4].d[1]*A(4,1)+mA[4].d[2]*A(4,2)+mA[4].d[3]*A(4,3)+mA[4].d[4]*A(4,4));


         // C=ABA^T elements C15,C55
#define A1555 _mm_setr_pd(A(0,4),A(4,4)) 
         temp2.v=_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A1112,temp.v[0]),_mm_mul_pd(A5152,temp.v[0])),
                  _mm_hadd_pd(_mm_mul_pd(A1314,temp.v[1]),_mm_mul_pd(A5354,temp.v[1]))),
               _mm_mul_pd(A1555,temp.v[2]));
         double C15=temp2.d[0];
         double C55=temp2.d[1];

         return DMatrix5x5(C11,C12,C13,C14,C15,C22,C23,C24,C25,C33,C34,C35,C44,C45,C55);
      }

      // Matrix multiplication. Requires the SSE3 instruction HADD (horizontal add)
      DMatrix5x5 operator*(const DMatrix5x5 &m2){
         ALIGNED_16_BLOCK_WITH_PTR(__m128d, 8, p)
            __m128d &A11A12=p[0];
         __m128d &A13A14=p[1];
         __m128d &A21A22=p[2];
         __m128d &A23A24=p[3];
         __m128d &A31A32=p[4];
         __m128d &A33A34=p[5];
         __m128d &A41A42=p[6];
         __m128d &A43A44=p[7];
         A11A12=_mm_setr_pd(mA[0].d[0],mA[1].d[0]);
         A13A14=_mm_setr_pd(mA[2].d[0],mA[3].d[0]);
         A21A22=_mm_setr_pd(mA[0].d[1],mA[1].d[1]);
         A23A24=_mm_setr_pd(mA[2].d[1],mA[3].d[1]);
         A31A32=_mm_setr_pd(mA[0].d[2],mA[1].d[2]);
         A33A34=_mm_setr_pd(mA[2].d[2],mA[3].d[2]);
         A41A42=_mm_setr_pd(mA[0].d[3],mA[1].d[3]);
         A43A44=_mm_setr_pd(mA[2].d[3],mA[3].d[3]);
         return
            DMatrix5x5(_mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A11A12,m2.GetV(0,0)),_mm_mul_pd(A21A22,m2.GetV(0,0))),
                        _mm_hadd_pd(_mm_mul_pd(A13A14,m2.GetV(1,0)),_mm_mul_pd(A23A24,m2.GetV(1,0)))),
                     _mm_mul_pd(mA[4].v[0],_mm_set1_pd(m2(4,0)))),
                  _mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A11A12,m2.GetV(0,1)),_mm_mul_pd(A21A22,m2.GetV(0,1))),
                        _mm_hadd_pd(_mm_mul_pd(A13A14,m2.GetV(1,1)),_mm_mul_pd(A23A24,m2.GetV(1,1)))),
                     _mm_mul_pd(mA[4].v[0],_mm_set1_pd(m2(4,1)))),
                  _mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A11A12,m2.GetV(0,2)),_mm_mul_pd(A21A22,m2.GetV(0,2))),
                        _mm_hadd_pd(_mm_mul_pd(A13A14,m2.GetV(1,2)),_mm_mul_pd(A23A24,m2.GetV(1,2)))),
                     _mm_mul_pd(mA[4].v[0],_mm_set1_pd(m2(4,2)))),
                  _mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A11A12,m2.GetV(0,3)),_mm_mul_pd(A21A22,m2.GetV(0,3))),
                        _mm_hadd_pd(_mm_mul_pd(A13A14,m2.GetV(1,3)),_mm_mul_pd(A23A24,m2.GetV(1,3)))),
                     _mm_mul_pd(mA[4].v[0],_mm_set1_pd(m2(4,3)))),
                  _mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A11A12,m2.GetV(0,4)),_mm_mul_pd(A21A22,m2.GetV(0,4))),
                        _mm_hadd_pd(_mm_mul_pd(A13A14,m2.GetV(1,4)),_mm_mul_pd(A23A24,m2.GetV(1,4)))),
                     _mm_mul_pd(mA[4].v[0],_mm_set1_pd(m2(4,4)))),

                  _mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A31A32,m2.GetV(0,0)),_mm_mul_pd(A41A42,m2.GetV(0,0))),
                        _mm_hadd_pd(_mm_mul_pd(A33A34,m2.GetV(1,0)),_mm_mul_pd(A43A44,m2.GetV(1,0)))),
                     _mm_mul_pd(mA[4].v[1],_mm_set1_pd(m2(4,0)))),
                  _mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A31A32,m2.GetV(0,1)),_mm_mul_pd(A41A42,m2.GetV(0,1))),
                        _mm_hadd_pd(_mm_mul_pd(A33A34,m2.GetV(1,1)),_mm_mul_pd(A43A44,m2.GetV(1,1)))),
                     _mm_mul_pd(mA[4].v[1],_mm_set1_pd(m2(4,1)))),	
                  _mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A31A32,m2.GetV(0,2)),_mm_mul_pd(A41A42,m2.GetV(0,2))),
                           _mm_hadd_pd(_mm_mul_pd(A33A34,m2.GetV(1,2)),_mm_mul_pd(A43A44,m2.GetV(1,2)))),
                        _mm_mul_pd(mA[4].v[1],_mm_set1_pd(m2(4,2)))),
                  _mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A31A32,m2.GetV(0,3)),_mm_mul_pd(A41A42,m2.GetV(0,3))),
                           _mm_hadd_pd(_mm_mul_pd(A33A34,m2.GetV(1,3)),_mm_mul_pd(A43A44,m2.GetV(1,3)))),
                        _mm_mul_pd(mA[4].v[1],_mm_set1_pd(m2(4,3)))),	
                  _mm_add_pd(_mm_add_pd(_mm_hadd_pd(_mm_mul_pd(A31A32,m2.GetV(0,4)),_mm_mul_pd(A41A42,m2.GetV(0,4))),
                           _mm_hadd_pd(_mm_mul_pd(A33A34,m2.GetV(1,4)),_mm_mul_pd(A43A44,m2.GetV(1,4)))),
                        _mm_mul_pd(mA[4].v[1],_mm_set1_pd(m2(4,4)))),	

                  _mm_set_sd(mA[0].d[4]*m2(0,0)+mA[1].d[4]*m2(1,0)+mA[2].d[4]*m2(2,0)+mA[3].d[4]*m2(3,0)+mA[4].d[4]*m2(4,0)),
                  _mm_set_sd(mA[0].d[4]*m2(0,1)+mA[1].d[4]*m2(1,1)+mA[2].d[4]*m2(2,1)+mA[3].d[4]*m2(3,1)+mA[4].d[4]*m2(4,1)),
                  _mm_set_sd(mA[0].d[4]*m2(0,2)+mA[1].d[4]*m2(1,2)+mA[2].d[4]*m2(2,2)+mA[3].d[4]*m2(3,2)+mA[4].d[4]*m2(4,2)),
                  _mm_set_sd(mA[0].d[4]*m2(0,3)+mA[1].d[4]*m2(1,3)+mA[2].d[4]*m2(2,3)+mA[3].d[4]*m2(3,3)+mA[4].d[4]*m2(4,3)),
                  _mm_set_sd(mA[0].d[4]*m2(0,4)+mA[1].d[4]*m2(1,4)+mA[2].d[4]*m2(2,4)+mA[3].d[4]*m2(3,4)+mA[4].d[4]*m2(4,4))
                     );
      }


#else 
      // Matrix multiplication: (5x5) x (5x5)
      DMatrix5x5 operator*(const DMatrix5x5 &m2){
         struct dvec1{
            __m128d v[5][5];
         };
         ALIGNED_16_BLOCK_WITH_PTR(struct dvec1, 1, p)
            struct dvec1 &temp=p[0];
         for (unsigned int i=0;i<5;i++){
            for (unsigned int j=0;j<5;j++){
               temp.v[i][j]=_mm_set1_pd(m2(i,j));
            }
         }
         // Preprocessor macro for multiplying two __m128d elements together
#define MUL(i,j,k) _mm_mul_pd(GetV((i),(j)),temp.v[(j)][(k)])

         return DMatrix5x5(_mm_add_pd(MUL(0,0,0),
                  _mm_add_pd(MUL(0,1,0),
                     _mm_add_pd(MUL(0,2,0),
                        _mm_add_pd(MUL(0,3,0),
                           MUL(0,4,0))))),
               _mm_add_pd(MUL(0,0,1),
                  _mm_add_pd(MUL(0,1,1),
                     _mm_add_pd(MUL(0,2,1),
                        _mm_add_pd(MUL(0,3,1),
                           MUL(0,4,1))))),
               _mm_add_pd(MUL(0,0,2),
                  _mm_add_pd(MUL(0,1,2),
                     _mm_add_pd(MUL(0,2,2),
                        _mm_add_pd(MUL(0,3,2),
                           MUL(0,4,2))))),
               _mm_add_pd(MUL(0,0,3),
                  _mm_add_pd(MUL(0,1,3),
                     _mm_add_pd(MUL(0,2,3),
                        _mm_add_pd(MUL(0,3,3),
                           MUL(0,4,3))))),  
               _mm_add_pd(MUL(0,0,4),
                  _mm_add_pd(MUL(0,1,4),
                     _mm_add_pd(MUL(0,2,4),
                        _mm_add_pd(MUL(0,3,4),
                           MUL(0,4,4))))),
               _mm_add_pd(MUL(1,0,0),
                     _mm_add_pd(MUL(1,1,0),
                        _mm_add_pd(MUL(1,2,0),
                           _mm_add_pd(MUL(1,3,0),
                              MUL(1,4,0))))),
               _mm_add_pd(MUL(1,0,1),
                     _mm_add_pd(MUL(1,1,1),
                        _mm_add_pd(MUL(1,2,1),
                           _mm_add_pd(MUL(1,3,1),
                              MUL(1,4,1))))),
               _mm_add_pd(MUL(1,0,2),
                     _mm_add_pd(MUL(1,1,2),
                        _mm_add_pd(MUL(1,2,2),
                           _mm_add_pd(MUL(1,3,2),
                              MUL(1,4,2))))),
               _mm_add_pd(MUL(1,0,3),
                     _mm_add_pd(MUL(1,1,3),
                        _mm_add_pd(MUL(1,2,3),
                           _mm_add_pd(MUL(1,3,3),
                              MUL(1,4,3))))),  
               _mm_add_pd(MUL(1,0,4),
                     _mm_add_pd(MUL(1,1,4),
                        _mm_add_pd(MUL(1,2,4),
                           _mm_add_pd(MUL(1,3,4),
                              MUL(1,4,4))))),
               _mm_add_pd(MUL(2,0,0),
                     _mm_add_pd(MUL(2,1,0),
                        _mm_add_pd(MUL(2,2,0),
                           _mm_add_pd(MUL(2,3,0),
                              MUL(2,4,0))))),
               _mm_add_pd(MUL(2,0,1),
                     _mm_add_pd(MUL(2,1,1),
                        _mm_add_pd(MUL(2,2,1),
                           _mm_add_pd(MUL(2,3,1),
                              MUL(2,4,1))))),
               _mm_add_pd(MUL(2,0,2),
                     _mm_add_pd(MUL(2,1,2),
                        _mm_add_pd(MUL(2,2,2),
                           _mm_add_pd(MUL(2,3,2),
                              MUL(2,4,2))))),
               _mm_add_pd(MUL(2,0,3),
                     _mm_add_pd(MUL(2,1,3),
                        _mm_add_pd(MUL(2,2,3),
                           _mm_add_pd(MUL(2,3,3),
                              MUL(2,4,3))))),  
               _mm_add_pd(MUL(2,0,4),
                     _mm_add_pd(MUL(2,1,4),
                        _mm_add_pd(MUL(2,2,4),
                           _mm_add_pd(MUL(2,3,4),
                              MUL(2,4,4))))));    
      }

      // The following code performs the matrix operation ABA^T, where B is a symmetric matrix
      DMatrix5x5 SandwichMultiply(const DMatrix5x5 &A){
         ALIGNED_16_BLOCK_WITH_PTR(__m128d, 5, p)
            __m128d &A1=p[0];
         __m128d &A2=p[1];
         __m128d &A3=p[2];
         __m128d &A4=p[3];
         __m128d &A5=p[4];
         A5=_mm_set1_pd(A(0,4));
         A4=_mm_set1_pd(A(0,3));
         A3=_mm_set1_pd(A(0,2));
         A2=_mm_set1_pd(A(0,1));
         A1=_mm_set1_pd(A(0,0));

         // BA^T column 1 
         union dvec2{
            __m128d v;
            double d[2];
         };
         ALIGNED_16_BLOCK_WITH_PTR(union dvec2, 1, p2)
            union dvec2 &temp2=p2[0];
         temp2.v=_mm_add_pd(_mm_mul_pd(A1,GetV(0,0)),
               _mm_add_pd(_mm_mul_pd(A2,GetV(0,1)),
                  _mm_add_pd(_mm_mul_pd(A3,GetV(0,2)),
                     _mm_add_pd(_mm_mul_pd(A4,GetV(0,3)),
                        _mm_mul_pd(A5,GetV(0,4))))));
         ALIGNED_16_BLOCK_WITH_PTR(__m128d, 5, p3)
            __m128d &BA1=p3[0];
         __m128d &BA2=p3[1];
         __m128d &BA3=p3[2];
         __m128d &BA4=p3[3];
         __m128d &BA5=p3[4];
         BA1=_mm_set1_pd(temp2.d[0]);
         BA2=_mm_set1_pd(temp2.d[1]);
         temp2.v=_mm_add_pd(_mm_mul_pd(A1,GetV(1,0)),
               _mm_add_pd(_mm_mul_pd(A2,GetV(1,1)),
                  _mm_add_pd(_mm_mul_pd(A3,GetV(1,2)),
                     _mm_add_pd(_mm_mul_pd(A4,GetV(1,3)),
                        _mm_mul_pd(A5,GetV(1,4))))));
         BA3=_mm_set1_pd(temp2.d[0]);
         BA4=_mm_set1_pd(temp2.d[1]);
         BA5=_mm_set1_pd(mA[0].d[4]*A(0,0)+mA[1].d[4]*A(0,1)+mA[2].d[4]*A(0,2)+mA[3].d[4]*A(0,3)+mA[4].d[4]*A(0,4));

         temp2.v=_mm_add_pd(_mm_mul_pd(A.GetV(0,0),BA1),
               _mm_add_pd(_mm_mul_pd(A.GetV(0,1),BA2),
                  _mm_add_pd(_mm_mul_pd(A.GetV(0,2),BA3),
                     _mm_add_pd(_mm_mul_pd(A.GetV(0,3),BA4),
                        _mm_mul_pd(A.GetV(0,4),BA5)))));
         double C11=temp2.d[0];
         double C12=temp2.d[1];

         temp2.v=_mm_add_pd(_mm_mul_pd(A.GetV(1,0),BA1),
               _mm_add_pd(_mm_mul_pd(A.GetV(1,1),BA2),
                  _mm_add_pd(_mm_mul_pd(A.GetV(1,2),BA3),
                     _mm_add_pd(_mm_mul_pd(A.GetV(1,3),BA4),
                        _mm_mul_pd(A.GetV(1,4),BA5)))));

         double C13=temp2.d[0];
         double C14=temp2.d[1];

         temp2.v=_mm_add_pd(_mm_mul_pd(A.GetV(2,0),BA1),
               _mm_add_pd(_mm_mul_pd(A.GetV(2,1),BA2),
                  _mm_add_pd(_mm_mul_pd(A.GetV(2,2),BA3),
                     _mm_add_pd(_mm_mul_pd(A.GetV(2,3),BA4),
                        _mm_mul_pd(A.GetV(2,4),BA5)))));

         double C15=temp2.d[0];


         // BA^T column 2 
         A5=_mm_set1_pd(A(1,4));
         A4=_mm_set1_pd(A(1,3));
         A3=_mm_set1_pd(A(1,2));
         A2=_mm_set1_pd(A(1,1));
         A1=_mm_set1_pd(A(1,0));
         temp2.v=_mm_add_pd(_mm_mul_pd(A1,GetV(0,0)),
               _mm_add_pd(_mm_mul_pd(A2,GetV(0,1)),
                  _mm_add_pd(_mm_mul_pd(A3,GetV(0,2)),
                     _mm_add_pd(_mm_mul_pd(A4,GetV(0,3)),
                        _mm_mul_pd(A5,GetV(0,4))))));
         BA1=_mm_set1_pd(temp2.d[0]);
         BA2=_mm_set1_pd(temp2.d[1]);
         temp2.v=_mm_add_pd(_mm_mul_pd(A1,GetV(1,0)),
               _mm_add_pd(_mm_mul_pd(A2,GetV(1,1)),
                  _mm_add_pd(_mm_mul_pd(A3,GetV(1,2)),
                     _mm_add_pd(_mm_mul_pd(A4,GetV(1,3)),
                        _mm_mul_pd(A5,GetV(1,4))))));
         BA3=_mm_set1_pd(temp2.d[0]);
         BA4=_mm_set1_pd(temp2.d[1]);
         BA5=_mm_set1_pd(mA[0].d[4]*A(1,0)+mA[1].d[4]*A(1,1)+mA[2].d[4]*A(1,2)+mA[3].d[4]*A(1,3)+mA[4].d[4]*A(1,4));

         temp2.v=_mm_add_pd(_mm_mul_pd(A.GetV(0,0),BA1),
               _mm_add_pd(_mm_mul_pd(A.GetV(0,1),BA2),
                  _mm_add_pd(_mm_mul_pd(A.GetV(0,2),BA3),
                     _mm_add_pd(_mm_mul_pd(A.GetV(0,3),BA4),
                        _mm_mul_pd(A.GetV(0,4),BA5)))));

         double C22=temp2.d[1];

         temp2.v=_mm_add_pd(_mm_mul_pd(A.GetV(1,0),BA1),
               _mm_add_pd(_mm_mul_pd(A.GetV(1,1),BA2),
                  _mm_add_pd(_mm_mul_pd(A.GetV(1,2),BA3),
                     _mm_add_pd(_mm_mul_pd(A.GetV(1,3),BA4),
                        _mm_mul_pd(A.GetV(1,4),BA5)))));

         double C23=temp2.d[0];
         double C24=temp2.d[1];

         temp2.v=_mm_add_pd(_mm_mul_pd(A.GetV(2,0),BA1),
               _mm_add_pd(_mm_mul_pd(A.GetV(2,1),BA2),
                  _mm_add_pd(_mm_mul_pd(A.GetV(2,2),BA3),
                     _mm_add_pd(_mm_mul_pd(A.GetV(2,3),BA4),
                        _mm_mul_pd(A.GetV(2,4),BA5)))));

         double C25=temp2.d[0];


         // BA^T column 3 
         A5=_mm_set1_pd(A(2,4));
         A4=_mm_set1_pd(A(2,3));
         A3=_mm_set1_pd(A(2,2));
         A2=_mm_set1_pd(A(2,1));
         A1=_mm_set1_pd(A(2,0));
         temp2.v=_mm_add_pd(_mm_mul_pd(A1,GetV(0,0)),
               _mm_add_pd(_mm_mul_pd(A2,GetV(0,1)),
                  _mm_add_pd(_mm_mul_pd(A3,GetV(0,2)),
                     _mm_add_pd(_mm_mul_pd(A4,GetV(0,3)),
                        _mm_mul_pd(A5,GetV(0,4))))));
         BA1=_mm_set1_pd(temp2.d[0]);
         BA2=_mm_set1_pd(temp2.d[1]);
         temp2.v=_mm_add_pd(_mm_mul_pd(A1,GetV(1,0)),
               _mm_add_pd(_mm_mul_pd(A2,GetV(1,1)),
                  _mm_add_pd(_mm_mul_pd(A3,GetV(1,2)),
                     _mm_add_pd(_mm_mul_pd(A4,GetV(1,3)),
                        _mm_mul_pd(A5,GetV(1,4))))));
         BA3=_mm_set1_pd(temp2.d[0]);
         BA4=_mm_set1_pd(temp2.d[1]);
         BA5=_mm_set1_pd(mA[0].d[4]*A(2,0)+mA[1].d[4]*A(2,1)+mA[2].d[4]*A(2,2)+mA[3].d[4]*A(2,3)+mA[4].d[4]*A(2,4));

         temp2.v=_mm_add_pd(_mm_mul_pd(A.GetV(1,0),BA1),
               _mm_add_pd(_mm_mul_pd(A.GetV(1,1),BA2),
                  _mm_add_pd(_mm_mul_pd(A.GetV(1,2),BA3),
                     _mm_add_pd(_mm_mul_pd(A.GetV(1,3),BA4),
                        _mm_mul_pd(A.GetV(1,4),BA5)))));

         double C33=temp2.d[0];
         double C34=temp2.d[1];

         temp2.v=_mm_add_pd(_mm_mul_pd(A.GetV(2,0),BA1),
               _mm_add_pd(_mm_mul_pd(A.GetV(2,1),BA2),
                  _mm_add_pd(_mm_mul_pd(A.GetV(2,2),BA3),
                     _mm_add_pd(_mm_mul_pd(A.GetV(2,3),BA4),
                        _mm_mul_pd(A.GetV(2,4),BA5)))));

         double C35=temp2.d[0];


         // BA^T column 4
         A5=_mm_set1_pd(A(3,4));
         A4=_mm_set1_pd(A(3,3));
         A3=_mm_set1_pd(A(3,2));
         A2=_mm_set1_pd(A(3,1));
         A1=_mm_set1_pd(A(3,0));
         temp2.v=_mm_add_pd(_mm_mul_pd(A1,GetV(0,0)),
               _mm_add_pd(_mm_mul_pd(A2,GetV(0,1)),
                  _mm_add_pd(_mm_mul_pd(A3,GetV(0,2)),
                     _mm_add_pd(_mm_mul_pd(A4,GetV(0,3)),
                        _mm_mul_pd(A5,GetV(0,4))))));
         BA1=_mm_set1_pd(temp2.d[0]);
         BA2=_mm_set1_pd(temp2.d[1]);
         temp2.v=_mm_add_pd(_mm_mul_pd(A1,GetV(1,0)),
               _mm_add_pd(_mm_mul_pd(A2,GetV(1,1)),
                  _mm_add_pd(_mm_mul_pd(A3,GetV(1,2)),
                     _mm_add_pd(_mm_mul_pd(A4,GetV(1,3)),
                        _mm_mul_pd(A5,GetV(1,4))))));
         BA3=_mm_set1_pd(temp2.d[0]);
         BA4=_mm_set1_pd(temp2.d[1]);
         BA5=_mm_set1_pd(mA[0].d[4]*A(3,0)+mA[1].d[4]*A(3,1)+mA[2].d[4]*A(3,2)+mA[3].d[4]*A(3,3)+mA[4].d[4]*A(3,4));

         temp2.v=_mm_add_pd(_mm_mul_pd(A.GetV(1,0),BA1),
               _mm_add_pd(_mm_mul_pd(A.GetV(1,1),BA2),
                  _mm_add_pd(_mm_mul_pd(A.GetV(1,2),BA3),
                     _mm_add_pd(_mm_mul_pd(A.GetV(1,3),BA4),
                        _mm_mul_pd(A.GetV(1,4),BA5)))));
         double C44=temp2.d[1];


         temp2.v=_mm_add_pd(_mm_mul_pd(A.GetV(2,0),BA1),
               _mm_add_pd(_mm_mul_pd(A.GetV(2,1),BA2),
                  _mm_add_pd(_mm_mul_pd(A.GetV(2,2),BA3),
                     _mm_add_pd(_mm_mul_pd(A.GetV(2,3),BA4),
                        _mm_mul_pd(A.GetV(2,4),BA5)))));

         double C45=temp2.d[0];

         // BA^T column 5
         A5=_mm_set1_pd(A(4,4));
         A4=_mm_set1_pd(A(4,3));
         A3=_mm_set1_pd(A(4,2));
         A2=_mm_set1_pd(A(4,1));
         A1=_mm_set1_pd(A(4,0));
         temp2.v=_mm_add_pd(_mm_mul_pd(A1,GetV(0,0)),
               _mm_add_pd(_mm_mul_pd(A2,GetV(0,1)),
                  _mm_add_pd(_mm_mul_pd(A3,GetV(0,2)),
                     _mm_add_pd(_mm_mul_pd(A4,GetV(0,3)),
                        _mm_mul_pd(A5,GetV(0,4))))));
         BA1=_mm_set1_pd(temp2.d[0]);
         BA2=_mm_set1_pd(temp2.d[1]);
         temp2.v=_mm_add_pd(_mm_mul_pd(A1,GetV(1,0)),
               _mm_add_pd(_mm_mul_pd(A2,GetV(1,1)),
                  _mm_add_pd(_mm_mul_pd(A3,GetV(1,2)),
                     _mm_add_pd(_mm_mul_pd(A4,GetV(1,3)),
                        _mm_mul_pd(A5,GetV(1,4))))));

         BA3=_mm_set1_pd(temp2.d[0]);
         BA4=_mm_set1_pd(temp2.d[1]);
         BA5=_mm_set1_pd(mA[0].d[4]*A(4,0)+mA[1].d[4]*A(4,1)+mA[2].d[4]*A(4,2)+mA[3].d[4]*A(4,3)+mA[4].d[4]*A(4,4));


         temp2.v=_mm_add_pd(_mm_mul_pd(A.GetV(2,0),BA1),
               _mm_add_pd(_mm_mul_pd(A.GetV(2,1),BA2),
                  _mm_add_pd(_mm_mul_pd(A.GetV(2,2),BA3),
                     _mm_add_pd(_mm_mul_pd(A.GetV(2,3),BA4),
                        _mm_mul_pd(A.GetV(2,4),BA5)))));
         double C55=temp2.d[0];

         return DMatrix5x5(C11,C12,C13,C14,C15,C22,C23,C24,C25,C33,C34,C35,C44,C45,C55);
      }

#endif


      // Find the transpose of this matrix
      DMatrix5x5 Transpose(){
#define SWAP(i,j,k,m) _mm_setr_pd(mA[(j)].d[(i)],mA[(m)].d[(k)])

         return DMatrix5x5(SWAP(0,0,0,1),SWAP(1,0,1,1),SWAP(2,0,2,1),SWAP(3,0,3,1),SWAP(4,0,4,1),
               SWAP(0,2,0,3),SWAP(1,2,1,3),SWAP(2,2,2,3),SWAP(3,2,3,3),SWAP(4,2,4,3),

               _mm_setr_pd(mA[4].d[0],0.),_mm_setr_pd(mA[4].d[1],0.),
               _mm_setr_pd(mA[4].d[2],0.),_mm_setr_pd(mA[4].d[3],0.),
               _mm_setr_pd(mA[4].d[4],0.));
      }

      // Matrix inversion by blocks   
      DMatrix5x5 Invert(){
         DMatrix2x2 A(GetV(0,0),GetV(0,1));
         DMatrix3x2 C(GetV(1,0),GetV(1,1),GetV(2,0),GetV(2,1));
         DMatrix3x2 CAinv=C*A.Invert();
         DMatrix3x3 D(GetV(1,2),GetV(1,3),GetV(1,4),GetV(2,2),GetV(2,3),GetV(2,4));
         DMatrix2x3 B(GetV(0,2),GetV(0,3),GetV(0,4));
         DMatrix2x3 BDinv=B*D.Invert();
         DMatrix2x2 AA=(A-BDinv*C).Invert();
         DMatrix3x3 DD=(D-CAinv*B).Invert();
         return DMatrix5x5(AA,-AA*BDinv,-DD*CAinv,DD);
      }
      // Matrix inversion by blocks for a symmetric matrix
      DMatrix5x5 InvertSym(){
         DMatrix2x2 A(GetV(0,0),GetV(0,1));
         DMatrix3x2 C(GetV(1,0),GetV(1,1),GetV(2,0),GetV(2,1));
         DMatrix3x2 CAinv=C*A.Invert();
         DMatrix3x3 D(GetV(1,2),GetV(1,3),GetV(1,4),GetV(2,2),GetV(2,3),GetV(2,4));
         DMatrix2x3 B(GetV(0,2),GetV(0,3),GetV(0,4));
         DMatrix2x3 BDinv=B*D.InvertSym();
         DMatrix2x2 AA=(A-BDinv*C).Invert();
         DMatrix3x3 DD=(D-CAinv*B).InvertSym();
         return DMatrix5x5(AA,-AA*BDinv,-DD*CAinv,DD);
      }



      // Zero out the matrix
      DMatrix5x5 Zero(){
         for (unsigned int i=0;i<5;i++){
            for (unsigned int j=0;j<3;j++){
               mA[i].v[j]=_mm_setzero_pd();
            }
         }
         return *this;
      }

      // Scale 5x5 matrix by a floating point number
      DMatrix5x5 operator*=(const double c){
         ALIGNED_16_BLOCK_WITH_PTR(__m128d, 1, p)
            __m128d &scale=p[0];
         scale=_mm_set1_pd(c);

         for (unsigned int i=0;i<5;i++){
            for (unsigned int j=0;j<3;j++){
               mA[i].v[j]=_mm_mul_pd(scale,mA[i].v[j]);
            }
         }
         return *this;
      }

      void Print(){
         cout << "DMatrix5x5:" <<endl;
         cout << "     |      0    |      1    |      2    |      3    |      4    |" <<endl;
         cout << "----------------------------------------------------------------------" <<endl;

         for (unsigned int i=0;i<5;i++){
            cout <<"   "<< i << " |";
            for (unsigned int j=0;j<5;j++){
               cout << setw(11)<<setprecision(4)<<mA[j].d[i] <<" "; 
            } 
            cout << endl;
         }      
      }

   private:
      union dvec{
         __m128d v[3];
         double d[6];
      };
      ALIGNED_16_BLOCK(union dvec, 5, mA)
};

// Scale 5x5 matrix by a floating point number
inline DMatrix5x5 operator*(const double c,const DMatrix5x5 &M){
   ALIGNED_16_BLOCK_WITH_PTR(__m128d, 1, p)
      __m128d &scale=p[0];
   scale=_mm_set1_pd(c);

#define SCALE(i,j) _mm_mul_pd(scale,M.GetV((i),(j)))
   return DMatrix5x5(SCALE(0,0),SCALE(0,1),SCALE(0,2),SCALE(0,3),SCALE(0,4),
         SCALE(1,0),SCALE(1,1),SCALE(1,2),SCALE(1,3),SCALE(1,4),
         SCALE(2,0),SCALE(2,1),SCALE(2,2),SCALE(2,3),SCALE(2,4));

}

#endif


