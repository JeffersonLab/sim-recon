// $Id$
//
//    File: DMatrix.h
// Created: Fri Apr  6 14:34:26 EDT 2007
// Creator: davidl (on Darwin swire-b241.jlab.org 8.9.0 powerpc)
//

/// Implement a linear algebra package
///
/// This file implements one of two definitions of the DMatrix
/// class depending upon whether USE_MATRIX_CLHEP is defined or
/// not. The DMatrix class itself does not do any of the matrix
/// operations and instead relies on either the ROOT package or
/// the CLHEP package to do all the work.
///
/// The idea here is to write code using DMatrix which povides
/// a layer which can be used to easily switch the underlying
/// package that actually implements to algorithms. This allows
/// codes to be compiled "ROOT-free" or "CLHEP-free" depending
/// on the end user's preference.
///
/// The end result is that authors can also write code using the
/// methods of either package's class. For instance, one can
/// use both the num_row() and GetNrows() methods to get the
/// number of rows. Both invert() and Invert(), (the CLHEP and
/// ROOT methods respectively), are available.
///
/// Since DMatrix inherits from a class of one of the underlying
/// packages, all the methods of that class are inherited. The
/// mapping must be done for the other package's equivalent class.
/// This makes it a bit dangerous in that if one develops
/// code using ROOT as the underlying package and calls ROOT
/// methods that are not implemented when CLHEP is used,
/// the code will not compile when switching to CLHEP.


#ifndef _DMatrix_
#define _DMatrix_

#include <iostream>

#ifndef USE_MATRIX_CLHEP
//=============================== ROOT ===============================
// If USE_MATRIX_CLHEP is not specified, then base DMatrix on
// the ROOT linear algebra package

#include <TMatrixD.h>

class DMatrix:public TMatrixD{
	public:
		DMatrix(){}
		DMatrix(int nrows, int ncols):TMatrixD(nrows, ncols){}
		DMatrix(const TMatrixD &m1):TMatrixD(m1){}
		virtual ~DMatrix(){}
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DMatrix";}
		
		// CLHEP methods
		inline double determinant(void) const{return Determinant();}
		inline int num_col(void) const{return GetNcols();}
		inline int num_row(void) const{return GetNrows();}
		inline void invert(int &ierr){Invert(); ierr=IsValid();}
		
		// DMatrix methods
		void Print(void);
};

inline DMatrix operator+(const DMatrix &m1, const DMatrix &m2){return static_cast<DMatrix>(static_cast<TMatrixD>(m1)+static_cast<TMatrixD>(m2));}
inline DMatrix operator-(const DMatrix &m1, const DMatrix &m2){return static_cast<DMatrix>(static_cast<TMatrixD>(m1)-static_cast<TMatrixD>(m2));}
inline DMatrix operator*(const DMatrix &m1, const DMatrix &m2){return static_cast<DMatrix>(static_cast<TMatrixD>(m1)*static_cast<TMatrixD>(m2));}
inline DMatrix operator/(const DMatrix &m1, const DMatrix &m2){return static_cast<DMatrix>(static_cast<TMatrixD>(m1)/static_cast<TMatrixD>(m2));}

//=============================== ROOT ===============================

#else // USE_MATRIX_CLHEP

//=============================== CLHEP ===============================
// If USE_MATRIX_CLHEP is specified, then base DMatrix on
// the CLHEP linear algebra package

#include <CLHEP/Matrix/Matrix.h>
using namespace CLHEP;

class DMatrix : public HepMatrix {
	public:
		DMatrix(){}
		DMatrix(int nrows, int ncols):HepMatrix(nrows, ncols){}
		DMatrix(const HepMatrix &m1):HepMatrix(m1){}
		virtual ~DMatrix(){}
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DMatrix";}
		
		// ROOT methods
		double E2Norm(void) const;
		inline int GetNcols(void) const{return num_col();}
		inline int GetNrows(void) const{return num_row();}
		inline double Determinant(void) const{return determinant();}
		inline void Invert(double *det){int ierr; invert(ierr); *det=ierr>0 ? -1:0.0;}
				
		// DMatrix methods
		void Print(void);

	protected:
		bool is_valid;
	
	private:

};

inline DMatrix operator+(const DMatrix &m1, const DMatrix &m2){return static_cast<DMatrix>(static_cast<HepMatrix>(m1)+static_cast<HepMatrix>(m2));}
inline DMatrix operator-(const DMatrix &m1, const DMatrix &m2){return static_cast<DMatrix>(static_cast<HepMatrix>(m1)-static_cast<HepMatrix>(m2));}
inline DMatrix operator*(const DMatrix &m1, const DMatrix &m2){return static_cast<DMatrix>(static_cast<HepMatrix>(m1)*static_cast<HepMatrix>(m2));}
inline DMatrix operator/(const DMatrix &m1, const DMatrix &m2){return static_cast<DMatrix>(static_cast<HepMatrix>(m1)/static_cast<HepMatrix>(m2));}

//=============================== CLHEP ===============================
#endif // USE_MATRIX_CLHEP


inline void DMatrix::Print(void)
{
	std::cout<<std::endl;
	std::cout<<"DMatrix("<<GetNrows()<<", "<<GetNcols()<<"):"<<std::endl;
	for(int row=0; row<GetNrows(); row++){
		for(int col=0; col<GetNcols(); col++){
			std::cout<<(*this)[row][col]<<"\t";
		}
		std::cout<<std::endl;
	}	
}

#endif // _DMatrix_

