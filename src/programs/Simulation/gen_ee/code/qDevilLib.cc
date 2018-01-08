#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include "TRandom.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "qDevilLib.h"

// define overloaded + (plus) operator
Complx Complx::operator+ (const Complx& c) const
{
  Complx result;
  result.real = (this->real + c.real);
  result.imag = (this->imag + c.imag);
  return result;
}
// define overloaded - (minus) operator
Complx Complx::operator- (const Complx& c) const
{
  Complx result;
  result.real = (this->real - c.real);
  result.imag = (this->imag - c.imag);
  return result;
}

// define overloaded * (mult) operator
Complx Complx::operator* (const Complx& c) const
{
  Complx result;
  result.real = (this->real * c.real - this->imag * c.imag);
  result.imag = (this->real * c.imag + c.real * this->imag);
  return result;
}

// define overloaded = (equal) operator
Complx Complx::operator= (const Complx& c)
{
  this->real = c.real;
  this->imag = c.imag;

  return *this;
}

void Complx::Show(void)
{
  cout<<real<<"+"<<imag<<"i"<<endl;
}
void Complx::Set(double rVal,double iVal){real = rVal;imag=iVal;};

//////////////////////
// define overloaded + (plus) operator
QedElement QedElement::operator+ (const QedElement& q) const
{
  QedElement result;
  //START SCALAR + SCALAR
  if (this->QedType == "scalar" && q.QedType == "scalar") {
    result.SetZeroNone();
    if (this->lIndexName == "none" && q.lIndexName == "none"){
      result.scalar = (this->scalar + q.scalar);
      result.QedType = "scalar";
    }
    if (this->lIndexName != "none" && q.lIndexName == "none"){
      result.lIndexName = this->lIndexName;
      result.lIndexPosition = this->lIndexPosition;
      result.scalar0 = (this->scalar0 + q.scalar);
      result.scalarX = (this->scalarX + q.scalar);
      result.scalarY = (this->scalarY + q.scalar);
      result.scalarZ = (this->scalarZ + q.scalar);
      result.scalar5 = (this->scalar5 + q.scalar);
      result.QedType = "scalar";
      result.lIndexName = this->lIndexName;
      result.lIndexName = this->lIndexPosition;
    }
    if (this->lIndexName == "none" && q.lIndexName != "none"){
      result.lIndexName = q.lIndexName;
      result.lIndexPosition = q.lIndexPosition;
      result.scalar0 = (this->scalar + q.scalar0);
      result.scalarX = (this->scalar + q.scalarX);
      result.scalarY = (this->scalar + q.scalarY);
      result.scalarZ = (this->scalar + q.scalarZ);
      result.scalar5 = (this->scalar + q.scalar5);
      result.QedType = "scalar";
      result.lIndexName = q.lIndexName;
      result.lIndexName = q.lIndexPosition;
    }
    if (this->lIndexName != "none" && q.lIndexName != "none"){
      result.lIndexName = this->lIndexName;
      result.lIndexPosition = this->lIndexPosition;
      result.scalar0 = (this->scalar0 + q.scalar0);
      result.scalarX = (this->scalarX + q.scalarX);
      result.scalarY = (this->scalarY + q.scalarY);
      result.scalarZ = (this->scalarZ + q.scalarZ);
      result.scalar5 = (this->scalar5 + q.scalar5);
      result.QedType = "scalar";
      result.lIndexName = this->lIndexName;
      result.lIndexName = this->lIndexPosition;
    }
  }  //FINISHED SCALAR + SCALAR
  //START SCALAR + MATRIX
  if (this->QedType == "scalar" && q.QedType == "matrix") {
    result.SetZeroNone();
    if (this->lIndexName == "none" && q.lIndexName == "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  if (i != j) {
	    result.matrix[i][j] = (q.matrix[i][j]);
	  } else {
	    result.matrix[i][j] = (this->scalar + q.matrix[i][j]);
	  }
	}
      }
      result.QedType = "matrix";
    }
    if (this->lIndexName != "none" && q.lIndexName == "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  if (i != j) {
	    result.matrix0[i][j] = (q.matrix[i][j]);
	    result.matrixX[i][j] = (q.matrix[i][j]);
	    result.matrixY[i][j] = (q.matrix[i][j]);
	    result.matrixZ[i][j] = (q.matrix[i][j]);
	    result.matrix5[i][j] = (q.matrix[i][j]);
	  } else {
	    result.matrix0[i][j] = (this->scalar0 + q.matrix[i][j]);
	    result.matrixX[i][j] = (this->scalarX + q.matrix[i][j]);
	    result.matrixY[i][j] = (this->scalarY + q.matrix[i][j]);
	    result.matrixZ[i][j] = (this->scalarZ + q.matrix[i][j]);
	    result.matrix5[i][j] = (this->scalar5 + q.matrix[i][j]);
	  }
	}
      }
      result.QedType = "matrix";
      result.lIndexName = this->lIndexName;
      result.lIndexName = this->lIndexPosition;
    }
    if (this->lIndexName == "none" && q.lIndexName != "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  if (i != j) {
	    result.matrix0[i][j] = (q.matrix0[i][j]);
	    result.matrixX[i][j] = (q.matrixX[i][j]);
	    result.matrixY[i][j] = (q.matrixY[i][j]);
	    result.matrixZ[i][j] = (q.matrixZ[i][j]);
	    result.matrix5[i][j] = (q.matrix5[i][j]);
	  } else {
	    result.matrix0[i][j] = (this->scalar + q.matrix0[i][j]);
	    result.matrixX[i][j] = (this->scalar + q.matrixX[i][j]);
	    result.matrixY[i][j] = (this->scalar + q.matrixY[i][j]);
	    result.matrixZ[i][j] = (this->scalar + q.matrixZ[i][j]);
	    result.matrix5[i][j] = (this->scalar + q.matrix5[i][j]);
	  }
	}
      }
      result.QedType = "matrix";
      result.lIndexName = q.lIndexName;
      result.lIndexName = q.lIndexPosition;
    }
    if (this->lIndexName != "none" && q.lIndexName != "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  if (i != j) {
	    result.matrix0[i][j] = (q.matrix0[i][j]);
	    result.matrixX[i][j] = (q.matrixX[i][j]);
	    result.matrixY[i][j] = (q.matrixY[i][j]);
	    result.matrixZ[i][j] = (q.matrixZ[i][j]);
	    result.matrix5[i][j] = (q.matrix5[i][j]);
	  } else {
	    result.matrix0[i][j] = (this->scalar0 + q.matrix0[i][j]);
	    result.matrixX[i][j] = (this->scalarX + q.matrixX[i][j]);
	    result.matrixY[i][j] = (this->scalarY + q.matrixY[i][j]);
	    result.matrixZ[i][j] = (this->scalarZ + q.matrixZ[i][j]);
	    result.matrix5[i][j] = (this->scalar5 + q.matrix5[i][j]);
	  }
	}
      }
      result.QedType = "matrix";
      result.lIndexName = this->lIndexName;
      result.lIndexName = this->lIndexPosition;
    }
  }  //FINISHED SCALAR + MATRIX
  //START MATRIX + SCALAR
  if (this->QedType == "matrix" && q.QedType == "scalar") {
    result.SetZeroNone();
    if (this->lIndexName == "none" && q.lIndexName == "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  if (i != j) {
	    result.matrix[i][j] = (this->matrix[i][j]);
	  } else {
	    result.matrix[i][j] = (q.scalar + this->matrix[i][j]);
	  }
	}
      }
      result.QedType = "matrix";
    }
    if (this->lIndexName != "none" && q.lIndexName == "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  if (i != j) {
	    result.matrix0[i][j] = (this->matrix0[i][j]);
	    result.matrixX[i][j] = (this->matrixX[i][j]);
	    result.matrixY[i][j] = (this->matrixY[i][j]);
	    result.matrixZ[i][j] = (this->matrixZ[i][j]);
	    result.matrix5[i][j] = (this->matrix5[i][j]);
	  } else {
	    result.matrix0[i][j] = (q.scalar + this->matrix0[i][j]);
	    result.matrixX[i][j] = (q.scalar + this->matrixX[i][j]);
	    result.matrixY[i][j] = (q.scalar + this->matrixY[i][j]);
	    result.matrixZ[i][j] = (q.scalar + this->matrixZ[i][j]);
	    result.matrix5[i][j] = (q.scalar + this->matrix5[i][j]);
	  }
	}
      }
      result.QedType = "matrix";
      result.lIndexName = this->lIndexName;
      result.lIndexName = this->lIndexPosition;
    }
    if (this->lIndexName == "none" && q.lIndexName != "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  if (i != j) {
	    result.matrix0[i][j] = (this->matrix[i][j]);
	    result.matrixX[i][j] = (this->matrix[i][j]);
	    result.matrixY[i][j] = (this->matrix[i][j]);
	    result.matrixZ[i][j] = (this->matrix[i][j]);
	    result.matrix5[i][j] = (this->matrix[i][j]);
	  } else {
	    result.matrix0[i][j] = (q.scalar0 + this->matrix[i][j]);
	    result.matrixX[i][j] = (q.scalarX + this->matrix[i][j]);
	    result.matrixY[i][j] = (q.scalarY + this->matrix[i][j]);
	    result.matrixZ[i][j] = (q.scalarZ + this->matrix[i][j]);
	    result.matrix5[i][j] = (q.scalar5 + this->matrix[i][j]);
	  }
	}
      }
      result.QedType = "matrix";
      result.lIndexName = q.lIndexName;
      result.lIndexName = q.lIndexPosition;
    }
    if (this->lIndexName != "none" && q.lIndexName != "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  if (i != j) {
	    result.matrix0[i][j] = (this->matrix0[i][j]);
	    result.matrixX[i][j] = (this->matrixX[i][j]);
	    result.matrixY[i][j] = (this->matrixY[i][j]);
	    result.matrixZ[i][j] = (this->matrixZ[i][j]);
	    result.matrix5[i][j] = (this->matrix5[i][j]);
	  } else {
	    result.matrix0[i][j] = (q.scalar0 + this->matrix0[i][j]);
	    result.matrixX[i][j] = (q.scalarX + this->matrixX[i][j]);
	    result.matrixY[i][j] = (q.scalarY + this->matrixY[i][j]);
	    result.matrixZ[i][j] = (q.scalarZ + this->matrixZ[i][j]);
	    result.matrix5[i][j] = (q.scalar5 + this->matrix5[i][j]);
	  }
	}
      }
      result.QedType = "matrix";
      result.lIndexName = this->lIndexName;
      result.lIndexName = this->lIndexPosition;
    }
  }  //FINISHED MATRIX + SCALAR
  //START !SCALAR + !SCALAR
  if (this->QedType != "scalar" && q.QedType != "scalar") {
    result.SetZeroNone();
    if (this->lIndexName == "none" && q.lIndexName == "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix[i][j] = (this->matrix[i][j] +  q.matrix[i][j]);
	}
      }
      result.QedType = "matrix";
    }
    if (this->lIndexName != "none" && q.lIndexName == "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix0[i][j] = (this->matrix0[i][j] +  q.matrix[i][j]);
	  result.matrixX[i][j] = (this->matrixX[i][j] +  q.matrix[i][j]);
	  result.matrixY[i][j] = (this->matrixY[i][j] +  q.matrix[i][j]);
	  result.matrixZ[i][j] = (this->matrixZ[i][j] +  q.matrix[i][j]);
	  result.matrix5[i][j] = (this->matrix5[i][j] +  q.matrix[i][j]);
	}
      }
      result.QedType = "matrix";
      result.lIndexName = this->lIndexName;
      result.lIndexName = this->lIndexPosition;
    }
    if (this->lIndexName == "none" && q.lIndexName != "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix0[i][j] = (this->matrix[i][j] +  q.matrix0[i][j]);
	  result.matrixX[i][j] = (this->matrix[i][j] +  q.matrixX[i][j]);
	  result.matrixY[i][j] = (this->matrix[i][j] +  q.matrixY[i][j]);
	  result.matrixZ[i][j] = (this->matrix[i][j] +  q.matrixZ[i][j]);
	  result.matrix5[i][j] = (this->matrix[i][j] +  q.matrix5[i][j]);
	}
      }
      result.QedType = "matrix";
      result.lIndexName = q.lIndexName;
      result.lIndexName = q.lIndexPosition;
    }
    if (this->lIndexName != "none" && q.lIndexName != "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix0[i][j] = (this->matrix0[i][j] +  q.matrix0[i][j]);
	  result.matrixX[i][j] = (this->matrixX[i][j] +  q.matrixX[i][j]);
	  result.matrixY[i][j] = (this->matrixY[i][j] +  q.matrixY[i][j]);
	  result.matrixZ[i][j] = (this->matrixZ[i][j] +  q.matrixZ[i][j]);
	  result.matrix5[i][j] = (this->matrix5[i][j] +  q.matrix5[i][j]);
	}
      }
      result.QedType = "matrix";
      result.lIndexName = this->lIndexName;
      result.lIndexName = this->lIndexPosition;
    }
  }//FINISHED !SCALAR + !SCALAR
  return result;
} //FINISHED DEFINING +

// define overloaded - (minus) operator
QedElement QedElement::operator- (const QedElement& q) const
{
  QedElement result;
  //START SCALAR - SCALAR
  Complx Zero(0.0,0.0);
  if (this->QedType == "scalar" && q.QedType == "scalar") {
    result.SetZeroNone();
    if (this->lIndexName == "none" && q.lIndexName == "none"){
      result.scalar = (this->scalar - q.scalar);
      result.QedType = "scalar";
    }
    if (this->lIndexName != "none" && q.lIndexName == "none"){
      result.lIndexName = this->lIndexName;
      result.lIndexPosition = this->lIndexPosition;
      result.scalar0 = (this->scalar0 - q.scalar);
      result.scalarX = (this->scalarX - q.scalar);
      result.scalarY = (this->scalarY - q.scalar);
      result.scalarZ = (this->scalarZ - q.scalar);
      result.scalar5 = (this->scalar5 - q.scalar);
      result.QedType = "scalar";
      result.lIndexName = this->lIndexName;
      result.lIndexName = this->lIndexPosition;
    }
    if (this->lIndexName == "none" && q.lIndexName != "none"){
      result.lIndexName = q.lIndexName;
      result.lIndexPosition = q.lIndexPosition;
      result.scalar0 = (this->scalar - q.scalar0);
      result.scalarX = (this->scalar - q.scalarX);
      result.scalarY = (this->scalar - q.scalarY);
      result.scalarZ = (this->scalar - q.scalarZ);
      result.scalar5 = (this->scalar - q.scalar5);
      result.QedType = "scalar";
      result.lIndexName = q.lIndexName;
      result.lIndexName = q.lIndexPosition;
    }
    if (this->lIndexName != "none" && q.lIndexName != "none"){
      result.lIndexName = this->lIndexName;
      result.lIndexPosition = this->lIndexPosition;
      result.scalar0 = (this->scalar0 - q.scalar0);
      result.scalarX = (this->scalarX - q.scalarX);
      result.scalarY = (this->scalarY - q.scalarY);
      result.scalarZ = (this->scalarZ - q.scalarZ);
      result.scalar5 = (this->scalar5 - q.scalar5);
      result.QedType = "scalar";
      result.lIndexName = this->lIndexName;
      result.lIndexName = this->lIndexPosition;
    }
  }  //FINISHED SCALAR - SCALAR
  //START SCALAR - MATRIX
  if (this->QedType == "scalar" && q.QedType == "matrix") {
    result.SetZeroNone();
    if (this->lIndexName == "none" && q.lIndexName == "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  if (i != j) {
	    result.matrix[i][j] = (Zero - q.matrix[i][j]);
	  } else {
	    result.matrix[i][j] = (this->scalar - q.matrix[i][j]);
	  }
	}
      }
      result.QedType = "matrix";
    }
    if (this->lIndexName != "none" && q.lIndexName == "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  if (i != j) {
	    result.matrix0[i][j] = (Zero - q.matrix[i][j]);
	    result.matrixX[i][j] = (Zero - q.matrix[i][j]);
	    result.matrixY[i][j] = (Zero - q.matrix[i][j]);
	    result.matrixZ[i][j] = (Zero - q.matrix[i][j]);
	    result.matrix5[i][j] = (Zero - q.matrix[i][j]);
	  } else {
	    result.matrix0[i][j] = (this->scalar0 - q.matrix[i][j]);
	    result.matrixX[i][j] = (this->scalarX - q.matrix[i][j]);
	    result.matrixY[i][j] = (this->scalarY - q.matrix[i][j]);
	    result.matrixZ[i][j] = (this->scalarZ - q.matrix[i][j]);
	    result.matrix5[i][j] = (this->scalar5 - q.matrix[i][j]);
	  }
	}
      }
      result.QedType = "matrix";
      result.lIndexName = this->lIndexName;
      result.lIndexName = this->lIndexPosition;
    }
    if (this->lIndexName == "none" && q.lIndexName != "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  if (i != j) {
	    result.matrix0[i][j] = (Zero - q.matrix0[i][j]);
	    result.matrixX[i][j] = (Zero - q.matrixX[i][j]);
	    result.matrixY[i][j] = (Zero - q.matrixY[i][j]);
	    result.matrixZ[i][j] = (Zero - q.matrixZ[i][j]);
	    result.matrix5[i][j] = (Zero - q.matrix5[i][j]);
	  } else {
	    result.matrix0[i][j] = (this->scalar - q.matrix0[i][j]);
	    result.matrixX[i][j] = (this->scalar - q.matrixX[i][j]);
	    result.matrixY[i][j] = (this->scalar - q.matrixY[i][j]);
	    result.matrixZ[i][j] = (this->scalar - q.matrixZ[i][j]);
	    result.matrix5[i][j] = (this->scalar - q.matrix5[i][j]);
	  }
	}
      }
      result.QedType = "matrix";
      result.lIndexName = q.lIndexName;
      result.lIndexName = q.lIndexPosition;
    }
    if (this->lIndexName != "none" && q.lIndexName != "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  if (i != j) {
	    result.matrix0[i][j] = (Zero - q.matrix0[i][j]);
	    result.matrixX[i][j] = (Zero - q.matrixX[i][j]);
	    result.matrixY[i][j] = (Zero - q.matrixY[i][j]);
	    result.matrixZ[i][j] = (Zero - q.matrixZ[i][j]);
	    result.matrix5[i][j] = (Zero - q.matrix5[i][j]);
	  } else {
	    result.matrix0[i][j] = (this->scalar0 - q.matrix0[i][j]);
	    result.matrixX[i][j] = (this->scalarX - q.matrixX[i][j]);
	    result.matrixY[i][j] = (this->scalarY - q.matrixY[i][j]);
	    result.matrixZ[i][j] = (this->scalarZ - q.matrixZ[i][j]);
	    result.matrix5[i][j] = (this->scalar5 - q.matrix5[i][j]);
	  }
	}
      }
      result.QedType = "matrix";
      result.lIndexName = this->lIndexName;
      result.lIndexName = this->lIndexPosition;
    }
  }  //FINISHED SCALAR - MATRIX
  //START MATRIX - SCALAR
  if (this->QedType == "matrix" && q.QedType == "scalar") {
    result.SetZeroNone();
    if (this->lIndexName == "none" && q.lIndexName == "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  if (i != j) {
	    result.matrix[i][j] = (this->matrix[i][j]);
	  } else {
	    result.matrix[i][j] = (this->matrix[i][j] - q.scalar);
	  }
	}
      }
      result.QedType = "matrix";
    }
    if (this->lIndexName != "none" && q.lIndexName == "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  if (i != j) {
	    result.matrix0[i][j] = (this->matrix0[i][j]);
	    result.matrixX[i][j] = (this->matrixX[i][j]);
	    result.matrixY[i][j] = (this->matrixY[i][j]);
	    result.matrixZ[i][j] = (this->matrixZ[i][j]);
	    result.matrix5[i][j] = (this->matrix5[i][j]);
	  } else {
	    result.matrix0[i][j] = (this->matrix0[i][j] - q.scalar);
	    result.matrixX[i][j] = (this->matrixX[i][j] - q.scalar);
	    result.matrixY[i][j] = (this->matrixY[i][j] - q.scalar);
	    result.matrixZ[i][j] = (this->matrixZ[i][j] - q.scalar);
	    result.matrix5[i][j] = (this->matrix5[i][j] - q.scalar);
	  }
	}
      }
      result.QedType = "matrix";
      result.lIndexName = this->lIndexName;
      result.lIndexName = this->lIndexPosition;
    }
    if (this->lIndexName == "none" && q.lIndexName != "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  if (i != j) {
	    result.matrix0[i][j] = (this->matrix[i][j]);
	    result.matrixX[i][j] = (this->matrix[i][j]);
	    result.matrixY[i][j] = (this->matrix[i][j]);
	    result.matrixZ[i][j] = (this->matrix[i][j]);
	    result.matrix5[i][j] = (this->matrix[i][j]);
	  } else {
	    result.matrix0[i][j] = (this->matrix[i][j] - q.scalar0);
	    result.matrixX[i][j] = (this->matrix[i][j] - q.scalarX);
	    result.matrixY[i][j] = (this->matrix[i][j] - q.scalarY);
	    result.matrixZ[i][j] = (this->matrix[i][j] - q.scalarZ);
	    result.matrix5[i][j] = (this->matrix[i][j] - q.scalar5);
	  }
	}
      }
      result.QedType = "matrix";
      result.lIndexName = q.lIndexName;
      result.lIndexName = q.lIndexPosition;
    }
    if (this->lIndexName != "none" && q.lIndexName != "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  if (i != j) {
	    result.matrix0[i][j] = (this->matrix0[i][j]);
	    result.matrixX[i][j] = (this->matrixX[i][j]);
	    result.matrixY[i][j] = (this->matrixY[i][j]);
	    result.matrixZ[i][j] = (this->matrixZ[i][j]);
	    result.matrix5[i][j] = (this->matrix5[i][j]);
	  } else {
	    result.matrix0[i][j] = (this->matrix0[i][j] - q.scalar);
	    result.matrixX[i][j] = (this->matrixX[i][j] - q.scalar);
	    result.matrixY[i][j] = (this->matrixY[i][j] - q.scalar);
	    result.matrixZ[i][j] = (this->matrixZ[i][j] - q.scalar);
	    result.matrix5[i][j] = (this->matrix5[i][j] - q.scalar);
	  }
	}
      }
      result.QedType = "matrix";
      result.lIndexName = this->lIndexName;
      result.lIndexName = this->lIndexPosition;
    }
  }  //FINISHED MATRIX - SCALAR
  //START !SCALAR - !SCALAR
  if (this->QedType != "scalar" && q.QedType != "scalar") {
    result.SetZeroNone();
    if (this->lIndexName == "none" && q.lIndexName == "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix[i][j] = (this->matrix[i][j] -  q.matrix[i][j]);
	}
      }
      result.QedType = "matrix";
    }
    if (this->lIndexName != "none" && q.lIndexName == "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix0[i][j] = (this->matrix0[i][j] -  q.matrix[i][j]);
	  result.matrixX[i][j] = (this->matrixX[i][j] -  q.matrix[i][j]);
	  result.matrixY[i][j] = (this->matrixY[i][j] -  q.matrix[i][j]);
	  result.matrixZ[i][j] = (this->matrixZ[i][j] -  q.matrix[i][j]);
	  result.matrix5[i][j] = (this->matrix5[i][j] -  q.matrix[i][j]);
	}
      }
      result.QedType = "matrix";
      result.lIndexName = this->lIndexName;
      result.lIndexName = this->lIndexPosition;
    }
    if (this->lIndexName == "none" && q.lIndexName != "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix0[i][j] = (this->matrix[i][j] -  q.matrix0[i][j]);
	  result.matrixX[i][j] = (this->matrix[i][j] -  q.matrixX[i][j]);
	  result.matrixY[i][j] = (this->matrix[i][j] -  q.matrixY[i][j]);
	  result.matrixZ[i][j] = (this->matrix[i][j] -  q.matrixZ[i][j]);
	  result.matrix5[i][j] = (this->matrix[i][j] -  q.matrix5[i][j]);
	}
      }
      result.QedType = "matrix";
      result.lIndexName = q.lIndexName;
      result.lIndexName = q.lIndexPosition;
    }
    if (this->lIndexName != "none" && q.lIndexName != "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix0[i][j] = (this->matrix0[i][j] -  q.matrix0[i][j]);
	  result.matrixX[i][j] = (this->matrixX[i][j] -  q.matrixX[i][j]);
	  result.matrixY[i][j] = (this->matrixY[i][j] -  q.matrixY[i][j]);
	  result.matrixZ[i][j] = (this->matrixZ[i][j] -  q.matrixZ[i][j]);
	  result.matrix5[i][j] = (this->matrix5[i][j] -  q.matrix5[i][j]);
	}
      }
      result.QedType = "matrix";
      result.lIndexName = this->lIndexName;
      result.lIndexName = this->lIndexPosition;
    }
  }//FINISHED !SCALAR - !SCALAR
  return result;
} //FINISHED DEFINING -



// define overloaded * (mult) operator
QedElement QedElement::operator* (const QedElement& q) const
{
  QedElement result;
  //START SCALAR * SCALAR
  if (this->QedType == "scalar" && q.QedType == "scalar") {
    result.SetZeroNone();
    result.QedType = "scalar";
    if (this->lIndexName == "none" && q.lIndexName == "none"){
      result.scalar = (this->scalar * q.scalar);
    }
    if (this->lIndexName != "none" && q.lIndexName == "none"){
      result.lIndexName = this->lIndexName;
      result.lIndexPosition = this->lIndexPosition;
      result.scalar0 = (this->scalar0 * q.scalar);
      result.scalarX = (this->scalarX * q.scalar);
      result.scalarY = (this->scalarY * q.scalar);
      result.scalarZ = (this->scalarZ * q.scalar);
      result.scalar5 = (this->scalar5 * q.scalar);
    }
    if (this->lIndexName == "none" && q.lIndexName != "none"){
      result.lIndexName = q.lIndexName;
      result.lIndexPosition = q.lIndexPosition;
      result.scalar0 = (this->scalar * q.scalar0);
      result.scalarX = (this->scalar * q.scalarX);
      result.scalarY = (this->scalar * q.scalarY);
      result.scalarZ = (this->scalar * q.scalarZ);
      result.scalar5 = (this->scalar * q.scalar5);
    }
    if (this->lIndexName != "none" && q.lIndexName != "none"){
      result.scalar = (
		       this->scalar0 * q.scalar0 -
		       this->scalarX * q.scalarX -
		       this->scalarY * q.scalarY -
		       this->scalarZ * q.scalarZ);
    }
  } //FINISHED SCALAR * SCALAR
  //START SCALAR * !SCALAR
  if (this->QedType == "scalar" && q.QedType != "scalar") {
    result.SetZeroNone();
    result.QedType = q.QedType;
    if (this->lIndexName == "none" && q.lIndexName == "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix[i][j] = (this->scalar * q.matrix[i][j]);
	}
      }
    }
    if (this->lIndexName != "none" && q.lIndexName == "none"){
      result.lIndexName = this->lIndexName;
      result.lIndexPosition = this->lIndexPosition;
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix0[i][j] = (this->scalar0 * q.matrix[i][j]);    
	  result.matrixX[i][j] = (this->scalarX * q.matrix[i][j]);    
	  result.matrixY[i][j] = (this->scalarY * q.matrix[i][j]);    
	  result.matrixZ[i][j] = (this->scalarZ * q.matrix[i][j]);    
	  result.matrix5[i][j] = (this->scalar5 * q.matrix[i][j]);    
	}
      }
    }
    if (this->lIndexName == "none" && q.lIndexName != "none"){
      result.lIndexName = q.lIndexName;
      result.lIndexPosition = q.lIndexPosition;
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix0[i][j] = (this->scalar * q.matrix0[i][j]);    
	  result.matrixX[i][j] = (this->scalar * q.matrixX[i][j]);    
	  result.matrixY[i][j] = (this->scalar * q.matrixY[i][j]);    
	  result.matrixZ[i][j] = (this->scalar * q.matrixZ[i][j]);    
	  result.matrix5[i][j] = (this->scalar * q.matrix5[i][j]);    
	}
      }
    }
    if (this->lIndexName != "none" && q.lIndexName != "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix0[i][j] = (this->scalar0 * q.matrix0[i][j] -    
				  this->scalarX * q.matrixX[i][j] -
				  this->scalarY * q.matrixY[i][j] -    
				  this->scalarZ * q.matrixZ[i][j]);
	}
      }
    }
  } //FINISHED SCALAR * !SCALAR
  //START !SCALAR * SCALAR
  if (this->QedType != "scalar" && q.QedType == "scalar") {
    result.SetZeroNone();
    result.QedType = q.QedType;
    if (this->lIndexName == "none" && q.lIndexName == "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix[i][j] = (this->matrix[i][j] * q.scalar);
	}
      }
    }
    if (this->lIndexName != "none" && q.lIndexName == "none"){
      result.lIndexName = this->lIndexName;
      result.lIndexPosition = this->lIndexPosition;
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix0[i][j] = (this->matrix0[i][j] * q.scalar);    
	  result.matrixX[i][j] = (this->matrixX[i][j] * q.scalar);    
	  result.matrixY[i][j] = (this->matrixY[i][j] * q.scalar);    
	  result.matrixZ[i][j] = (this->matrixZ[i][j] * q.scalar);    
	  result.matrix5[i][j] = (this->matrix5[i][j] * q.scalar);    
	}
      }
    }
    if (this->lIndexName == "none" && q.lIndexName != "none"){
      result.lIndexName = q.lIndexName;
      result.lIndexPosition = q.lIndexPosition;
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix0[i][j] = (this->matrix[i][j] * q.scalar0);    
	  result.matrixX[i][j] = (this->matrix[i][j] * q.scalarX);    
	  result.matrixY[i][j] = (this->matrix[i][j] * q.scalarY);    
	  result.matrixZ[i][j] = (this->matrix[i][j] * q.scalarZ);    
	  result.matrix5[i][j] = (this->matrix[i][j] * q.scalar5);    
	}
      }
    }
    if (this->lIndexName != "none" && q.lIndexName != "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix[i][j] = (this->matrix0[i][j] * q.scalar0 -    
				  this->matrixX[i][j] * q.scalarX -
				  this->matrixY[i][j] * q.scalarY -    
				  this->matrixZ[i][j] * q.scalarZ);     
	}
      }
    }
  } //FINISHED !SCALAR * SCALAR
  //START !SCALAR * !SCALAR
  if (this->QedType != "scalar" && q.QedType != "scalar") {
    Complx zero(0.0,0.0);
    result.SetZeroNone();
    if (this->lIndexName == "none" && q.lIndexName == "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix[i][j] = zero;
	  for (int k=0; k<4; k++) {
	    result.matrix[i][j] = result.matrix[i][j] + (this->matrix[i][k] * q.matrix[k][j]);
	  }
	}
      }
    }
    if (this->lIndexName != "none" && q.lIndexName == "none"){
      result.lIndexName = this->lIndexName;
      result.lIndexPosition = this->lIndexPosition;
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix0[i][j] = zero;
	  result.matrixX[i][j] = zero;
	  result.matrixY[i][j] = zero;
	  result.matrixZ[i][j] = zero;
	  result.matrix5[i][j] = zero;
	  for (int k=0; k<4; k++) {
	    result.matrix0[i][j] = result.matrix0[i][j] + (this->matrix0[i][k] * q.matrix[k][j]);
	    result.matrixX[i][j] = result.matrixX[i][j] + (this->matrixX[i][k] * q.matrix[k][j]);
	    result.matrixY[i][j] = result.matrixY[i][j] + (this->matrixY[i][k] * q.matrix[k][j]);
	    result.matrixZ[i][j] = result.matrixZ[i][j] + (this->matrixZ[i][k] * q.matrix[k][j]);
	    result.matrix5[i][j] = result.matrix5[i][j] + (this->matrix5[i][k] * q.matrix[k][j]);
	  }
	}
      }
    }
    if (this->lIndexName == "none" && q.lIndexName != "none"){
      result.lIndexName = q.lIndexName;
      result.lIndexPosition = q.lIndexPosition;
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix0[i][j] = zero;
	  result.matrixX[i][j] = zero;
	  result.matrixY[i][j] = zero;
	  result.matrixZ[i][j] = zero;
	  result.matrix5[i][j] = zero;
	  for (int k=0; k<4; k++) {
	    result.matrix0[i][j] = result.matrix0[i][j] + (this->matrix[i][k] * q.matrix0[k][j]);
	    result.matrixX[i][j] = result.matrixX[i][j] + (this->matrix[i][k] * q.matrixX[k][j]);
	    result.matrixY[i][j] = result.matrixY[i][j] + (this->matrix[i][k] * q.matrixY[k][j]);
	    result.matrixZ[i][j] = result.matrixZ[i][j] + (this->matrix[i][k] * q.matrixZ[k][j]);
	    result.matrix5[i][j] = result.matrix5[i][j] + (this->matrix[i][k] * q.matrix5[k][j]);
	  }
	}
      }
    }
    if (this->lIndexName != "none" && q.lIndexName != "none"){
      for (int i=0; i<4; i++) {
	for (int j=0; j<4; j++) {
	  result.matrix[i][j] = zero;
	  for (int k=0; k<4; k++) {
	    result.matrix[i][j] = result.matrix[i][j] + (this->matrix0[i][k] * q.matrix0[k][j] -
				   this->matrixX[i][k] * q.matrixX[k][j] -
				   this->matrixY[i][k] * q.matrixY[k][j] -
				   this->matrixZ[i][k] * q.matrixZ[k][j]);
	  }
	}
      }
    }    
    if (this->QedType == "matrix" && q.QedType == "matrix") {
      result.QedType = "matrix";
    }
    if (this->QedType == "matrix" && q.QedType == "vector") {
      result.QedType = "vector";
    }
    if (this->QedType == "matrix" && q.QedType == "vectorT") {
      result.QedType = "matrix";
    }
    if (this->QedType == "vector" && q.QedType == "matrix") {
      result.QedType = "matrix";
    }
    if (this->QedType == "vectorT" && q.QedType == "matrix") {
      result.QedType = "vectorT";
    }
    if (this->QedType == "vector" && q.QedType == "vectorT") {
      result.QedType = "matrix";
    }
    if (this->QedType == "vectorT" && q.QedType == "vector") {
      result.QedType = "scalar";
      result.scalar = result.matrix[0][0];
      result.scalar0 = result.matrix0[0][0];
      result.scalarX = result.matrixX[0][0];
      result.scalarY = result.matrixY[0][0];
      result.scalarZ = result.matrixZ[0][0];
      result.scalar5 = result.matrix5[0][0];
    }
  }//FINISHED !SCALAR * !SCALAR

  return result;

}
// define overloaded = (equal) operator
QedElement QedElement::operator= (const QedElement& q) 
{
  QedElement result;
  this->scalar = q.scalar;
  this->scalar0 = q.scalar0;
  this->scalarX = q.scalarX;
  this->scalarY = q.scalarY;
  this->scalarZ = q.scalarZ;
  this->scalar5 = q.scalar5;

  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      this->matrix[i][j] = q.matrix[i][j];
      this->matrix0[i][j] = q.matrix0[i][j];
      this->matrixX[i][j] = q.matrixX[i][j];
      this->matrixY[i][j] = q.matrixY[i][j];
      this->matrixZ[i][j] = q.matrixZ[i][j];
      this->matrix5[i][j] = q.matrix5[i][j];
    }
  }
  this->QedType = q.QedType;
  this->lIndexName = q.lIndexName;
  this->lIndexPosition = q.lIndexPosition;
  return *this;
}
void QedElement::SetScalar(Complx sVal){
  SetZeroNone();
  QedType = "scalar";
  scalar = sVal;
}

void QedElement::SetZeroNone(){
  QedType = "none";
  lIndexName = "none";
  lIndexPosition = 0;
  scalar.Set(0,0);
  scalar0.Set(0,0);
  scalarX.Set(0,0);
  scalarY.Set(0,0);
  scalarZ.Set(0,0);
  scalar5.Set(0,0);
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      matrix[i][j].Set(0,0);
      matrix0[i][j].Set(0,0);
      matrixX[i][j].Set(0,0);
      matrixY[i][j].Set(0,0);
      matrixZ[i][j].Set(0,0);
      matrix5[i][j].Set(0,0);
    }
  }
}
void QedElement::SetGamma(string lIndexNameVal){
  SetZeroNone();
  lIndexName = lIndexNameVal;
  lIndexPosition = 1;
  QedType = "matrix";
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      matrix0[i][j] = gamma0[i][j];
      matrixX[i][j] = gammaX[i][j];
      matrixY[i][j] = gammaY[i][j];
      matrixZ[i][j] = gammaZ[i][j];
      matrix5[i][j] = gamma5[i][j];
    }
  }
}
void QedElement::Star(){
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      matrix[i][j] = matrix[i][j].Star();
      matrix0[i][j] = matrix0[i][j].Star();
      matrixX[i][j] = matrixX[i][j].Star();
      matrixY[i][j] = matrixY[i][j].Star();
      matrixZ[i][j] = matrixZ[i][j].Star();
      matrix5[i][j] = matrix5[i][j].Star();
    }
  }
  scalar = scalar.Star();
}
void QedElement::Transpose(){
  Complx tmpMatrix[4][4],
    tmpMatrix0[4][4],
    tmpMatrixX[4][4],
    tmpMatrixY[4][4],
    tmpMatrixZ[4][4],
    tmpMatrix5[4][4];
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      tmpMatrix[j][i] = matrix[i][j];
      tmpMatrix0[j][i] = matrix0[i][j];
      tmpMatrixX[j][i] = matrixX[i][j];
      tmpMatrixY[j][i] = matrixY[i][j];
      tmpMatrixZ[j][i] = matrixZ[i][j];
      tmpMatrix5[j][i] = matrix5[i][j];
    }
  }
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      matrix[i][j] = tmpMatrix[i][j];
      matrix0[i][j] = tmpMatrix0[i][j];
      matrixX[i][j] = tmpMatrixX[i][j];
      matrixY[i][j] = tmpMatrixY[i][j];
      matrixZ[i][j] = tmpMatrixZ[i][j];
      matrix5[i][j] = tmpMatrix5[i][j];
    }
  }

}
void QedElement::SetU(double px,double py,
	      double pz,double mass,
	      int spin){
  SetZeroNone();
  QedType = "vector";

  double energy,momentum,norm;
  Complx a00,a10,a20,a30;
  double tmpValR,tmpValI;
  momentum = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
  energy = sqrt(pow(momentum,2) + pow(mass,2));
  norm = sqrt(energy + mass);
  if (spin == 1) {
    a00.Set(norm,0.0);
    a10.Set(0.0,0.0);

    tmpValR = norm*pz/(energy + mass);
    a20.Set(tmpValR,0.0);

    tmpValR = norm*px/(energy + mass);
    tmpValI = norm*py/(energy + mass);
    a30.Set(tmpValR,tmpValI);

  }
  if (spin == -1) {
    a00.Set(0.0,0.0);
    a10.Set(norm,0.0);

    tmpValR = norm*px/(energy + mass);
    tmpValI = -norm*py/(energy + mass);
    a20.Set(tmpValR,tmpValI);
    
    tmpValR = -norm*pz/(energy + mass);
    a30.Set(tmpValR,0.0);

  }
  matrix[0][0] = a00;
  matrix[1][0] = a10;
  matrix[2][0] = a20;
  matrix[3][0] = a30;
}
void QedElement::SetV(double px,double py,
	      double pz,double mass,
	      int spin){
  SetZeroNone();

  double energy,momentum,norm;
  double tmpValR,tmpValI;
  Complx a00,a10,a20,a30;
  QedType = "vector";
  momentum = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
  energy = sqrt(pow(momentum,2) + pow(mass,2));
  norm = sqrt(energy + mass);
  if (spin == 1) {
    tmpValR = norm*px/(energy + mass);
    tmpValI = -norm*py/(energy + mass);
    a00.Set(tmpValR,tmpValI);   

    tmpValR = -norm*pz/(energy + mass);
    a10.Set(tmpValR,0.0);

    a20.Set(0.0,0.0);
    a30.Set(norm,0.0);
  }
  if (spin == -1) {
    tmpValR = -norm*pz/(energy + mass);
    a00.Set(tmpValR,0.0);

    tmpValR = -norm*px/(energy + mass);
    tmpValI = -norm*py/(energy + mass);
    a10.Set(tmpValR,tmpValI);

    a20.Set(-norm,0.0);
    a30.Set(0.0,0.0);
  }
  matrix[0][0] = a00;
  matrix[1][0] = a10;
  matrix[2][0] = a20;
  matrix[3][0] = a30;
}
void QedElement::SetUBar(double px,double py,
	      double pz,double mass,
	      int spin){
  Complx tmp[4][4];
  SetU(px,py,pz,mass,spin);
  Star();
  Transpose();
  //Multiply Udagger by gamma0
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      tmp[i][j] = 0.0;
      for (int k=0; k<4; k++) {
	tmp[i][j] = tmp[i][j] + (matrix[i][k] * gamma0[k][j]);
      }
      matrix[i][j] = tmp[i][j];
    }
  }
    
  QedType = "vectorT";
}
void QedElement::SetVBar(double px,double py,
	      double pz,double mass,
	      int spin){
  Complx tmp[4][4];
  SetV(px,py,pz,mass,spin);
  Star();
  Transpose();
  //Multiply Udagger by gamma0
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      tmp[i][j] = 0.0;
      for (int k=0; k<4; k++) {
	tmp[i][j] = tmp[i][j] + (matrix[i][k] * gamma0[k][j]);
      }
      matrix[i][j] = tmp[i][j];
    }
  }
  QedType = "vectorT";
}
void QedElement::SetMomentumSlash(double px,double py,
				  double pz,double mass){
  SetZeroNone();
  double energy,momentum;
  Complx cE,cPx,cPy,cPz;
  momentum = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
  energy = sqrt(pow(momentum,2) + pow(mass,2));
  cE.Set(energy,0.0);
  cPx.Set(-px,0.0);
  cPy.Set(-py,0.0);
  cPz.Set(-pz,0.0);
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
	matrix[i][j] = cE * gamma0[i][j] + 
	  cPx* gammaX[i][j] +
	  cPy* gammaY[i][j] +
	  cPz* gammaZ[i][j];
    }
  }
  QedType = "matrix";
}
void QedElement::SetEpsilonSlash(double cx,double cy,
				  double cz){
  SetZeroNone();
  Complx cCx,cCy,cCz;
  cCx.Set(-cx,0.0);
  cCy.Set(-cy,0.0);
  cCz.Set(-cz,0.0);
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
	matrix[i][j] = cCx*gammaX[i][j] +
	  cCy*gammaY[i][j] +
	  cCz*gammaZ[i][j];
    }
  }
  QedType = "matrix";
}
void QedElement::ShowMatrix(int matNumber){
  if (matNumber == -1) {
    double rVal,iVal;
    for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
	rVal = matrix[i][j].r();
	iVal = matrix[i][j].i();
	cout<<"("<<rVal<<"\t"<<iVal<<"i) ";
      }
      cout<<" "<<endl;
    }
  }
  if (matNumber == 0) {
    double rVal,iVal;
    for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
	rVal = matrix0[i][j].r();
	iVal = matrix0[i][j].i();
	cout<<"("<<rVal<<"\t"<<iVal<<"i) ";
      }
      cout<<" "<<endl;
    }
  }
  if (matNumber == 1) {
    double rVal,iVal;
    for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
	rVal = matrixX[i][j].r();
	iVal = matrixX[i][j].i();
	cout<<"("<<rVal<<"\t"<<iVal<<"i) ";
      }
      cout<<" "<<endl;
    }
  }
  if (matNumber == 2) {
    double rVal,iVal;
    for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
	rVal = matrixY[i][j].r();
	iVal = matrixY[i][j].i();
	cout<<"("<<rVal<<"\t"<<iVal<<"i) ";
      }
      cout<<" "<<endl;
    }
  }
  if (matNumber == 3) {
    double rVal,iVal;
    for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
	rVal = matrixZ[i][j].r();
	iVal = matrixZ[i][j].i();
	cout<<"("<<rVal<<"\t"<<iVal<<"i) ";
      }
      cout<<" "<<endl;
    }
  }
  if (matNumber == 4) {
    double rVal,iVal;
    for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
	rVal = matrix5[i][j].r();
	iVal = matrix5[i][j].i();
	cout<<"("<<rVal<<"\t"<<iVal<<"i) ";
      }
      cout<<" "<<endl;
    }
  }

  if (matNumber == 5) {
    double rVal,iVal;
    for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
	rVal = gamma0[i][j].r();
	iVal = gamma0[i][j].i();
	cout<<"("<<rVal<<"\t"<<iVal<<"i) ";
      }
      cout<<" "<<endl;
    }
  }

}
void QedElement::ShowAll(){
  cout<<"QedType = "<<QedType<<endl;
  cout<<"lIndexName = "<<lIndexName<<endl;
  cout<<"scalar = ";scalar.Show();
  cout<<"scalar0 = ";scalar0.Show();
  cout<<"scalarX = ";scalarX.Show();
  cout<<"scalarY = ";scalarY.Show();
  cout<<"scalarZ = ";scalarZ.Show();
  cout<<"scalar5 = ";scalar5.Show();
  cout<<"matrix = "<<endl;
  ShowMatrix(-1);
  cout<<"matrix0 = "<<endl;
  ShowMatrix(0);
  cout<<"matrixX = "<<endl;
  ShowMatrix(1);
  cout<<"matrixY = "<<endl;
  ShowMatrix(2);
  cout<<"matrixZ = "<<endl;
  ShowMatrix(3);
  cout<<"matrix5 = "<<endl;
  ShowMatrix(4);
}

double ampSqPT(int rxnType,int polDirection, TLorentzVector target, TLorentzVector beam,
                 TLorentzVector recoil,TLorentzVector q1,TLorentzVector q2){

  double sumFactor = 0.5;
  double massElectron = 0.51099907e-3; //Mass in GeV                                 
  double massTarget = target.Mag();

  //IF UNPOLARIZED SET SUM FACTOR
  if (polDirection == 0) {
    sumFactor = 0.25;
  }

  //Define the kinematics                                                       
  //1 -> initial photon                                                         
  //2 -> initial target particle                                                       
  //3 -> final electron                                                         
  //4 -> final positron                                                         
  //5 -> final target particle (recoil)                                                

  //Using LAB FRAME                                                             
  double p1  = beam.E();
  double E1  = p1;
  double px1 = 0.0;
  double py1 = 0.0;
  double pz1 = p1;

  double E2  = target.E();
  double px2 = target.Px();
  double py2 = target.Py();
  double pz2 = target.Pz();
  
  double E3  = q1.E();
  double px3 = q1.Px();
  double py3 = q1.Py();
  double pz3 = q1.Pz();
  
  double E4  = q2.E();
  double px4 = q2.Px();
  double py4 = q2.Py();
  double pz4 = q2.Pz();
  
  double E5  = recoil.E();
  double px5 = recoil.Px();
  double py5 = recoil.Py();
  double pz5 = recoil.Pz();

  TLorentzVector moTransfer = recoil - target;

  if (p1 <= 2*massElectron) return 0.0;

  double q21,q22,q23,q24,q25,q26,q27,q28;
  double p1p2,p1p3,p1p4,p1p5,p2p3,p2p5,p3p4,p4p5;
  double kinFac1,kinFac2,kinFac3,kinFac4;
  double kinFac5,kinFac6,kinFac7,kinFac8;
  double m34,m35,m1,m2,m3,m4,m5;
  m1 = 0.0;m2 = massTarget;m3 = massElectron;m4 = massElectron;m5 = massTarget;

  int spin1,spin2,spin3,spin4,spin5;
  QedElement u2;
  QedElement u3Bar;
  QedElement v4;
  QedElement u5Bar;
  QedElement ep1Slash;
  QedElement ep2Slash;
  QedElement p1Slash;
  QedElement p2Slash;
  QedElement p3Slash;
  QedElement p4Slash;
  QedElement p5Slash;
  QedElement gammaMu1;
  QedElement gammaMu2;
  gammaMu1.SetGamma("mu");
  gammaMu2.SetGamma("mu");
  QedElement matrixElementTmp1;
  QedElement matrixElementTmp2;
  QedElement matrixElementTmp3;
  QedElement matrixElementTmp4;
  QedElement matrixElementTmp5;
  QedElement matrixElementTmp6;
  QedElement matrixElementTmp7;
  QedElement matrixElementTmp8;
  QedElement matrixElementTmp;
  Complx matrixElement(0.0,0.0);

  QedElement mValElectron;
  Complx mCmplxElectron(massElectron,0.0);
  mValElectron.SetScalar(mCmplxElectron);

  QedElement mValTarget;
  Complx mCmplxTarget(massTarget,0.0);
  mValTarget.SetScalar(mCmplxTarget);

  Complx cTmp(0.0);
  double mAmplSq = 0.0;
  double mAmplSqSum = 0.0;

  //Create the invariant mass square of 3,4 and 3,5 systems                     
  m34 = sqrt(pow(E3+E4,2)-pow(px3+px4,2)-pow(py3+py4,2)-pow(pz3+pz4,2));
  m35 = sqrt(pow(E3+E5,2)-pow(px3+px5,2)-pow(py3+py5,2)-pow(pz3+pz5,2));

  //Define some dot products                                                    
  p1p2 = E1*E2 - px1*px2 - py1*py2 - pz1*pz2;
  p1p3 = E1*E3 - px1*px3 - py1*py3 - pz1*pz3;
  p1p4 = E1*E4 - px1*px4 - py1*py4 - pz1*pz4;
  p1p5 = E1*E5 - px1*px5 - py1*py5 - pz1*pz5;
  p2p3 = E2*E3 - px2*px3 - py2*py3 - pz2*pz3;
  p2p5 = E2*E5 - px2*px5 - py2*py5 - pz2*pz5;
  p3p4 = E3*E4 - px3*px4 - py3*py4 - pz3*pz4;
  p4p5 = E4*E5 - px4*px5 - py4*py5 - pz4*pz5;

  //Define q2 for each diagram                                                  
  q21 = pow(massTarget,2) + pow(massTarget,2) - 2*p2p5;
  q22 = q21;
  q23 = pow(massElectron,2) + pow(massTarget,2) - 2*p2p3;
  q24 = q23;

  q25 = 2*pow(massElectron,2) + 2*p3p4;
  q26 = q25;
  q27 = pow(massTarget,2) + pow(massElectron,2) + 2*p4p5;
  q28 = q27;

  //Define kinematic factors for each diagram                                   
  kinFac1 = -1.0/(2.0*p1p3);
  kinFac2 = -1.0/(2.0*p1p4);
  //Swap legs 3 and 5                                                           
  kinFac3 = -1.0/(2.0*p1p5);
  kinFac4 = -1.0/(2.0*p1p4);

  kinFac5 = 1.0/(2.0*p1p2);
  kinFac6 = -1.0/(2.0*p1p5);
  //Swap legs 3 and 5                                                           
  kinFac7 = 1.0/(2.0*p1p2);
  kinFac8 = -1.0/(2.0*p1p3);

  int Z = 1;
  double V1,V2,V3,V4,V5,V6,V7,V8;
  V1 = Z*1.0/q21;
  V2 = Z*1.0/q22;
  V3 = Z*1.0/q23;
  V4 = Z*1.0/q24;
  V5 = Z*1.0/q25;
  V6 = Z*1.0/q26;
  V7 = Z*1.0/q27;
  V8 = Z*1.0/q28;

  //Set up multiplicative factors for each diagram                              
  QedElement kf1,kf2,kf3,kf4,kf5,kf6,kf7,kf8;
  Complx kf1Val(kinFac1*V1,0.0);
  Complx kf2Val(kinFac2*V2,0.0);
  Complx kf3Val(kinFac3*V3,0.0);
  Complx kf4Val(kinFac4*V4,0.0);
  Complx kf5Val(kinFac5*V5,0.0);
  Complx kf6Val(kinFac6*V6,0.0);
  Complx kf7Val(kinFac7*V7,0.0);
  Complx kf8Val(kinFac8*V8,0.0);

  kf1.SetScalar(kf1Val);
  kf2.SetScalar(kf2Val);
  kf3.SetScalar(kf3Val);
  kf4.SetScalar(kf4Val);
  kf5.SetScalar(kf5Val);
  kf6.SetScalar(kf6Val);
  kf7.SetScalar(kf7Val);
  kf8.SetScalar(kf8Val);

  mAmplSq = 0.0;

  //FILL SOME SLASH TERMS
  p1Slash.SetMomentumSlash(px1,py1,pz1,0.0);
  p2Slash.SetMomentumSlash(px2,py2,pz2,massTarget);
  p3Slash.SetMomentumSlash(px3,py3,pz3,massElectron);
  p4Slash.SetMomentumSlash(px4,py4,pz4,massElectron);
  p5Slash.SetMomentumSlash(px5,py5,pz5,massTarget);

  int phoSpin1 = -1;
  int phoSpin2 = 1;
  if (polDirection == 2){
    phoSpin1 = -1; //Y-DIRECTION
    phoSpin2 = -1; //Y-DIRECTION 
  }
  if (polDirection == 1){
    phoSpin1 = 1; //X-DIRECTION
    phoSpin2 = 1; //X-DIRECTION
  }

  //Sum over spin assignments
  for (spin1=phoSpin1;spin1<=phoSpin2;spin1+=2){//SUMMING OVER INITIAL PHOTON SPINS
    for (spin2=-1;spin2<=+1;spin2+=2){ //SUMMING OVER INITIAL TARGET SPINS
      for (spin3=-1;spin3<=+1;spin3+=2){//SUMMING OVER ELECTRON SPINS
        for (spin4=-1;spin4<=+1;spin4+=2){//SUMMING OVER POSITRON SPINS
          for (spin5=-1;spin5<=+1;spin5+=2){//SUMMING OVER RECOIL SPINS
            u2.SetU(px2,py2,pz2,massTarget,spin2);
            u3Bar.SetUBar(px3,py3,pz3,massElectron,spin3);
            v4.SetV(px4,py4,pz4,massElectron,spin4);
            u5Bar.SetUBar(px5,py5,pz5,massTarget,spin5);

            if (spin1 == -1) ep1Slash.SetEpsilonSlash(0,1,0); //PHOTON POLARIZED IN Y-DIRECTION
            if (spin1 ==  1) ep1Slash.SetEpsilonSlash(1,0,0); //PHOTON POLARIZED IN X-DIRECTION

            //create the matrix elements
            matrixElementTmp1 = u3Bar*ep1Slash*((p3Slash - p1Slash) + mValElectron)*gammaMu1*v4*
              (u5Bar*gammaMu2*u2);
            matrixElementTmp2 = u3Bar*gammaMu1*((p1Slash - p4Slash) + mValElectron)*ep1Slash*v4*
              (u5Bar*gammaMu2*u2);

	    if (rxnType == 3) {
	      //Swap the 5 and 3 electron legs                                    
	      matrixElementTmp3 = u5Bar*ep1Slash*((p5Slash - p1Slash) + mValElectron)*gammaMu1*v4*
		(u3Bar*gammaMu2*u2);
	      matrixElementTmp4 = u5Bar*gammaMu1*((p1Slash - p4Slash) + mValElectron)*ep1Slash*v4*
		(u3Bar*gammaMu2*u2);
	    }

            //Compton like diagrams
            matrixElementTmp5 = u5Bar*gammaMu1*(p1Slash + p2Slash + mValTarget)*ep1Slash*u2*
	      (u3Bar*gammaMu2*v4);
            matrixElementTmp6 = u5Bar*ep1Slash*(p5Slash - p1Slash + mValTarget)*gammaMu1*u2*
	      (u3Bar*gammaMu2*v4);

	    if (rxnType == 3) {	    
	      //Swap the 5 and 3 electron legs                                                   
	      matrixElementTmp7 = u3Bar*gammaMu1*(p1Slash + p2Slash + mValTarget)*ep1Slash*u2*
		(u5Bar*gammaMu2*v4);
	      matrixElementTmp8 = u3Bar*ep1Slash*(p3Slash - p1Slash + mValTarget)*gammaMu1*u2*
		(u5Bar*gammaMu2*v4);
	    }

            //Full caluculation                                 
	    if (rxnType == 2) {	    
	      matrixElementTmp =
              kf1*matrixElementTmp1 + kf2*matrixElementTmp2;

          //kf1*matrixElementTmp1 + kf2*matrixElementTmp2
          //+ kf5*matrixElementTmp5 + kf6*matrixElementTmp6;
	    }
	    if (rxnType == 3) {	    
	      matrixElementTmp =
		kf1*matrixElementTmp1 + kf2*matrixElementTmp2
		+ kf5*matrixElementTmp3 + kf6*matrixElementTmp4
		+ kf5*matrixElementTmp5 + kf6*matrixElementTmp6
		+ kf5*matrixElementTmp7 + kf6*matrixElementTmp7;
	    }

            //Convert from qedElement type to Complx type                                      
            matrixElement = matrixElementTmp.GetScalar();
            //Need to get squared absolute value                                               
            cTmp = matrixElement*matrixElement.Star();
            //Need to convert from Complx to double and add up amplitudes:                     
            mAmplSq += cTmp.r();

          }
        }
      }
    }
  } //Finished with the spin sum                                          
  
  mAmplSqSum = mAmplSq*sumFactor;

  if (mAmplSqSum >= 0 || mAmplSqSum <=0){
    //do nothing
  } else {
    cout<<"!!!!!Bad Event!!!!!"<<endl;
    mAmplSqSum = 0.0;
  }

  return mAmplSqSum;
}


//////////////////////////////////////////////

int PhasePT::Gen(TRandom3 *random){
  Double_t event_E0;
  Double_t event_Epos;
  Double_t event_phi12;
  Double_t event_Mpair;
  Double_t event_qR2;
  Double_t event_phiR;
  Double_t event_thetaR;
  Double_t PIval = 2*atan2(1,0);
  Double_t mElectron = 0.51099907e-3;
  
  event_weight = 1.0;
  event_E0 = beam.E();
  double kin = event_E0;
  
  // generate E+ uniform on [0,E0]                                                          
  event_Epos = random->Uniform(event_E0);
  event_weight *= event_E0;

  // generate phi12 uniform on [0,2pi]                                                      
  event_phi12 = random->Uniform(2*PIval);
  event_weight *= 2*PIval;

  // generate phiR uniform on [0,2pi]                                                       
  event_phiR = random->Uniform(2*PIval);
  event_weight *= 2*PIval;
  
  // generate Mpair with weight (1/M) / (Mcut^2 + M^2)                                      
  Double_t Mmin=2*mElectron;
  Double_t Mcut=m12Cut;
  Double_t um0 = 1 + pow(Mcut/Mmin,2);
  Double_t um = pow(um0,random->Uniform(1));
  event_Mpair = Mcut/sqrt(um-1);
  event_weight *= event_Mpair*(pow(Mcut,2)+pow(event_Mpair,2))
    *log(um0)/(2*pow(Mcut,2));
  
  // generate qR^2 with weight (1/qR^2) / sqrt(qRcut^2 + qR^2)                              
  Double_t qRmin = pow(event_Mpair,2)/(2*event_E0);
  Double_t qRcut = rCut;                         
  Double_t uq0 = qRmin/(qRcut+sqrt(pow(qRcut,2)+pow(qRmin,2)));
  Double_t uq = pow(uq0,random->Uniform(1));
  event_qR2 = pow(2*qRcut*uq/pow(1-pow(uq,2),2),2);
  event_weight *= event_qR2*sqrt(1+event_qR2/pow(qRcut,2))
    *(-2*log(uq0));
  
  // overall measure Jacobian factor                                                       
  event_weight *= event_Mpair/(2*event_E0);
  
  // compute recoil polar angle thetaR                                                     
  Double_t E3 = sqrt(event_qR2 + target.Mag2());
  Double_t costhetaR = (pow(event_Mpair,2)/2 + (kin + target.Mag())*(E3-target.Mag())
			)/(kin*sqrt(event_qR2));
  
  int killVal = 0;
  if (fabs(costhetaR) > 1) {
    //cout << "no kinematic solution because |costhetaR| > 1" << endl;                    
    event_thetaR = 99;
    killVal = 1;
  }
  else {
    event_thetaR = acos(costhetaR);
  }
  

  //RECOIL:                                                                                
  double moR  = sqrt(event_qR2);
  double eR   = sqrt(pow(moR,2) + target.Mag2());
  double moRx = moR*cos(event_phiR)*sin(event_thetaR);
  double moRy = moR*sin(event_phiR)*sin(event_thetaR);
  double moRz = moR*cos(event_thetaR);
  
  //POSITRON:                                                                              
  double ePos    = event_Epos;
  Double_t k12star2 = pow(event_Mpair/2,2) - pow(mElectron,2); //mo^2 of e+ in Mpair rest frame
  if (k12star2 < 0) {
    //cout << "no kinematic solution because k12star2 < 0" << endl;                        
    if (killVal == 0) killVal = 2;
    //return 0;                                                                            
  }

  Double_t k12star = sqrt(k12star2); //mo of e+ in Mpair rest frame                        
  Double_t E12 = event_E0 + target.Mag() - eR; //Energy of pair in lab frame                  
  Double_t q12mag=sqrt(pow(E12,2) - pow(event_Mpair,2)); //momentum of pair in lab frame   
  Double_t costhetastar=(ePos-E12/2)*event_Mpair/(k12star*q12mag); //????                  
  if (fabs(costhetastar) > 1) {
    //cout << "no kinematic solution because |costhetastar| > 1" << endl;                  
    if (killVal == 0) killVal = 3;
    if (killVal == 1) killVal = 4;
    //return 0;                                                                            
  }
  
  Double_t sinthetastar = sqrt(1-pow(costhetastar,2));
  Double_t moK12starX   = k12star*sinthetastar*cos(event_phi12);
  Double_t moK12starY   = k12star*sinthetastar*sin(event_phi12);
  Double_t moK12starZ   = k12star*costhetastar;
  
  TLorentzVector q1star(-moK12starX,-moK12starY,-moK12starZ,event_Mpair/2.0);
  TLorentzVector q2star(moK12starX,moK12starY,moK12starZ,event_Mpair/2.0);

  if (killVal == 1110){
    cout<<"q1star P,E,M = "<<q1star.P()<<", "<<q1star.E()<<", "<<q1star.Mag()<<endl;
    cout<<"sinthetastar = "<<sinthetastar<<endl;
    cout<<"costhetastar = "<<costhetastar<<endl;
    cout<<"moK12star X,Y,Z = "<<moK12starX<<", "<<moK12starY<<", "<<moK12starZ<<endl;
    cout<<"ePos, E12, event_Mpair, k12star, q12mag = "<<ePos<<", "<<E12<<", "
	<<event_Mpair<<", "<<k12star<<", "<<q12mag<<endl;
  }

  //DEFINE THE LAB FRAME 4-VECTORS: beam,target,recoil,q1,q2                               
  recoil.SetXYZT(moRx,moRy,moRz,eR);
  TLorentzVector cms(0.0,0.0,0.0,0.0);
  TLorentzVector m12Vec(0.0,0.0,0.0,0.0);
  cms = beam + target;
  m12Vec = cms - recoil;
  q1 = q1star;
  q2 = q2star;
  q1.Boost(m12Vec.BoostVector());
  q2.Boost(m12Vec.BoostVector());
  q12 = q1 + q2;
  
  q23 = q2 + recoil;

  if (killVal == 1) return -1;
  if (killVal == 2) return -2;
  if (killVal == 3) return -3;
  if (killVal == 4) return -4;

  if (m12Vec.Mag() < 0.0) return -5;
  
  double sVal = cms.Mag2();
  double e2Star = q12.Mag()/2.0;
  double e3Star = (sVal - q12.Mag2() - target.Mag2())/(2*q12.Mag());
  double p2Star = sqrt(pow(e2Star,2) - pow(mElectron,2));
  double p3Star = sqrt(pow(e3Star,2) - target.Mag2());
  double m23Max = pow(e2Star + e3Star,2) - pow(p2Star - p3Star,2);
  double m23Min = pow(e2Star + e3Star,2) - pow(p2Star + p3Star,2);

  if (q23.Mag2() > m23Max) return -6;         
  if (q23.Mag2() < m23Min) return -7;

  return 1;
}



