

// This file contains C-callable routine signatures for the routines
// in calibDB.cc. It is specifically intended to give access to the 
// calibration DB to C and FORTRAN routines. This is needed since the
// calibration DB API is in C++.
//
// This file should contain no C++ code.

typedef struct {
  char str[128];
}mystr_t;


void initcalibdb_(char *bfield_type, char *bfield_map,int *runno);
void gufld_db_(float *r, float *B);
int GetCalib(const char* namepath, unsigned int *Nvals, float* vals);
void GetLorentzDeflections(float *lorentz_x, float *lorentz_z, 
			   float **lorentz_nx, float **lorentz_nz, 
			   const unsigned int Nxpoints, 
			   const unsigned int Nzpoints);
int GetConstants(const char* namepath, int *Nvals, float* vals, mystr_t* strings);
int GetArrayConstants(const char* namepath, int *Nvals, float* vals, mystr_t* strings);
int GetColumn(const char* namepath, int *Nvals, float* vals, char *key_cstr);
