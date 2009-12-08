

// This file contains C-callable routine signatures for the routines
// in calibDB.cc. It is specifically intended to give access to the 
// calibration DB to C and FORTRAN routines. This is needed since the
// calibration DB API is in C++.
//
// This file should contain no C++ code.

void initcalibdb_(char *bfield_type, char *bfield_map);
void gufld2_(float *r, float *B);
int GetCalib(const char* namepath, unsigned int *Nvals, float* vals);
void GetLorentzDefelections(float *lorentz_x, float *lorentz_z, float **lorentz_nx, float **lorentz_nz
		, const unsigned int Nxpoints, const unsigned int Nzpoints);
