#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <assert.h>
#include <algorithm>


#ifndef GPDS_H
#define GPDS_H

//#ifdef __cplusplus
//extern"C" {
//#endif

using namespace std;

class GPDs
{
 public:
  GPDs(const char*, int, int, int, double, double, double);
  double GetReH() const;
  double GetImH() const;
  double GetReE() const;
  double GetImE() const;
  double GetReHtild() const;
  double GetImHtild() const;
  double GetDterm() const;
  void Set_q2_t_eta(double, double, double);
  ~GPDs();

 private:
  int n_q2, n_t, n_eta;
  double tM, Q2, eta;
  double q21, t1, eta1;
  double q22, t2, eta2;
  double *arr_t;
  double *arr_eta;
  double *arr_q2;
  double ReH, ImH, ReE, ImE, ImHtild, ReHtild, Dterm;
  char *fname;
  bool ReadFile();
  void DefineValues();
  
  double *ImH_2;
  double *ReH_2;
  double *ImE_2;
  double *ReE_2;
  double *ImHtild_2;
  double *ReHtild_2;
  double *Dterm_2;
  
  bool fstat; // Status of the file '1' if file read succsessfully '0' if not
  
  vector<double> v_q2;
  vector<double> v_t;
  vector<double> v_eta;
};

//#ifdef __cplusplus
//}
//#endif

#endif
