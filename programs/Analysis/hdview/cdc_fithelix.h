


#ifndef _CDC_FITHELIX_H_
#define _CDC_FITHELIX_H_

#include <TVector3.h>
#include "derror.h"


derror_t cdc_fithelix(TVector3 *v, int Npoints);
derror_t cdc_firstguess(TVector3 *points, int Npoints, float &theta, float &phi, float &p, float &q);
void HelixChisq(Int_t& npar, Double_t *x, Double_t &chisq, Double_t *par, Int_t iflag);

#endif //_CDC_FITHELIX_H_
