
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"





__global__ void
GPUpolCoef_kernel(GPU_AMP_PROTO , int polBeam, GDouble polFrac)
{
  int pol=(polBeam==1 ? +1 : -1); // y and x-pol. respectively

  //(1+g) for x-pol, (1-g) for y-pol
  WCUComplex ans = { sqrt((1.0-pol*polFrac)*0.5), 0 };
  
  pcDevAmp[GPU_THIS_EVENT] = ans;

}

void
GPUpolCoef_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                   int polBeam, GDouble polFrac)
{
  GPUpolCoef_kernel<<< dimGrid, dimBlock >>>
    ( GPU_AMP_ARGS, polBeam, polFrac );

}
