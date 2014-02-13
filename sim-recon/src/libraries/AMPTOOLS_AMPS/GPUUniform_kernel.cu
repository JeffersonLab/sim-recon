
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"


__global__ void
GPUUniform_kernel(GPU_AMP_PROTO)
{

  WCUComplex ans = { 1, 0};  
  pcDevAmp[GPU_THIS_EVENT] = ans;

}

void
GPUUniform_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO)
{
  GPUUniform_kernel<<< dimGrid, dimBlock >>>(GPU_AMP_ARGS);
}
