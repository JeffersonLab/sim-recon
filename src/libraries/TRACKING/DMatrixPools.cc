#include "TRACKING/DReferenceTrajectory.h"

//Declare thread_local resource pools
thread_local std::shared_ptr<DResourcePool<TMatrixFSym>> DReferenceTrajectory::dResourcePool_TMatrixFSym = std::make_shared<DResourcePool<TMatrixFSym>>(10, 10, 50);
