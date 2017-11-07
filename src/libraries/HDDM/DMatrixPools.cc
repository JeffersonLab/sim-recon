#include "HDDM/DEventSourceHDDM.h"
#include "HDDM/DEventSourceREST.h"

//Declare thread_local resource pools
thread_local std::shared_ptr<DResourcePool<TMatrixFSym>> DEventSourceHDDM::dResourcePool_TMatrixFSym = std::make_shared<DResourcePool<TMatrixFSym>>(10, 10, 50);
thread_local std::shared_ptr<DResourcePool<TMatrixFSym>> DEventSourceREST::dResourcePool_TMatrixFSym = std::make_shared<DResourcePool<TMatrixFSym>>(10, 10, 50);
