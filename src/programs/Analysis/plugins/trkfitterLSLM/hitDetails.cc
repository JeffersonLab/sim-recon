#include "hitDetails.h"

int FDCHitDetails::instanceCount = 0;
int CDCHitDetails::instanceCount = 0;

FDCHitDetails::FDCHitDetails(): rCorr(3) {
  instanceCount++;
  return;
}

FDCHitDetails::~FDCHitDetails() {
  instanceCount--;
  return;
}

int FDCHitDetails::getInstanceCount() {
  return instanceCount;
}

CDCHitDetails::CDCHitDetails() {
  instanceCount++;
  return;
}

CDCHitDetails::~CDCHitDetails() {
  instanceCount--;
  return;
}

int CDCHitDetails::getInstanceCount() {
  return instanceCount;
}
