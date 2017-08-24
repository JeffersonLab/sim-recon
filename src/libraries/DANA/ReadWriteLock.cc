/*
 * ReadWriteLock.cc
 *
 *  Created on: Apr 25, 2017
 *      Author: Hovanes Egiyan
 */

#include "ReadWriteLock.h"

//using namespace MutexLockBase;

// Define global functions for locking and unlocking as lambdas

bool MutexLockBase::functionsAreInitialized __attribute__((used)) = false;

std::function<int(pthread_mutex_t&)> MutexLockBase::plainLockLock __attribute__((used));
std::function<int(pthread_mutex_t&)> MutexLockBase::plainLockUnlock __attribute__((used));
std::function<int(pthread_rwlock_t&)> MutexLockBase::readLockLock __attribute__((used));
std::function<int(pthread_rwlock_t&)> MutexLockBase::readLockUnock __attribute__((used));
std::function<int(pthread_rwlock_t&)> MutexLockBase::writeLockLock __attribute__((used));
std::function<int(pthread_rwlock_t&)> MutexLockBase::writeLockUnlock __attribute__((used));

