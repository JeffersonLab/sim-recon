/*
 * ReadWriteLock.h
 *
 * Define classes to guard the old C-style mutexes that are used in Hall D
 *
 *  Created on: Apr 25, 2017
 *      Author: Hovanes Egiyan
 */

#ifndef READWRITELOCK_H_
#define READWRITELOCK_H_

#include <pthread.h>

#include <functional>
#include <iostream>

class MutexLockBase {
public:
	static bool functionsAreInitialized;

	MutexLockBase() {
		if( ! functionsAreInitialized ) {
			initFunctions();
		}
	}
	virtual ~MutexLockBase() {}

	static int initFunctions() {
		MutexLockBase::plainLockLock =
				[=](pthread_mutex_t& mutex) -> int {return pthread_mutex_lock(&mutex);};
		MutexLockBase::plainLockUnlock =
				[=](pthread_mutex_t& mutex) -> int {return pthread_mutex_unlock(&mutex);};

		MutexLockBase::readLockLock =
				[=](pthread_rwlock_t& mutex) -> int {return pthread_rwlock_rdlock(&mutex);};
		MutexLockBase::readLockUnock =
				[=](pthread_rwlock_t& mutex) -> int {return pthread_rwlock_unlock(&mutex);};

		MutexLockBase::writeLockLock =
				[=](pthread_rwlock_t& mutex) -> int {return pthread_rwlock_rdlock(&mutex);};
		MutexLockBase::writeLockUnlock =
				[=](pthread_rwlock_t& mutex) -> int {return pthread_rwlock_unlock(&mutex);};
		MutexLockBase::functionsAreInitialized = true;
		return 0;
	}

	static std::function<int(pthread_mutex_t&)> plainLockLock  ;
	static std::function<int(pthread_mutex_t&)> plainLockUnlock  ;

	static std::function<int(pthread_rwlock_t&)> readLockLock ;
	static std::function<int(pthread_rwlock_t&)> readLockUnock  ;

	static std::function<int(pthread_rwlock_t&)> writeLockLock  ;
	static std::function<int(pthread_rwlock_t&)> writeLockUnlock  ;
};

// Define the template class
// This class will be the main class used for locking/unlocking
template<typename T>
class MutexLock : public MutexLockBase {
protected:
	T& m_;
	std::function<int(T&)> lockFunc;
	std::function<int(T&)> unlockFunc;
public:

	//! Construct and lock mutex
	MutexLock(T& m, std::function<int(T&)>& lock,
			std::function<int(T&)>& unlock) : MutexLockBase(),
			m_(m),  lockFunc(lock), unlockFunc(unlock){
		lockFunc(m_);
	}
	//! Destruct and unlock mutex
	virtual ~MutexLock() {
		unlockFunc(m_);
	}
	//! Explicitly lock mutex, need to make sure you do not get into racing condition
	virtual int lock() {
		return lockFunc(m_);
	}
	//! Explicitly unlock existing mutex
	virtual int unlock() {
		return unlockFunc(m_);
	}
	virtual T& getLock() {
		return m_;
	}
};


// class for simple mutex
class PlainLock: public MutexLock<pthread_mutex_t> {
public:
	PlainLock(pthread_mutex_t& m) :
			MutexLock<pthread_mutex_t>(m, plainLockLock,
					plainLockUnlock) {
	}
};

// class for read-write mutex for read locking
class ReadLock: public MutexLock<pthread_rwlock_t> {
public:
	ReadLock(pthread_rwlock_t& m) :
			MutexLock<pthread_rwlock_t>(m, readLockLock,
					readLockUnock) {
	}
};

// class for read-write mutex for write locking
class WriteLock: public MutexLock<pthread_rwlock_t> {
public:
	WriteLock(pthread_rwlock_t& m) :
			MutexLock<pthread_rwlock_t>(m, writeLockLock,
					writeLockUnlock) {
	}
};

#endif /* READWRITELOCK_H_ */
