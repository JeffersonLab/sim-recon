#ifndef DResourcePool_h
#define DResourcePool_h

#include <atomic>
#include <typeinfo>
#include <vector>
#include <mutex>
#include <iostream>
#include <type_traits>
#include <memory>

#include "DResettable.h"

using namespace std;

/****************************************************** OVERVIEW ******************************************************
 *
 * This class can be used to pool resources between pools for a given object type, even on different threads.
 * The way this works is that, in addition to a "local" resource pool for each instantiated DResourcePool object,
 * there is a private static resource pool, which is thus shared amongst all pools of that type.
 *
 * So, if an object is requested and the local pool is empty, it will then try to retrieve one from the shared pool.
 * If none exist, it makes a new one.
 * Also, if you a recycle an object and the local pool is "full," it will save it to the shared pool instead.
 * The reason there is both a "local" pool and a shared pool is to minimize the locking needed.
 * Control variables can be set to define exactly how this behaves.
 *
 * A pool counter keeps track of how many DResourcePool's there are for a given type.
 * Once it drops to zero, all of the remaining objects are deleted.
 *
 ***************************************** DEFINITION AND USE: NON-SHARED_PTR'S ***************************************
 *
 * If you do not intend to retrieve the resources from this pool as shared_ptr's, then just define a pool (in e.g. a factory) as:
 * DResourcePool<MyClass> dMyClassPool;
 *
 * You can then retrieve and recycle resources via the Get_Resource() and Recycle() functions.
 * Be sure to recycle the memory once you are done with it (e.g. beginning of factory evnt() call) or else it will leak.
 *
 ******************************************* DEFINITION AND USE: SHARED_PTR'S *****************************************
 *
 * You can retrieve shared_ptr's of the objects by calling the Get_SharedResource() method.
 * The advantage of using shared_ptr's is that they automatically keep track of when they are out of scope.
 * These shared_ptr's have been created with a DSharedPtrRecycler functor:
 * Once the shared_ptr goes out of scope, the contained resource is automatically recycled back to the pool.
 *
 * Note that there is a tricky situation that arises when a shared_ptr outlives the life of the DResourcePool.
 * This can happen if (e.g.) shared_ptr's are stored as members of a factory alongside a DResourcePool, and the pool is destroyed first.
 * Or if they are created by a thread_local pool, and the objects outlive the thread. 
 *
 * To combat this, we have the DSharedPtrRecycler hold a weak_ptr to the DResourcePool (where the object gets recycled to).
 * That way, if the pool has been deleted, the weak_ptr will have expired and we can manually delete the object instead of trying to recycle it.
 * This only works if the pool itself has been created within a shared_ptr in the first place, so it is recommended that these pools be defined as:
 *
 * YOU MUST DO THIS!!
 * auto dMyClassPool = std::make_shared<DResourcePool<MyClass>>();
 *
 *************************************************** FACTORY OBJECTS **************************************************
 *
 * If you want the _data objects for the factory to be managed by a resource pool instead of JANA, in the factory init() call:
 * SetFactoryFlag(NOT_OBJECT_OWNER);
 * With this flag set, JANA will not delete the objects in _data, but it will clear them (_data.clear()) prior to the evnt() method.
 * So that means you must also put them in a separate factory vector so that you have a handle on them to recycle them at the beginning of the factory evnt().
 *
 *************************************************** CLASS COMPONENTS **************************************************
 *
 * Note that the components of the particle classes (DKinematicData, DChargedTrackHypothesis, DNeutralParticleHypothesis) contain shared_ptr's.
 * That's because the components (kinematics, timing info, etc.) are often identical between objects, and instead of duplicating the memory, it's cheaper to just shared
 * For example, the kinematics for the pre-kinfit DChargedTrackHypothesis are identical to those from the DTrackTimeBased.
 * Also, the tracking information for the post-kinfit DChargedTrackHypothesis is identical to the pre-kinfit information.
 *
 * The resource pools for these shared_ptr's are defined to be private thread_local within the classes themselves.
 * That way each thread has an instance of the pool, while still sharing a common pool underneath.
 *
 *********************************************** BEWARE: CLASS COMPONENTS **********************************************
 *
 * Problem: Creating a resource on one thread, and recycling it on another thread: race condition in the resource pool LOCAL pool.
 *
 * How does this happen?  You create a charged-hypo on one thread, which creates a DKinematicInfo object on that thread.
 * You then recycle the hypo, and another thread picks it up.  The other thread calls Hypo::Reset(), which clears the DKinematicInfo.
 * Since it's stored in a shared_ptr, it goes back to the thread it was created.
 *
 * Solution 1: Disable intra-thread sharing of objects that CONTAIN shared_ptrs: Everything that is or inherits from DKinematicData
 * Do this by setting the max size of the shared pool to zero for these types.
 *
 * Solution 2: Have those objects inherit from DResettable (below), and define the member functions.
 *
 **********************************************************************************************************************/
 
// Apple compiler does not currently support alignas. Make this an empty definition if it is not already defined.
#ifndef alignas
#define alignas(A) 
#endif

template <typename DType> class DResourcePool : public std::enable_shared_from_this<DResourcePool<DType>>
{
	//TYPE TRAIT REQUIREMENTS
	//If these statements are false, this won't compile
	static_assert(!std::is_pointer<DType>::value, "The template type for DResourcePool must not be a pointer (the stored type IS a pointer though).");
	static_assert(!std::is_const<DType>::value, "The template type for DResourcePool must not be const.");
	static_assert(!std::is_volatile<DType>::value, "The template type for DResourcePool must not be volatile.");

	public:
		DResourcePool(void);
		DResourcePool(size_t locGetBatchSize, size_t locNumToAllocateAtOnce, size_t locMaxLocalPoolSize);
		DResourcePool(size_t locGetBatchSize, size_t locNumToAllocateAtOnce, size_t locMaxLocalPoolSize, size_t locMaxSharedPoolSize, size_t locDebugLevel);
		~DResourcePool(void);

		void Set_ControlParams(size_t locGetBatchSize, size_t locNumToAllocateAtOnce, size_t locMaxLocalPoolSize);
		void Set_ControlParams(size_t locGetBatchSize, size_t locNumToAllocateAtOnce, size_t locMaxLocalPoolSize, size_t locMaxSharedPoolSize, size_t locDebugLevel);
		DType* Get_Resource(void);
		shared_ptr<DType> Get_SharedResource(void);

		//RECYCLE CONST OBJECTS //these methods just const_cast and call the non-const versions
		void Recycle(const DType* locResource){Recycle(const_cast<DType*>(locResource));}
		void Recycle(vector<const DType*>& locResources); //move-clears the input vector

		//RECYCLE NON-CONST OBJECTS
		void Recycle(DType* locResource);
		void Recycle(vector<DType*>& locResources); //move-clears the input vector

		size_t Get_SharedPoolSize(void) const;
		size_t Get_PoolSize(void) const{return dResourcePool_Local.size();}
		size_t Get_NumObjectsAllThreads(void) const{return dObjectCounter;}

		static constexpr unsigned int Get_CacheLineSize(void)
		{
			/// Returns the cache line size for the processor of the target platform.
			/*! The cache line size is useful for creating a buffer to make sure that a variable accessed by multiple threads does not share the cache line it is on.
			    This is useful for variables that may be written-to by one of the threads, because the thread will acquire locked access to the entire cache line.
			    This blocks other threads from operating on the other data stored on the cache line. Note that it is also important to align the shared data as well.
			    See http://www.drdobbs.com/parallel/eliminate-false-sharing/217500206?pgno=4 for more details. */

			//cache line size is 64 for ifarm1402, gcc won't allow larger than 128
			//the cache line size is in /sys/devices/system/cpu/cpu0/cache/index0/coherency_line_size
			return 64; //units are in bytes
		}

	private:

		//Enable this version if type inherits from DResettable //void: is return type
		template <typename RType> typename std::enable_if<std::is_base_of<DResettable, RType>::value, void>::type Release_Resources(DResettable* locResource){locResource->Release();}
		template <typename RType> typename std::enable_if<std::is_base_of<DResettable, RType>::value, void>::type Reset(DResettable* locResource){locResource->Reset();}

		//Enable this version if type does NOT inherit from DResettable //void: is return type
		template <typename RType> typename std::enable_if<!std::is_base_of<DResettable, RType>::value, void>::type Release_Resources(RType* locResource){};
		template <typename RType> typename std::enable_if<!std::is_base_of<DResettable, RType>::value, void>::type Reset(RType* locResource){};

		//Assume that access to the shared pool won't happen very often: will mostly access the thread-local pool (this object)
		void Get_Resources_StaticPool(void);
		void Recycle_Resources_StaticPool(vector<DType*>& locResources);

		alignas(Get_CacheLineSize()) size_t dDebugLevel = 0;
		alignas(Get_CacheLineSize()) size_t dGetBatchSize = 100;
		alignas(Get_CacheLineSize()) size_t dNumToAllocateAtOnce = 20;
		alignas(Get_CacheLineSize()) size_t dMaxLocalPoolSize = 2000;
		alignas(Get_CacheLineSize()) vector<DType*> dResourcePool_Local;

		//static class members have external linkage: same instance shared between every translation unit (would be globally, put only private access)
		alignas(Get_CacheLineSize()) static mutex dSharedPoolMutex;
		alignas(Get_CacheLineSize()) static vector<DType*> dResourcePool_Shared;
		alignas(Get_CacheLineSize()) static size_t dMaxSharedPoolSize;
		alignas(Get_CacheLineSize()) static size_t dPoolCounter; //must be accessed within a lock due to how it's used in destructor: freeing all resources
		alignas(Get_CacheLineSize()) static atomic<size_t> dObjectCounter; //can be accessed without a lock!
};

/********************************************************************************* DSharedPtrRecycler *********************************************************************************/

template <typename DType> class DSharedPtrRecycler
{
	public:
		DSharedPtrRecycler(void) = delete;
		DSharedPtrRecycler(const std::shared_ptr<DResourcePool<DType>>& locResourcePool) : dResourcePool(locResourcePool) {};
		void operator()(const DType* locResource) const{(*this)(const_cast<DType*>(locResource));}
		void operator()(DType* locResource) const;

	private:
		std::weak_ptr<DResourcePool<DType>> dResourcePool;
};

template <typename DType> void DSharedPtrRecycler<DType>::operator()(DType* locResource) const
{
	auto locSharedPtr = dResourcePool.lock();
	if(locSharedPtr == nullptr)
		delete locResource;
	else
		locSharedPtr->Recycle(locResource);
}

/************************************************************************* STATIC MEMBER DEFINITIONS, STRUCTORS *************************************************************************/

//STATIC MEMBER DEFINITIONS
//Since these are part of a template, these statics will only be defined once, no matter how much this header is included
template <typename DType> mutex DResourcePool<DType>::dSharedPoolMutex;
template <typename DType> vector<DType*> DResourcePool<DType>::dResourcePool_Shared = {};
template <typename DType> size_t DResourcePool<DType>::dMaxSharedPoolSize{1000};
template <typename DType> size_t DResourcePool<DType>::dPoolCounter{0};
template <typename DType> atomic<size_t> DResourcePool<DType>::dObjectCounter{0};

//CONSTRUCTORS
template <typename DType> DResourcePool<DType>::DResourcePool(size_t locGetBatchSize, size_t locNumToAllocateAtOnce, size_t locMaxLocalPoolSize, size_t locMaxSharedPoolSize, size_t locDebugLevel) : DResourcePool()
{
	Set_ControlParams(locGetBatchSize, locNumToAllocateAtOnce, locMaxLocalPoolSize, locMaxSharedPoolSize, locDebugLevel);
}

template <typename DType> DResourcePool<DType>::DResourcePool(size_t locGetBatchSize, size_t locNumToAllocateAtOnce, size_t locMaxLocalPoolSize) : DResourcePool()
{
	Set_ControlParams(locGetBatchSize, locNumToAllocateAtOnce, locMaxLocalPoolSize);
}

template <typename DType> DResourcePool<DType>::DResourcePool(void)
{
	dResourcePool_Local.reserve(dMaxLocalPoolSize);
	{
		std::lock_guard<std::mutex> locLock(dSharedPoolMutex); //LOCK
		++dPoolCounter;
		if(dDebugLevel > 0)
			cout << "CONSTRUCTOR THREAD COUNTER " << typeid(DType).name() << ": " << dPoolCounter << endl;
		if(dPoolCounter == 1)
			dResourcePool_Shared.reserve(dMaxSharedPoolSize);
	} //UNLOCK
}

//DESTRUCTOR
template <typename DType> DResourcePool<DType>::~DResourcePool(void)
{
	//Move all objects into the shared pool
	Recycle_Resources_StaticPool(dResourcePool_Local);

	//if this was the last thread, delete all of the remaining resources
	//first move them outside of the vector, then release the lock
	vector<DType*> locResources;
	{
		std::lock_guard<std::mutex> locLock(dSharedPoolMutex); //LOCK
		--dPoolCounter;
		if(dDebugLevel > 0)
			cout << "DESTRUCTOR THREAD COUNTER " << typeid(DType).name() << ": " << dPoolCounter << endl;
		if(dPoolCounter > 0)
			return; //not the last thread

		//last thread: move all resources out of the shared pool
		if(dDebugLevel > 0)
			cout << "DESTRUCTOR MOVING FROM SHARED POOL " << typeid(DType).name() << ": " << std::distance(dResourcePool_Shared.begin(), dResourcePool_Shared.end()) << endl;
		std::move(dResourcePool_Shared.begin(), dResourcePool_Shared.end(), std::back_inserter(locResources));
		dResourcePool_Shared.clear();
	} //UNLOCK

	//delete the resources
	if(dDebugLevel > 0)
		cout << "DESTRUCTOR DELETING " << typeid(DType).name() << ": " << locResources.size() << endl;
	for(auto locResource : locResources)
		delete locResource;
	dObjectCounter -= locResources.size(); //I sure hope this is zero!
	if(dDebugLevel > 0)
		cout << "All objects (ought) to be destroyed, theoretical # remaining: " << dObjectCounter;
}

/************************************************************************* NON-SHARED-POOL-ACCESSING MEMBER FUNCTIONS *************************************************************************/

template <typename DType> DType* DResourcePool<DType>::Get_Resource(void)
{
	if(dDebugLevel >= 10)
		cout << "GET RESOURCE " << typeid(DType).name() << endl;
	if(dResourcePool_Local.empty())
		Get_Resources_StaticPool();
	if(dResourcePool_Local.empty())
	{
		//perhaps instead use custom allocator
		auto locPotentialNewSize = dResourcePool_Local.size() + dNumToAllocateAtOnce - 1;
		auto locNumToAllocate = (locPotentialNewSize <= dMaxLocalPoolSize) ? dNumToAllocateAtOnce : (dMaxLocalPoolSize - dResourcePool_Local.size() + 1);
		for(size_t loc_i = 0; loc_i < locNumToAllocate - 1; ++loc_i)
			dResourcePool_Local.push_back(new DType);
		dObjectCounter += locNumToAllocate;
		return new DType();
	}

	auto locResource = dResourcePool_Local.back();
	dResourcePool_Local.pop_back();
	Reset<DType>(locResource);
	return locResource;
}

template <typename DType> void DResourcePool<DType>::Set_ControlParams(size_t locGetBatchSize, size_t locNumToAllocateAtOnce, size_t locMaxLocalPoolSize, size_t locMaxSharedPoolSize, size_t locDebugLevel)
{
	dDebugLevel = locDebugLevel;
	dGetBatchSize = locGetBatchSize;
	dNumToAllocateAtOnce = (locNumToAllocateAtOnce > 0) ? locNumToAllocateAtOnce : 1;
	dMaxLocalPoolSize = locMaxLocalPoolSize;
	{
		std::lock_guard<std::mutex> locLock(dSharedPoolMutex); //LOCK
		dMaxSharedPoolSize = locMaxSharedPoolSize;
		dResourcePool_Shared.reserve(dMaxSharedPoolSize);
	} //UNLOCK
}

template <typename DType> void DResourcePool<DType>::Set_ControlParams(size_t locGetBatchSize, size_t locNumToAllocateAtOnce, size_t locMaxLocalPoolSize)
{
	dGetBatchSize = locGetBatchSize;
	dNumToAllocateAtOnce = (locNumToAllocateAtOnce > 0) ? locNumToAllocateAtOnce : 1;
	dMaxLocalPoolSize = locMaxLocalPoolSize;
}

template <typename DType> shared_ptr<DType> DResourcePool<DType>::Get_SharedResource(void)
{
	return shared_ptr<DType>(Get_Resource(), DSharedPtrRecycler<DType>(this->shared_from_this()));
}

template <typename DType> void DResourcePool<DType>::Recycle(vector<const DType*>& locResources)
{
	vector<DType*> locNonConstResources;
	locNonConstResources.reserve(locResources.size());

	auto Deconstifier = [](const DType* locConstPointer) -> DType* {return const_cast<DType*>(locConstPointer);};
	std::transform(locResources.begin(), locResources.end(), std::back_inserter(locNonConstResources), Deconstifier);
	locResources.clear();

	Recycle(locNonConstResources);
}

template <typename DType> void DResourcePool<DType>::Recycle(vector<DType*>& locResources)
{
	for(auto& locResource : locResources)
		Release_Resources<DType>(locResource);

	size_t locFirstToMoveIndex = 0;
	auto locPotentialNewPoolSize = dResourcePool_Local.size() + locResources.size();
	if(locPotentialNewPoolSize > dMaxLocalPoolSize) //we won't move all of the resources into the local pool, as it would be too large: only move a subset
		locFirstToMoveIndex = locResources.size() - (dMaxLocalPoolSize - dResourcePool_Local.size());

	std::move(locResources.begin() + locFirstToMoveIndex, locResources.end(), std::back_inserter(dResourcePool_Local));
	locResources.resize(locFirstToMoveIndex);
	if(!locResources.empty())
		Recycle_Resources_StaticPool(locResources);
}

template <typename DType> void DResourcePool<DType>::Recycle(DType* locResource)
{
	if(locResource == nullptr)
		return;
	vector<DType*> locResourceVector{locResource};
	Recycle(locResourceVector);
}

/************************************************************************* SHARED-POOL-ACCESSING MEMBER FUNCTIONS *************************************************************************/

template <typename DType> void DResourcePool<DType>::Get_Resources_StaticPool(void)
{
	auto locPotentialNewLocalSize = dResourcePool_Local.size() + dGetBatchSize;
	auto locGetBatchSize = (locPotentialNewLocalSize <= dMaxLocalPoolSize) ? dGetBatchSize : dMaxLocalPoolSize - dResourcePool_Local.size();
	{
		std::lock_guard<std::mutex> locLock(dSharedPoolMutex); //LOCK
		if(dResourcePool_Shared.empty())
			return;
		auto locFirstToMoveIndex = (locGetBatchSize >= dResourcePool_Shared.size()) ? 0 : dResourcePool_Shared.size() - locGetBatchSize;
		if(dDebugLevel > 0)
			cout << "MOVING FROM SHARED POOL " << typeid(DType).name() << ": " << dResourcePool_Shared.size() - locFirstToMoveIndex << endl;
		std::move(dResourcePool_Shared.begin() + locFirstToMoveIndex, dResourcePool_Shared.end(), std::back_inserter(dResourcePool_Local));
		dResourcePool_Shared.resize(locFirstToMoveIndex);
	} //UNLOCK
}

template <typename DType> void DResourcePool<DType>::Recycle_Resources_StaticPool(vector<DType*>& locResources)
{
	size_t locFirstToMoveIndex = 0;
	{
		std::lock_guard<std::mutex> locLock(dSharedPoolMutex); //LOCK

		auto locPotentialNewPoolSize = dResourcePool_Shared.size() + locResources.size();
		if(locPotentialNewPoolSize > dMaxSharedPoolSize) //we won't move all of the resources into the shared pool, as it would be too large: only move a subset
			locFirstToMoveIndex = locResources.size() - (dMaxSharedPoolSize - dResourcePool_Shared.size());

		if(dDebugLevel > 0)
			cout << "MOVING TO SHARED POOL " << typeid(DType).name() << ": " << dResourcePool_Local.size() - locFirstToMoveIndex << endl;

		std::move(locResources.begin() + locFirstToMoveIndex, locResources.end(), std::back_inserter(dResourcePool_Shared));
	} //UNLOCK

	if(dDebugLevel > 0)
		cout << "DELETING " << typeid(DType).name() << ": " << locFirstToMoveIndex << endl;

	//any resources that were not moved into the shared pool are deleted instead (too many)
	auto Deleter = [](DType* locResource) -> void {delete locResource;};
	std::for_each(locResources.begin(), locResources.begin() + locFirstToMoveIndex, Deleter);

	dObjectCounter -= locFirstToMoveIndex;
	locResources.clear();
}

template <typename DType> size_t DResourcePool<DType>::Get_SharedPoolSize(void) const
{
	std::lock_guard<std::mutex> locLock(dSharedPoolMutex); //LOCK
	return dResourcePool_Shared.size();
} //UNLOCK

#endif // DResourcePool_h
