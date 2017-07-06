#ifndef DResourcePool_h
#define DResourcePool_h

#include <atomic>
#include <typeinfo>
#include <vector>
#include <mutex>
#include <iostream>
#include <type_traits>
#include <memory>

using namespace std;

//typical implementation: DResourcePool<MyClass> dMyClassPool; //stores pointers of type MyClass*
template <typename DType> class DResourcePool : public std::enable_shared_from_this<DResourcePool<DType>>
{
	//TYPE TRAIT REQUIREMENTS
	//If these statements are false, this won't compile
	static_assert(!std::is_pointer<DType>::value, "The template type for DResourcePool must not be a pointer (the stored type IS a pointer though).");
	static_assert(!std::is_const<DType>::value, "The template type for DResourcePool must not be const.");
	static_assert(!std::is_volatile<DType>::value, "The template type for DResourcePool must not be volatile.");

	public:
		DResourcePool(void);
		~DResourcePool(void);

		void Set_ControlParams(size_t locGetBatchSize = 100, size_t locNumToAllocateAtOnce = 20, size_t locRecycleBatchSize = 1000, size_t locWhenToRecyclePoolSize = 2000, size_t locMaxSharedPoolSize = 10000);
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

		//Assume that access to the shared pool won't happen very often: will mostly access the thread-local pool (this object)
		void Get_Resources_StaticPool(void);
		void Recycle_Resources_StaticPool(void);
		void Recycle_Resources_StaticPool(typename vector<DType*>::iterator locRemoveIterator);

		alignas(Get_CacheLineSize()) size_t dGetBatchSize = 100;
		alignas(Get_CacheLineSize()) size_t dNumToAllocateAtOnce = 20;
		alignas(Get_CacheLineSize()) size_t dRecycleBatchSize = 1000;
		alignas(Get_CacheLineSize()) size_t dWhenToRecyclePoolSize = 2000; //what size the pool should be before recycling dRecycleBatchSize objects
		alignas(Get_CacheLineSize()) size_t dMaxSharedPoolSize = 10000;
		alignas(Get_CacheLineSize()) vector<DType*> dResourcePool_Local;

		//static class members have external linkage: same instance shared between every translation unit (would be globally, put only private access)
		alignas(Get_CacheLineSize()) static mutex dSharedPoolMutex;
		alignas(Get_CacheLineSize()) static vector<DType*> dResourcePool_Shared;
		alignas(Get_CacheLineSize()) static size_t dThreadCounter; //must be accessed within a lock due to how it's used in destructor: freeing all resources

		alignas(Get_CacheLineSize()) size_t dContainerResourceMaxCapacity = 1000;
		alignas(Get_CacheLineSize()) size_t dContainerResourceReduceCapacityTo = 100;
};

template <typename DType> class DSharedPtrRecycler
{
	public:
		DSharedPtrRecycler(void) = delete;
		DSharedPtrRecycler(const std::shared_ptr<DResourcePool<DType>>& locResourcePool) : dResourcePool(locResourcePool) {};
		void operator()(DType* locResource) const
		{
			cout << "SHARED_PTR RECYCLE " << typeid(DType).name() << ": " << locResource << endl;
			auto locSharedPtr = dResourcePool.lock();
			if(locSharedPtr == nullptr)
				delete locResource;
			else
				locSharedPtr->Recycle(locResource);
		}
		void operator()(const DType* locResource) const
		{
			cout << "SHARED_PTR RECYCLE " << typeid(DType).name() << ": " << locResource << endl;
			auto locSharedPtr = dResourcePool.lock();
			if(locSharedPtr == nullptr)
				delete locResource;
			else
				locSharedPtr->Recycle(locResource);
		}
	private:
		std::weak_ptr<DResourcePool<DType>> dResourcePool;
};

/************************************************************************* STATIC MEMBER DEFINITIONS, STRUCTORS *************************************************************************/

//STATIC MEMBER DEFINITIONS
//Since these are part of a template, these statics will only be defined once, no matter how much this header is included
template <typename DType> mutex DResourcePool<DType>::dSharedPoolMutex;
template <typename DType> vector<DType*> DResourcePool<DType>::dResourcePool_Shared = {};
template <typename DType> size_t DResourcePool<DType>::dThreadCounter{0};

//CONSTRUCTOR
template <typename DType> DResourcePool<DType>::DResourcePool(void)
{
	dResourcePool_Local.reserve(dGetBatchSize);
	{
		std::lock_guard<std::mutex> locLock(dSharedPoolMutex); //LOCK
		++dThreadCounter;
		cout << "CONSTRUCTOR THREAD COUNTER " << typeid(DType).name() << ": " << dThreadCounter << endl;
	}
}

//DESTRUCTOR
template <typename DType> DResourcePool<DType>::~DResourcePool(void)
{
	//Move all objects into the shared pool
	Recycle_Resources_StaticPool(dResourcePool_Local.begin());

	//if this was the last thread, delete all of the remaining resources
	//first move them outside of the vector, then release the lock
	vector<DType*> locResources;
	{
		std::lock_guard<std::mutex> locLock(dSharedPoolMutex); //LOCK
		--dThreadCounter;
		cout << "DESTRUCTOR THREAD COUNTER " << typeid(DType).name() << ": " << dThreadCounter << endl;
		if(dThreadCounter > 0)
			return; //not the last thread

		//last thread: move all resources out of the shared pool
		cout << "DESTRUCTOR MOVING FROM SHARED POOL " << typeid(DType).name() << ": " << std::distance(dResourcePool_Shared.begin(), dResourcePool_Shared.end()) << endl;
		std::move(dResourcePool_Shared.begin(), dResourcePool_Shared.end(), std::back_inserter(locResources));
		dResourcePool_Shared.clear();
	}

	//delete the resources
	cout << "DESTRUCTOR DELETING " << typeid(DType).name() << ": " << locResources.size() << endl;
	for(auto locResource : locResources)
		delete locResource;
}

/************************************************************************* NON-SHARED-POOL-ACCESSING MEMBER FUNCTIONS *************************************************************************/

template <typename DType> DType* DResourcePool<DType>::Get_Resource(void)
{
	cout << "GET RESOURCE " << typeid(DType).name() << endl;
	if(dResourcePool_Local.empty())
		Get_Resources_StaticPool();
	if(dResourcePool_Local.empty())
	{
		//perhaps instead use custom allocator
		for(size_t loc_i = 0; loc_i < dNumToAllocateAtOnce - 1; ++loc_i)
			dResourcePool_Local.push_back(new DType);
		return new DType();
	}

	auto locResource = dResourcePool_Local.back();
	dResourcePool_Local.pop_back();
	return locResource;
}

template <typename DType> void DResourcePool<DType>::Set_ControlParams(size_t locGetBatchSize, size_t locNumToAllocateAtOnce, size_t locRecycleBatchSize, size_t locWhenToRecyclePoolSize, size_t locMaxSharedPoolSize)
{
	dGetBatchSize = locGetBatchSize;
	dNumToAllocateAtOnce = locNumToAllocateAtOnce;
	dRecycleBatchSize = locRecycleBatchSize;
	dWhenToRecyclePoolSize = locWhenToRecyclePoolSize;
	dMaxSharedPoolSize = locMaxSharedPoolSize;
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
	dResourcePool_Local.reserve(dResourcePool_Local.size() + locResources.size());
	std::move(locResources.begin(), locResources.end(), std::back_inserter(dResourcePool_Local));
	locResources.clear();
	if(dResourcePool_Local.size() > dWhenToRecyclePoolSize)
		Recycle_Resources_StaticPool();
}

//http://stackoverflow.com/questions/8314827/how-can-i-specialize-a-template-member-function-for-stdvectort
template <typename DType> void DResourcePool<DType>::Recycle(DType* locResource)
{
	cout << "RECYCLE " << typeid(DType).name() << ": " << locResource << endl;
	if(locResource == nullptr)
		return;
	dResourcePool_Local.push_back(locResource);
	if(dResourcePool_Local.size() > dWhenToRecyclePoolSize)
		Recycle_Resources_StaticPool();
	cout << "DONE RECYCLING" << endl;
}

/************************************************************************* SHARED-POOL-ACCESSING MEMBER FUNCTIONS *************************************************************************/

template <typename DType> void DResourcePool<DType>::Get_Resources_StaticPool(void)
{
	dResourcePool_Local.reserve(dResourcePool_Local.size() + dGetBatchSize);
	{
		std::lock_guard<std::mutex> locLock(dSharedPoolMutex); //LOCK
		if(dResourcePool_Shared.empty())
			return;
		auto locFirstMoveIterator = (dGetBatchSize >= dResourcePool_Shared.size()) ? dResourcePool_Shared.begin() : dResourcePool_Shared.end() - dGetBatchSize;
		cout << "MOVING FROM SHARED POOL " << typeid(DType).name() << ": " << std::distance(locFirstMoveIterator, dResourcePool_Shared.end()) << endl;
		std::move(locFirstMoveIterator, dResourcePool_Shared.end(), std::back_inserter(dResourcePool_Local));
		dResourcePool_Shared.erase(locFirstMoveIterator, dResourcePool_Shared.end());
	}
}

template <typename DType> void DResourcePool<DType>::Recycle_Resources_StaticPool(void)
{
	//we will remove dGetBatchSize resources from the local resource pool (or all if size < batch size)
	auto locRemoveIterator = (dGetBatchSize >= dResourcePool_Local.size()) ? dResourcePool_Local.begin() : dResourcePool_Local.end() - dGetBatchSize;
	Recycle_Resources_StaticPool(locRemoveIterator);
}

template <typename DType> void DResourcePool<DType>::Recycle_Resources_StaticPool(typename vector<DType*>::iterator locRemoveIterator)
{
	auto locMoveIterator = locRemoveIterator; //we will move resources into the shared pool, starting at this spot
	{
		std::lock_guard<std::mutex> locLock(dSharedPoolMutex); //LOCK
		auto locNewPoolSize = dResourcePool_Shared.size() + dRecycleBatchSize;
		if(locNewPoolSize > dMaxSharedPoolSize)
		{
			//we won't move all of the resources into the shared pool, as it would be too large: only move a subset
			locMoveIterator = dResourcePool_Local.end() - (dMaxSharedPoolSize - dResourcePool_Shared.size());
			dResourcePool_Shared.reserve(dMaxSharedPoolSize);
		}
		else
			dResourcePool_Shared.reserve(locNewPoolSize);
		cout << "MOVING TO SHARED POOL " << typeid(DType).name() << ": " << std::distance(locMoveIterator, dResourcePool_Local.end()) << endl;
		std::move(locMoveIterator, dResourcePool_Local.end(), std::back_inserter(dResourcePool_Shared));
	}


	cout << "DELETING " << typeid(DType).name() << ": " << std::distance(locRemoveIterator, locMoveIterator) << endl;

	//any resources that were not moved into the shared pool are deleted instead (too many)
	auto Deleter = [](DType* locResource) -> void {delete locResource;};
	if(locMoveIterator != locRemoveIterator)
		std::for_each(locRemoveIterator, locMoveIterator, Deleter);

	dResourcePool_Local.erase(locRemoveIterator, dResourcePool_Local.end());
}

template <typename DType> size_t DResourcePool<DType>::Get_SharedPoolSize(void) const
{
	size_t locSharedPoolSize = 0;
	{
		std::lock_guard<std::mutex> locLock(dSharedPoolMutex);
		locSharedPoolSize = dResourcePool_Shared.size();
	}
	return locSharedPoolSize;
}

#endif // DResourcePool_h
