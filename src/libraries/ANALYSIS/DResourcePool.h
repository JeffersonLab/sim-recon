#ifndef DResourcePool_h
#define DResourcePool_h

#include <atomic>
#include <vector>
#include <type_traits>
#include <memory>

using namespace std;

//typical implementation: DResourcePool<MyClass> dMyClassPool; //stores pointers of type MyClass*
template <typename DType> class DResourcePool
{
	//TYPE TRAIT REQUIREMENTS
	//If these statements are false, this won't compile
	static_assert(!std::is_pointer<DType>::value, "The template type for DResourcePool must not be a pointer (the stored type IS a pointer though).");
	static_assert(!std::is_const<DType>::value, "The template type for DResourcePool must not be const.");
	static_assert(!std::is_volatile<DType>::value, "The template type for DResourcePool must not be volatile.");
	static_assert(!std::is_default_constructible<DType>::value, "The template type for DResourcePool must have a default constructor.");

	public:
		DResourcePool(void);
		~DResourcePool(void);

		DType* Get_Resource(void);
		shared_ptr<DType> Get_SharedResource(void);

		void Recycle(DType* locResource);
		void Recycle(const DType* locResource){Recycle(const_cast<DType*>(locResource));}

		//These methods clear the input vectors
		void Recycle(vector<DType*>& locResources);
		void Recycle(vector<const DType*>& locResources);

		size_t Get_SharedPoolSize(void) const;
		size_t Get_PoolSize(void) const{return dResourcePool_Local.size();}

		static constexpr unsigned int Get_CacheLineSize(void)
		{
			/// Returns the cache line size for the processor of the target platform.
			/*! The cache line size is useful for creating a buffer to make sure that a variable accessed by multiple threads does not share the cache line it is on.
			    This is useful for variables that may be written-to by one of the threads, because the thread will acquire locked access to the entire cache line.
			    This blocks other threads from operating on the other data stored on the cache line. Note that it is also important to align the shared data as well.
			    See http://www.drdobbs.com/parallel/eliminate-false-sharing/217500206?pgno=4 for more details. */

			//cache line size is 64 for ifarm1402, so this is conservatively larger
			//the cache line size is in /sys/devices/system/cpu/cpu0/cache/index0/coherency_line_size
			return 256; //units are in bytes
		}

	private:

		//Assume that access to the shared pool won't happen very often: will mostly access the thread-local pool (this object)
		void Get_Resources_StaticPool(void);
		void Recycle_Resources_StaticPool(void);

		size_t alignas(Get_CacheLineSize()) dGetBatchSize = 100;
		size_t alignas(Get_CacheLineSize()) dNumToAllocateAtOnce = 20;
		size_t alignas(Get_CacheLineSize()) dRecycleBatchSize = 1000;
		size_t alignas(Get_CacheLineSize()) dWhenToRecyclePoolSize = 2000; //what size the pool should be before recycling dRecycleBatchSize objects
		size_t alignas(Get_CacheLineSize()) dMaxSharedPoolSize = 10000;
		vector<DType*> alignas(Get_CacheLineSize()) dResourcePool_Local;

		//static class members have external linkage: same instance shared between every translation unit (would be globally, put only private access)
		static vector<DType*> alignas(Get_CacheLineSize()) dResourcePool_Shared;
		static atomic<bool> alignas(Get_CacheLineSize()) dSharedPoolLock;
};

//CONSIDER SPECIALIZING FOR CONTAINER TYPES
//So that they clear() the container and reduce the capacity: std::swap() with temporary
template <typename DType> class DSharedPtrRecycler
{
	public:
		DSharedPtrRecycler(DResourcePool<DType>* locResourcePool) : dResourcePool(locResourcePool) {};
		void operator()(DType* locResource) const {dResourcePool->Recycle(locResource);}
		void operator()(const DType* locResource) const {dResourcePool->Recycle(locResource);}
	private:
		DResourcePool<DType>* dResourcePool;
};

/************************************************************************* STATIC MEMBER DEFINITIONS, STRUCTORS *************************************************************************/

//STATIC MEMBER DEFINITIONS
//Since these are part of a template, these statics will only be defined once, no matter how much this header is included
template <typename DType> vector<DType*> DResourcePool<DType>::dResourcePool_Shared = {};
template <typename DType> atomic<bool> DResourcePool<DType>::dSharedPoolLock = false;

//CONSTRUCTOR
template <typename DType> DResourcePool<DType>::DResourcePool(void)
{
	dResourcePool_Local.reserve(dGetBatchSize);
}

//DESTRUCTOR
template <typename DType> DResourcePool<DType>::~DResourcePool(void)
{
	//Move all objects into the shared pool
	Recycle_SharedResources(dResourcePool_Local);
}

/************************************************************************* NON-SHARED-POOL-ACCESSING MEMBER FUNCTIONS *************************************************************************/

template <typename DType> DType* DResourcePool<DType>::Get_Resource(void)
{
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

template <typename DType> shared_ptr<DType> DResourcePool<DType>::Get_SharedResource(void)
{
	return shared_ptr<DType>(Get_Resource(), DSharedPtrRecycler<DType>(this));
}

template <typename DType> void DResourcePool<DType>::Recycle(DType* locResource)
{
	dResourcePool_Local.push_back(locResource);
	if(dResourcePool_Local.size() > dWhenToRecyclePoolSize)
		Recycle_Resources_StaticPool();
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

/************************************************************************* SHARED-POOL-ACCESSING MEMBER FUNCTIONS *************************************************************************/

template <typename DType> void DResourcePool<DType>::Get_Resources_StaticPool(void)
{
	dResourcePool_Local.reserve(dResourcePool_Local.size() + dGetBatchSize);
	while(dSharedPoolLock.exchange(true)){} //LOCK
	{
		auto locFirstMoveIterator = (dGetBatchSize >= dResourcePool_Shared.size()) ? dResourcePool_Shared.begin() : dResourcePool_Shared.end() - dGetBatchSize;
		std::move(locFirstMoveIterator, dResourcePool_Shared.end(), std::back_inserter(dResourcePool_Local));
		dResourcePool_Shared.erase(locFirstMoveIterator, dResourcePool_Shared.end());
	}
	dSharedPoolLock = false; //UNLOCK
}

template <typename DType> void DResourcePool<DType>::Recycle_Resources_StaticPool(void)
{
	//we will remove dGetBatchSize resources from the local resource pool (or all if size < batch size)
	auto locRemoveIterator = (dGetBatchSize >= dResourcePool_Local.size()) ? dResourcePool_Local.begin() : dResourcePool_Local.end() - dGetBatchSize;
	auto locMoveIterator = locRemoveIterator; //we will move resources into the shared pool, starting at this spot
	while(dSharedPoolLock.exchange(true)){} //LOCK
	{
		auto locNewPoolSize = dResourcePool_Shared.size() + dRecycleBatchSize;
		if(locNewPoolSize > dMaxSharedPoolSize)
		{
			//we won't move all of the resources into the shared pool, as it would be too large: only move a subset
			locMoveIterator = dResourcePool_Local.end() - (dMaxSharedPoolSize - dResourcePool_Shared.size());
			dResourcePool_Shared.reserve(dMaxSharedPoolSize);
		}
		else
			dResourcePool_Shared.reserve(locNewPoolSize);
		std::move(locMoveIterator, dResourcePool_Local.end(), std::back_inserter(dResourcePool_Shared));
	}
	dSharedPoolLock = false; //UNLOCK

	//any resources that were not moved into the shared pool are deleted instead (too many)
	auto Deleter = [](DType* locResource) -> void {delete locResource;};
	if(locMoveIterator != locRemoveIterator)
		std::for_each(locRemoveIterator, locMoveIterator, Deleter);

	dResourcePool_Local.erase(locRemoveIterator, dResourcePool_Local.end());
}

template <typename DType> size_t DResourcePool<DType>::Get_SharedPoolSize(void) const
{
	size_t locSharedPoolSize = 0;
	while(dSharedPoolLock.exchange(true)){} //LOCK
	{
		locSharedPoolSize = dResourcePool_Shared.size();
	}
	dSharedPoolLock = false; //UNLOCK
	return locSharedPoolSize;
}

#endif // DResourcePool_h
