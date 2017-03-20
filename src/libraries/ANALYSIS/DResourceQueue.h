#include <atomic>

using namespace std;

//This is implemented as a singly-linked-list with 2 locks: One for acquiring, one for recycling
//http://www.drdobbs.com/parallel/writing-a-generalized-concurrent-queue/211601363
template <typename DType> class DResourceQueue
{
	public:
		DResourceQueue(void);
		~DResourceQueue(void);

		void Recycle(DType* locObject);
		DType* Acquire(void);
		long long Get_ApproximateQueueSize(void){return dApproximateQueueSize;}

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

		class alignas(Get_CacheLineSize()) DQueueNode
		{
			DQueueNode( DType* locValue ) : dObject(locValue), dNextNode(nullptr) { }
			DType* dObject;
			atomic<DQueueNode*> dNextNode;
		};

		// list nodes
		DQueueNode* alignas(Get_CacheLineSize()) dFirstNode;
		DQueueNode* alignas(Get_CacheLineSize()) dLastNode;

		// approximate list size //not unsigned in case it drops below 0
		atomic<long long> alignas(Get_CacheLineSize()) dApproximateQueueSize;

		// locks
		atomic<bool> alignas(Get_CacheLineSize()) dAcquisitionLock;
		atomic<bool> alignas(Get_CacheLineSize()) dRecycleLock;
};

inline template <typename DType> DResourceQueue<DType>::DResourceQueue(void)
{
	dFirstNode = dLastNode = new DQueueNode(nullptr);
	dRecycleLock = dAcquisitionLock = false;
}

inline template <typename DType> DResourceQueue<DType>::~DResourceQueue(void)
{
    // delete the list
	while(dFirstNode != nullptr)
	{
		DQueueNode* locTempNode = dFirstNode;
		dFirstNode = locTempNode->dNextNode;
		delete locTempNode->dObject; // no-op if nullptr
		delete locTempNode;
	}
}

inline template <typename DType> void DResourceQueue::Recycle(DType* locObject)
{
	if(locObject == nullptr)
		return;

	DQueueNode* locTempNode = new DQueueNode(locObject);

	//acquire lock
	while(dRecycleLock.exchange(true)){}

	//add to the queue
	dLastNode->dNextNode = locTempNode;
	dLastNode = locTempNode;

	//release lock
	dRecycleLock = false;

	//increment atomic counter
	++dApproximateQueueSize;
}

inline template <typename DType> DType* DResourceQueue::Acquire(void)
{
	//acquire lock
	while(dAcquisitionLock.exchange(true)){}

	//get pointers to nodes
	DQueueNode* locFirstNode = dFirstNode;
	DQueueNode* locNextNode = locFirstNode->dNextNode;

	//check if queue is empty
	if(locNextNode == nullptr)
	{
		dAcquisitionLock = false;
		return nullptr;
	}

	//take out the object
	DType* locObject = locNextNode->dObject;
	locNextNode->dObject = nullptr;

	//update first node
	dFirstNode = locNextNode;

	//release lock
	dAcquisitionLock = false;

	//decrement atomic counter
	--dApproximateQueueSize;

	//cleanup
	delete locValue;
	delete locFirstNode;

	//return object
	return locObject;
}
