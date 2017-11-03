#ifndef DResettable_h
#define DResettable_h

//primarily for DResourcePool objects
class DResettable
{
	public:
		virtual ~DResettable(void){}
		virtual void Release(void) = 0; //Release all (pointers to) resources, called when recycled to pool
		virtual void Reset(void) = 0; //Re-initialize the object, called when retrieved from pool
};

#endif // DResettable_h
