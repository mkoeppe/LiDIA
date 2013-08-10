//
// overload global new and delete using the allocator defined in gmm
//

#include	<LiDIA/gmm.h>

void *operator new(size_t SIZE)
{ return gmm::allocate_uncollectable(SIZE); }

void operator delete(void *PTR)
{ gmm::release(PTR); }

#ifdef HAVE_ARRAY_NEW

void *operator new[](size_t SIZE)
{ return gmm::allocate_uncollectable(SIZE); }

void operator delete[](void *PTR)
{ gmm::release(PTR); }

#endif
