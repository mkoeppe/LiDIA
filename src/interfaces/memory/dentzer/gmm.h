#ifndef LIDIA_MM_H
#define LIDIA_MM_H

#include	<new.h>
#include	<sys/types.h>

/////////////////////////////////////////////////////////////
// Implementation of gmm.h using libI's modified mm
/////////////////////////////////////////////////////////////

#include	"LiDIA/mm.h"

#define USE_DENTZER 1

#define GMM_INITIALIZATION                                      \
class LiDIA_mm                                                  \
{                                                               \
  private:                                                      \
    static int initialized; \
                                                                \
  public:                                                       \
                                                                \
  LiDIA_mm()                                                    \
  {                                                             \
    if (LiDIA_mm::initialized == 0)                             \
      {                                                         \
        mm_initialize(); \
        LiDIA_mm::initialized = 1; \
      }                                                         \
  }                                                             \
                                                                \
  ~LiDIA_mm()                                                   \
  { }                                                           \
                                                                \
}; \
                                                                \
                                                                \
static LiDIA_mm LiDIA_mm_init_dummy;

#define GMM_IMPLEMENTATION                                      \
                                                                \
inline void *gmm::allocate(size_t NMEMB, size_t SIZE)           \
{ return (void *) mm_calloc_without_init(NMEMB, SIZE); }        \
                                                                \
inline void *gmm::allocate(size_t SIZE)                         \
{ return (void *) mm_malloc_without_init(SIZE); }               \
                                                                \
inline void *gmm::allocate_atomic(size_t SIZE)                  \
{ return (void *) mm_malloc_without_init(SIZE); }               \
                                                                \
inline void *gmm::resize(void *PTR, size_t NSIZE, size_t OSIZE) \
{ OSIZE = OSIZE; return (void *) mm_realloc(PTR, NSIZE); }       \
                                                                \
inline void *gmm::allocate_uncollectable(size_t SIZE)           \
{ return (void *) mm_malloc_without_init(SIZE); }               \
                                                                \
inline void gmm::release(void *PTR)                             \
{ mm_free(PTR); }                                               \
                                                                \
inline void gmm::collect()                                      \
{ mm_collect(); }



/////////////////////////////////////////////////////////////
// INTERFACE BEGIN
/////////////////////////////////////////////////////////////

GMM_INITIALIZATION

//
// Check  if  the  compiler  supports overloading of operator
// new[]
//
#if ! defined(HAVE_ARRAY_NEW) \
 && (__BORLANDC__ >= 0x450 || (__GNUC__ >= 2 && __GNUC_MINOR__ >= 6))
#define HAVE_ARRAY_NEW
#endif

//
// memory_manager_mode  is  used  to determine if objects are
// allocated   in  a  collection  allowing  or  a  collection
// disallowing  manner.  This  mode  has  no  meaning  if the
// implementation   memory   manager   provides   no  garbage
// collection.
//
enum memory_manager_mode
{
	GC, AtGC, NoGC
};

//
// LiDIA's memory manager
//
// INTERFACE
//
// 1-allocate ..................... allocate collected memory
// 2-allocate_uncollectable ..... allocate uncollected memory
// 3-allocate_atomic ................. allocate atomic memory
// 4-resize ......................... expand or shrink memory
// 5-release .................................... free memory
// 6-new .......................... get memory for C++ object
// 7-new_with_mode ..... get memory for C++ object using mode
// 8-delete ...................... free memory for C++ object
// 9-collect ....... (probably) initiate a garbage collection
//
class gmm
{
 public:


	//
	// allocate() allocates memory for an array of NMEMB elements
	// of  SIZE bytes each and returns a pointer to the allocated
	// memory.  The memory is set to zero.
	//
	// For  allocate(),  the  value  returned is a pointer to the
	// allocated  memory,  which is suitably aligned for any kind
	// of variable, or NULL if the request fails.
	//
	inline static void *allocate(size_t NMEMB, size_t SIZE);


	//
	// allocate()  allocates  SIZE bytes and returns a pointer to
	// the allocated memory.  The memory ist not cleared.
	//
	// For  allocate(),  the  value  returned is a pointer to the
	// allocated  memory,  which is suitably aligned for any kind
	// of variable, or NULL if the request fails.
	//
	inline static void *allocate(size_t SIZE);


	//
	// allocate_uncollectable()  allocates SIZE bytes and returns
	// a  pointer  to  the  allocated  memory. The memory ist not
	// cleared  and  will  not  be  collected by the implementing
	// memory manager.
	//
	// For  allocate_uncollectable(),  the  value  returned  is a
	// pointer   to  the  allocated  memory,  which  is  suitably
	// aligned  for  any kind of variable, or NULL if the request
	// fails.
	//
	inline static void *allocate_uncollectable(size_t SIZE);

	//
	// allocate_atomic()  allocates  SIZE  bytes  and  returns  a
	// pointer   to   the   allocated  memory.  allocate_atomic()
	// assumes that there are no relevant pointers in the object.
	// The memory ist not cleared.
	//
	// For allocate_atomic(), the value  returned is a pointer to
	// the  allocated memory,  which  is suitably aligned for any
	// kind of variable, or NULL if the request fails.
	//
	inline static void *allocate_atomic(size_t SIZE);

	// resize() changes  the size of the memory block pointed  to
	// by  PTR  to NSIZE bytes. The contents will be unchanged to
	// the  minimum  of the OSIZE an NSIZE sizes; newly allocated
	// memory will be uninitialized.  If PTR is NULL, the call is
	// equivalent  to  allocate(SIZE);  if size is equal to zero,
	// the  call  is  equivalent  to  release(PTR). Unless PTR is
	// NULL,  it  must  have  been returned by an earlier call to
	// allocate() or resize().
	//
	// resize()  returns a pointer to the newly allocated memory,
	// which is suitably aligned for any kind of variable and may
	// be  different from PTR, or NULL if the request fails or if
	// NSIZE was equal to 0.
	//
	inline static void *resize(void *PTR, size_t NSIZE, size_t OSIZE);


	//
	// release() frees the memory  space pointed to by PTR, which
	// must have been returned by a previous  call to allocate().
	// If PTR is NULL, no operation is performed.
	//
	// release() returns no value.
	//
	inline static void release(void *PTR);


	//
	// collect()   initiates   a   garbage   collection   if  the
	// implementing  memory  manager supports garbage collection.
	//
	inline static void collect();

	//
	// operator new overloading
	//
	void *operator new(size_t SIZE) {
			return gmm::allocate(SIZE);
		}


	//
	// operator new overloading with a predefined mode
	//
	inline void *operator new(size_t SIZE, memory_manager_mode MODE) {
			if (MODE == GC)
				return gmm::allocate(SIZE);
			else if (MODE == AtGC)
				return gmm::allocate_atomic(SIZE);
			else
				return gmm::allocate_uncollectable(SIZE);
		}


	//
	// operator delete overloading
	//
	void operator delete(void *PTR) {
			gmm::release(PTR);
		}

#ifdef HAVE_ARRAY_NEW

	inline void *operator new[] (size_t SIZE) {
			return gmm::operator new(SIZE);
		}

	inline void *operator new[] (size_t SIZE, memory_manager_mode MODE) {
			return gmm::operator new(SIZE, MODE);
		}

	inline void operator delete[] (void *PTR) {
			gmm::operator delete(PTR);
		}

#endif

};

GMM_IMPLEMENTATION

/////////////////////////////////////////////////////////////
// INTERFACE END
/////////////////////////////////////////////////////////////

#endif

