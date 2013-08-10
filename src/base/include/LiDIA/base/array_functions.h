// -*- C++ -*-
//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//	$Id$
//
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_ARRAY_FUNCTIONS_H_GUARD_
#define LIDIA_ARRAY_FUNCTIONS_H_GUARD_


#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// copy functions
//

template< class T >
inline void
copy_data(T *v, const T *w, lidia_size_t len)
{
	for (--len; len >= 0; len--)
		v[len] = w[len];
}



template< class T >
inline void
copy_data2(T **V, const T **W, lidia_size_t len1, lidia_size_t len2)
{
	for (--len1; len1 >= 0; len1--)
		copy_data(V[len1], W[len1], len2);
}



//
// memory allocation and copying
//

template< class T >
inline T *
allocate_memory_and_copy_data(const T *&w, lidia_size_t len)
{
	register T *v = new T[len];
	copy_data(v, w, len);
	return v;
}



template< class T >
inline T **
allocate_memory_and_copy_data2(const T **W, lidia_size_t len1, lidia_size_t len2)
{
	register T **V = new T *[len1];
	for (--len1; len1 >= 0; len1--)
		V[len1] = allocate_memory_and_copy_data(W[len1], len2);
	return V;
}



//
// swap function
//

template< class T >
inline void
swap_data(T *v, T *w, lidia_size_t len)
{
	for (--len; len >= 0; len--)
		swap(v[len], w[len]);
}



template< class T >
inline void
swap_data2(T **V, T **W, lidia_size_t len1, lidia_size_t len2)
{
	for (--len1; len1 >= 0; len1--)
		swap_data(V[len1], W[len1], len2);
}



//
// memory allocation and swaping
//

template< class T >
inline T *
allocate_memory_and_swap_data(T *w, lidia_size_t len)
{
	register T *v = new T[len];
	for (--len; len >= 0; len--)
		swap(v[len], w[len]);
	return v;
}



template< class T >
inline T **
allocate_memory_and_swap_data2(T **W, lidia_size_t len1, lidia_size_t len2)
{
	register T **V = new T *[len1];

	for (--len1; len1 >= 0; len1--)
		V[len1] = allocate_memory_and_swap_data(W[len1], len2);
	return V;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_ARRAY_FUNCTIONS_H_GUARD_
