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


#ifndef LIDIA_LIDIA_VECTOR_H_GUARD_
#define LIDIA_LIDIA_VECTOR_H_GUARD_



#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif
#ifndef LIDIA_SORT_VECTOR_H_GUARD_
# include	"LiDIA/sort_vector.h"
#endif
#ifndef LIDIA_FILE_VECTOR_H_GUARD_
# include	"LiDIA/file_vector.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T >
class lidia_vector: virtual public file_vector< T >,
		    virtual public sort_vector< T >,
		    virtual public math_vector< T >
{
	//
	// constructors
	//

public:

	lidia_vector():
		base_vector< T > () {}

	explicit lidia_vector(const vector_flags &md):
		base_vector< T > (md) {}
	explicit lidia_vector(lidia_size_t all):
		base_vector< T > (all) {}

	lidia_vector(lidia_size_t all, const vector_flags &md):
		base_vector< T > (all, md) {}
	lidia_vector(lidia_size_t all, lidia_size_t len):
		base_vector< T > (all, len) {}
	lidia_vector(lidia_size_t all, lidia_size_t len, const vector_flags &md):
		base_vector< T > (all, len, md) {}
	lidia_vector(const base_vector< T > &v):
		base_vector< T > (v) {}
	lidia_vector(const base_vector< T > &v, const vector_flags &md):
		base_vector< T > (v, md) {}
	lidia_vector(const T *v, lidia_size_t len):
		base_vector< T > (v, len) {}
	lidia_vector(const T *v, lidia_size_t len, const vector_flags &md):
		base_vector< T > (v, len, md) {}

	//
	// destructor
	//

public:

	~lidia_vector() {}


	//
	// assignment
	//

	lidia_vector< T > &
	operator = (const base_vector< T > & v)
	{
		debug_handler_l ("lidia_vector", "operator = ", LDBL_VECTOR +2);

		if (& v != this)
			base_vector< T >::operator = (v);
		return *this;
	}
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_LIDIA_VECTOR_H_GUARD_
