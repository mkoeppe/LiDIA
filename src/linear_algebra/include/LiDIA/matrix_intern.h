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


#ifndef LIDIA_MATRIX_INTERN_H_GUARD_
#define LIDIA_MATRIX_INTERN_H_GUARD_


#ifndef LIDIA_BASE_MATRIX_H_GUARD_
# include	"LiDIA/base_matrix.h"
#endif
#ifndef LIDIA_MATH_MATRIX_H_GUARD_
# include	"LiDIA/math_matrix.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T >
class matrix : public base_matrix< T >
{
	//
	// constructors
	//

public:

	matrix()
		: base_matrix< T > () {}
	matrix(const matrix_flags &flags)
		: base_matrix< T > (flags) {}
	matrix(lidia_size_t a, lidia_size_t b)
		: base_matrix< T > (a, b) {}
	matrix(lidia_size_t a, lidia_size_t b, const matrix_flags &flags)
		: base_matrix< T > (a, b, flags) {}
	matrix(const base_matrix< T > &A)
		: base_matrix< T > (A) {}
	matrix(const base_matrix< T > &A, const matrix_flags &flags)
		: base_matrix< T > (A, flags) {}
	matrix(const dense_base_matrix< T > &A)
		: base_matrix< T > (A) {}
	matrix(const sparse_base_matrix< T > &A)
		: base_matrix< T > (A) {}
	matrix(const base_vector< T > &v)
		: base_matrix< T > (v) {}
	matrix(const base_vector< T > &v, const matrix_flags &flags)
		: base_matrix< T > (v, flags) {}
	matrix(lidia_size_t a, lidia_size_t b, const T **v)
		: base_matrix< T > (a, b, v) {}
	matrix(lidia_size_t a, lidia_size_t b, const T **v, const matrix_flags &flags)
		: base_matrix< T > (a, b, v, flags) {}

	//
	// destructor
	//

public:

	~matrix() {}

	//
	// assignment
	//

public:

	matrix< T > & operator = (const matrix< T > &B)
	{
		assign(B);
		return *this;
	}
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_MATRIX_INTERN_H_GUARD_
