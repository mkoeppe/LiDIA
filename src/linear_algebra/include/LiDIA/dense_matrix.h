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


#ifndef LIDIA_DENSE_MATRIX_H_GUARD_
#define LIDIA_DENSE_MATRIX_H_GUARD_


#ifndef LIDIA_DENSE_BASE_MATRIX_H_GUARD_
# include	"LiDIA/dense_base_matrix.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T >
class dense_matrix : public dense_base_matrix< T >
{

	//
	// constructors
	//

public:

	dense_matrix() : dense_base_matrix< T > () {}
	dense_matrix(lidia_size_t a, lidia_size_t b) : dense_base_matrix< T > (a, b) {}
	dense_matrix(lidia_size_t a, lidia_size_t b, const T **v) : dense_base_matrix< T > (a, b, v) {}
	dense_matrix(const matrix_representation< T > &A) : dense_base_matrix< T > (A) {}
	dense_matrix(const dense_base_matrix< T > &A) : dense_base_matrix< T > (A) {}

	//
	// destructor
	//

public:

	~dense_matrix() {}
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_DENSE_MATRIX_H_GUARD_
