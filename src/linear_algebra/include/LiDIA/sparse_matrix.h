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


#ifndef LIDIA_SPARSE_MATRIX_H_GUARD_
#define LIDIA_SPARSE_MATRIX_H_GUARD_


#ifndef LIDIA_SPARSE_BASE_MATRIX_H_GUARD_
# include	"LiDIA/sparse_base_matrix.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T >
class sparse_matrix : public sparse_base_matrix< T >
{

	//
	// constructors
	//

public:

	sparse_matrix() : sparse_base_matrix< T > () {}
	sparse_matrix(lidia_size_t a, lidia_size_t b) : sparse_base_matrix< T > (a, b) {}
	sparse_matrix(lidia_size_t a, lidia_size_t b, const T **v) : sparse_base_matrix< T > (a, b, v) {}
	sparse_matrix(const matrix_representation< T > &A) : sparse_base_matrix< T > (A) {}
	sparse_matrix(const sparse_base_matrix< T > &A) : sparse_base_matrix< T > (A) {}

	//
	// destructor
	//

public:

	~sparse_matrix() {}
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_SPARSE_MATRIX_H_GUARD_
