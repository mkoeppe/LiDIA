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


#ifndef LIDIA_BASE_MATRIX_ALGORITHMS_H_GUARD_
#define LIDIA_BASE_MATRIX_ALGORITHMS_H_GUARD_



#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_MATRIX_REPRESENTATION_H_GUARD_
# include	"LiDIA/matrix/matrix_representation.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#define base_matrix_algorithms BMA



template< class T, class ARG, class ARG1 >
class base_matrix_algorithms
{

private:

	const ARG modul; // member matrix
	const ARG1 modul1; // argument matrix 1

public:

	//
	// constructor
	//

	base_matrix_algorithms () {}

	//
	// destructor
	//

	~base_matrix_algorithms () {}

	//
	// constructor kernel
	//

	void constructor(MR< T > &) const;
	void constructor(MR< T > &, const base_vector< T > &) const;
	void constructor(MR< T > &, const MR< T > &) const;

	//
	// exchange functions / swap functions
	//

	void swap_columns(MR< T > &A, lidia_size_t i, lidia_size_t j) const;
	void swap_rows(MR< T > &A, lidia_size_t i, lidia_size_t j) const;

	//
	// insert_at
	//

	void insert_at(MR< T > &, lidia_size_t, lidia_size_t,
		       const MR< T > &, lidia_size_t, lidia_size_t,
		       lidia_size_t, lidia_size_t) const;

	//
	// structur functions
	//

	void resize(MR< T > &A, lidia_size_t r, lidia_size_t c) const;
	void kill(MR< T > &A) const;

	//
	// assignments
	//

	void assign(MR< T > &A, const MR< T > &M) const;

	//
	// diagonal function
	//

	void diag(MR< T > &A, const T &a, const T &b) const;

	//
	// transpose function
	//

	void trans(MR< T > &R, const MR< T > &A) const;
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#undef base_matrix_algorithms



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/matrix/base_matrix_algorithms.cc"
#endif



#endif	// LIDIA_BASE_MATRIX_ALGORITHMS_H_GUARD_
