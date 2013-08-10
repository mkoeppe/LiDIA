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


#ifndef LIDIA_DENSE_RING_MATRIX_KERNEL_H_GUARD_
#define LIDIA_DENSE_RING_MATRIX_KERNEL_H_GUARD_



#ifndef LIDIA_DENSE_BASE_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/dense_base_matrix_kernel.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#define dense_ring_matrix_kernel DRMK



template< class T >
class dense_ring_matrix_kernel : public DBMK< T >
{

public:

	//
	// constructor
	//

	dense_ring_matrix_kernel () {}

	//
	// destructor
	//

	~dense_ring_matrix_kernel () {}

	//
	// BEGIN: arithmetic procedures
	//

	//
	// addition
	//

	void add(MR< T > &RES, const MR< T > &M, const MR< T > &N) const;
	void add(MR< T > &RES, const MR< T > &M, const T &a) const;
	void add(MR< T > &RES, const T &a, const MR< T > &M) const;

	//
	// subtraction
	//

	void subtract(MR< T > &RES, const MR< T > &M, const MR< T > &N) const;
	void subtract(MR< T > &RES, const MR< T > &M, const T &a) const;
	void subtract(MR< T > &RES, const T &a, const MR< T > &M) const;

	//
	// multiplication
	//

	void multiply(MR< T > &RES, const MR< T > &A, const MR< T > &B) const;
	void multiply(MR< T > &RES, const MR< T > &A, const T &k) const;
	void multiply(MR< T > &RES, const T &k, const MR< T > &A) const;

	void compwise_multiply(MR< T > &RES, const MR< T > &A, const MR< T > &B) const;

	void multiply_right(const MR< T > &RES, T *&c, const T *v) const;
	void multiply_left(const MR< T > &RES, T *&c, const T *v) const;

	//
	// negation
	//

	void negate(MR< T > &RES, const MR< T > &B) const;

	//
	// END: arithmetical procedures
	//

	//
	// comparisons
	//

	bool equal(const MR< T > &RES, const MR< T > &N) const;

	//
	// trace
	//

	void trace(const MR< T > &RES, T &tr) const;
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/matrix/dense_ring_matrix_kernel.cc"
#endif



#undef dense_ring_matrix_kernel



#endif	// LIDIA_DENSE_RING_MATRIX_KERNEL_H_GUARD_
