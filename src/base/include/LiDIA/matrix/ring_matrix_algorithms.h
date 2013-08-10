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


#ifndef LIDIA_RING_MATRIX_ALGORITHMS_H_GUARD_
#define LIDIA_RING_MATRIX_ALGORITHMS_H_GUARD_



#ifndef LIDIA_MATRIX_REPRESENTATION_H_GUARD_
# include	"LiDIA/matrix/matrix_representation.h"
#endif
#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#define ring_matrix_algorithms RMA



template< class T , class ARG1, class ARG2, class ARG3 >
class ring_matrix_algorithms
{
private:

	const ARG1 modul1;
	const ARG2 modul2;
	const ARG3 modul3;

public:

	//
	// constructor
	//

	ring_matrix_algorithms () {}

	//
	// destructor
	//

	~ring_matrix_algorithms () {}

	//////////////////////////////////
	// BEGIN: arithmetic procedures //
	//////////////////////////////////

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

	void multiply_right(const MR< T > &RES, math_vector< T > &res, const math_vector< T > &v) const;
	void multiply_right(const MR< T > &RES, T *&c, const T *v) const;

	void multiply_left(const MR< T > &RES, math_vector< T > &res, const math_vector< T > &v) const;
	void multiply_left(const MR< T > &RES, T *&c, const T *v) const;

	//
	// negation
	//

	void negate(MR< T > &RES, const MR< T > &B) const;

	//////////////////////////////////
	// END: arithmetical procedures //
	//////////////////////////////////

	//
	// comparisons
	//

	bool equal(const MR< T > &RES, const MR< T > &N) const;

	//
	// trace
	//

	void trace(const MR< T > &RES, T &tr) const;

	//
	// size reduction
	//

	void size_reduction(MR< T > &) const;
};



#undef ring_matrix_algorithms



#define modular_ring_matrix_algorithms MRMA



template< class T , class ARG1, class ARG2, class ARG3 >
class modular_ring_matrix_algorithms
{
private:

	const ARG1 modul1;
	const ARG2 modul2;
	const ARG3 modul3;

public:

	//
	// constructor
	//

	modular_ring_matrix_algorithms () {}

	//
	// destructor
	//

	~modular_ring_matrix_algorithms () {}

	//////////////////////////////////
	// BEGIN: arithmetic procedures //
	//////////////////////////////////

	//
	// addition
	//

	void add(MR< T > &RES, const MR< T > &M, const MR< T > &N, const T &mod) const;
	void add(MR< T > &RES, const MR< T > &M, const T &a, const T &mod) const;
	void add(MR< T > &RES, const T &a, const MR< T > &M, const T &mod) const;

	//
	// subtraction
	//

	void subtract(MR< T > &RES, const MR< T > &M, const MR< T > &N, const T &mod) const;
	void subtract(MR< T > &RES, const MR< T > &M, const T &a, const T &mod) const;
	void subtract(MR< T > &RES, const T &a, const MR< T > &M, const T &mod) const;

	//
	// multiplication
	//

	void multiply(MR< T > &RES, const MR< T > &A, const MR< T > &B, const T &mod) const;
	void multiply(MR< T > &RES, const MR< T > &A, const T &k, const T &mod) const;
	void multiply(MR< T > &RES, const T &k, const MR< T > &A, const T &mod) const;

	void compwise_multiply(MR< T > &RES, const MR< T > &A, const MR< T > &B, const T &mod) const;

	void multiply_right(const MR< T > &RES, math_vector< T > &res, const math_vector< T > &v, const T &mod) const;
	void multiply_right(const MR< T > &RES, T *&c, const T *v, const T &mod) const;

	void multiply_left(const MR< T > &RES, math_vector< T > &res, const math_vector< T > &v, const T &mod) const;
	void multiply_left(const MR< T > &RES, T *&c, const T *v, const T &mod) const;

	//
	// negation
	//

	void negate(MR< T > &RES, const MR< T > &B, const T &mod) const;

	//////////////////////////////////
	// END: arithmetical procedures //
	//////////////////////////////////

	//
	// trace
	//

	void trace(const MR< T > &RES, T &tr, const T &mod) const;
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/matrix/ring_matrix_algorithms.cc"
#endif



#undef modular_ring_matrix_algorithms



#endif	// LIDIA_RING_MATRIX_ALGORITHMS_H_GUARD_
