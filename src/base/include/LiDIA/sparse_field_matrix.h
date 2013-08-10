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
//	$Id: sparse_field_matrix.h,v 2.7 2004/06/15 10:19:18 lidiaadm Exp $
//
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_SPARSE_FIELD_MATRIX_H_GUARD_
#define LIDIA_SPARSE_FIELD_MATRIX_H_GUARD_



#ifndef LIDIA_SPARSE_RING_MATRIX_H_GUARD_
# include	"LiDIA/sparse_ring_matrix.h"
#endif
#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif
#ifndef LIDIA_SPARSE_FIELD_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/sparse_field_matrix_kernel.h"
#endif
#ifndef LIDIA_FIELD_MATRIX_ALGORITHMS_H_GUARD_
# include	"LiDIA/matrix/field_matrix_algorithms.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#define SFMKex SFMK< T >



template< class T >
class sparse_field_matrix : public sparse_ring_matrix< T >
{

private:

	const SFMKex S_field_modul;

	const FMA< T, SFMKex, SFMKex, SFMKex > SSS_field_modul;

	//
	// constructors
	//

public:

	sparse_field_matrix()
		: sparse_ring_matrix< T > () {}
	sparse_field_matrix(lidia_size_t a, lidia_size_t b)
		: sparse_ring_matrix< T > (a, b) {}
	sparse_field_matrix(const base_vector< T > &v)
		: sparse_ring_matrix< T > (v) {}
	sparse_field_matrix(const sparse_base_matrix< T > &A)
		: sparse_ring_matrix< T > (A) {}
	sparse_field_matrix(const sparse_field_matrix< T > &A)
		: sparse_ring_matrix< T > (A) {}
	sparse_field_matrix(lidia_size_t a, lidia_size_t b, const T **v)
		: sparse_ring_matrix< T > (a, b, v) {}

	//
	// destructor
	//

public:

	~sparse_field_matrix() {}

	//////////////////////////////////
	// BEGIN: arithmetic procedures //
	//////////////////////////////////

	//
	// assignment
	//

public:

	using sparse_ring_matrix< T >::assign;

	sparse_field_matrix< T > & operator = (const sparse_field_matrix< T > &B)
	{
//Changed by G.A.
		//assign(B);
		sparse_ring_matrix< T >::assign(B);
		return *this;
	}

	//
	// division
	//

	void divide(const sparse_field_matrix< T > &, const T &);

	void compwise_divide(const sparse_field_matrix< T > &, const sparse_field_matrix< T > &);

	/////////////////////////////////
	// END: arithmetic procedures  //
	// BEGIN: arithmetic operators //
	/////////////////////////////////

	///////////////////////////////
	// END: arithmetic operators //
	///////////////////////////////
};



//////////////////////////////////
// BEGIN: arithmetic procedures //
//////////////////////////////////

//
// division
//

template< class T >
inline void
divide(sparse_field_matrix< T > &RES, const sparse_field_matrix< T > &M, const T &a)
{
	RES.divide(M, a);
}



template< class T >
inline void
compwise_divide(sparse_field_matrix< T > &RES, const sparse_field_matrix< T > &M,
		const sparse_field_matrix< T > &N)
{
	RES.compwise_divide(M, N);
}



/////////////////////////////////
// END: arithmetic procedures  //
// BEGIN: arithmetic operators //
/////////////////////////////////

//
// division
//

template< class T >
inline sparse_field_matrix< T >
operator / (const sparse_field_matrix< T > &A, const T &m)
{
	sparse_field_matrix< T > RES(A.rows, A.columns);
	RES.divide(A, m);
	return RES;
}



template< class T >
inline sparse_field_matrix< T > &
operator /= (sparse_field_matrix< T > &A, const T &m)
{
	A.divide(A, m);
	return A;
}



///////////////////////////////
// END: arithmetic operators //
///////////////////////////////



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#undef SFMKex



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/sparse_field_matrix.cc"
#endif



#endif	// LIDIA_SPARSE_FIELD_MATRIX_H_GUARD_
