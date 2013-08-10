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


#ifndef LIDIA_DENSE_FIELD_MATRIX_H_GUARD_
#define LIDIA_DENSE_FIELD_MATRIX_H_GUARD_



#ifndef LIDIA_DENSE_RING_MATRIX_H_GUARD_
# include	"LiDIA/dense_ring_matrix.h"
#endif
#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif
#ifndef LIDIA_DENSE_RING_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/dense_ring_matrix_kernel.h"
#endif
#ifndef LIDIA_DENSE_FIELD_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/dense_field_matrix_kernel.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#define DFMKex DFMK< T >

template< class T >
class dense_field_matrix : public dense_ring_matrix< T >
{

private:

	const DFMKex dense_modul;

	//
	// constructors
	//

public:

	dense_field_matrix()
		: dense_ring_matrix< T > () {}
	dense_field_matrix(lidia_size_t a, lidia_size_t b)
		: dense_ring_matrix< T > (a, b) {}
	dense_field_matrix(const base_vector< T > &v)
		: dense_ring_matrix< T > (v) {}
	dense_field_matrix(const dense_base_matrix< T > &A)
		: dense_ring_matrix< T > (A) {}
	dense_field_matrix(const dense_field_matrix< T > &A)
		: dense_ring_matrix< T > (A) {}
	dense_field_matrix(lidia_size_t a, lidia_size_t b, const T **v)
		: dense_ring_matrix< T > (a, b, v) {}

	//
	// destructor
	//

public:

	~dense_field_matrix() {}

	//////////////////////////////////
	// BEGIN: arithmetic procedures //
	//////////////////////////////////

	//
	// assignment
	//

public:

	using dense_ring_matrix< T >::assign;

	dense_field_matrix< T > & operator = (const dense_field_matrix< T > &B)
	{
		this->assign(B);
		return *this;
	}

        void assign(const dense_field_matrix< T > &B)
        {
                dense_ring_matrix< T >::assign(B);
        }

	//
	// division
	//

	void divide(const dense_field_matrix< T > &, const T &);

	void compwise_divide(const dense_field_matrix< T > &, const dense_field_matrix< T > &);

	/////////////////////////////////
	// END: arithmetic procedures  //
	// BEGIN: arithmetic operators //
	/////////////////////////////////

	//
	// division
	//

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
divide(dense_field_matrix< T > &RES, const dense_field_matrix< T > &M, const T &a)
{
	RES.divide(M, a);
}



template< class T >
inline void
compwise_divide(dense_field_matrix< T > &RES, const dense_field_matrix< T > &M,
		const dense_field_matrix< T > &N)
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
inline dense_field_matrix< T >
operator / (const dense_field_matrix< T > &A, const T &m)
{
	dense_field_matrix< T > RES(A.rows, A.columns);
	RES.divide(A, m);
	return RES;
}



template< class T >
inline dense_field_matrix< T > &
operator /= (dense_field_matrix< T > &A, const T &m)
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



#undef DFMKex



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/dense_field_matrix.cc"
#endif


#endif	// LIDIA_DENSE_FIELD_MATRIX_H_GUARD_
