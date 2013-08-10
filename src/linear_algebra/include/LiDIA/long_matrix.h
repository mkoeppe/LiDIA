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


#ifndef LIDIA_LONG_MATRIX_H_GUARD_
#define LIDIA_LONG_MATRIX_H_GUARD_


#ifndef LIDIA_MATRIX_INTERN_H_GUARD_
# include	"LiDIA/matrix_intern.h"
#endif
#ifndef LIDIA_RING_MATRIX_H_GUARD_
# include	"LiDIA/ring_matrix.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template<>
class matrix< long > : public ring_matrix< long >
{

	//
	// constructors
	//

public:

	matrix()
		: ring_matrix< long > () {}
	matrix(const matrix_flags &flags)
		: ring_matrix< long > (flags) {}
	matrix(lidia_size_t a, lidia_size_t b)
		: ring_matrix< long > (a, b) {}
	matrix(lidia_size_t a, lidia_size_t b, const matrix_flags &flags)
		: ring_matrix< long > (a, b, flags) {}
	matrix(const base_matrix< long > &A)
		: ring_matrix< long > (A) {}
	matrix(const base_matrix< long > &A, const matrix_flags &flags)
		: ring_matrix< long > (A, flags) {}
	matrix(const dense_base_matrix< long > &A)
		: ring_matrix< long > (A) {}
	matrix(const sparse_base_matrix< long > &A)
		: ring_matrix< long > (A) {}
	matrix(const base_vector< long > &v)
		: ring_matrix< long > (v) {}
	matrix(const base_vector< long > &v, const matrix_flags &flags)
		: ring_matrix< long > (v, flags) {}
	matrix(lidia_size_t a, lidia_size_t b, const long **v)
		: ring_matrix< long > (a, b, v) {}
	matrix(lidia_size_t a, lidia_size_t b, const long **v, const matrix_flags &flags)
		: ring_matrix< long > (a, b, v, flags) {}

	//
	// destructor
	//

public:

	~matrix() {}

	//
	// assignment
	//

	matrix< long > & operator = (const base_matrix< long > & B)
	{
		base_matrix< long >::assign(B);
		return *this;
	}

	//
	// pseudo-division
	//

public:

	friend void divide(matrix< long > &, const matrix< long > &, const long &);
	friend void compwise_divide(matrix< long > &, const matrix< long > &,
				    const matrix< long > &);

protected:

	void divide(const matrix< long > &, const long &);
	void compwise_divide(const matrix< long > &, const matrix< long > &);

	//
	// pseudo - division
	//

public:

	friend matrix< long > operator / (const matrix< long > &, const long &);
	friend matrix< long > & operator /= (matrix< long > &, const long &);

	//
	// remainder
	//

public:

	friend matrix< long > operator % (const matrix< long > &, const long &);

	friend matrix< long > &operator %= (matrix< long > &, const long &);

	friend void remainder(matrix< long > &, const matrix< long > &, const long &);

};



//////////////////////////////////
// BEGIN: arithmetic procedures //
//////////////////////////////////

//
// division
//

inline void
divide (matrix< long > &RES, const matrix< long > &M, const long &a)
{
	RES.divide(M, a);
}



inline void
compwise_divide (matrix< long > &RES, const matrix< long > &M,
		 const matrix< long > &N)
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

inline matrix< long >
operator / (const matrix< long > &A, const long &m)
{
	matrix< long > RES(A.rows, A.columns);
	RES.divide(A, m);
	return RES;
}



inline matrix< long > &
operator /= (matrix< long > &A, const long &m)
{
	A.divide(A, m);
	return A;
}



//
// remainder
//

inline matrix< long >
operator % (const matrix< long > &A, const long &mod)
{
	matrix< long > RES(A.rows, A.columns);
	remainder(RES, A, mod);
	return RES;
}



inline matrix< long > &
operator %= (matrix< long > &A, const long &mod)
{
	remainder(A, A, mod);
	return A;
}

///////////////////////////////
// END: arithmetic operators //
///////////////////////////////

typedef matrix< long > long_matrix;



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_LONG_MATRIX_H_GUARD_
