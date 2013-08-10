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
//	$Id: field_matrix.h,v 2.6 2004/06/15 10:19:17 lidiaadm Exp $
//
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_FIELD_MATRIX_H_GUARD_
#define LIDIA_FIELD_MATRIX_H_GUARD_



#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif
#ifndef LIDIA_RING_MATRIX_H_GUARD_
# include	"LiDIA/ring_matrix.h"
#endif
#ifndef LIDIA_DENSE_FIELD_MATRIX_H_GUARD_
# include	"LiDIA/dense_field_matrix.h"
#endif
#ifndef LIDIA_SPARSE_FIELD_MATRIX_H_GUARD_
# include	"LiDIA/sparse_field_matrix.h"
#endif
#ifndef LIDIA_DENSE_FIELD_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/dense_field_matrix_kernel.h"
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



#define DFMKex DFMK< T >
#define SFMKex DFMK< T >

template< class T >
class field_matrix : public ring_matrix< T >
{
	//
	// modul definitions
	//

private:

	const DFMKex D_field_modul;
	const SFMKex S_field_modul;

	const FMA< T, DFMKex, DFMKex, DFMKex > DDD_field_modul;
	const FMA< T, DFMKex, DFMKex, SFMKex > DDS_field_modul;
	const FMA< T, DFMKex, SFMKex, DFMKex > DSD_field_modul;
	const FMA< T, SFMKex, DFMKex, DFMKex > SDD_field_modul;
	const FMA< T, DFMKex, SFMKex, SFMKex > DSS_field_modul;
	const FMA< T, SFMKex, DFMKex, SFMKex > SDS_field_modul;
	const FMA< T, SFMKex, SFMKex, DFMKex > SSD_field_modul;
	const FMA< T, SFMKex, SFMKex, SFMKex > SSS_field_modul;

  //
  // constructors
  //

public:

	field_matrix()
		: ring_matrix< T > () {}
	field_matrix(const matrix_flags &flags)
		: ring_matrix< T > (flags) {}
	field_matrix(lidia_size_t a, lidia_size_t b)
		: ring_matrix< T > (a, b) {}
	field_matrix(lidia_size_t a, lidia_size_t b, const matrix_flags &flags)
		: ring_matrix< T > (a, b, flags) {}
	field_matrix(const base_vector< T > &v)
		: ring_matrix< T > (v) {}
	field_matrix(const base_vector< T > &v, const matrix_flags &flags)
		: ring_matrix< T > (v, flags) {}
	field_matrix(const base_matrix< T > &A)
		: ring_matrix< T > (A) {}
	field_matrix(const base_matrix< T > &A, const matrix_flags &flags)
		: ring_matrix< T > (A, flags) {}
	field_matrix(const field_matrix< T > &A)
		: ring_matrix< T > (A) {}
	field_matrix(const field_matrix< T > &A, const matrix_flags &flags)
		: ring_matrix< T > (A, flags) {}
	field_matrix(const dense_base_matrix< T > &A)
		: ring_matrix< T > (A) {}
	field_matrix(const sparse_base_matrix< T > &A)
		: ring_matrix< T > (A) {}
	field_matrix(lidia_size_t a, lidia_size_t b, const T **v)
		: ring_matrix< T > (a, b, v) {}
	field_matrix(lidia_size_t a, lidia_size_t b, const T **v, const matrix_flags &flags)
		: ring_matrix< T > (a, b, v, flags) {}

	//
	// destructor
	//

public:

	~field_matrix() {}

	//
	// BEGIN: arithmetic procedures
	//

  //
  // assignment
  //

public:

	using ring_matrix< T >::assign;

	field_matrix< T > & operator = (const field_matrix< T > &B)
	{
// Changed by G.A.		
		//assign(B);
		ring_matrix< T >::assign(B);
		return *this;
	}

	//
	// division
	//

#if 0
public:

	friend void divide(field_matrix< T > &,
			   const field_matrix< T > &, const T &);

protected:
#endif

	void divide(const field_matrix< T > &, const T &);

#if 0
public:

	friend void compwise_divide(field_matrix< T > &,
				    const field_matrix< T > &,
				    const field_matrix< T > &);

protected:
#endif

	void compwise_divide(const field_matrix< T > &, const field_matrix< T > &);

	//
	// division
	//

#if 0
public:

	friend field_matrix< T > operator / (const field_matrix< T > &,
					     const T &);

	friend field_matrix< T > & operator /= (field_matrix< T > &,
						const T &);

#endif


	//
	// transpose function for field_matrix
	//

	field_matrix< T > trans() const;
	void trans(const field_matrix< T > & A);

};



//////////////////////////////////
// BEGIN: arithmetic procedures //
//////////////////////////////////

//
// division
//

template< class T >
inline void
divide(field_matrix< T > &RES, const field_matrix< T > &M, const T &a)
{
	RES.divide(M, a);
}



template< class T >
inline void
compwise_divide(field_matrix< T > &RES, const field_matrix< T > &M,
		const field_matrix< T > &N)
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
inline field_matrix< T >
operator / (const field_matrix< T > &A, const T &m)
{
	field_matrix< T > RES(A.rows, A.columns, A.bitfield);
	RES.divide(A, m);
	return RES;
}



template< class T >
inline field_matrix< T > &
operator /= (field_matrix< T > &A, const T &m)
{
	A.divide(A, m);
	return A;
}



//
// transpose function for field matrix
//

template< class T >
inline void
field_matrix< T >::trans(const field_matrix< T > &A)
{
	base_matrix< T >::trans(A);
}



template< class T >
inline field_matrix< T >
field_matrix< T >::trans() const
{
	return base_matrix< T >::trans();
}



template< class T >
inline field_matrix< T >
trans(const field_matrix< T > &A)
{
	return A.trans();
}



///////////////////////////////
// END: arithmetic operators //
///////////////////////////////



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/field_matrix.cc"
#endif



#undef DFMKex
#undef SFMKex



#endif	// LIDIA_FIELD_MATRIX_H_GUARD_
