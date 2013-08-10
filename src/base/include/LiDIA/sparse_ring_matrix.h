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
//	$Id: sparse_ring_matrix.h,v 2.6 2004/06/15 10:19:19 lidiaadm Exp $
//
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_SPARSE_RING_MATRIX_H_GUARD_
#define LIDIA_SPARSE_RING_MATRIX_H_GUARD_



#ifndef LIDIA_SPARSE_BASE_MATRIX_H_GUARD_
# include	"LiDIA/sparse_base_matrix.h"
#endif
#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif
#ifndef LIDIA_SPARSE_RING_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/sparse_ring_matrix_kernel.h"
#endif
#ifndef LIDIA_RING_MATRIX_ALGORITHMS_H_GUARD_
# include	"LiDIA/matrix/ring_matrix_algorithms.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#define SRMKex SRMK< T >



template< class T, class R > class solver_kernel;



template< class T >
class sparse_ring_matrix : public sparse_base_matrix< T >
{

	//
	// friend classes
	//

	friend class solver_kernel< T, sparse_ring_matrix < T > >;

private:

	const SRMKex S_ring_modul;

	const RMA< T, SRMKex, SRMKex, SRMKex > SSS_ring_modul;

	//
	// constructors
	//

public:

	sparse_ring_matrix()
		: sparse_base_matrix< T > () {}
	sparse_ring_matrix(lidia_size_t a, lidia_size_t b)
		: sparse_base_matrix< T > (a, b) {}
	sparse_ring_matrix(const base_vector< T > &v)
		: sparse_base_matrix< T > (v) {}
	sparse_ring_matrix(const sparse_base_matrix< T > &A)
		: sparse_base_matrix< T > (A) {}
	sparse_ring_matrix(const sparse_ring_matrix< T > &A)
		: sparse_base_matrix< T > (A) {}
	sparse_ring_matrix(lidia_size_t a, lidia_size_t b, const T **v)
		: sparse_base_matrix< T > (a, b, v) {}

	//
	// destructor
	//

public:

	~sparse_ring_matrix() {}

	//
	// assignment
	//

public:

	using sparse_base_matrix< T >::assign;

	sparse_ring_matrix< T > & operator = (const sparse_ring_matrix< T > &B)
	{
//Changed by G.A.		
		//assign(B);
		sparse_base_matrix< T >::assign(B);
		return *this;
	}

	//////////////////////////////////
	// BEGIN: arithmetic procedures //
	//////////////////////////////////

	//
	// addition
	//

public:

#if 0
	friend void add(sparse_ring_matrix< T > &,
			const sparse_ring_matrix< T > &,
			const sparse_ring_matrix< T > &);

	friend void add(sparse_ring_matrix< T > &,
			const sparse_ring_matrix< T > &, const T &);

	friend void add(sparse_ring_matrix< T > &, const T &,
			const sparse_ring_matrix< T > &);

//protected:
#endif

	void add(const sparse_ring_matrix< T > &, const sparse_ring_matrix< T > &);
	void add(const sparse_ring_matrix< T > &, const T &);
	void add(const T &, const sparse_ring_matrix< T > &);

	//
	// subtraction
	//

#if 0
public:

	friend void subtract(sparse_ring_matrix< T > &,
			     const sparse_ring_matrix< T > &,
			     const sparse_ring_matrix< T > &);

	friend void subtract(sparse_ring_matrix< T > &,
			     const sparse_ring_matrix< T > &, const T &);

	friend void subtract(sparse_ring_matrix< T > &,
			     const T &, const sparse_ring_matrix< T > &);

//protected:
#endif

	void subtract(const sparse_ring_matrix< T > &, const sparse_ring_matrix< T > &);
	void subtract(const sparse_ring_matrix< T > &, const T &);
	void subtract(const T &, const sparse_ring_matrix< T > &);

	//
	// multiplication
	//

#if 0
public:

	friend void multiply(sparse_ring_matrix< T > &,
			     const sparse_ring_matrix< T > &,
			     const sparse_ring_matrix< T > &);

	friend void multiply(sparse_ring_matrix< T > &,
			     const sparse_ring_matrix< T > &, const T &);

	friend void multiply(sparse_ring_matrix< T > &,
			     const T &, const sparse_ring_matrix< T > &);

//protected:
#endif

	void multiply(const sparse_ring_matrix< T > &, const sparse_ring_matrix< T > &);
	void multiply(const sparse_ring_matrix< T > &, const T &);
	void multiply(const T &, const sparse_ring_matrix< T > &);

#if 0
public:

	friend void compwise_multiply(sparse_ring_matrix< T > &,
				      const sparse_ring_matrix< T > &,
				      const sparse_ring_matrix< T > &);

//protected:
#endif

	void compwise_multiply(const sparse_ring_matrix< T > &,
			       const sparse_ring_matrix< T > &);

#if 0
public:

	friend void multiply(T *&, const sparse_ring_matrix< T > &, const T *);

	friend void multiply(T *&, const T *, const sparse_ring_matrix< T > &);

	friend void multiply(math_vector< T > &, const sparse_ring_matrix< T > &,
			     const math_vector< T > &);

	friend void multiply(math_vector< T > &, const math_vector< T > &,
			     const sparse_ring_matrix< T > &);

//protected:
#endif

	void multiply_right(T *&, const T *) const;
	void multiply_right(math_vector< T > &, const math_vector< T > &) const;
	void multiply_left(T *&, const T *) const;
	void multiply_left(math_vector< T > &, const math_vector< T > &) const;

	//
	// negation
	//

public:

#if 0
	friend void negate(sparse_ring_matrix< T > &, const sparse_ring_matrix< T > &);

//protected:
#endif

	void negate(const sparse_ring_matrix< T > &);

	/////////////////////////////////
	// END: arithmetic procedures  //
	// BEGIN: arithmetic operators //
	/////////////////////////////////

#if 0
	//
	// addition
	//

public:

	friend sparse_ring_matrix< T > operator + (const sparse_ring_matrix< T > &,
						   const sparse_ring_matrix< T > &);

	friend sparse_ring_matrix< T > operator + (const sparse_ring_matrix< T > &,
						   const T &);

	friend sparse_ring_matrix< T > operator + (const T &,
						   const sparse_ring_matrix< T > &);

	friend sparse_ring_matrix< T > & operator += (sparse_ring_matrix< T > &,
						      const sparse_ring_matrix< T > &);

	friend sparse_ring_matrix< T > & operator += (sparse_ring_matrix< T > &,
						      const T &);

	//
	// subtraction
	//

public:

	friend sparse_ring_matrix< T > operator - (const sparse_ring_matrix< T > &,
						   const sparse_ring_matrix< T > &);

	friend sparse_ring_matrix< T > operator - (const sparse_ring_matrix< T > &,
						   const T &);

	friend sparse_ring_matrix< T > operator - (const T &,
						   const sparse_ring_matrix< T > &);

	friend sparse_ring_matrix< T > & operator -= (sparse_ring_matrix< T > &,
						      const sparse_ring_matrix< T > &);

	friend sparse_ring_matrix< T > & operator -= (sparse_ring_matrix< T > &,
						      const T &);

	//
	// multiplication
	//

public:

	friend sparse_ring_matrix< T > operator * (const sparse_ring_matrix< T > &,
						   const sparse_ring_matrix< T > &);

	friend sparse_ring_matrix< T > operator * (const sparse_ring_matrix< T > &,
						   const T &);

	friend sparse_ring_matrix< T > operator * (const T &,
						   const sparse_ring_matrix< T > &);

	friend T * operator * (const sparse_ring_matrix< T > &, const T *);

	friend T * operator * (const T *, const sparse_ring_matrix< T > &);

public:

	friend math_vector< T > operator * (const sparse_ring_matrix< T > &,
					    const math_vector< T > &);

	friend math_vector< T > operator * (const math_vector< T > &,
					    const sparse_ring_matrix< T > &);

	friend sparse_ring_matrix< T > & operator *= (sparse_ring_matrix< T > &,
						      const sparse_ring_matrix< T > &);

	friend sparse_ring_matrix< T > & operator *= (sparse_ring_matrix< T > &,
						      const T &);

	//
	// negation
	//


public:

	friend sparse_ring_matrix< T > operator - (const sparse_ring_matrix< T > &);

#endif

	///////////////////////////////
	// END: arithmetic operators //
	///////////////////////////////

	//
	// comparisons
	//

	bool equal(const sparse_ring_matrix< T > &) const;
	bool unequal(const sparse_ring_matrix< T > &A) const
	{
		return (!equal(A));
	}

	//
	// trace
	//

	void trace(T &) const;

	T trace() const
	{
		T tr;
		trace(tr);
		return tr;
	}

	//
	// size reduction
	//

	void size_reduction();
};



//////////////////////////////////
// BEGIN: arithmetic procedures //
//////////////////////////////////

//
// addition
//

template< class T >
inline void
add(sparse_ring_matrix< T > &RES, const sparse_ring_matrix< T > &M,
    const sparse_ring_matrix< T > &N)
{
	RES.add(M, N);
}



template< class T >
inline void
add(sparse_ring_matrix< T > &RES, const sparse_ring_matrix< T > &M, const T &a)
{
	RES.add(M, a);
}



template< class T >
inline void
add(sparse_ring_matrix< T > &RES, const T &a, const sparse_ring_matrix< T > &M)
{
	RES.add(a, M);
}



//
// subtraction
//

template< class T >
inline void
subtract(sparse_ring_matrix< T > &RES, const sparse_ring_matrix< T > &M,
	 const sparse_ring_matrix< T > &N)
{
	RES.subtract(M, N);
}



template< class T >
inline void
subtract(sparse_ring_matrix< T > &RES, const sparse_ring_matrix< T > &M, const T &a)
{
	RES.subtract(M, a);
}



template< class T >
inline void
subtract(sparse_ring_matrix< T > &RES, const T &a, const sparse_ring_matrix< T > &M)
{
	RES.subtract(a, M);
}



//
// multiplication
//

template< class T >
inline void
multiply(sparse_ring_matrix< T > &RES, const sparse_ring_matrix< T > &M,
	 const sparse_ring_matrix< T > &N)
{
	RES.multiply(M, N);
}



template< class T >
inline void
multiply(sparse_ring_matrix< T > &RES, const sparse_ring_matrix< T > &M, const T &a)
{
	RES.multiply(M, a);
}



template< class T >
inline void
multiply(sparse_ring_matrix< T > &RES, const T &a, const sparse_ring_matrix< T > &M)
{
	RES.multiply(a, M);
}



template< class T >
inline void
compwise_multiply(sparse_ring_matrix< T > &RES, const sparse_ring_matrix< T > &M,
		  const sparse_ring_matrix< T > &N)
{
	RES.compwise_multiply(M, N);
}



template< class T >
inline void
multiply(T *&res, const sparse_ring_matrix< T > &A, const T *v)
{
	A.multiply_right(res, v);
}



template< class T >
inline void
multiply(T *&res, const T *v, const sparse_ring_matrix< T > &A)
{
	A.multiply_left(res, v);
}



template< class T >
inline void
multiply(math_vector< T > &res, const sparse_ring_matrix< T > &A,
	 const math_vector< T > &v)
{
	A.multiply_right(res, v);
}



template< class T >
inline void
multiply(math_vector< T > &res, const math_vector< T > &v,
	 const sparse_ring_matrix< T > &A)
{
	A.multiply_left(res, v);
}



//
// negation
//

template< class T >
inline void
negate(sparse_ring_matrix< T > &A, const sparse_ring_matrix< T > &B)
{
	A.negate(B);
}



/////////////////////////////////
// END: arithmetic procedures  //
// BEGIN: arithmetic operators //
/////////////////////////////////

//
// addition
//

template< class T >
inline sparse_ring_matrix< T >
operator + (const sparse_ring_matrix< T > &A, const sparse_ring_matrix< T > &B)
{
	sparse_ring_matrix< T > RES(A.rows, A.columns);
	RES.add(A, B);
	return RES;
}



template< class T >
inline sparse_ring_matrix< T >
operator + (const sparse_ring_matrix< T > &A, const T &a)
{
	sparse_ring_matrix< T > RES(A.rows, A.columns);
	RES.add(A, a);
	return RES;
}



template< class T >
inline sparse_ring_matrix< T >
operator + (const T &a, const sparse_ring_matrix< T > &A)
{
	sparse_ring_matrix< T > RES(A.rows, A.columns);
	RES.add(a, A);
	return RES;
}



template< class T >
inline sparse_ring_matrix< T > &
operator += (sparse_ring_matrix< T > &A, const sparse_ring_matrix< T > &B)
{
	A.add(A, B);
	return A;
}



template< class T >
inline sparse_ring_matrix< T > &
operator += (sparse_ring_matrix< T > &A, const T &a)
{
	A.add(A, a);
	return A;
}



//
// subtraction
//

template< class T >
inline sparse_ring_matrix< T >
operator - (const sparse_ring_matrix< T > &A, const sparse_ring_matrix< T > &B)
{
	sparse_ring_matrix< T > RES(A.rows, A.columns);
	RES.subtract(A, B);
	return RES;
}



template< class T >
inline sparse_ring_matrix< T >
operator - (const sparse_ring_matrix< T > &A, const T &a)
{
	sparse_ring_matrix< T >
		RES(A.rows, A.columns);
	RES.subtract(A, a);
	return RES;
}



template< class T >
inline sparse_ring_matrix< T >
operator - (const T &a, const sparse_ring_matrix< T > &A)
{
	sparse_ring_matrix< T > RES(A.rows, A.columns);
	RES.subtract(a, A);
	return RES;
}



template< class T >
inline sparse_ring_matrix< T > &
operator -= (sparse_ring_matrix< T > &A, const sparse_ring_matrix< T > &B)
{
	A.subtract(A, B);
	return A;
}



template< class T >
inline sparse_ring_matrix< T > &
operator -= (sparse_ring_matrix< T > &A, const T &a)
{
	A.subtract(A, a);
	return A;
}



//
// multiplication
//

template< class T >
inline sparse_ring_matrix< T >
operator * (const sparse_ring_matrix< T > &A, const sparse_ring_matrix< T > &B)
{
	sparse_ring_matrix< T > RES(A.rows, B.columns);
	RES.multiply(A, B);
	return RES;
}



template< class T >
inline sparse_ring_matrix< T >
operator * (const sparse_ring_matrix< T > &A, const T &m)
{
	sparse_ring_matrix< T > RES(A.rows, A.columns);
	RES.multiply(A, m);
	return RES;
}



template< class T >
inline sparse_ring_matrix< T >
operator * (const T &m, const sparse_ring_matrix< T > &A)
{
	sparse_ring_matrix< T > RES(A.rows, A.columns);
	RES.multiply(m, A);
	return RES;
}



template< class T >
inline T *
operator * (const sparse_ring_matrix< T > &A, const T *v)
{
	T *b = new T[A.rows];
	memory_handler(b, "sparse_ring_matrix", "operator * :: "
		       "Error in memory allocation (b)");
	A.multiply_right(b, v);
	return b;
}



template< class T >
inline T *
operator * (const T *v, const sparse_ring_matrix< T > &A)
{
	T *b = new T[A.columns];
	memory_handler(b, "sparse_ring_matrix", "operator * :: "
		       "Error in memory allocation (b)");
	A.multiply_left(b, v);
	return b;
}



template< class T >
inline math_vector< T >
operator * (const sparse_ring_matrix< T > &A, const math_vector< T > &v)
{
	math_vector< T > b;
	A.multiply_right(b, v);
	return b;
}



template< class T >
inline math_vector< T >
operator * (const math_vector< T > &v, const sparse_ring_matrix< T > &A)
{
	math_vector< T > b;
	A.multiply_left(b, v);
	return b;
}



template< class T >
inline sparse_ring_matrix< T > &
operator *= (sparse_ring_matrix< T > &A, const sparse_ring_matrix< T > &B)
{
	A.multiply(A, B);
	return A;
}



template< class T >
inline sparse_ring_matrix< T > &
operator *= (sparse_ring_matrix< T > &A, const T &m)
{
	A.multiply(A, m);
	return A;
}



//
// negation
//

template< class T >
inline sparse_ring_matrix< T >
operator - (const sparse_ring_matrix< T > &A)
{
	sparse_ring_matrix< T > B(A.rows, A.columns);
	B.negate(A);
	return B;
}



///////////////////////////////
// END: arithmetic operators //
///////////////////////////////

//
// comparisons
//

template< class T >
inline bool
operator == (const sparse_ring_matrix< T > &A, const sparse_ring_matrix< T > &B)
{
	return A.equal(B);
}



template< class T >
inline bool
equal(const sparse_ring_matrix< T > &A, const sparse_ring_matrix< T > &B)
{
	return A.equal(B);
}



template< class T >
inline bool
operator != (const sparse_ring_matrix< T > &A, const sparse_ring_matrix< T > &B)
{
	return !A.equal(B);
}



template< class T >
inline bool
unequal(const sparse_ring_matrix< T > &A, const sparse_ring_matrix< T > &B)
{
	return !A.equal(B);
}



//
// trace
//

template< class T >
inline T trace(const sparse_ring_matrix< T > &A)
{
	T tr;
	A.trace(tr);
	return tr;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#undef SRMKex



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/sparse_ring_matrix.cc"
#endif



#endif	// LIDIA_SPARSE_RING_MATRIX_H_GUARD_
