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


#ifndef LIDIA_DENSE_RING_MATRIX_H_GUARD_
#define LIDIA_DENSE_RING_MATRIX_H_GUARD_



#ifndef LIDIA_DENSE_BASE_MATRIX_H_GUARD_
# include	"LiDIA/dense_base_matrix.h"
#endif
#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif
#ifndef LIDIA_DENSE_RING_MATRIX_KERNEL_H_GUARD_
# include	"LiDIA/matrix/dense_ring_matrix_kernel.h"
#endif
#ifndef LIDIA_RING_MATRIX_ALGORITHMS_H_GUARD_
# include	"LiDIA/matrix/ring_matrix_algorithms.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#define DRMKex DRMK< T >

template< class T >
class dense_ring_matrix : public dense_base_matrix< T >
{

private:

	const DRMKex D_ring_modul;

	const RMA< T, DRMKex, DRMKex, DRMKex > DDD_ring_modul;

	//
	// constructors
	//

public:

	dense_ring_matrix()
		: dense_base_matrix< T > () {}
	dense_ring_matrix(lidia_size_t a, lidia_size_t b)
		: dense_base_matrix< T > (a, b) {}
	dense_ring_matrix(const base_vector< T > &v)
		: dense_base_matrix< T > (v) {}
	dense_ring_matrix(const dense_base_matrix< T > &A)
		: dense_base_matrix< T > (A) {}
	dense_ring_matrix(const dense_ring_matrix< T > &A)
		: dense_base_matrix< T > (A) {}
	dense_ring_matrix(lidia_size_t a, lidia_size_t b, const T **v)
		: dense_base_matrix< T > (a, b, v) {}

	//
	// destructor
	//

public:

	~dense_ring_matrix() {}

	//
	// assignment
	//

public:

	using dense_base_matrix< T >::assign;

	dense_ring_matrix< T > & operator = (const dense_ring_matrix< T > &B)
	{
		this->assign(B);
		return *this;
	}

	//////////////////////////////////
	// BEGIN: arithmetic procedures //
	//////////////////////////////////

	//
	// addition
	//

#if 0
public:

	friend void add(dense_ring_matrix< T > &,
			const dense_ring_matrix< T > &,
			const dense_ring_matrix< T > &);
	friend void add(dense_ring_matrix< T > &,
			const dense_ring_matrix< T > &, const T &);
	friend void add(dense_ring_matrix< T > &,
			const T &, const dense_ring_matrix< T > &);

protected:
#endif

	void add(const dense_ring_matrix< T > &, const dense_ring_matrix< T > &);
	void add(const dense_ring_matrix< T > &, const T &);
	void add(const T &, const dense_ring_matrix< T > &);

	//
	// subtraction
	//

#if 0
public:

	friend void subtract(dense_ring_matrix< T > &,
			     const dense_ring_matrix< T > &,
			     const dense_ring_matrix< T > &);
	friend void subtract(dense_ring_matrix< T > &,
			     const dense_ring_matrix< T > &, const T &);
	friend void subtract(dense_ring_matrix< T > &,
			     const T &, const dense_ring_matrix< T > &);

protected:
#endif

	void subtract(const dense_ring_matrix< T > &, const dense_ring_matrix< T > &);
	void subtract(const dense_ring_matrix< T > &, const T &);
	void subtract(const T &, const dense_ring_matrix< T > &);

	//
	// multiplication
	//

#if 0
public:

	friend void multiply(dense_ring_matrix< T > &,
			     const dense_ring_matrix< T > &,
			     const dense_ring_matrix< T > &);
	friend void multiply(dense_ring_matrix< T > &,
			     const dense_ring_matrix< T > &, const T &);
	friend void multiply(dense_ring_matrix< T > &,
			     const T &, const dense_ring_matrix< T > &);

protected:
#endif

	void multiply(const dense_ring_matrix< T > &, const dense_ring_matrix< T > &);
	void multiply(const dense_ring_matrix< T > &, const T &);
	void multiply(const T &, const dense_ring_matrix< T > &);

#if 0
public:

	friend void compwise_multiply(dense_ring_matrix< T > &,
				      const dense_ring_matrix< T > &,
				      const dense_ring_matrix< T > &);

protected:
#endif

	void compwise_multiply(const dense_ring_matrix< T > &, const dense_ring_matrix< T > &);

#if 0
public:

	friend void multiply(T *&, const dense_ring_matrix< T > &, const T *);
	friend void multiply(T *&, const T *, const dense_ring_matrix< T > &);
	friend void multiply(math_vector< T > &, const dense_ring_matrix< T > &,
			     const math_vector< T > &);

	friend void multiply(math_vector< T > &, const math_vector< T > &,
			     const dense_ring_matrix< T > &);

protected:
#endif

	void multiply_right(T *&, const T *) const;
	void multiply_right(math_vector< T > &, const math_vector< T > &) const;
	void multiply_left(T *&, const T *) const;
	void multiply_left(math_vector< T > &, const math_vector< T > &) const;

	//
	// negation
	//

#if 0
public:

	friend void negate(dense_ring_matrix< T > &, const dense_ring_matrix< T > &);

protected:
#endif

	void negate(const dense_ring_matrix< T > &);

	/////////////////////////////////
	// END: arithmetic procedures  //
	// BEGIN: arithmetic operators //
	/////////////////////////////////

#if 0
	//
	// addition
	//

public:

	friend dense_ring_matrix< T > operator + (const dense_ring_matrix< T > &,
						  const dense_ring_matrix< T > &);

	friend dense_ring_matrix< T > operator + (const dense_ring_matrix< T > &, const T &);

	friend dense_ring_matrix< T > operator + (const T &, const dense_ring_matrix< T > &);

	friend dense_ring_matrix< T > & operator += (dense_ring_matrix< T > &,
						     const dense_ring_matrix< T > &);

	friend dense_ring_matrix< T > & operator += (dense_ring_matrix< T > &, const T &);

	//
	// subtraction
	//

public:

	friend dense_ring_matrix< T > operator - (const dense_ring_matrix< T > &,
						  const dense_ring_matrix< T > &);

	friend dense_ring_matrix< T > operator - (const dense_ring_matrix< T > &, const T &);

	friend dense_ring_matrix< T > operator - (const T &, const dense_ring_matrix< T > &);


	friend dense_ring_matrix< T > & operator -= (dense_ring_matrix< T > &,
						     const dense_ring_matrix< T > &);

	friend dense_ring_matrix< T > & operator -= (dense_ring_matrix< T > &, const T &);

	//
	// multiplication
	//

public:

	friend dense_ring_matrix< T > operator * (const dense_ring_matrix< T > &,
						  const dense_ring_matrix< T > &);

	friend dense_ring_matrix< T > operator * (const dense_ring_matrix< T > &, const T &);

	friend dense_ring_matrix< T > operator * (const T &, const dense_ring_matrix< T > &);

	friend T * operator * (const dense_ring_matrix< T > &, const T *);

	friend T * operator * (const T *, const dense_ring_matrix< T > &);

public:

	friend math_vector< T > operator * (const dense_ring_matrix< T > &,
					    const math_vector< T > &);

	friend math_vector< T > operator * (const math_vector< T > &,
					    const dense_ring_matrix< T > &);

	friend dense_ring_matrix< T > & operator *= (dense_ring_matrix< T > &,
						     const dense_ring_matrix< T > &);

	friend dense_ring_matrix< T > & operator *= (dense_ring_matrix< T > &, const T &);

	//
	// negation
	//

public:

	friend dense_ring_matrix< T > operator - (const dense_ring_matrix< T > &);

#endif

	///////////////////////////////
	// END: arithmetic operators //
	///////////////////////////////

	//
	// comparisons
	//

	bool equal(const dense_ring_matrix< T > &) const;
	bool unequal(const dense_ring_matrix< T > &A) const
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
add(dense_ring_matrix< T > &RES, const dense_ring_matrix< T > &M,
    const dense_ring_matrix< T > &N)
{
	RES.add(M, N);
}



template< class T >
inline void
add(dense_ring_matrix< T > &RES, const dense_ring_matrix< T > &M, const T &a)
{
	RES.add(M, a);
}



template< class T >
inline void
add(dense_ring_matrix< T > &RES, const T &a, const dense_ring_matrix< T > &M)
{
	RES.add(a, M);
}



//
// subtraction
//

template< class T >
inline void
subtract(dense_ring_matrix< T > &RES, const dense_ring_matrix< T > &M,
	 const dense_ring_matrix< T > &N)
{
	RES.subtract(M, N);
}



template< class T >
inline void
subtract(dense_ring_matrix< T > &RES, const dense_ring_matrix< T > &M, const T &a)
{
	RES.subtract(M, a);
}



template< class T >
inline void
subtract(dense_ring_matrix< T > &RES, const T &a, const dense_ring_matrix< T > &M)
{
	RES.subtract(a, M);
}



//
// multiplication
//

template< class T >
inline void
multiply(dense_ring_matrix< T > &RES, const dense_ring_matrix< T > &M,
	 const dense_ring_matrix< T > &N)
{
	RES.multiply(M, N);
}



template< class T >
inline void
multiply(dense_ring_matrix< T > &RES, const dense_ring_matrix< T > &M, const T &a)
{
	RES.multiply(M, a);
}



template< class T >
inline void
multiply(dense_ring_matrix< T > &RES, const T &a, const dense_ring_matrix< T > &M)
{
	RES.multiply(a, M);
}



template< class T >
inline void
compwise_multiply(dense_ring_matrix< T > &RES, const dense_ring_matrix< T > &M,
		  const dense_ring_matrix< T > &N)
{
	RES.compwise_multiply(M, N);
}



template< class T >
inline void
multiply(T *&res, const dense_ring_matrix< T > &A, const T *v)
{
	A.multiply_right(res, v);
}



template< class T >
inline void
multiply(T *&res, const T *v, const dense_ring_matrix< T > &A)
{
	A.multiply_left(res, v);
}



template< class T >
inline void
multiply(math_vector< T > &res, const dense_ring_matrix< T > &A,
	 const math_vector< T > &v)
{
	A.multiply_right(res, v);
}



template< class T >
inline void
multiply(math_vector< T > &res, const math_vector< T > &v,
	 const dense_ring_matrix< T > &A)
{
	A.multiply_left(res, v);
}



//
// negation
//

template< class T >
inline void
negate(dense_ring_matrix< T > &A, const dense_ring_matrix< T > &B)
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
inline dense_ring_matrix< T >
operator + (const dense_ring_matrix< T > &A, const dense_ring_matrix< T > &B)
{
	dense_ring_matrix< T > RES(A.rows, A.columns);
	RES.add(A, B);
	return RES;
}



template< class T >
inline dense_ring_matrix< T >
operator + (const dense_ring_matrix< T > &A, const T &a)
{
	dense_ring_matrix< T > RES(A.rows, A.columns);
	RES.add(A, a);
	return RES;
}



template< class T >
inline dense_ring_matrix< T >
operator + (const T &a, const dense_ring_matrix< T > &A)
{
	dense_ring_matrix< T > RES(A.rows, A.columns);
	RES.add(a, A);
	return RES;
}



template< class T >
inline dense_ring_matrix< T > &
operator += (dense_ring_matrix< T > &A, const dense_ring_matrix< T > &B)
{
	A.add(A, B);
	return A;
}



template< class T >
inline dense_ring_matrix< T > &
operator += (dense_ring_matrix< T > &A, const T &a)
{
	A.add(A, a);
	return A;
}



//
// subtraction
//

template< class T >
inline dense_ring_matrix< T >
operator - (const dense_ring_matrix< T > &A, const dense_ring_matrix< T > &B)
{
	dense_ring_matrix< T > RES(A.rows, A.columns);
	RES.subtract(A, B);
	return RES;
}



template< class T >
inline dense_ring_matrix< T >
operator - (const dense_ring_matrix< T > &A, const T &a)
{
	dense_ring_matrix< T >
		RES(A.rows, A.columns);
	RES.subtract(A, a);
	return RES;
}



template< class T >
inline dense_ring_matrix< T >
operator - (const T &a, const dense_ring_matrix< T > &A)
{
	dense_ring_matrix< T > RES(A.rows, A.columns);
	RES.subtract(a, A);
	return RES;
}



template< class T >
inline dense_ring_matrix< T > &
operator -= (dense_ring_matrix< T > &A, const dense_ring_matrix< T > &B)
{
	A.subtract(A, B);
	return A;
}



template< class T >
inline dense_ring_matrix< T > &
operator -= (dense_ring_matrix< T > &A, const T &a)
{
	A.subtract(A, a);
	return A;
}



//
// multiplication
//

template< class T >
inline dense_ring_matrix< T >
operator * (const dense_ring_matrix< T > &A, const dense_ring_matrix< T > &B)
{
	dense_ring_matrix< T > RES(A.rows, B.columns);
	RES.multiply(A, B);
	return RES;
}



template< class T >
inline dense_ring_matrix< T >
operator * (const dense_ring_matrix< T > &A, const T &m)
{
	dense_ring_matrix< T > RES(A.rows, A.columns);
	RES.multiply(A, m);
	return RES;
}



template< class T >
inline dense_ring_matrix< T >
operator * (const T &m, const dense_ring_matrix< T > &A)
{
	dense_ring_matrix< T > RES(A.rows, A.columns);
	RES.multiply(m, A);
	return RES;
}



template< class T >
inline T *
operator * (const dense_ring_matrix< T > &A, const T *v)
{
	T *b = new T[A.rows];
	memory_handler(b, "dense_ring_matrix", "operator * :: "
		       "Error in memory allocation (b)");
	A.multiply_right(b, v);
	return b;
}



template< class T >
inline T *
operator * (const T *v, const dense_ring_matrix< T > &A)
{
	T *b = new T[A.columns];
	memory_handler(b, "dense_ring_matrix", "operator * :: "
		       "Error in memory allocation (b)");
	A.multiply_left(b, v);
	return b;
}



template< class T >
inline math_vector< T >
operator * (const dense_ring_matrix< T > &A, const math_vector< T > &v)
{
	math_vector< T > b;
	A.multiply_right(b, v);
	return b;
}



template< class T >
inline math_vector< T >
operator * (const math_vector< T > &v, const dense_ring_matrix< T > &A)
{
	math_vector< T > b;
	A.multiply_left(b, v);
	return b;
}



template< class T >
inline dense_ring_matrix< T > &
operator *= (dense_ring_matrix< T > &A, const dense_ring_matrix< T > &B)
{
	A.multiply(A, B);
	return A;
}



template< class T >
inline dense_ring_matrix< T > &
operator *= (dense_ring_matrix< T > &A, const T &m)
{
	A.multiply(A, m);
	return A;
}



//
// negation
//

template< class T >
inline dense_ring_matrix< T >
operator - (const dense_ring_matrix< T > &A)
{
	dense_ring_matrix< T > B(A.rows, A.columns);
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
operator == (const dense_ring_matrix< T > &A, const dense_ring_matrix< T > &B)
{
	return A.equal(B);
}



template< class T >
inline bool
equal(const dense_ring_matrix< T > &A, const dense_ring_matrix< T > &B)
{
	return A.equal(B);
}



template< class T >
inline bool
operator != (const dense_ring_matrix< T > &A, const dense_ring_matrix< T > &B)
{
	return !A.equal(B);
}



template< class T >
inline bool
unequal(const dense_ring_matrix< T > &A, const dense_ring_matrix< T > &B)
{
	return !A.equal(B);
}



//
// trace
//

template< class T >
inline T
trace(const dense_ring_matrix< T > &A)
{
	T tr;
	A.trace(tr);
	return tr;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/dense_ring_matrix.cc"
#endif

#undef DRMKex



#endif	// LIDIA_DENSE_RING_MATRIX_H_GUARD_
