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
//	$Id: ring_matrix.h,v 2.6 2004/06/15 10:19:18 lidiaadm Exp $
//
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_RING_MATRIX_H_GUARD_
#define LIDIA_RING_MATRIX_H_GUARD_



#ifndef LIDIA_BASE_MATRIX_H_GUARD_
# include	"LiDIA/base_matrix.h"
#endif
#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif
#ifndef LIDIA_DENSE_RING_MATRIX_H_GUARD_
# include	"LiDIA/dense_ring_matrix.h"
#endif
#ifndef LIDIA_SPARSE_RING_MATRIX_H_GUARD_
# include	"LiDIA/sparse_ring_matrix.h"
#endif
#ifndef LIDIA_RING_MATRIX_ALGORITHMS_H_GUARD_
# include	"LiDIA/matrix/ring_matrix_algorithms.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



#define DRMKex DRMK< T >
#define SRMKex SRMK< T >



template< class T >
class ring_matrix : public base_matrix< T >
{
	//
	// modul definitions
	//

	const DRMKex D_ring_modul;
	const SRMKex S_ring_modul;

	const RMA< T, DRMKex, DRMKex, DRMKex > DDD_ring_modul;
	const RMA< T, DRMKex, DRMKex, SRMKex > DDS_ring_modul;
	const RMA< T, DRMKex, SRMKex, DRMKex > DSD_ring_modul;
	const RMA< T, SRMKex, DRMKex, DRMKex > SDD_ring_modul;
	const RMA< T, DRMKex, SRMKex, SRMKex > DSS_ring_modul;
	const RMA< T, SRMKex, DRMKex, SRMKex > SDS_ring_modul;
	const RMA< T, SRMKex, SRMKex, DRMKex > SSD_ring_modul;
	const RMA< T, SRMKex, SRMKex, SRMKex > SSS_ring_modul;

	//
	// constructors
	//

public:

	ring_matrix()
		: base_matrix< T > () {}
	ring_matrix(const matrix_flags &flags)
		: base_matrix< T > (flags) {}
	ring_matrix(lidia_size_t a, lidia_size_t b)
		: base_matrix< T > (a, b) {}
	ring_matrix(lidia_size_t a, lidia_size_t b, const matrix_flags &flags)
		: base_matrix< T > (a, b, flags) {}
	explicit ring_matrix(const base_vector< T > &v)
		: base_matrix< T > (v) {}
	ring_matrix(const base_vector< T > &v, const matrix_flags &flags)
		: base_matrix< T > (v, flags) {}
	ring_matrix(const base_matrix< T > &A)
		: base_matrix< T > (A) {}
	ring_matrix(const base_matrix< T > &A, const matrix_flags &flags)
		: base_matrix< T > (A, flags) {}
	ring_matrix(const ring_matrix< T > &A)
		: base_matrix< T > (A) {}
	ring_matrix(const ring_matrix< T > &A, const matrix_flags &flags)
		: base_matrix< T > (A, flags) {}
	ring_matrix(const dense_base_matrix< T > &A)
		: base_matrix< T > (A) {}
	ring_matrix(const sparse_base_matrix< T > &A)
		: base_matrix< T > (A) {}
	ring_matrix(lidia_size_t a, lidia_size_t b, const T **v)
		: base_matrix< T > (a, b, v) {}
	ring_matrix(lidia_size_t a, lidia_size_t b, const T **v, const matrix_flags &flags)
		: base_matrix< T > (a, b, v, flags) {}

	//
	// destructor
	//

public:

	~ring_matrix() {}

	//
	// assignment
	//

public:

	using base_matrix< T >::assign;

	ring_matrix< T > & operator = (const ring_matrix< T > &B)
	{
// Changed by G.A.	
		base_matrix< T >::assign(B);
		return *this;
	}

	//////////////////////////////////
	// BEGIN: arithmetic procedures //
	//////////////////////////////////

	//
	// addition
	//

	void add(const ring_matrix< T > &, const ring_matrix< T > &);
	void add(const ring_matrix< T > &, const T &);
	void add(const T &, const ring_matrix< T > &);

	//
	// subtraction
	//

	void subtract(const ring_matrix< T > &, const ring_matrix< T > &);
	void subtract(const ring_matrix< T > &, const T &);
	void subtract(const T &, const ring_matrix< T > &);

	//
	// multiplication
	//

	void multiply(const ring_matrix< T > &, const ring_matrix< T > &);
	void multiply(const ring_matrix< T > &, const T &);
	void multiply(const T &, const ring_matrix< T > &);

	void compwise_multiply(const ring_matrix< T > &,
			       const ring_matrix< T > &);

	void multiply_right(T *&, const T *) const;
	void multiply_right(math_vector< T > &, const math_vector< T > &) const;
	void multiply_left(T *&, const T *) const;
	void multiply_left(math_vector< T > &, const math_vector< T > &) const;

	//
	// negation
	//

	void negate(const ring_matrix< T > &);

	/////////////////////////////////
	// END: arithmetic procedures  //
	// BEGIN: arithmetic operators //
	/////////////////////////////////

	//
	// addition
	//

	///////////////////////////////
	// END: arithmetic operators //
	///////////////////////////////

	//
	// comparisons
	//

	bool equal(const ring_matrix< T > &) const;
	bool unequal(const ring_matrix< T > &A) const
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
	// transpose function for ring_matrix
	//

	ring_matrix< T > trans() const;
	void trans(const ring_matrix< T > & A);


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
add(ring_matrix< T > &RES, const ring_matrix< T > &M,
    const ring_matrix< T > &N)
{
	RES.add(M, N);
}



template< class T >
inline void
add(ring_matrix< T > &RES, const ring_matrix< T > &M, const T &a)
{
	RES.add(M, a);
}



template< class T >
inline void
add(ring_matrix< T > &RES, const T &a, const ring_matrix< T > &M)
{
	RES.add(a, M);
}



//
// subtraction
//

template< class T >
inline void
subtract(ring_matrix< T > &RES, const ring_matrix< T > &M,
	 const ring_matrix< T > &N)
{
	RES.subtract(M, N);
}



template< class T >
inline void
subtract(ring_matrix< T > &RES, const ring_matrix< T > &M, const T &a)
{
	RES.subtract(M, a);
}



template< class T >
inline void
subtract(ring_matrix< T > &RES, const T &a, const ring_matrix< T > &M)
{
	RES.subtract(a, M);
}



//
// multiplication
//

template< class T >
inline void
multiply(ring_matrix< T > &RES, const ring_matrix< T > &M,
	 const ring_matrix< T > &N)
{
	RES.multiply(M, N);
}



template< class T >
inline void
multiply(ring_matrix< T > &RES, const ring_matrix< T > &M, const T &a)
{
	RES.multiply(M, a);
}



template< class T >
inline void
multiply(ring_matrix< T > &RES, const T &a, const ring_matrix< T > &M)
{
	RES.multiply(a, M);
}



template< class T >
inline void
compwise_multiply(ring_matrix< T > &RES, const ring_matrix< T > &M,
		  const ring_matrix< T > &N)
{
	RES.compwise_multiply(M, N);
}



template< class T >
inline void
multiply(T *&res, const ring_matrix< T > &A, const T *v)
{
	A.multiply_right(res, v);
}



template< class T >
inline void
multiply(T *&res, const T *v, const ring_matrix< T > &A)
{
	A.multiply_left(res, v);
}



template< class T >
inline void
multiply(math_vector< T > &res, const ring_matrix< T > &A,
	 const math_vector< T > &v)
{
	A.multiply_right(res, v);
}



template< class T >
inline void
multiply(math_vector< T > &res, const math_vector< T > &v,
	 const ring_matrix< T > &A)
{
	A.multiply_left(res, v);
}



//
// negation
//

template< class T >
inline void
negate(ring_matrix< T > &A, const ring_matrix< T > &B)
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
inline ring_matrix< T >
operator + (const ring_matrix< T > &A, const ring_matrix< T > &B)
{
	ring_matrix< T > RES(A.rows, A.columns, A.bitfield);
	RES.add(A, B);
	return RES;
}



template< class T >
inline ring_matrix< T >
operator + (const ring_matrix< T > &A, const T &a)
{
	ring_matrix< T > RES(A.rows, A.columns, A.bitfield);
	RES.add(A, a);
	return RES;
}



template< class T >
inline ring_matrix< T >
operator + (const T &a, const ring_matrix< T > &A)
{
	ring_matrix< T > RES(A.rows, A.columns, A.bitfield);
	RES.add(a, A);
	return RES;
}



template< class T >
inline ring_matrix< T > &
operator += (ring_matrix< T > &A, const ring_matrix< T > &B)
{
	A.add(A, B);
	return A;
}



template< class T >
inline ring_matrix< T > &
operator += (ring_matrix< T > &A, const T &a)
{
	A.add(A, a);
	return A;
}



//
// subtraction
//

template< class T >
inline ring_matrix< T >
operator - (const ring_matrix< T > &A, const ring_matrix< T > &B)
{
	ring_matrix< T > RES(A.rows, A.columns, A.bitfield);
	RES.subtract(A, B);
	return RES;
}



template< class T >
inline ring_matrix< T >
operator - (const ring_matrix< T > &A, const T &a)
{
	ring_matrix< T > RES(A.rows, A.columns, A.bitfield);
	RES.subtract(A, a);
	return RES;
}



template< class T >
inline ring_matrix< T >
operator - (const T &a, const ring_matrix< T > &A)
{
	ring_matrix< T > RES(A.rows, A.columns, A.bitfield);
	RES.subtract(a, A);
	return RES;
}



template< class T >
inline ring_matrix< T > &
operator -= (ring_matrix< T > &A, const ring_matrix< T > &B)
{
	A.subtract(A, B);
	return A;
}



template< class T >
inline ring_matrix< T > &
operator -= (ring_matrix< T > &A, const T &a)
{
	A.subtract(A, a);
	return A;
}



//
// multiplication
//

template< class T >
inline ring_matrix< T >
operator * (const ring_matrix< T > &A, const ring_matrix< T > &B)
{
	ring_matrix< T > RES(A.rows, B.columns, A.bitfield);
	RES.multiply(A, B);
	return RES;
}



template< class T >
inline ring_matrix< T >
operator * (const ring_matrix< T > &A, const T &m)
{
	ring_matrix< T > RES(A.rows, A.columns, A.bitfield);
	RES.multiply(A, m);
	return RES;
}



template< class T >
inline ring_matrix< T >
operator * (const T &m, const ring_matrix< T > &A)
{
	ring_matrix< T > RES(A.rows, A.columns, A.bitfield);
	RES.multiply(m, A);
	return RES;
}



template< class T >
inline T *
operator * (const ring_matrix< T > &A, const T *v)
{
	T *b = new T[A.rows];
	memory_handler(b, "ring_matrix", "operator * :: "
		       "Error in memory allocation (b)");
	A.multiply_right(b, v);
	return b;
}



template< class T >
inline T *
operator * (const T *v, const ring_matrix< T > &A)
{
	T *b = new T[A.columns];
	memory_handler(b, "ring_matrix", "operator * :: "
		       "Error in memory allocation (b)");
	A.multiply_left(b, v);
	return b;
}



template< class T >
inline math_vector< T >
operator * (const ring_matrix< T > &A, const math_vector< T > &v)
{
	math_vector< T > b;
	A.multiply_right(b, v);
	return b;
}



template< class T >
inline math_vector< T >
operator * (const math_vector< T > &v, const ring_matrix< T > &A)
{
	math_vector< T > b;
	A.multiply_left(b, v);
	return b;
}



template< class T >
inline ring_matrix< T > &
operator *= (ring_matrix< T > &A, const ring_matrix< T > &B)
{
	A.multiply(A, B);
	return A;
}



template< class T >
inline ring_matrix< T > &
operator *= (ring_matrix< T > &A, const T &m)
{
	A.multiply(A, m);
	return A;
}



//
// negation
//

template< class T >
inline ring_matrix< T >
operator - (const ring_matrix< T > &A)
{
	ring_matrix< T > B(A.rows, A.columns, A.bitfield);
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
operator == (const ring_matrix< T > &A, const ring_matrix< T > &B)
{
	return A.equal(B);
}



template< class T >
inline bool
equal(const ring_matrix< T > &A, const ring_matrix< T > &B)
{
	return A.equal(B);
}



template< class T >
inline bool
operator != (const ring_matrix< T > &A, const ring_matrix< T > &B)
{
	return !A.equal(B);
}



template< class T >
inline bool
unequal(const ring_matrix< T > &A, const ring_matrix< T > &B)
{
	return !A.equal(B);
}



//
// trace
//

template< class T >
inline T
trace(const ring_matrix< T > &A)
{
	T tr;
	A.trace(tr);
	return tr;
}



//
// transpose function for ring matrix
//

template< class T >
inline void
ring_matrix< T >::trans(const ring_matrix< T > &A)
{
	base_matrix< T >::trans(A);
}



template< class T >
inline ring_matrix< T >
ring_matrix< T >::trans() const
{
	return ring_matrix< T >(base_matrix< T >::trans());
}



template< class T >
inline void
trans(ring_matrix< T > &AT, const base_matrix< T > &A)
{
	AT.trans(A);
}



template< class T >
inline ring_matrix< T >
trans(const ring_matrix< T > &A)
{
	return A.trans();
}



//
// specialization for type bigmod
//

template<>
class ring_matrix< bigmod > : public base_matrix< bigmod >
{
	//
	// modul definitions
	//

	const DRMK< bigint > D_ring_modul_bigint;
	const SRMK< bigint > S_ring_modul_bigint;

	const DRMK< long > D_ring_modul_long;
	const SRMK< long > S_ring_modul_long;

	const RMA< bigint, DRMK< bigint >, DRMK< bigint >, DRMK< bigint > > DDD_ring_modul_bigint;
	const RMA< bigint, DRMK< bigint >, DRMK< bigint >, SRMK< bigint > > DDS_ring_modul_bigint;
	const RMA< bigint, DRMK< bigint >, SRMK< bigint >, DRMK< bigint > > DSD_ring_modul_bigint;
	const RMA< bigint, SRMK< bigint >, DRMK< bigint >, DRMK< bigint > > SDD_ring_modul_bigint;
	const RMA< bigint, DRMK< bigint >, SRMK< bigint >, SRMK< bigint > > DSS_ring_modul_bigint;
	const RMA< bigint, SRMK< bigint >, DRMK< bigint >, SRMK< bigint > > SDS_ring_modul_bigint;
	const RMA< bigint, SRMK< bigint >, SRMK< bigint >, DRMK< bigint > > SSD_ring_modul_bigint;
	const RMA< bigint, SRMK< bigint >, SRMK< bigint >, SRMK< bigint > > SSS_ring_modul_bigint;

	const RMA< long, DRMK< long >, DRMK< long >, DRMK< long > > DDD_ring_modul_long;
	const RMA< long, DRMK< long >, DRMK< long >, SRMK< long > > DDS_ring_modul_long;
	const RMA< long, DRMK< long >, SRMK< long >, DRMK< long > > DSD_ring_modul_long;
	const RMA< long, SRMK< long >, DRMK< long >, DRMK< long > > SDD_ring_modul_long;
	const RMA< long, DRMK< long >, SRMK< long >, SRMK< long > > DSS_ring_modul_long;
	const RMA< long, SRMK< long >, DRMK< long >, SRMK< long > > SDS_ring_modul_long;
	const RMA< long, SRMK< long >, SRMK< long >, DRMK< long > > SSD_ring_modul_long;
	const RMA< long, SRMK< long >, SRMK< long >, SRMK< long > > SSS_ring_modul_long;

	//
	// constructors
	//

public:

	ring_matrix< bigmod > ()
		: base_matrix< bigmod > (0) {}
	ring_matrix< bigmod > (const bigint &mod)
		: base_matrix< bigmod > (mod) {}
	ring_matrix< bigmod > (const bigint &mod, const matrix_flags &flags)
		: base_matrix< bigmod > (mod, flags) {}
	ring_matrix< bigmod > (lidia_size_t a, lidia_size_t b, const bigint &mod)
		: base_matrix< bigmod > (a, b, mod) {}
	ring_matrix< bigmod > (lidia_size_t a, lidia_size_t b, const bigint &mod, const matrix_flags &flags)
		: base_matrix< bigmod > (a, b, mod, flags) {}
	ring_matrix< bigmod > (const ring_matrix< bigmod > &A)
		: base_matrix< bigmod > (A) {}
	ring_matrix< bigmod > (const ring_matrix< bigmod > &A, const matrix_flags &flags)
		: base_matrix< bigmod > (A, flags) {}

	ring_matrix< bigmod > (const ring_matrix< bigint > &A, const bigint &mod)
		: base_matrix< bigmod > (A, mod) {}
	ring_matrix< bigmod > (const ring_matrix< long > &A, const long &mod)
		: base_matrix< bigmod > (A, mod) {}

	//
	// destructor
	//

public:

	~ring_matrix() {}

	//
	// assignment
	//

public:

	ring_matrix< bigmod > & operator = (const ring_matrix< bigmod > &B)
	{
		assign(B);
		return *this;
	}

	//////////////////////////////////
	// BEGIN: arithmetic procedures //
	//////////////////////////////////

  //
  // addition
  //

public:

	friend void add(ring_matrix< bigmod > &,
			const ring_matrix< bigmod > &,
			const ring_matrix< bigmod > &);

	friend void add(ring_matrix< bigmod > &,
			const ring_matrix< bigmod > &, const bigmod &);

	friend void add(ring_matrix< bigmod > &,
			const bigmod &, const ring_matrix< bigmod > &);

protected:

	void add(const ring_matrix< bigmod > &, const ring_matrix< bigmod > &);
	void add(const ring_matrix< bigmod > &, const bigmod &);
	void add(const bigmod &, const ring_matrix< bigmod > &);

	//
	// subtraction
	//

public:

	friend void subtract(ring_matrix< bigmod > &,
			     const ring_matrix< bigmod > &,
			     const ring_matrix< bigmod > &);

	friend void subtract(ring_matrix< bigmod > &,
			     const ring_matrix< bigmod > &, const bigmod &);

	friend void subtract(ring_matrix< bigmod > &,
			     const bigmod &, const ring_matrix< bigmod > &);

protected:

	void subtract(const ring_matrix< bigmod > &, const ring_matrix< bigmod > &);
	void subtract(const ring_matrix< bigmod > &, const bigmod &);
	void subtract(const bigmod &, const ring_matrix< bigmod > &);

	//
	// multiplication
	//

public:

	friend void multiply(ring_matrix< bigmod > &,
			     const ring_matrix< bigmod > &,
			     const ring_matrix< bigmod > &);

	friend void multiply(ring_matrix< bigmod > &,
			     const ring_matrix< bigmod > &, const bigmod &);

	friend void multiply(ring_matrix< bigmod > &,
			     const bigmod &, const ring_matrix< bigmod > &);

protected:

	void multiply(const ring_matrix< bigmod > &, const ring_matrix< bigmod > &);
	void multiply(const ring_matrix< bigmod > &, const bigmod &);
	void multiply(const bigmod &, const ring_matrix< bigmod > &);

public:

	friend void compwise_multiply(ring_matrix< bigmod > &,
				      const ring_matrix< bigmod > &,
				      const ring_matrix< bigmod > &);

protected:

	void compwise_multiply(const ring_matrix< bigmod > &,
			       const ring_matrix< bigmod > &);

public:

	friend void multiply(bigmod *&, const ring_matrix< bigmod > &, const bigmod *);

	friend void multiply(bigmod *&, const bigmod *, const ring_matrix< bigmod > &);

	friend void multiply(math_vector< bigmod > &, const ring_matrix< bigmod > &,
			     const math_vector< bigmod > &);

	friend void multiply(math_vector< bigmod > &, const math_vector< bigmod > &,
			     const ring_matrix< bigmod > &);

protected:

	void multiply_right(bigmod *&, const bigmod *) const;
	void multiply_right(math_vector< bigmod > &, const math_vector< bigmod > &) const;
	void multiply_left(bigmod *&, const bigmod *) const;
	void multiply_left(math_vector< bigmod > &, const math_vector< bigmod > &) const;

	//
	// negation
	//

public:

	friend void negate(ring_matrix< bigmod > &,
			   const ring_matrix< bigmod > &);

protected:

	void negate(const ring_matrix< bigmod > &);

	/////////////////////////////////
	// END: arithmetic procedures  //
	// BEGIN: arithmetic operators //
	/////////////////////////////////

	//
	// addition
	//

public:

	friend ring_matrix< bigmod > operator + (const ring_matrix< bigmod > &, const ring_matrix< bigmod > &);

	friend ring_matrix< bigmod > operator + (const ring_matrix< bigmod > &, const bigmod &);

	friend ring_matrix< bigmod > operator + (const bigmod &, const ring_matrix< bigmod > &);

	friend ring_matrix< bigmod > & operator += (ring_matrix< bigmod > &, const ring_matrix< bigmod > &);

	friend ring_matrix< bigmod > & operator += (ring_matrix< bigmod > &, const bigmod &);

	//
	// subtraction
	//

public:

	friend ring_matrix< bigmod > operator - (const ring_matrix< bigmod > &,
						 const ring_matrix< bigmod > &);

	friend ring_matrix< bigmod > operator - (const ring_matrix< bigmod > &,
						 const bigmod &);

	friend ring_matrix< bigmod > operator - (const bigmod &,
						 const ring_matrix< bigmod > &);

	friend ring_matrix< bigmod > & operator -= (ring_matrix< bigmod > &,
						    const ring_matrix< bigmod > &);

	friend ring_matrix< bigmod > & operator -= (ring_matrix< bigmod > &,
						    const bigmod &);

	//
	// multiplication
	//

public:

	friend ring_matrix< bigmod > operator * (const ring_matrix< bigmod > &,
						 const ring_matrix< bigmod > &);

	friend ring_matrix< bigmod > operator * (const ring_matrix< bigmod > &,
						 const bigmod &);

	friend ring_matrix< bigmod > operator * (const bigmod &,
						 const ring_matrix< bigmod > &);

	friend bigmod * operator * (const ring_matrix< bigmod > &, const bigmod *);

	friend bigmod * operator * (const bigmod *, const ring_matrix< bigmod > &);

public:

	friend math_vector< bigmod > operator * (const ring_matrix< bigmod > &,
						 const math_vector< bigmod > &);

	friend math_vector< bigmod > operator * (const math_vector< bigmod > &,
						 const ring_matrix< bigmod > &);

	friend ring_matrix< bigmod > & operator *= (ring_matrix< bigmod > &,
						    const ring_matrix< bigmod > &);

	friend ring_matrix< bigmod > & operator *= (ring_matrix< bigmod > &,
						    const bigmod &);
	//
	// negation
	//

public:

	friend ring_matrix< bigmod > operator - (const ring_matrix< bigmod > &);

	///////////////////////////////
	// END: arithmetic operators //
	///////////////////////////////

	//
	// comparisons
	//

	friend bool operator == (const ring_matrix< bigmod > &,
				 const ring_matrix< bigmod > &);

	bool equal(const ring_matrix< bigmod > &) const;
	friend bool equal(const ring_matrix< bigmod > &,
			  const ring_matrix< bigmod > &);

	friend bool operator != (const ring_matrix< bigmod > &,
				 const ring_matrix< bigmod > &);

	bool unequal(const ring_matrix< bigmod > &A) const
	{
		return (!equal(A));
	}

	friend bool unequal(const ring_matrix< bigmod > &,
			    const ring_matrix< bigmod > &);

	//
	// trace
	//

	void trace(bigmod &) const;

	bigmod trace() const
	{
		bigmod tr;
		trace(tr);
		return tr;
	}

	friend bigmod trace(const ring_matrix< bigmod > &);


	//
	// transpose funstion for ring_matrix< bigmod >
	//

	ring_matrix< bigmod > trans () const;
	void trans (const ring_matrix< bigmod > & B);



	//
	// size reduction
	//

	void size_reduction();
};



//
// transpose function for ring_matrix< bigmod >
//

inline void
ring_matrix< bigmod >::trans(const ring_matrix< bigmod > & B)
{
	base_matrix< bigmod >::trans(B);
}



inline ring_matrix< bigmod >
ring_matrix< bigmod >::trans() const
{
	ring_matrix< bigmod > AT;

	AT.base_matrix< bigmod >::trans(*this);
	return AT;
}



inline ring_matrix< bigmod >
trans(const ring_matrix< bigmod > &A)
{
	return A.trans();
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#undef DRMKex
#undef SRMKex



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/ring_matrix.cc"
#endif



#endif	// LIDIA_RING_MATRIX_H_GUARD_
