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


#ifndef LIDIA_MATH_VECTOR_H_GUARD_
#define LIDIA_MATH_VECTOR_H_GUARD_



#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_MODULAR_OPERATIONS_INL_GUARD_
# include	"LiDIA/modular_operations.inl"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T >
class math_vector : public virtual base_vector< T >
{

	//
	// constructors
	//

public:

	math_vector():
		base_vector< T > () {}

	explicit math_vector(const vector_flags &md):
		base_vector< T > (md) {}
	explicit math_vector(lidia_size_t all):
		base_vector< T > (all) {}

	math_vector(lidia_size_t all, const vector_flags &md):
		base_vector< T > (all, md) {}
	math_vector(lidia_size_t all, lidia_size_t len):
		base_vector< T > (all, len) {}
	math_vector(lidia_size_t all, lidia_size_t len, const vector_flags &md):
		base_vector< T > (all, len, md) {}
	math_vector(const base_vector< T > &v):
		base_vector< T > (v) {}
	math_vector(const base_vector< T > & v, const vector_flags &md):
		base_vector< T > (v, md) {}
	math_vector(const T *v, lidia_size_t len):
		base_vector< T > (v, len) {}
	math_vector(const T *v, lidia_size_t len, const vector_flags &md):
		base_vector< T > (v, len, md) {}

	//
	// destructor
	//

public:

	~math_vector() {}


	//
	// math procedures
	//

	//
	// addition
	//

	void add(const math_vector< T > &, const math_vector< T > &);
	void add(const math_vector< T > &, const T &);
	void add(const T &, const math_vector< T > &);

	//
	// subtraction
	//

	void subtract(const math_vector< T > &, const math_vector< T > &);
	void subtract(const math_vector< T > &, const T &);
	void subtract(const T &, const math_vector< T > &);

	//
	// multiplication
	//

	void compwise_multiply(const math_vector< T > &, const math_vector< T > &);
	void multiply(T &, const math_vector< T > &) const;
	void right_multiply(const math_vector< T > &, const T &);
	void left_multiply(const T &, const math_vector< T > &);

	//
	// divide
	//

	void compwise_divide(const math_vector< T > &, const math_vector< T > &);

	void divide(const math_vector< T > &, const T &);
	void divide(const T &, const math_vector< T > &);

	//
	// negation
	//

	void negate(const math_vector< T > &);

	//
	// miscellaneous
	//

public:

	T sum_of_squares() const;

	//
	// comparisons
	//

public:

	bool equal(const math_vector< T > &) const;
	bool unequal(const math_vector< T > &v) const
	{
		return (!equal(v));
	}
};



//
// math procedures
//

//
// addition
//

template< class T >
inline void
add(math_vector< T > &r, const math_vector< T > &v, const math_vector< T > &w)
{
	r.add(v, w);
}



template< class T >
inline void
add(math_vector< T > &r, const math_vector< T > &v, const T &t)
{
	r.add(v, t);
}



template< class T >
inline void
add(math_vector< T > &r, const T &t, const math_vector< T > &v)
{
	r.add(t, v);
}



//
// subtraction
//

template< class T >
inline void
subtract(math_vector< T > &r, const math_vector< T > &v, const math_vector< T > &w)
{
	r.subtract(v, w);
}



template< class T >
inline void
subtract(math_vector< T > &r, const math_vector< T > &v, const T &t)
{
	r.subtract(v, t);
}



template< class T >
inline void
subtract(math_vector< T > &r, const T &t, const math_vector< T > &v)
{
	r.subtract(t, v);
}



//
// multiplication
//

template< class T >
inline void
compwise_multiply(math_vector< T > &r, const math_vector< T > &v, const math_vector< T > &w)
{
	r.compwise_multiply(v, w);
}



template< class T >
inline void
multiply(T &res, const math_vector< T > &v, const math_vector< T > &w)
{
	v.multiply(res, w);
}



template< class T >
inline void
multiply(math_vector< T > &r, const math_vector< T > &v, const T &t)
{
	r.right_multiply(v, t);
}



template< class T >
inline void
multiply(math_vector< T > &r, const T &t, const math_vector< T > &v)
{
	r.left_multiply(t, v);
}



//
// divide
//

template< class T >
inline void
compwise_divide(math_vector< T > &r, const math_vector< T > &v, const math_vector< T > &w)
{
	r.compwise_divide(v, w);
}



template< class T >
inline void
divide(math_vector< T > &r, const math_vector< T > &v, const T &t)
{
	r.divide(v, t);
}



template< class T >
inline void
divide(math_vector< T > &r, const T &t, const math_vector< T > &v)
{
	r.divide(t, v);
}



//
// negation
//

template< class T >
inline void
negate(math_vector< T > &r, const math_vector< T > &v)
{
	r.negate(v);
}



//
// math operations
//

//
// addition
//

template< class T >
inline math_vector< T >
operator + (const math_vector< T > &v, const math_vector< T > &w)
{
	math_vector< T > sum(v.size(), vector_flags::fixed);

	sum.add(v, w);
	return sum;
}



template< class T >
inline math_vector< T >
operator + (const math_vector< T > &v, const T &t)
{
	math_vector< T > sum(v.size(), vector_flags::fixed);

	sum.add(v, t);
	return sum;
}



template< class T >
inline math_vector< T >
operator + (const T &t, const math_vector< T > &v)
{
	math_vector< T > sum(v.size(), vector_flags::fixed);

	sum.add(t, v);
	return sum;
}



template< class T >
inline math_vector< T > &
operator += (math_vector< T > &sum, const math_vector< T > &v)
{
	sum.add(sum, v);
	return sum;
}



template< class T >
inline math_vector< T > &
operator += (math_vector< T > &sum, const T &t)
{
	sum.add(sum, t);
	return sum;
}



//
// subtraction
//

template< class T >
inline math_vector< T >
operator - (const math_vector< T > &v, const math_vector< T > &w)
{
	math_vector< T > dif (v.size(), vector_flags::fixed);

	dif.subtract(v, w);
	return dif;
}



template< class T >
inline math_vector< T >
operator - (const math_vector< T > &v, const T &t)
{
	math_vector< T > dif (v.size(), vector_flags::fixed);

	dif.subtract(v, t);
	return dif;
}



template< class T >
inline math_vector< T >
operator - (const T &t, const math_vector< T > & v)
{
	math_vector< T > dif (v.size(), vector_flags::fixed);

	dif.subtract(t, v);
	return dif;
}



template< class T >
inline math_vector< T > &
operator -= (math_vector< T > &dif, const math_vector< T > &v)
{
	dif.subtract(dif, v);
	return dif;
}



template< class T >
inline math_vector< T > &
operator -= (math_vector< T > &dif, const T &t)
{
	dif.subtract(dif, t);
	return dif;
}



//
// multiplication
//

template< class T >
inline T
operator * (const math_vector< T > &v, const math_vector< T > &w)
{
	T pro;

	v.multiply(pro, w);
	return pro;
}



template< class T >
inline math_vector< T >
operator * (const math_vector< T > &v, const T &t)
{
	math_vector< T > pro(v.size(), vector_flags::fixed);

	pro.right_multiply(v, t);
	return pro;
}



template< class T >
inline math_vector< T >
operator * (const T &t, const math_vector< T > &v)
{
	math_vector< T > pro(v.size(), vector_flags::fixed);

	pro.left_multiply(t, v);
	return pro;
}



template< class T >
inline math_vector< T > &
operator *= (math_vector< T > &pro, const T &t)
{
	pro.right_multiply(pro, t);
	return pro;
}



//
// division
//

template< class T >
inline math_vector< T >
operator / (const math_vector< T > &v, const T &t)
{
	math_vector< T > quo(v.size(), vector_flags::fixed);

	quo.divide(v, t);
	return quo;
}



template< class T >
inline math_vector< T >
operator / (const T &t, const math_vector< T > &v)
{
	math_vector< T > quo(v.size(), vector_flags::fixed);

	quo.divide(t, v);
	return quo;
}



template< class T >
inline math_vector< T > &
operator /= (math_vector< T > &quo, const T &t)
{
	quo.divide(quo, t);
	return quo;
}



//
// negation
//

template< class T >
inline math_vector< T >
operator - (const math_vector< T > &v)
{
	math_vector< T > neg(v.size(), vector_flags::fixed);

	neg.negate(v);
	return neg;
}



//
// miscellaneous
//

template< class T >
inline T
sum_of_squares(const math_vector< T > &v)
{
	return v.sum_of_squares();
}



//
// comparisons
//

template< class T >
inline bool
operator == (const math_vector< T > &w, const math_vector< T > &v)
{
	return w.equal(v);
}



template< class T >
inline bool
equal(const math_vector< T > &w, const math_vector< T > &v)
{
	return w.equal(v);
}



template< class T >
inline bool
operator != (const math_vector< T > &w, const math_vector< T > &v)
{
	return !w.equal(v);
}



template< class T >
inline bool
unequal(const math_vector< T > &w, const math_vector< T > &v)
{
	return !w.equal(v);
}



//
// special
//

inline math_vector< bigint >
operator % (const math_vector< bigint > &v, bigint &mod)
{
	lidia_size_t i, l = v.size();

	math_vector< bigint > w(l, l);
	for (i = 0; i < l; i++) {
		if (v[i] >= mod || v[i] <= -mod)
			w[i] = v[i] % mod;
		else
			w[i] = v[i];
	}
	return w;
}



inline math_vector< long >
operator % (const math_vector< long > &v, long mod)
{
	lidia_size_t i, l = v.size();

	math_vector< long > w(l, l);
	for (i = 0; i < l; i++) {
		if (v[i] >= mod || v[i] <= -mod)
			best_remainder(w[i], v[i], mod);
		else
			w[i] = v[i];
	}
	return w;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/math_vector.cc"
#endif



#endif	// LIDIA_MATH_VECTOR_H_GUARD_
