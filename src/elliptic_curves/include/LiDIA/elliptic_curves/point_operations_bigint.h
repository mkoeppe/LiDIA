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
//	Author	: Markus Maurer (MM), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_POINT_OPERATIONS_BIGINT_H_GUARD_
#define LIDIA_POINT_OPERATIONS_BIGINT_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_ELLIPTIC_CURVE_FLAGS_H_GUARD_
# include	"LiDIA/elliptic_curve_flags.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T > class point_operations;
template< class T > class point;

template<>
class point_operations< bigint >
{
	//
	// types of function pointers
	//
	typedef void (*generic_add_ptr)(point< bigint > &, const point< bigint > &, const point< bigint > &);
	typedef void (*generic_negate_ptr)(point< bigint > &, const point< bigint > &);
	typedef void (*generic_mult_by_2_ptr)(point< bigint > &, const point< bigint > &);

	//
	// friends
	//
#if defined(_MSC_VER)
	friend template point< bigint >;
#else
	friend class point< bigint >;
#endif

	friend void multiply_by_2(point< bigint > &,
				  const point< bigint > &);

	friend void negate(point< bigint > &,
			   const point< bigint > &);

	friend void add(point< bigint > &,
			const point< bigint > &,
			const point< bigint > &);

	friend void subtract(point< bigint > &,
			     const point< bigint > &,
			     const point< bigint > &);

private:
	generic_add_ptr       _add;
	generic_negate_ptr    _negate;
	generic_mult_by_2_ptr _mult_by_2;

	//
	// constructor / destructor
	//
public:
	point_operations();

public:
	~point_operations();

	//
	// assigment
	//
public:
	void set(elliptic_curve_flags::curve_parametrization cp,
		 elliptic_curve_flags::curve_model cm,
		 const bigint& characteristic);

public:
	point_operations< bigint > & operator = (const point_operations< bigint > & p);

};

//
// R = P + Q,
//
// Condition P != +- Q not verified.
//
void add_swnf_projective(point< bigint > & R,
			 const point< bigint > & P,
			 const point< bigint > & Q);

void add_lwnf_projective(point< bigint > &,
			 const point< bigint > &,
			 const point< bigint > &);

//
// R = -P
//
void negate_swnf_projective(point< bigint > &R,
			    const point< bigint > &P);

void negate_lwnf_projective(point< bigint > &,
			    const point< bigint > &);

//
// R = 2 * P
//
void mult_by_2_swnf_projective(point< bigint > &R,
			       const point< bigint > &P);

void mult_by_2_lwnf_projective(point< bigint > &,
			       const point< bigint > &);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_POINT_OPERATIONS_BIGINT_H_GUARD_
