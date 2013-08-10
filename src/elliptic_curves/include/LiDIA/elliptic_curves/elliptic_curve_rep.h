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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_ELLIPTIC_CURVE_REP_H_GUARD_
#define LIDIA_ELLIPTIC_CURVE_REP_H_GUARD_


#ifndef LIDIA_BASE_ELLIPTIC_CURVE_REP_H_GUARD_
# include	"LiDIA/elliptic_curves/base_elliptic_curve_rep.h"
#endif
#ifndef LIDIA_GF_ELEMENT_H_GUARD_
# include	"LiDIA/gf_element.h"
#endif
#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_RATIONAL_FACTORIZATION_H_GUARD_
# include	"LiDIA/rational_factorization.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T > class point;
template< class T > class elliptic_curve;

template< class T >
class elliptic_curve_rep : public base_elliptic_curve_rep< T >
{
	//
	// constructor / destructor
	//
public:
	elliptic_curve_rep();
	elliptic_curve_rep(const elliptic_curve_rep< T > &);

	elliptic_curve_rep(
		const T & x4, const T & x6,
		elliptic_curve_flags::curve_model m =
		elliptic_curve_flags::AFFINE);

	elliptic_curve_rep(
		const T & x1, const T & x2,
		const T & x3, const T & x4, const T & x6,
		elliptic_curve_flags::curve_model m =
		elliptic_curve_flags::AFFINE);

	~elliptic_curve_rep();

  //
  // assignment
  //
	elliptic_curve_rep< T > & operator = (const elliptic_curve_rep< T > &);
	void assign(const elliptic_curve_rep< T > &);

	T j_invariant() const;
};



template<>
class elliptic_curve_rep< gf_element > : public base_elliptic_curve_rep< gf_element >
{
	friend class elliptic_curve< gf_element >;

	//
	// constructor / destructor
	//
public:
	elliptic_curve_rep();
	elliptic_curve_rep(const elliptic_curve_rep< gf_element > &);

	elliptic_curve_rep(
		const gf_element & x4, const gf_element & x6,
		elliptic_curve_flags::curve_model m =
		elliptic_curve_flags::AFFINE);

	elliptic_curve_rep(
		const gf_element & x1, const gf_element & x2,
		const gf_element & x3, const gf_element & x4,
		const gf_element & x6,
		elliptic_curve_flags::curve_model m =
		elliptic_curve_flags::AFFINE);

	~elliptic_curve_rep();

	//
	// assignment
	//
	elliptic_curve_rep< gf_element > & operator = (const elliptic_curve_rep< gf_element > &);
	void assign(const elliptic_curve_rep< gf_element > &);

	//
	// invariants
	//
	gf_element j_invariant() const;

	//
	// properties
	//
	unsigned int degree_of_definition() const;
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/elliptic_curve_rep.cc"
#endif



#endif	// LIDIA_ELLIPTIC_CURVE_REP_H_GUARD_
