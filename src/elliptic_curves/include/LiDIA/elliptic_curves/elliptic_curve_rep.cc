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


#ifndef LIDIA_ELLIPTIC_CURVE_REP_CC_GUARD_
#define LIDIA_ELLIPTIC_CURVE_REP_CC_GUARD_


#ifndef LIDIA_ELLIPTIC_CURVE_REP_H_GUARD_
# include	"LiDIA/elliptic_curves/elliptic_curve_rep.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// constructor / destructor
//
template< class T > elliptic_curve_rep< T >
::elliptic_curve_rep()

	: base_elliptic_curve_rep< T > ()
{
	debug_handler ("elliptic_curve_rep< T >",
		       "elliptic_curve_rep()");
}



template< class T > elliptic_curve_rep< T >
::elliptic_curve_rep(const elliptic_curve_rep< T > & e)

	: base_elliptic_curve_rep< T > ()
{
	debug_handler ("elliptic_curve_rep< T >",
		       "elliptic_curve_rep(const elliptic_curve_rep< T > &)");

	this->assign(e);
}



template< class T > elliptic_curve_rep< T >
::elliptic_curve_rep(const T & x4, const T & x6,
		     elliptic_curve_flags::curve_model m)

	: base_elliptic_curve_rep< T > ()
{
	debug_handler ("elliptic_curve_rep< T >",
		       "elliptic_curve_rep< T > (const T& x4, const T& x6, "
		       "elliptic_curve_flags::curve_model)");

	this->set_coefficients(x4, x6, m);
}



template< class T > elliptic_curve_rep< T >
::elliptic_curve_rep(const T & x1, const T & x2, const T & x3,
		     const T & x4, const T & x6,
		     elliptic_curve_flags::curve_model m)

	: base_elliptic_curve_rep< T > ()
{
	debug_handler ("elliptic_curve_rep< T >",
		       "elliptic_curve_rep< T > (const T&x1, const T&x2, "
		       "const T&x3, const T&x4, const T&x6, "
		       "elliptic_curve_flags::curve_model)");

	this->set_coefficients(x1, x2, x3, x4, x6, m);
}



template< class T > elliptic_curve_rep< T >
::~elliptic_curve_rep()
{
	debug_handler ("elliptic_curve_rep< T >",
		       "~elliptic_curve_rep()");
}



//
// assignment
//

template< class T >
elliptic_curve_rep< T > & elliptic_curve_rep< T >
::operator = (const elliptic_curve_rep< T > & e)
{
	debug_handler ("base_elliptic_curve< T >",
		       "operator = (const elliptic_curve_rep< T >");

	if (this != &e)
		base_elliptic_curve_rep< T >::assign(e);
	return *this;
}



template< class T >
void elliptic_curve_rep< T >
::assign (const elliptic_curve_rep< T > & e)
{
	debug_handler ("base_elliptic_curve< T >",
		       "assign(const elliptic_curve_rep< T >");

	if (this != &e)
		base_elliptic_curve_rep< T >::assign(e);
}



//
// invariants
//

template< class T >
T elliptic_curve_rep< T >
::j_invariant() const
{
	debug_handler ("elliptic_curve_rep< T >",
		       "j_invariant() const");

	return (this->c4 * this->c4 * this->c4) / this->delta;
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_ELLIPTIC_CURVE_REP_CC_GUARD_
