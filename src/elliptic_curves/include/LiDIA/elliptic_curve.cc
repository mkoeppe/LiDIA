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


#ifndef LIDIA_ELLIPTIC_CURVE_CC_GUARD_
#define LIDIA_ELLIPTIC_CURVE_CC_GUARD_


#ifndef LIDIA_BASE_ELLIPTIC_CURVE_H_GUARD_
# include	"LiDIA/elliptic_curves/base_elliptic_curve.h"
#endif
#ifndef LIDIA_ELLIPTIC_CURVE_H_GUARD_
# include	"LiDIA/elliptic_curve.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// constructor / destructor
//
template< class T >
elliptic_curve< T >::elliptic_curve ()
	: base_elliptic_curve< T > ()
{
	debug_handler("elliptic_curve< T >",
		      "elliptic_curve()");
}



template< class T >
elliptic_curve< T >::elliptic_curve (const elliptic_curve< T > & E)
	: base_elliptic_curve< T > (E.e)
{
	debug_handler("elliptic_curve< T >",
		      "elliptic_curve(const elliptic_curve< T > &)");
}



template< class T >
elliptic_curve< T >::elliptic_curve (const T & x4, const T & x6,
				     elliptic_curve_flags::curve_model m)
	: base_elliptic_curve< T > (x4, x6, m)
{
	debug_handler("elliptic_curve< T >",
		      "elliptic_curve(const T&x4, const T&x6, "
		      "elliptic_curve_flags::curve_model m)");
}



template< class T >
elliptic_curve< T >::elliptic_curve (const T & x1, const T & x2,
				     const T & x3, const T & x4, const T & x6,
				     elliptic_curve_flags::curve_model m)
	: base_elliptic_curve< T > (x1, x2, x3, x4, x6, m)
{
	debug_handler("elliptic_curve< T >",
		      "elliptic_curve(const T&x1, const T&x2, const T&x3, "
		      "const T&x4, const T&x6, "
		      "elliptic_curve_flags::curve_model m)");
}



template< class T >
elliptic_curve< T >::~elliptic_curve ()
{
	debug_handler("elliptic_curve< T >",
		      "~elliptic_curve()");

	// Cleaning done by ~base_elliptic_curve < T > .
}



//
// assignment
//

template< class T >
elliptic_curve< T > &
elliptic_curve< T >::operator = (const elliptic_curve< T > & E)
{
	debug_handler ("elliptic_curve< T >",
		       "operator = (const elliptic_curve< T >");

	this->assign(E);
	return *this;
}



template< class T >
void
elliptic_curve< T >::assign (const elliptic_curve< T > & E)
{
	debug_handler ("elliptic_curve< T >",
		       "assign(const elliptic_curve< T >");

	if (this != &E) {
		if (E.e == NULL)
			lidia_error_handler("elliptic_curve< T >::"
					    "assign(const elliptic_curve< T > & E)",
					    "E.e == NULL");
		else if (this->e == NULL) {
			this->e = E.e;
			this->e->inc_ref_counter();
		}
		else if (this->e != E.e) {
			if (this->e->get_ref_counter() == 1)
				delete this->e;
			else
				this->e->dec_ref_counter();

			this->e = E.e;
			this->e->inc_ref_counter();
		}
	}
}



//
// invariants
//

template< class T >
T
elliptic_curve< T >::j_invariant () const
{
	debug_handler ("elliptic_curve< T >",
		       "j_invariant() const");

	if (this->e == NULL)
		lidia_error_handler("elliptic_curve< T >::"
				    "j_invariant() const",
				    "e == NULL");
	return this->e->j_invariant();
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_ELLIPTIC_CURVE_CC_GUARD_
