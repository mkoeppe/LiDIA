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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"
#include	"LiDIA/gf_element.h"
#include	"LiDIA/elliptic_curves/elliptic_curve_rep.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// constructor / destructor
//
elliptic_curve_rep< gf_element >
::elliptic_curve_rep()

	: base_elliptic_curve_rep< gf_element > ()
{
	debug_handler ("elliptic_curve_rep< gf_element >",
		       "elliptic_curve_rep()");
}



elliptic_curve_rep< gf_element >
::elliptic_curve_rep(const elliptic_curve_rep< gf_element > & e)

	: base_elliptic_curve_rep< gf_element > ()
{
	debug_handler ("elliptic_curve_rep< gf_element >",
		       "elliptic_curve_rep(const elliptic_curve_rep< gf_element > &)");

	this->assign(e);
}



elliptic_curve_rep< gf_element >
::elliptic_curve_rep(const gf_element & x4, const gf_element & x6,
		     elliptic_curve_flags::curve_model m)

	: base_elliptic_curve_rep< gf_element > ()
{
	debug_handler ("elliptic_curve_rep< gf_element >",
		       "elliptic_curve_rep< gf_element > ("
		       "const gf_element& x4, const gf_element& x6, "
		       "elliptic_curve_flags::curve_model)");

	this->set_coefficients(x4, x6, m);
}



elliptic_curve_rep< gf_element >
::elliptic_curve_rep(const gf_element & x1, const gf_element & x2,
		     const gf_element & x3, const gf_element & x4,
		     const gf_element & x6,
		     elliptic_curve_flags::curve_model m)

	: base_elliptic_curve_rep< gf_element > ()
{
	debug_handler ("elliptic_curve_rep< gf_element >",
		       "elliptic_curve_rep< gf_element > ("
		       "const gf_element&x1, const gf_element&x2, "
		       "const gf_element&x3, const gf_element&x4, "
		       "const gf_element&x6, "
		       "elliptic_curve_flags::curve_model)");

	this->set_coefficients(x1, x2, x3, x4, x6, m);
}



elliptic_curve_rep< gf_element >
::~elliptic_curve_rep()
{
	debug_handler ("elliptic_curve_rep< gf_element >",
		       "~elliptic_curve_rep()");
}



//
// assignment
//

elliptic_curve_rep< gf_element > &
elliptic_curve_rep< gf_element >
::operator = (const elliptic_curve_rep< gf_element > & e)
{
	debug_handler ("base_elliptic_curve< gf_element >",
		       "operator = (const elliptic_curve_rep< gf_element >");

	if (this != &e) {
		base_elliptic_curve_rep< gf_element >::assign(e);
	}
	return *this;
}



void elliptic_curve_rep< gf_element >
::assign (const elliptic_curve_rep< gf_element > & e)
{
	debug_handler ("base_elliptic_curve< gf_element >",
		       "assign(const elliptic_curve_rep< gf_element >");

	if (this != &e) {
		base_elliptic_curve_rep< gf_element >::assign(e);
	}
}



//
// invariants
//

gf_element elliptic_curve_rep< gf_element >
::j_invariant() const
{
	debug_handler ("elliptic_curve_rep< gf_element >",
		       "j_invariant() const");

	return (c4*c4*c4)/delta;
}



//
// properties
//

unsigned int elliptic_curve_rep< gf_element >
::degree_of_definition() const
{
	debug_handler ("elliptic_curve_rep< gf_element >",
		       "degree_of_definition() const");

	unsigned int d = 1, h;

	if (!a1.is_zero()) {
		h = a1.relative_degree();
		d = d*(h/gcd(static_cast<long>(d), static_cast<long>(h)));
	}
	if (!a2.is_zero()) {
		h = a2.relative_degree();
		d = d*(h/gcd(static_cast<long>(d), static_cast<long>(h)));
	}
	if (!a3.is_zero()) {
		h = a3.relative_degree();
		d = d*(h/gcd(static_cast<long>(d), static_cast<long>(h)));
	}
	if (!a4.is_zero()) {
		h = a4.relative_degree();
		d = d*(h/gcd(static_cast<long>(d), static_cast<long>(h)));
	}
	if (!a6.is_zero()) {
		h = a6.relative_degree();
		d = d*(h/gcd(static_cast<long>(d), static_cast<long>(h)));
	}
	return d;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
