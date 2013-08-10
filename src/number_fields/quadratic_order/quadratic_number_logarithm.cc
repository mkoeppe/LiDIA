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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/arith.inl"
#include	"LiDIA/xbigfloat.h"
#include	"LiDIA/quadratic_ideal.h"
#include	"LiDIA/quadratic_order.h"
#include	"LiDIA/quadratic_number_logarithm.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
//  constructors / destructor
//

quadratic_number_logarithm
::quadratic_number_logarithm ()
{
	debug_handler ("quadratic_number_logarithm",
		       "quadratic_number_logarithm()");
	t = -1;
}



quadratic_number_logarithm
::quadratic_number_logarithm (
	const quadratic_number_logarithm & L)
{
	debug_handler ("quadratic_number_logarithm",
		       "quadratic_number_logarithm(L)");

	this->assign(L);
}



quadratic_number_logarithm
::~quadratic_number_logarithm ()
{
	debug_handler ("quadratic_number_logarithm",
		       "~quadratic_number_logarithm()");
}



//
//  access
//

bigint quadratic_number_logarithm
::get_discriminant () const
{
	debug_handler ("quadratic_number_logarithm",
		       "get_discriminant()const");

	if (t == -1)
		lidia_error_handler ("quadratic_number_logarithm::get_discriminant()",
				     "Not initialized.");

	return A.which_order().discriminant();
}



const quadratic_order & quadratic_number_logarithm
::which_order () const
{
	debug_handler ("quadratic_number_logarithm",
		       "which_order()const");

	if (t == -1)
		lidia_error_handler ("quadratic_number_logarithm::which_order()",
				     "Not initialized.");

	return A.which_order();
}



const quadratic_order & quadratic_number_logarithm
::get_order () const
{
	debug_handler ("quadratic_number_logarithm",
		       "get_order()const");

	if (t == -1)
		lidia_error_handler ("quadratic_number_logarithm::get_order()",
				     "Not initialized.");

	return A.get_order();
}



const quadratic_ideal & quadratic_number_logarithm
::get_A () const
{
	debug_handler ("quadratic_number_logarithm",
		       "get_A()const");

	if (t == -1)
		lidia_error_handler ("quadratic_number_logarithm::get_A()",
				     "Not initialized.");

	return A;
}



const quadratic_ideal & quadratic_number_logarithm
::get_B () const
{
	debug_handler ("quadratic_number_logarithm",
		       "get_B()const");

	if (t == -1)
		lidia_error_handler ("quadratic_number_logarithm::get_B()",
				     "Not initialized.");

	return B;
}



const xbigfloat & quadratic_number_logarithm
::get_log_approximation () const
{
	debug_handler ("quadratic_number_logarithm",
		       "get_log_approximation()const");

	if (t == -1)
		lidia_error_handler ("quadratic_number_logarithm::"
				     "get_log_approximation()",
				     "Not initialized.");

	return a;
}



long quadratic_number_logarithm
::get_log_accuracy () const
{
	debug_handler ("quadratic_number_logarithm",
		       "get_log_accuracy()const");

	if (t == -1)
		lidia_error_handler ("quadratic_number_logarithm::"
				     "get_log_accuracy()",
				     "Not initialized.");

	return k;
}



int quadratic_number_logarithm
::get_log_type () const
{
	debug_handler ("quadratic_number_logarithm",
		       "get_log_type()const");

	if (t == -1)
		lidia_error_handler ("quadratic_number_logarithm::"
				     "get_log_type()",
				     "Not initialized.");

	return t;
}



int quadratic_number_logarithm
::get_sign () const
{
	debug_handler ("quadratic_number_logarithm",
		       "get_sign()const");

	if (t == -1)
		lidia_error_handler ("quadratic_number_logarithm::"
				     "get_sign()",
				     "Not initialized.");

	return s;
}



//
//  assignments
//

quadratic_number_logarithm &
quadratic_number_logarithm::operator = (const quadratic_number_logarithm & x)
{
	debug_handler ("quadratic_number_logarithm",
		       "operator = (L)");
	this->assign(x);
	return *this;
}



void quadratic_number_logarithm
::assign (const quadratic_number_logarithm & x)
{
	debug_handler ("quadratic_number_logarithm",
		       "assign(L)");
	if (this != &x) {
		A.assign(x.A);
		B.assign(x.B);
		a.assign(x.a);
		k = x.k;
		s = x.s;
		t = x.t;
	}
}



void quadratic_number_logarithm
::assign (const quadratic_ideal & pA,
	  const quadratic_ideal & pB,
	  const xbigfloat       & pa,
	  long                    pk,
	  int                     ps,
	  int                     pt)
{
	debug_handler ("quadratic_number_logarithm",
		       "assign(A, B, a, k, s, t)");

	if (pt == 0 || pt == 1) {
		if (pk >= 4) {
			A.assign(pA);
			B.assign(pB);
			a.assign(pa);
			k = pk;
			s = ps;
			t = pt;
		}
		else
			lidia_error_handler("quadratic_number_logarithm::"
					    "assign(A, B, a, k, s, t)",
					    "k < 4");
	}
	else
		lidia_error_handler("quadratic_number_logarithm::"
				    "assign(A, B, a, k, s, t)",
				    "Invalid value for t");
}



void swap (quadratic_number_logarithm & x,
	   quadratic_number_logarithm & y)
{
	debug_handler ("quadratic_number_logarithm",
		       "swap(x, y)");

	if (&x != &y) {
		swap(x.A, y.A);
		swap(x.B, y.B);
		swap(x.a, y.a);
		swap(x.k, y.k);
		swap(x.s, y.s);
		swap(x.t, y.t);
	}
}



//
//  arithmetic
//

void quadratic_number_logarithm
::invert ()
{
	debug_handler ("quadratic_number_logarithm",
		       "invert()");

	if (A.is_reduced()) {
		swap(A, B);
		a.negate();
	}
	else {
		// A/alpha = B, so B/(1/alpha) = A.
		// Reduce A by gamma, then
		//
		// (B/gamma) / (1/alpha) = A/gamma,
		// and 1/alpha in Min(B/gamma).
		//
		quadratic_number_standard gamma;
		A.reduce(gamma);
		divide(B, B, gamma);
		a.negate();
	}
}



//
//  input / output
//

std::istream &
operator >> (std::istream & in, quadratic_number_logarithm & x)
{
	debug_handler ("quadratic_number_logarithm",
		       "operator >> (std::istream, x)");

	std::cout << "quadratic_number_logarithm::operator >> not yet implemented.";
	std::cout << std::endl;
	return in;
}



std::ostream &
operator << (std::ostream & out, const quadratic_number_logarithm & x)
{
	debug_handler ("quadratic_number_logarithm",
		       "operator << (std::ostream, x)");

	std::cout << "(" << x.A << ", " << x.B << ", " << x.a << ", " << x.k;
	std::cout << ", " << x.s << ", " << x.t << ")";
	return out;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
