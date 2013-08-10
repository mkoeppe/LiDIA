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
//	Author	: Andrea Rau, Robert Carls
//	Changes	: See CVS log
//
//==============================================================================================



#include	"LiDIA/mv_term.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//---------------------------------------------------------------------
// checks whether a mv_term is a linear term, i.e. exactly one variable
// is used or var is of the form 2^i

bool mv_term::is_linear() const
{
	if (coeff.is_zero() || var.is_zero())
		return false;

	if (var.is_one())
		return true;

	bool found_first = false;
	int i;

	for (i = 0; i < var.bit_length(); i++) {
		if (var.bit(i))
			if (found_first == false)
				found_first = true;
			else
				return false;
	}
	return true;
}



bool mv_term::is_linear(lidia_size_t & k) const
{
	if (coeff.is_zero() || var.is_zero()) {
		k = -1;
		return false;
	}

	if (var.is_one()) {
		k = 0;
		return true;
	}

	int i;
	k = -1;

	for (i = 0; i < var.bit_length(); i++) {
		if (var.bit(i))
			if (k == -1)
				k = i;
			else {
				k = -1;
				return false;
			}
	}
	return true;
}



//
// functions for manipulation of terms
//

//---------------------------------------------------------------------
// in term t the Variable X_i is substituted by term s

void substitute(mv_term & t, const mv_term & s, unsigned int i)
{
	if (t.var.bit(i) == 0)
		return;

	bigint t_var(t.get_var());
	gf2n t_coeff(t.get_coeff());
	bigint ivar(bigint(1) << i);

	bitwise_xor(t_var, t_var, ivar); // eliminate coefficient i from term
	bitwise_or(t_var, t_var, s.get_var()); // substitute term s in t (as multiply)
	multiply(t_coeff, t_coeff, s.get_coeff());
	t.assign(t_coeff, t_var);
}



//---------------------------------------------------------------------
// in term t the Variable X_i which have a bit 0 in 'mask' are substituted
// by the corresponding bits in c, all others are not changed (note that
// mask is an inverted mask).

void substitute(mv_term & t, const bigint & cc, const bigint & bits)
{
	bigint t_var, ivar, c;

	bitwise_and(t_var, t.var, bits);

	if (t_var.is_zero())
		return;

	bitwise_and(c, cc, bits);
	bitwise_xor(ivar, t_var, c);
	bitwise_and(ivar, ivar, t_var);

	if (!ivar.is_zero())
		t.assign_zero();
	else {
		shift_left(ivar, bigint(1), t.var.bit_length());
		dec(ivar);
		bitwise_xor(ivar, ivar, c);
		bitwise_and(t.var, t.var, ivar);
	}
}



//
// output
//

std::ostream & operator << (std::ostream & out, const mv_term & a)
{
	gf2n coeff_a(a.get_coeff());
	bigint var_a(a.get_var());

	if (a.is_const() || a.is_zero())
		return (out << coeff_a);
	else {
		out << coeff_a << " *";

		for (unsigned int i = 0; i < static_cast<unsigned int>(var_a.bit_length()); i++) {
			if (var_a.bit(i))
				out << " X_" << i;
		}
		return out;
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
