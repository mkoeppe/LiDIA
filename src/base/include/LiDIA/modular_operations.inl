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
//	$Id: modular_operations.inl,v 2.3 2001/06/20 13:04:28 hamdy Exp $
//
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef __LIDIA_MODULAR_OPERATIONS_INL
#define __LIDIA_MODULAR_OPERATIONS_INL



#ifndef __LIDIA_INTERFACE_LIB_H
# include	<LiDIA/base/interface_lib.h>
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// base operations
//

inline void
pos_div_rem(bigint & q, bigint & r, const bigint & a, const bigint & b)
{
	debug_handler_l("modular_operations", "in inline pos_div_rem()", LDBL_MATRIX);

	div_rem(q, r, a, b);

	if (r.is_lt_zero())
		if (b.is_lt_zero()) {
			subtract(r, r, b);
			inc(q);
		}
		else {
			add(r, r, b);
			dec(q);
		}
}



inline void
pos_div_rem(long & q, long & r, long a, long b)
{
	debug_handler_l("modular_operations", "in inline pos_div_rem()", LDBL_MATRIX);

	q = a / b;
	r = a - q * b;

	if (r < 0)
		if (b < 0) {
			r -= b;
			q++;
		}
		else {
			r += b;
			q--;
		}
}



inline void
pos_div_rem(short int & q, short int & r, short int a, short int b)
{
	debug_handler_l("modular_operations", "in inline pos_div_rem()", LDBL_MATRIX);

	q = a / b;
	r = a - q * b;

	if (r < 0)
		if (b < 0) {
			r -= b;
			q++;
		}
		else {
			r += b;
			q--;
		}
}



inline void
best_remainder(bigint & a, const bigint & b, const bigint & mod)
{
	debug_handler_l("modular_operations", "in inline best_remainder()",
			LDBL_MATRIX);

	if (mod.is_one() || b.is_zero())
		a.assign_zero();
	else {
		bigint mod2;
		shift_right(mod2, mod, 1);
		if (b > -mod2 && b <= mod2)
			a = b;
		else
			remainder(a, b, mod);
		if (a <= -mod2)
			add(a, a, mod);
		if (a > mod2)
			subtract(a, a, mod);
	}
}



inline void
best_remainder(long & a, const bigint & b, long mod)
{
	debug_handler_l("modular_operations", "in inline best_remainder()", LDBL_MATRIX);

	if (mod == 1 || b.is_zero())
		a = 0;
	else {
		register long mod2 = mod >> 1;
		if (b > bigint(-mod2) && b <= bigint(mod2))
			b.longify(a);
		else
			remainder(a, b, mod);
		if (a <= -mod2)
			a += mod;
		if (a > mod2)
			a -= mod;
	}
}



inline void
add_mod(bigint & c, const bigint & a, const bigint & b, const bigint & mod)
{
	debug_handler_l("modular_operations", "in inline add_mod()", LDBL_MATRIX);

	if (mod.is_one())
		c.assign_zero();
	else {
		bigint mod2;
		shift_right(mod2, mod, 1);
		add(c, a, b);
		if (c <= -mod2)
			add(c, c, mod);
		if (c > mod2)
			subtract(c, c, mod);
	}
}



inline void
add_mod(long & d, long a, long b, long mod)
{
	debug_handler_l("modular_operations", "in inline add_mod()", LDBL_MATRIX);

	if (mod == 1)
		d = 0;
	else {
		register long mod2 = mod/2;
		double c = static_cast<double>(a) + static_cast<double>(b);
		if (c <= static_cast<double>(-mod2))
			c += static_cast<double>(mod);
		if (c > static_cast<double>(mod2))
			c -= static_cast<double>(mod);
		d = static_cast<long>(c);
	}
}



inline void
sub_mod(bigint & c, const bigint & a, const bigint & b, const bigint & mod)
{
	debug_handler("modular_operations", "in inline sub_mod()");

	if (mod.is_one())
		c.assign_zero();
	else {
		bigint mod2;
		shift_right(mod2, mod, 1);
		subtract(c, a, b);
		if (c <= -mod2)
			add(c, c, mod);
		if (c > mod2)
			subtract(c, c, mod);
	}
}



inline void
sub_mod(long & d, long a, long b, long mod)
{
	debug_handler_l("modular_operations", "in inline sub_mod()", LDBL_MATRIX);

	if (mod == 1)
		d = 0;
	else {
		register long mod2 = mod/2;
		double c = static_cast<double>(a) - static_cast<double>(b);
		if (c <= static_cast<double>(-mod2))
			c += static_cast<double>(mod);
		if (c > static_cast<double>(mod2))
			c -= static_cast<double>(mod);
		d = static_cast<long>(c);
	}
}



inline void
mult_mod(bigint & c, const bigint & a, const bigint & b, const bigint & mod)
{
	debug_handler_l("modular_operations", "in inline mult_mod()", LDBL_MATRIX);

	if (mod.is_one())
		c.assign_zero();
	else {
		multiply(c, a, b);
		best_remainder(c, c, mod);
	}
}



inline void
mult_mod(long & d, long a, long b, long mod)
{
	debug_handler_l("modular_operations", "in inline mult_mod()", LDBL_MATRIX);

	if (mod == 1)
		d = 0;
	else {
		register long mod2 = mod/2;
		double ab = static_cast<double>(a) * static_cast<double>(b);
		register long q = static_cast<long>(ab / static_cast<double>(mod));
		register long res = static_cast<long>(ab - (static_cast<double>(q) * static_cast<double>(mod)));
		if (res > mod2)
			res -= mod;
		if (res <= -mod2)
			res += mod;
		d = res;
	}
}



inline void
div_mod(bigint & c, const bigint & a, const bigint & b, const bigint & mod)
{
	debug_handler_l("modular_operations", "in inline div_mod()", LDBL_MATRIX);

	bigint u;
	bigint d = xgcd_left(u, b, mod);
	if (!d.is_one())
		lidia_error_handler("modular_operations", "div_mod - Version bigint :: "
				    "Inverse undefined");
	mult_mod(c, a, u, mod);
}



inline void
div_mod(long & c, long a, long b, long mod)
{
	debug_handler_l("modular_operations", "in inline div_mod()", LDBL_MATRIX);

	long u, v;
	register long d = xgcd(u, v, b, mod);
	if (d != 1)
		lidia_error_handler("modular_operations", "div_mod - Version long :: "
				    "Inverse undefined");
	mult_mod(c, a, u, mod);
}



inline void
inv_mod(bigint & c, const bigint & a, const bigint & mod)
{
	debug_handler_l("modular_operations", "in inline inv_mod()", LDBL_MATRIX);

	bigint d, mod2;
	shift_right(mod2, mod, 1);
	d = xgcd_left(c, a, mod);

	if (!d.is_one())
		lidia_error_handler("modular_operations", "inv_mod - Version bigint :: "
				    "Inverse undefined");
	if (c <= -mod2)
		add(c, c, mod);
	else
		if (c > mod2)
			subtract(c, c, mod);
}



inline void
inv_mod(long & e, long a, long mod)
{
	debug_handler_l("modular_operations", "in inline inv_mod()", LDBL_MATRIX);

	long t, mod2 = mod/2;
	register long d = xgcd(e, t, a, mod);
	if (d != 1)
		lidia_error_handler("modular_operations", "inv_mod - Version long :: "
				    "Inverse undefined");
	if (e <= -mod2)
		e += mod;
	else
		if (e > mod2)
			e -= mod;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// __LIDIA_MODULAR_OPERATIONS_INL
