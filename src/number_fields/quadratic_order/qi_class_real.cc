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
//	Author	: Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/qi_class_real.h"
#include	"LiDIA/number_fields/qo_list.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



int qi_class_real::info = 0;
int qi_class_real::do_verify = 0;



bigint
qi_class_real_key (const qi_class_real & G)
{
	return G.get_a();
}



inline xbigfloat
approx_Ln_reduce (bigint & a, bigint & b, bigint & Delta, long k)
{
	long b1, m, t;
	xbigfloat n, a1, a2, a3, sigma_n, l;

	// b1 = ceil(b(Delta)/2) + b(b)
	b1 = b_value(Delta);
	if (b1&1)
		b1 = (b1 >> 1) + 1;
	else
		b1 = b1 >> 1;
	b1 += b_value(b);

	// m = b1+max(b(a),b1)+1;
	m = b_value(a);
	if (m > b1)
		m += b1+1;
	else
		m = 2*b1+1;

	// abs. k+7 approximations to a+b*sqrt(D) and a-b*sqrt(D)
	//
	t = k+6;

	// relative t+m+4 approx. to b
	truncate(a1, b, t+m+5);

	// relative t+m+4 approx. to sqrt(Delta)
	sqrt(a2, Delta, t+m+4);

	// relative t+m+2 approx. to b*sqrt(Delta)
	multiply(a3, a1, a2);

	// relative t+m+1 approx. to b*sqrt(Delta)
	a3.truncate(t+m+4);

	// relative t+1 approx. to a+b*sqrt(Delta)
	add(n, a, a3);

	// relative t+1 approx. to a-b*sqrt(Delta)
	subtract(sigma_n, a, a3);


	// relative k+3 approx. to this/sigma(this)
	//

	// relative k+4 approx. to this/sigma(this)
	divide(n, n, sigma_n, k+6);

	// relative k+3 approx. to this/sigma(this)
	n.truncate(k+6);
	n.absolute_value();

	// absolute k+1 approx. to 1/2 ln|this/sigma(this)|
	log(l, n, k+2);
	shift_right(l, l, 1);

	// absolute k approx. to Ln(this)
	l.truncate(k+1+b_value(l));

	return l;
}



inline xbigfloat
approx_Ln (bigint & b, xbigfloat & rd, long xprec)
{
	xbigfloat dmod, dmod2;

	add(dmod, b, rd);
	subtract(dmod2, b, rd);

	divide(dmod, dmod, dmod2, xprec+6);
	dmod.truncate(xprec+6);

	dmod.absolute_value();
	log(dmod, dmod, xprec+2);
	shift_right(dmod, dmod, 1);

	dmod.truncate(xprec+1+b_value(dmod));

	return dmod;
}



//
// qi_class_real::normalize_real()
//
// Task:
//      normalizes the real ideal.
//

inline void
qi_class_real::normalize_real ()
{
	debug_handler("qi_class_real", "normalize_real");

	bigint a2, q, r, temp;
	bool is_normal;

	shift_left(a2, a, 1);
	if (a.compare(rootD) <= 0) {
		subtract(temp, rootD, a2);
		is_normal = ((temp.compare(b) < 0) && (b.compare(rootD) <= 0));
	}
	else
		is_normal = ((b.compare(-a) > 0) && (b.compare(a) <= 0));

	if (!is_normal) {
		if (a.compare(rootD) <= 0) {
			// q = floor((rootD - b) / 2a)
			subtract(temp, rootD, b);
		}
		else {
			// q = floor((a - b) / 2a)
			subtract(temp, a, b);
		}
		div_rem(q, r, temp, a2);
		if ((temp.is_lt_zero()) && (!r.is_zero()))
			dec(q);

		multiply(temp, a2, q);
		add(b, b, temp);
	}
}



//
// qi_class_real::reduce()
//
// Task:
//      reduces the ideal and updates its distance.
//

inline void
qi_class_real::reduce ()
{
	debug_handler("qi_class_real", "reduce");

	bigint q, r, a2, na, oa, nb, OB, BB, NB, temp;
	int s;
	xbigfloat dmod, temp2;

	normalize_real();

	shift_left(a2, abs(a), 1);
	subtract(temp, rootD, a2);
	if (temp.is_lt_zero())
		inc(temp);

	if (!((abs(temp) < b) && (b.compare(rootD) <= 0))) {
		// oa = (D - b^2) / 4a
		square(oa, b);
		subtract(oa, Delta, oa);
		divide(oa, oa, a);
		shift_right(oa, oa, 2);

		OB.assign_one();
		BB.assign_zero();

		while (!((abs(temp) < b) && (b.compare(rootD) <= 0))) {
			s = a.is_gt_zero() ? 0 : 1;

			// (rootD-b) = (2a)*q + r
			shift_left(a2, a, 1);
			add(temp, rootD, b);
			if (s)
				inc(temp);
			div_rem(q, r, temp, a2);

			multiply(temp, q, BB);
			add(NB, temp, OB);
			OB.assign(BB);
			BB.assign(NB);

			// nb = rootD + s - r;
			subtract(nb, rootD, r);
			if (s)
				inc(nb);

			// na = oa -q*((nb - b) >> 1);
			subtract(temp, nb, b);
			temp.divide_by_2();
			multiply(temp, temp, q);
			subtract(na, oa, temp);

			b.assign(nb);
			oa.assign(a);
			a.assign(na);

			shift_left(a2, abs(a), 1);
			subtract(temp, rootD, a2);
			if (temp.is_lt_zero())
				inc(temp);
		}

		multiply(a2, OB, a);
		a2.multiply_by_2();
		multiply(temp, BB, b);
		add(temp, temp, a2);

		// d += ln((2*OB*a + BB*b + BB*rd) / 2)
		add(d, d, approx_Ln_reduce(temp, BB, Delta, xprec));

		a.absolute_value();
	}
}



//
// constructor
//

qi_class_real::qi_class_real ()
{
	debug_handler("qi_class_real", "qi_class_real()");

	d.assign_zero();
}



//
// constructor
//
// Task:
//    set to the ideal belonging to the current order given by
//
//    [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    If the parameters do not form an ideal of the current quadratic order,
//    or the current quadratic order is zero, the ideal will be set to the
//    zero ideal.
//

qi_class_real::qi_class_real (const bigint & a2, const bigint & b2,
			      const xbigfloat & dist) : qi_class(a2, b2)
{
	debug_handler("qi_class_real", "qi_class_real(bigint, bigint, bigfloat)");

	if (Delta.is_le_zero()) {
		a.assign_zero();
		b.assign_zero();
	}
	d.assign(dist);
}



//
// constructor
//
// Task:
//    set to the ideal belonging to the current order given by
//
//    [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    If the parameters do not form an ideal of the current quadratic order,
//    or the current quadratic order is zero, the ideal will be set to the
//    zero ideal.
//

qi_class_real::qi_class_real (const long a2, const long b2,
			      const xbigfloat & dist) : qi_class(a2, b2)
{
	debug_handler("qi_class_real", "qi_class_real(long, long, bigfloat)");

	if (Delta.is_le_zero()) {
		a.assign_zero();
		b.assign_zero();
	}
	d.assign(dist);
}



//
// constructor
//
// Task:
//    set to the ideal belonging to the current order given by
//
//    [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    If the parameters do not form an ideal of the current quadratic order,
//    or the current quadratic order is zero, the ideal will be set to the
//    zero ideal.  The distance is set to zero.
//

qi_class_real::qi_class_real (const bigint & a2,
			      const bigint & b2) : qi_class(a2, b2)
{
	debug_handler("qi_class_real", "qi_class_real(bigint, bigint)");

	if (Delta.is_le_zero()) {
		a.assign_zero();
		b.assign_zero();
	}
	d.assign_zero();
}



//
// constructor
//
// Task:
//    set to the ideal belonging to the current order given by
//
//    [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    If the parameters do not form an ideal of the current quadratic order,
//    or the current quadratic order is zero, the ideal will be set to the
//    zero ideal.  The distance is set to zero.
//

qi_class_real::qi_class_real (const long a2, const long b2) : qi_class(a2, b2)
{
	debug_handler("qi_class_real", "qi_class_real(long, long)");

	if (Delta.is_le_zero()) {
		a.assign_zero();
		b.assign_zero();
	}
	d.assign_zero();
}



//
// constructor:
//    if qf is indefinite and primitive, initialize with qf reduced.
//

qi_class_real::qi_class_real (const quadratic_form & qf)
{
	debug_handler("qi_class_real", "qi_class_real(quadratic_form)");

	if ((!qf.is_indefinite()) || (!qf.is_primitive())) {
		a.assign_zero();
		b.assign_zero();
		d.assign_zero();
	}
	else {
		set_current_order(qf.which_order());
		a.assign(abs(qf.get_a()));
		b.assign(qf.get_b());
		d.assign_zero();
		reduce();
	}
}



//
// constructor:
//    if A is invertible, initialize with A reduced.
//

qi_class_real::qi_class_real (const quadratic_ideal & A)
{
	debug_handler("qi_class_real", "qi_class_real(quadratic_ideal)");

	if (!current_order)
		current_order = quadratic_order::qo_l().add_last_to_list(current_order);

	if ((!A.is_invertible()) || (A.discriminant().is_le_zero())) {
		a.assign_zero();
		b.assign_zero();
		d.assign_zero();
	}
	else {
		set_current_order((quadratic_order&)A.which_order());
		a.assign(A.get_a());
		b.assign(A.get_b());
		d.assign_zero();
		reduce();
	}
}



//
// constructor:
//    initialize with a copy of A.
//

qi_class_real::qi_class_real (const qi_class & A) : qi_class(A)
{
	debug_handler("qi_class_real", "qi_class_real(qi_class)");

	if (A.discriminant().is_le_zero()) {
		a.assign_zero();
		b.assign_zero();
	}

	d.assign_zero();
}



//
// constructor:
//    initialize with a copy of A.
//

qi_class_real::qi_class_real (const qi_class_real & A)
{
	debug_handler("qi_class_real", "qi_class_real(qi_class_real)");

	if (!current_order)
		current_order = quadratic_order::qo_l().add_last_to_list(current_order);

	a.assign(A.a);
	b.assign(A.b);
	d.assign(A.d);
}



//
// destructor
//

qi_class_real::~qi_class_real()
{
	debug_handler("qi_class_real", "~qi_class_real");
}



//
// qi_class_real::verbose()
//
// Task:
//      sets the verbosity of commands.  Currently, the following levels are
//      supported:
//         0 - nothing
//         > 0 - run times
//

void
qi_class_real::verbose (int state)
{
	debug_handler("qi_class_real", "verbose");

	if (state <= 0)
		info = 0;
	else
		info = state;
}



//
// qi_class_real::verification()
//
// Task:
//      sets the level of verifications.  Currently, the following levels are
//      supported:
//         0 - nothing
//         1 - computation is verified (true or false)
//         > 1 - computation is verified and intermediate data is output
//

void
qi_class_real::verification (int level)
{
	debug_handler("qi_class_real", "verification");

	if (level <= 0)
		do_verify = 0;
	else
		do_verify = level;
}



//
// qi_class_real::assign_zero()
//
// Task:
//      set to the zero ideal
//

void
qi_class_real::assign_zero ()
{
	debug_handler("qi_class_real", "assign_zero");

	a.assign_zero();
	b.assign_zero();
	d.assign_zero();
}



//
// qi_class_real::assign_one()
//
// Task:
//      set to the unit ideal of the current quadratic_order
//

void
qi_class_real::assign_one ()
{
	debug_handler("qi_class_real", "assign_one");

	if (Delta.is_le_zero()) {
		lidia_error_handler("qi_class_real", "assign_one - current quadratic "
				    "order is not defined or is not real");
		return;
	}

	a.assign_one();
	b.assign(Delta);
	normalize_real();
	d.assign_zero();
}



//
// qi_class_real::assign_principal()
//
// Task:
//    set to the reduced principal ideal generated by the element
//
//    (x + (D + sqrt(Delta))/2 y)
//

void
qi_class_real::assign_principal (const bigint & x, const bigint & y)
{
	debug_handler("qi_class_real", "assign_principal()");

	bigint x2, y2, n, m, k, l, temp, temp2;

	if (Delta.is_le_zero()) {
		lidia_error_handler("qi_class_real", "assign_principal - current quadratic "
				    "order is not defined or is not real");
		return;
	}

	// use (x + y \sqrt(Delta))/2 form
	y2.assign(y);
	shift_left(x2, x, 1);
	if (Delta.is_odd())
		add(x2, x2, y);

	// compute norm of (x2 + y2 sqrt(Delta))/2
	square(n, x2);
	square(temp, y2);
	multiply(temp, temp, Delta);
	subtract(n, n, temp);
	shift_right(n, n, 2);
	n.absolute_value();

	// solve m = k y2 + l (x2 + y2 Delta) / 2
	multiply(temp2, y2, Delta);
	add(temp2, temp2, x2);
	temp2.divide_by_2();
	m.assign(xgcd(k, l, y2, temp2));

	// a = abs(n / m^2)
	square(temp, m);
	divide(a, n, temp);

	// (b = kx2 + l(x2+y2)D/2 ) / m
	add(temp2, x2, y2);
	multiply(temp2, temp2, l);
	multiply(b, temp2, Delta);
	b.divide_by_2();
	multiply(temp, k, x2);
	add(b, b, temp);
	divide(b, b, m);
	shift_left(temp, a, 1);
	remainder(b, b, temp);

	d.assign_zero();
	reduce();
}



//
// qi_class_real::assign(bigint, bigint)
//
// Task:
//    set to the ideal belonging to the current order given by
//
//    q2 [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    If the parameters do not form an ideal, false will be returned and
//    the ideal will not be modified.
//

bool
qi_class_real::assign (const bigint & a2, const bigint & b2)
{
	debug_handler("qi_class_real", "assign(bigint, bigint)");

	bigint c, temp;

	d.assign_zero();

	if (Delta.is_le_zero()) {
		lidia_error_handler("qi_class_real", "assign(bigint, bigint) - current "
				    "quadratic order is not defined or is not real");
		return false;
	}

	if (a2.is_zero())
		if (b2.is_zero()) {
			a.assign_zero();
			b.assign_zero();
			return true;
		}
		else
			return false;
	else {
		square(c, b2);
		subtract(c, c, Delta);
		shift_left(temp, a2, 2);
		remainder(c, c, temp);
		if (c.is_zero()) {
			a.assign(a2);
			b.assign(b2);
			reduce();

			if (is_invertible())
				return true;
		}

		a.assign_zero();
		b.assign_zero();
		d.assign_zero();
		return false;
	}
}



//
// qi_class_real::assign(long, long)
//
// Task:
//    set to the ideal belonging to the current order given by
//
//    q2 [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    If the parameters do not form an ideal, false will be returned and
//    the ideal will not be modified.
//

bool
qi_class_real::assign (const long a2, const long b2)
{
	debug_handler("qi_class_real", "assign(long, long)");

	return assign(bigint(a2), bigint(b2));
}



//
// qi_class_real::assign(quadratic_form)
//
// Task:
//      set to the reduction of qf.  If qf is not indefinite or not primitive,
//      the error handler will be evoked.
//

void
qi_class_real::assign (const quadratic_form & qf)
{
	debug_handler("qi_class_real", "assign(quadratic_form)");

	if (!qf.is_indefinite()) {
		lidia_error_handler("qi_class_real", "assign(quaddratic_form) - the "
				    "quadratic form must be indefinite");
		return;
	}

	if (!qf.is_primitive()) {
		lidia_error_handler("qi_class_real", "assign(quadratic_form) - the "
				    "quadratic form is not primitive");
		return;
	}

	set_current_order(qf.which_order());
	a.assign(abs(qf.get_a()));
	b.assign(qf.get_b());
	d.assign_zero();
	reduce();
}



//
// qi_class_real::assign(quadratic_ideal)
//
// Task:
//      set to the reduction of B.  If B is not invertible, the error handler
//      will be evoked.
//

void
qi_class_real::assign (const quadratic_ideal & B)
{
	debug_handler("qi_class_real", "assign(quadratic_ideal)");

	if (B.discriminant().is_le_zero()) {
		lidia_error_handler("qi_class_real", "assign(quadratic_ideal) - the "
				    "quadratic ideal must belong to a real order");
		return;
	}

	if (!B.is_invertible()) {
		lidia_error_handler("qi_class_real", "assign(quadratic_ideal) - the "
				    "quadratic ideal is not invertible");
		return;
	}

	set_current_order((quadratic_order&)B.which_order());
	a.assign(B.get_a());
	b.assign(B.get_b());
	d.assign_zero();
	reduce();
}



//
// qi_class_real::assign(qi_class)
//
// Task:
//      set to a copy of B.
//

void
qi_class_real::assign (const qi_class & B)
{
	debug_handler("qi_class_real", "assign(qi_class)");

	if (B.discriminant().is_le_zero()) {
		lidia_error_handler("qi_class_real", "assign(qi_class) - the qi_class must "
				    "belong to a real order");
		return;
	}

	a.assign(B.get_a());
	b.assign(B.get_b());
	d.assign_zero();
}



//
// qi_class_real::assign_principal(bigfloat)
//
// Task:
//    set to the reduced principal ideal generated by the element
//
//    (x + (D + sqrt(Delta))/2 y)
//
//    and initializes the distance to dist.
//

void
qi_class_real::assign_principal (const bigint & x, const bigint & y,
				 const xbigfloat & dist)
{
	debug_handler("qi_class_real", "assign_principal(bigfloat)");

	assign_principal(x, y);
	d.assign(dist);
}



//
// qi_class_real::assign(bigint, bigint, xbigfloat)
//
// Task:
//    set to the ideal belonging to the current order given by
//
//    q2 [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    If the parameters do not form an ideal, false will be returned and
//    the ideal will not be modified.  The distance is initialized to dist.
//

bool
qi_class_real::assign (const bigint & a2, const bigint & b2,
		       const xbigfloat & dist)
{
	debug_handler("qi_class_real", "assign(bigint, bigint, xbigfloat)");

	bool OK;

	OK = assign(a2, b2);
	if (OK)
		d.assign(dist);

	return OK;
}



//
// qi_class_real::assign(long, long, xbigfloat)
//
// Task:
//    set to the ideal belonging to the current order given by
//
//    q2 [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    If the parameters do not form an ideal, false will be returned and
//    the ideal will not be modified.  The distance is initialized to dist.
//

bool
qi_class_real::assign (const long a2, const long b2, const xbigfloat & dist)
{
	debug_handler("qi_class_real", "assign(long, long, xbigfloat)");

	bool OK;

	OK = assign(a2, b2);
	if (OK)
		d.assign(dist);

	return OK;
}



//
// qi_class_real::assign(quadratic_form, xbigfloat)
//
// Task:
//      set to the reduction of qf.  If qf is not indefinite or not primitive,
//      the error handler will be evoked.  The distance is initialized to dist.
//

void
qi_class_real::assign (const quadratic_form & qf, const xbigfloat & dist)
{
	debug_handler("qi_class_real", "assign(quadratic_form, xbigfloat)");

	assign(qf);
	d.assign(dist);
}



//
// qi_class_real::assign(quadratic_ideal, xbigfloat)
//
// Task:
//      set to the reduction of B.  If B is not invertible, the error handler
//      will be evoked.  The distance is initialized to dist.
//

void
qi_class_real::assign (const quadratic_ideal & B, const xbigfloat & dist)
{
	debug_handler("qi_class_real", "assign(quadratic_ideal, xbigfloat)");

	assign(B);
	d.assign(dist);
}



//
// qi_class_real::assign(qi_class, xbigfloat)
//
// Task:
//      set to a copy of B and initializes with distance dist.
//

void
qi_class_real::assign (const qi_class & B, const xbigfloat & dist)
{
	debug_handler("qi_class_real", "assign(qi_class, xbigfloat)");

	assign(B);
	d.assign(dist);
}



//
// qi_class_real::assign(qi_class_real)
//
// Task:
//      set to a copy of B.
//

void
qi_class_real::assign (const qi_class_real & B)
{
	debug_handler("qi_class_real", "assign(qi_class_real)");

	a.assign(B.a);
	b.assign(B.b);
	d.assign(B.d);
}



//
// operator =
//
// Task:
//      make a copy of an existing qi_class_real
//

qi_class_real & qi_class_real::operator = (const qi_class_real & A)
{
	debug_handler("qi_class_real", "operator = ");

	assign(A);
	return *this;
}



//
// get_distance
//
// Task:
//      return the distance
//

bigfloat
qi_class_real::get_distance () const
{
	bigfloat bfd;
	long e_tmp;

	bfd.assign(d.get_mantissa());
	e_tmp = d.get_exponent()-d.get_mantissa().bit_length();
	if (e_tmp > 0)
		shift_left(bfd, bfd, e_tmp);
	else
		shift_right(bfd, bfd, -e_tmp);

	return bfd;
}



//
// get_distance_x
//
// Task:
//      return the distance as an xbigfloat
//

xbigfloat
qi_class_real::get_distance_x () const
{
	return d;
}



//
// multiply()
//
// Task:
//      computes the product of A and B (reduced ideal product) and the
//      new distance.
//

void
multiply (qi_class_real & C, const qi_class_real & A, const qi_class_real & B)
{
	debug_handler("qi_class_real", "multiply");

	multiply_real(C, A, B);
}



//
// multiply_real()
//
// Task:
//      multiplies two real ideal equivalence classes
//

void
multiply_real (qi_class_real & C, const qi_class_real & A,
	       const qi_class_real & B)
{
	debug_handler("qi_class_real", "multiply_real");

	bigint newa, newb, gpr, v, g, w, ab2, temp;

	// solve gpr = v A.a + w B.a
	gpr.assign(xgcd_left(v, A.a, B.a));

	// C.b = v A.a (B.b - A.b)
	multiply(newb, v, A.a);
	subtract(temp, B.b, A.b);
	multiply(newb, newb, temp);

	// C.a = A.a B.a
	multiply(newa, A.a, B.a);

	if (!gpr.is_one()) {
		add(ab2, A.b, B.b);
		ab2.divide_by_2();
		g.assign(xgcd(v, w, gpr, ab2));

		// C.b = (C.b*v + w(D - A.b^2)/2) / d
		multiply(newb, newb, v);

		square(temp, A.b);
		subtract(temp, qi_class_real::discriminant(), temp);
		temp.divide_by_2();
		multiply(temp, temp, w);

		add(newb, newb, temp);
		divide(newb, newb, g);

		// C.a = C.a / (g^2)
		square(temp, g);
		divide(newa, newa, temp);

		gpr.assign(g);
	}

	C.a.assign(newa);
	add(newb, newb, A.b);
	shift_left(ab2, C.a, 1);
	remainder(C.b, newb, ab2);

	add(C.d, A.d, B.d);
//  subtract(C.d,C.d,log(bigfloat(gpr)));

	C.reduce();
}



//
// qi_class_real::invert()
//
// Task:
//      inverts the qi_class_real and updates its distance
//

void
qi_class_real::invert ()
{
	debug_handler("qi_class_real", "invert");

	b.negate();
	normalize_real();
	d.negate();
//  add(d,d,log(bigfloat(a)));
}



//
// inverse(qi_class_real, qi_class_real)
//
// Task:
//      computes the inverse of B
//

void
inverse (qi_class_real & A, const qi_class_real & B)
{
	debug_handler("qi_class_real", "inverse(qi_class_real, qi_class_real)");

	A.assign(B);
	A.b.negate();
	A.normalize_real();
	A.d.negate();
//  add(A.d,A.d,log(bigfloat(A.a)));
}



//
// inverse(qi_class_real)
//
// Task:
//      returns the inverse of B
//

qi_class_real
inverse (const qi_class_real & B)
{
	debug_handler("qi_class_real", "inverse(qi_class_real)");

	qi_class_real A;

	A.assign(B);
	A.b.negate();
	A.normalize_real();
	A.d.negate();
//  add(A.d,A.d,log(bigfloat(A.a)));

	return A;
}



//
// divide()
//
// Task:
//      computes A*B^-1.
//

void
divide (qi_class_real & C, const qi_class_real & A, const qi_class_real & B)
{
	debug_handler("qi_class_real", "divide");

	qi_class_real temp;

	if (B.is_zero()) {
		lidia_error_handler("qi_class_real", "divide() - zero divisor");
		return;
	}

	temp = inverse(B);
	multiply_real(C, A, temp);
}



//
// square
//
// Task:
//      multiplies A with itself and computes the new distance.
//

void
square (qi_class_real & C, const qi_class_real & A)
{
	debug_handler("qi_class_real", "square");

	square_real(C, A);
}



//
// square_real()
//
// Task:
//      multiplies the real ideal equivalence class with itself and computes
//      the new distance.
//

void
square_real (qi_class_real & C, const qi_class_real & A)
{
	debug_handler("qi_class_real", "square_real");

	bigint newb, g, w, temp;

	// solve g = v A.a + w A.b
	g.assign(xgcd_right(w, A.a, A.b));

	// C.b = A.b + w (D - A.b^2) / (2d)
	square(newb, A.b);
	subtract(newb, qi_class_real::discriminant(), newb);
	shift_left(temp, g, 1);
	divide(newb, newb, temp);
	multiply(newb, newb, w);
	add(C.b, newb, A.b);

	// C.a = (A.a/g)^2
	divide(temp, A.a, g);
	square(C.a, temp);

	shift_left(temp, C.a, 1);
	remainder(C.b, C.b, temp);

	C.d.assign(A.d);
	C.d.multiply_by_2();
//  subtract(C.d,C.d,log(bigfloat(g)));

	C.reduce();
}



//
// power(bigint)
//
// Task:
//      computes A^i
//

void
power (qi_class_real & C, const qi_class_real & A, const bigint & i)
{
	debug_handler("qi_class_real", "power(bigint)");

	power_real(C, A, i);
}



//
// power(long)
//
// Task:
//      computes A^i

void
power (qi_class_real & C, const qi_class_real & A, const long i)
{
	debug_handler("qi_class_real", "power(long)");

	power_real(C, A, i);
}



//
// power_real(bigint)
//
// Task:
//      computes A^i using binary exponentiation
//

void
power_real (qi_class_real & C, const qi_class_real & A, const bigint & i)
{
	debug_handler("qi_class_real", "power_real(bigint)");

	qi_class_real B;
	bigint j;

	B.assign(A);
	j.assign(i);
	if (j.is_lt_zero()) {
		B.invert();
		j.absolute_value();
	}
	C.assign_one();
	while (j.is_gt_zero()) {
		if (j.is_odd())
			multiply_real(C, C, B);
		j.divide_by_2();
		if (j.is_gt_zero())
			square_real(B, B);
	}
}



//
// power_real(long)
//
// Task:
//      computes A^i using binary exponentiation
//

void
power_real (qi_class_real & C, const qi_class_real & A, const long i)
{
	debug_handler("qi_class_real", "power_real(long)");

	qi_class_real B;
	register long j;

	B.assign(A);
	j = i;
	if (j < 0) {
		B.invert();
		j = -i;
	}
	C.assign_one();
	while (j > 0) {
		if ((j & 1) == 1)
			multiply_real(C, C, B);
		j >>= 1;
		if (j > 0)
			square_real(B, B);
	}
}



//
// nearest()
//
// Task:
//      computes the ideal equivalent to S with distance from S as close as
//      possible to E.
//

qi_class_real
nearest (const qi_class_real & S, const bigfloat & E)
{
	debug_handler("qi_class_real", "nearest");

	long j, k;
	bigfloat EK, temp, temp2;
	qi_class_real A, B;

	if (E.is_zero())
		return S;

	// compute k such that 2^k < E < 2^(k+1)
	divide(temp, log(abs(E)), log(bigfloat(2.0)));
	floor(temp2, temp);
	temp2.longify(k);

	// EK = E / 2^k
	shift_right(EK, E, k);

	// compute B with distance < EK
	B.assign_one();
	if (EK.is_gt_zero())
		B.adjust_pos(EK);
	else
		B.adjust_neg(EK);

	for (j = 0; j < k; ++j) {
		EK.multiply_by_2();
		square(B, B);
		if (EK.is_gt_zero())
			B.adjust_pos(EK);
		else
			B.adjust_neg(EK);
	}

	if (S.is_one())
		A.assign(B);
	else {
		A.assign(S, 0.0);
		multiply_real(A, A, B);
	}
	if (E.is_gt_zero())
		A.adjust_pos(E);
	else
		A.adjust_neg(E);

	temp = abs(A.get_distance() - E);
	apply_rho(B, A);
	temp2 = abs(B.get_distance() - E);
	if (temp2 < temp)
		A.assign(B);
	else {
		apply_inverse_rho(B, A);
		temp2 = abs(B.get_distance() - E);
		if (temp2 < temp)
			A.assign(B);
	}

	return A;
}



//
// operator -
//
// Task:
//      returns A^-1
//

qi_class_real operator - (const qi_class_real & A)
{
	debug_handler("qi_class_real", "operator -");

	qi_class_real B;

	B = inverse(A);
	return B;
}



//
// operator *
//
// Task:
//      multiplies A and B
//

qi_class_real operator * (const qi_class_real & A, const qi_class_real & B)
{
	debug_handler("qi_class_real", "operator *");

	qi_class_real C;

	multiply_real(C, A, B);
	return C;
}



//
// operator /
//
// Task:
//      multiplies A and B^-1
//

qi_class_real operator / (const qi_class_real & A, const qi_class_real & B)
{
	debug_handler("qi_class_real", "operator /");

	qi_class_real C;

	divide(C, A, B);
	return C;
}



//
// operator *=
//
// Task:
//      *this = *this * A
//

qi_class_real & qi_class_real::operator *= (const qi_class_real & A)
{
	debug_handler("qi_class_real", "operator *= ");

	multiply_real(*this, *this, A);
	return *this;
}



//
// operator /=
//
// Task:
//      *this = *this / A
//

qi_class_real & qi_class_real::operator /= (const qi_class_real & A)
{
	debug_handler("qi_class_real", "operator /= ");

	multiply_real(*this, *this, inverse(A));
	return *this;
}



//
// swap
//
// Task:
//      swaps A and B
//

void
swap (qi_class_real & A, qi_class_real & B)
{
	debug_handler("qi_class_real", "swap");

	qi_class_real C;

	C.assign(A);
	A.assign(B);
	B.assign(C);
}



//
// generate_prime_ideal()
//
// Task:
//      computes a reduced representative of the equivalence class containing
//      the ideal lying over the prime p.  If such an ideal doesn't exist,
//      false is returned.
//

bool
generate_prime_ideal (qi_class_real & A, const bigint & p)
{
	debug_handler("qi_class_real", "generate_prime_ideal");

	bigint D, b, temp;
	int kro, Dp, ip;

	D.assign(qi_class_real::discriminant());
	if (D.is_zero()) {
		lidia_error_handler("qi_class_real", "generate_prime_ideal - current "
				    "quadratic order is not defined or is not real");
		return false;
	}

	kro = kronecker(D, p);
	if (kro < 0)
		return false;
	else {
		if (p == 2) {
			if (kro == 0) {
				Dp = static_cast<int>(remainder(D, 16L));
				if (Dp == 8)
					b.assign_zero();
				else if (Dp == 12)
					b.assign(p);
				else
					return false;
			}
			else
				b.assign(1);
		}
		else {
			if (kro == 0) {
				remainder(temp, D, p*p);
				if (temp.is_zero())
					return false;
				else
					b.assign_zero();
			}
			else {
				if (is_int(p)) {
					p.intify(ip);
					Dp = static_cast<int>(remainder(D, static_cast<long>(ip)));
					if (Dp < 0)
						Dp += ip;

					b.assign(ressol(Dp, ip));
				}
				else {
					remainder(temp, D, p);
					if (temp.is_lt_zero())
						add(temp, temp, p);
					ressol(b, temp, p);
				}

				if (b.is_lt_zero())
					add(b, b, p);
			}

			if ((D.is_odd()) != (b.is_odd()))
				subtract(b, p, b);
		}

		A.assign(p, b);
		return true;
	}
}



//
// qi_class_real::rho()
//
// Task:
//      applies the reduction operator once, resulting in another reduced
//      representative of the equivalence class.  The distance is adjusted
//      appropriately.
//

void
qi_class_real::rho ()
{
	debug_handler("qi_class_real", "rho");

	bigint q, r, a2, nb, na, oa, temp;
	int s;

	// oa = (D - b^2) / 4a
	square(oa, b);
	subtract(oa, Delta, oa);
	divide(oa, oa, a);
	shift_right(oa, oa, 2);

	s = a.is_gt_zero() ? 0 : 1;

	// (rootD-b) = (2a)*q + r
	shift_left(a2, a, 1);
	add(temp, rootD, b);
	if (s)
		inc(temp);
	div_rem(q, r, temp, a2);

	// nb = rootD + s - r;
	subtract(nb, rootD, r);
	if (s)
		inc(nb);

	// na = oa -q*((nb - b) >> 1);
	subtract(temp, nb, b);
	temp.divide_by_2();
	multiply(temp, temp, q);
	subtract(na, oa, temp);

	// d += ln( (nb + rd) / 2a )
	add(d, d, approx_Ln(nb, rd, xprec));

	b.assign(nb);
	a.assign(na);
}



//
// apply_rho(qi_class_real, qi_class_real)
//
// Task:
//      applies the reduction operator once, resulting in another reduced
//      representative of the equivalence class.  The distance is adjusted
//      appropriately.
//

void
apply_rho (qi_class_real & A, const qi_class_real & B)
{
	debug_handler("qi_class_real", "apply_rho(qi_class_real, qi_class_real)");

	A.assign(B);
	A.rho();
}



//
// qi_class_real::apply_rho(qi_class_real)
//
// Task:
//      applies the reduction operator once, resulting in another reduced
//      representative of the equivalence class.  The distance is adjusted
//      appropriately.
//

qi_class_real
apply_rho (const qi_class_real & A)
{
	debug_handler("qi_class_real", "apply_rho(qi_class_real)");

	qi_class_real B;

	B.assign(A);
	B.rho();
	return B;
}



//
// qi_class_real::inverse_rho()
//
// Task:
//      applies the inverse reduction operator once, resulting in another
//      reduced representative of the equivalence class.  The distance is
//      adjusted appropriately.
//

void
qi_class_real::inverse_rho ()
{
	debug_handler("qi_class_real", "inverse_rho");

	bigint q, a2, r, temp;

	// a = (D - b^2) / 4a
	square(temp, b);
	subtract(temp, Delta, temp);
	divide(temp, temp, a);
	shift_right(a, temp, 2);

	// d -= ln( (b + rd) / 2a )
	subtract(d, d, approx_Ln(b, rd, xprec));

	// q = floor((rootD + b) / 2a)
	shift_left(a2, a, 1);
	add(temp, rootD, b);
	div_rem(q, r, temp, a2);
	if ((temp.is_lt_zero()) && (!r.is_zero()))
		dec(q);

	// b = 2a(nq) - b
	multiply(temp, a2, q);
	subtract(b, temp, b);
}



//
// apply_inverse_rho(qi_class_real, qi_class_real)
//
// Task:
//      applies the inverse reduction operator once, resulting in another
//      reduced representative of the equivalence class.  The distance is
//      adjusted appropriately.
//

void
apply_inverse_rho (qi_class_real & A, const qi_class_real & B)
{
	debug_handler("qi_class_real", "apply_inverse_rho(qi_class_real, "
		      "qi_class_real)");

	A.assign(B);
	A.inverse_rho();
}



//
// apply_inverse_rho(qi_class_real)
//
// Task:
//      applies the inverse reduction operator once, resulting in another
//      reduced representative of the equivalence class.  The distance is
//      adjusted appropriately.
//

qi_class_real
apply_inverse_rho (const qi_class_real & A)
{
	debug_handler("qi_class_real", "apply_inverse_rho(qi_class_real)");

	qi_class_real B;

	B.assign(A);
	B.inverse_rho();
	return B;
}



//
// qi_class_real::convert_distance()
//
// Task:
//      given that the distance from A to the given ideal using Lenstra's
//      formulation is dist, returns the value of the distance from A to the
//      given ideal using standard formulation.
//

bigfloat
qi_class_real::convert_distance (const qi_class_real & A,
				 const bigfloat & dist) const
{
	debug_handler("qi_class_real", "convert_distance");

	qi_class_real AA, C, D, U;
	bigfloat newdist;
	bool done;

	AA.assign(A, 0);

	C = nearest(A, dist);
	D.assign(C);
	if (info > 3) {
		std::cout << "IN CONVERT:" << std::endl;
		std::cout << "*this = " << *this << std::endl;
		std::cout << "A = " << AA << std::endl;
		std::cout << "C = " << C << std::endl;
	}
	done = (is_equal(inverse(C)) && (dist > current_order->get_qo()->min_log)) ||
		(is_equal(inverse(D)) && (dist > current_order->get_qo()->min_log));
	while ((!is_equal(C)) && (!is_equal(D)) && (!done)) {
		C.rho();
		D.inverse_rho();
		if (info > 3) {
			std::cout << "\nC = " << C << std::endl;
			std::cout << "D = " << D << std::endl;
		}
		done = (is_equal(inverse(C)) && (dist > current_order->get_qo()->min_log)) ||
			(is_equal(inverse(D)) && (dist > current_order->get_qo()->min_log));
	}

	if (is_equal(C))
		newdist.assign(C.get_distance());
	else if (is_equal(D))
		newdist.assign(D.get_distance());
	else if (is_equal(inverse(C))) {
		subtract(newdist, current_order->get_qo()->regulator(), C.get_distance());

		U.assign_one();
		U = nearest(U, newdist);
		multiply(U, U, AA);
		newdist = U.get_distance();

		while (newdist > current_order->get_qo()->regulator())
			subtract(newdist, newdist, current_order->get_qo()->regulator());
		while (newdist < 0)
			add(newdist, newdist, current_order->get_qo()->regulator());
		newdist = convert_distance(A, newdist);
	}
	else {
		subtract(newdist, current_order->get_qo()->regulator(), D.get_distance());

		U.assign_one();
		U = nearest(U, newdist);
		multiply(U, U, AA);
		newdist = U.get_distance();

		while (newdist > current_order->get_qo()->regulator())
			subtract(newdist, newdist, current_order->get_qo()->regulator());
		while (newdist < 0)
			add(newdist, newdist, current_order->get_qo()->regulator());
		newdist = convert_distance(A, newdist);
	}

	while (newdist > current_order->get_qo()->regulator())
		subtract(newdist, newdist, current_order->get_qo()->regulator());
	while (newdist < 0)
		add(newdist, newdist, current_order->get_qo()->regulator());

	return newdist;
}



//
// qi_class_real::adjust_pos()
//
// Task:
//      computes the reduced representative in the same class as the ideal
//      with distance <= dist and as close to dist as possible.
//

void
qi_class_real::adjust_pos (const bigfloat & dist)
{
	debug_handler("qi_class_real", "adjust_pos");

	bigint q, r, a2, nb, na, oa, temp;
	int s;

	// oa = (D - b^2) / 4a
	square(oa, b);
	subtract(oa, Delta, oa);
	divide(oa, oa, a);
	shift_right(oa, oa, 2);

	while (d <= dist) {
		s = a.is_gt_zero() ? 0 : 1;

		// (rootD-b) = (2a)*q + r
		shift_left(a2, a, 1);
		add(temp, rootD, b);
		if (s)
			inc(temp);
		div_rem(q, r, temp, a2);

		// nb = rootD + s - r;
		subtract(nb, rootD, r);
		if (s)
			inc(nb);

		// na = oa -q*((nb - b) >> 1);
		subtract(temp, nb, b);
		temp.divide_by_2();
		multiply(temp, temp, q);
		subtract(na, oa, temp);

		// d += ln( (nb + rd) / 2a )
		add(d, d, approx_Ln(nb, rd, xprec));

		b.assign(nb);
		oa.assign(a);
		a.assign(na);
	}

	while (d > dist) {
		a.assign(oa);

		shift_left(a2, a, 1);

		// d -= ln( (b + rd) / 2a )
		subtract(d, d, approx_Ln(b, rd, xprec));

		add(temp, rootD, b);
		div_rem(q, r, temp, a2);
		if ((temp.is_lt_zero()) && (!r.is_zero()))
			dec(q);

		multiply(temp, a2, q);
		subtract(b, temp, b);

		if (d > dist) {
			square(oa, b);
			subtract(oa, Delta, oa);
			divide(oa, oa, a);
			shift_right(oa, oa, 2);
		}
	}
}



//
// qi_class_real::adjust_neg()
//
// Task:
//      computes the reduced representative in the same class as the ideal
//      with distance >= dist and as close to dist as possible.  The distance
//      of the ideal and dist are assumed to be negative.
//

void
qi_class_real::adjust_neg (const bigfloat & dist)
{
	debug_handler("qi_class_real", "adjust_neg");

	bigint q, r, a2, nb, na, oa, temp;
	int s;

	// oa = (D - b^2) / 4a
	square(oa, b);
	subtract(oa, Delta, oa);
	divide(oa, oa, a);
	shift_right(oa, oa, 2);

	while (d >= dist) {
		a.assign(oa);

		shift_left(a2, a, 1);

		// d -= ln( (b + rd) / 2a )
		subtract(d, d, approx_Ln(b, rd, xprec));

		add(temp, rootD, b);
		div_rem(q, r, temp, a2);
		if ((temp.is_lt_zero()) && (!r.is_zero()))
			dec(q);

		multiply(temp, a2, q);
		subtract(b, temp, b);

		square(oa, b);
		subtract(oa, Delta, oa);
		divide(oa, oa, a);
		shift_right(oa, oa, 2);
	}

	while (d < dist) {
		s = a.is_gt_zero() ? 0 : 1;

		// (rootD-b) = (2a)*q + r
		shift_left(a2, a, 1);
		add(temp, rootD, b);
		if (s)
			inc(temp);
		div_rem(q, r, temp, a2);

		// nb = rootD + s - r;
		subtract(nb, rootD, r);
		if (s)
			inc(nb);

		// na = oa -q*((nb - b) >> 1);
		subtract(temp, nb, b);
		temp.divide_by_2();
		multiply(temp, temp, q);
		subtract(na, oa, temp);

		// d += ln( (nb + rd) / 2a )
		add(d, d, approx_Ln(nb, rd, xprec));

		b.assign(nb);
		oa.assign(a);
		a.assign(na);
	}
}



//
// qi_class_real::adjust_abs()
//
// Task:
//      computes the reduced representative in the same class as the ideal
//      with |distance| >= dist and as close to dist as possible.
//

void
qi_class_real::adjust_abs (const bigfloat & dist)
{
	debug_handler("qi_class_real", "adjust_abs");

	bigint q, r, a2, nb, na, oa, temp;
	int s;
	xbigfloat dmod;

	// oa = (D - b^2) / 4a
	square(oa, b);
	subtract(oa, Delta, oa);
	divide(oa, oa, a);
	shift_right(oa, oa, 2);

	dmod.assign(d);
	dmod.absolute_value();
	while (dist >= dmod) {
		a.assign(oa);

		shift_left(a2, a, 1);

		// d -= ln( (b + rd) / 2a )
		subtract(d, d, approx_Ln(b, rd, xprec));

		add(temp, rootD, b);
		div_rem(q, r, temp, a2);
		if ((temp.is_lt_zero()) && (!r.is_zero()))
			dec(q);

		multiply(temp, a2, q);
		subtract(b, temp, b);

		square(oa, b);
		subtract(oa, Delta, oa);
		divide(oa, oa, a);
		shift_right(oa, oa, 2);

		dmod.assign(d);
		dmod.absolute_value();
	}

	while (dist < dmod) {
		s = a.is_gt_zero() ? 0 : 1;

		// (rootD-b) = (2a)*q + r
		shift_left(a2, a, 1);
		add(temp, rootD, b);
		if (s)
			inc(temp);
		div_rem(q, r, temp, a2);

		// nb = rootD + s - r;
		subtract(nb, rootD, r);
		if (s)
			inc(nb);

		// na = oa -q*((nb - b) >> 1);
		subtract(temp, nb, b);
		temp.divide_by_2();
		multiply(temp, temp, q);
		subtract(na, oa, temp);

		// d += ln( (nb + rd) / 2a )
		add(d, d, approx_Ln(nb, rd, xprec));

		b.assign(nb);
		oa.assign(a);
		a.assign(na);

		dmod.assign(d);
		dmod.absolute_value();
	}
}



//
// qi_class_real::make_list(bigint, bigfloat, hash_table)
//
// Task:
//      computes all reduced representatives of the qi_class_real with distance
//      from it less than dist and adds them to the hash table.  If we find an
//      ideal equal to 1, we exit and return its distance (the regulator).
//

bigfloat
qi_class_real::make_list (bigint & oa, const bigfloat & dist,
			  hash_table< qi_class_real > & HT)
{
	debug_handler("qi_class_real", "make_list(bigint, bigfloat, hash_table)");

	bigint q, r, a2, nb, na, temp;
	int s;
	bigfloat R;

	R.assign_zero();
	while (d < dist) {
		s = a.is_gt_zero() ? 0 : 1;

		// (rootD-b) = (2a)*q + r
		shift_left(a2, a, 1);
		add(temp, rootD, b);
		if (s)
			inc(temp);
		div_rem(q, r, temp, a2);

		// nb = rootD + s - r;
		subtract(nb, rootD, r);
		if (s)
			inc(nb);

		// na = oa -q*((nb - b) >> 1);
		subtract(temp, nb, b);
		temp.divide_by_2();
		multiply(temp, temp, q);
		subtract(na, oa, temp);

		// d += ln( (nb + rd) / 2a )
		add(d, d, approx_Ln(nb, rd, xprec));

		b.assign(nb);
		oa.assign(a);
		a.assign(na);

		HT.hash(*this);
		if (is_one()) {
			R.assign(d.get_mantissa());
			long e_tmp = d.get_exponent()-d.get_mantissa().bit_length();
			if (e_tmp > 0)
				shift_left(R, R, e_tmp);
			else
				shift_right(R, R, -e_tmp);
			break;
		}
	}

	return R;
}



//
// qi_class_real::make_list(bigint, bigfloat, qi_class_real, qi_class_real,
//   hash_table)
//
// Task:
//      computes all reduced representatives of the qi_class_real with distance
//      from it less than dist and adds them to the hash table.  If an ideal
//      equal to G or B is found, then we exit and return the appropriate
//      flag (1 if equal to G, 0 if equal to B, -1 if neither).
//

int
qi_class_real::make_list (bigint & oa, const bigfloat & dist,
			  const qi_class_real & G, const qi_class_real & U,
			  hash_table< qi_class_real > & HT)
{
	debug_handler("qi_class_real", "make_list(bigint, bigfloat, qi_class_real, "
		      "qi_class_real, hash_table");

	bigint q, r, a2, nb, na, temp;
	int s, flag = -1;

	while (d < dist) {
		s = a.is_gt_zero() ? 0 : 1;

		// (rootD-b) = (2a)*q + r
		shift_left(a2, a, 1);
		add(temp, rootD, b);
		if (s)
			inc(temp);
		div_rem(q, r, temp, a2);

		// nb = rootD + s - r;
		subtract(nb, rootD, r);
		if (s)
			inc(nb);

		// na = oa -q*((nb - b) >> 1);
		subtract(temp, nb, b);
		temp.divide_by_2();
		multiply(temp, temp, q);
		subtract(na, oa, temp);

		// d += ln( (nb + rd) / 2a )
		add(d, d, approx_Ln(nb, rd, xprec));

		b.assign(nb);
		oa.assign(a);
		a.assign(na);

		HT.hash(*this);
		if (is_equal(G)) {
			flag = 1;
			break;
		}
		if (is_equal(U)) {
			flag = 0;
			break;
		}
	}

	return flag;
}



//
// qi_class_real::make_list(bigfloat, long, hash_table)
//
// Task:
//      computes all reduced representative of the qi_class_real with distance
//      from it less than dist and adds them to the hash table.
//

void
qi_class_real::make_list (const bigfloat & dist, long x,
			  hash_table< ideal_node > & HT)
{
	debug_handler("qi_class_real", "make_list(bigfloat, long, hash_table)");

	bigint q, r, a2, nb, na, oa, temp;
	int s;

	// oa = (D - b^2) / 4a
	square(oa, b);
	subtract(oa, Delta, oa);
	divide(oa, oa, a);
	shift_right(oa, oa, 2);

	while (d < dist) {
		s = a.is_gt_zero() ? 0 : 1;

		// (rootD-b) = (2a)*q + r
		shift_left(a2, a, 1);
		add(temp, rootD, b);
		if (s)
			inc(temp);
		div_rem(q, r, temp, a2);

		// nb = rootD + s - r;
		subtract(nb, rootD, r);
		if (s)
			inc(nb);

		// na = oa -q*((nb - b) >> 1);
		subtract(temp, nb, b);
		temp.divide_by_2();
		multiply(temp, temp, q);
		subtract(na, oa, temp);

		// d += ln( (nb + rd) / 2a )
		add(d, d, approx_Ln(nb, rd, xprec));

		b.assign(nb);
		oa.assign(a);
		a.assign(na);

		HT.hash(ideal_node(qi_class(*this), x));
	}

	s = a.is_gt_zero() ? 0 : 1;

	// (rootD-b) = (2a)*q + r
	shift_left(a2, a, 1);
	add(temp, rootD, b);
	if (s)
		inc(temp);
	div_rem(q, r, temp, a2);

	// nb = rootD + s - r;
	subtract(nb, rootD, r);
	if (s)
		inc(nb);

	// na = oa -q*((nb - b) >> 1);
	subtract(temp, nb, b);
	temp.divide_by_2();
	multiply(temp, temp, q);
	subtract(na, oa, temp);

	// d += ln( (nb + rd) / 2a )
	add(d, d, approx_Ln(nb, rd, xprec));

	b.assign(nb);
	a.assign(na);
	HT.hash(ideal_node(qi_class(*this), x));
}



//
// qi_class_real::make_list(bigfloat, long, indexed_hash_table)
//
// Task:
//      computes all reduced representative of the qi_class_real with distance
//      from it less than dist and adds them to the hash table.
//

void
qi_class_real::make_list (const bigfloat & dist, long x,
			  indexed_hash_table< ideal_node > & HT)
{
	debug_handler("qi_class_real", "make_list(bigfloat, long, indexed_hash_table)");

	bigint q, r, a2, nb, na, oa, temp;
	int s;

	// oa = (D - b^2) / 4a
	square(oa, b);
	subtract(oa, Delta, oa);
	divide(oa, oa, a);
	shift_right(oa, oa, 2);

	while (d < dist) {
		s = a.is_gt_zero() ? 0 : 1;

		// (rootD-b) = (2a)*q + r
		shift_left(a2, a, 1);
		add(temp, rootD, b);
		if (s)
			inc(temp);
		div_rem(q, r, temp, a2);

		// nb = rootD + s - r;
		subtract(nb, rootD, r);
		if (s)
			inc(nb);

		// na = oa -q*((nb - b) >> 1);
		subtract(temp, nb, b);
		temp.divide_by_2();
		multiply(temp, temp, q);
		subtract(na, oa, temp);

		// d += ln( (nb + rd) / 2a )
		add(d, d, approx_Ln(nb, rd, xprec));

		b.assign(nb);
		oa.assign(a);
		a.assign(na);

		HT.hash(ideal_node(qi_class(*this), x));
	}

	s = a.is_gt_zero() ? 0 : 1;

	// (rootD-b) = (2a)*q + r
	shift_left(a2, a, 1);
	add(temp, rootD, b);
	if (s)
		inc(temp);
	div_rem(q, r, temp, a2);

	// nb = rootD + s - r;
	subtract(nb, rootD, r);
	if (s)
		inc(nb);

	// na = oa -q*((nb - b) >> 1);
	subtract(temp, nb, b);
	temp.divide_by_2();
	multiply(temp, temp, q);
	subtract(na, oa, temp);

	// d += ln( (nb + rd) / 2a )
	add(d, d, approx_Ln(nb, rd, xprec));

	b.assign(nb);
	a.assign(na);
	HT.hash(ideal_node(qi_class(*this), x));
}



//
// qi_class_real::is_principal()
//
// Task:
//      tests if the equivalence class is the principal class.  The alogirthm
//      is chosen based on the size of the discriminant.
//

bool
qi_class_real::is_principal () const
{
	debug_handler("qi_class_real", "is_principal()");

	bool is_prin;
	bigfloat dist;
	int num;
	quadratic_order *QO = current_order->get_qo();

	num = decimal_length(Delta);
	if ((!QO->is_R_computed() && (num < RSHANKSB)) ||
	    (QO->prin_list.no_of_elements() > 0))
		is_prin = is_principal_buch(dist);
	else
		is_prin = is_principal_subexp();

	return is_prin;
}



//
// qi_class_real::is_principal_buch()
//
// Task:
//      tests if the equivalence class is the principal class using a variation
//      of baby-step giant-step due to Buchmann.
//

bool
qi_class_real::is_principal_buch () const
{
	debug_handler("qi_class_real", "is_principal_buch()");

	bigfloat dist;
	quadratic_order *QO = current_order->get_qo();

	if (Delta.is_le_zero()) {
		lidia_error_handler("qi_class_real", "is_principal_buch() - no "
				    "current quadratic order has been defined, or it is not real");
		return false;
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class_real", "is_principal_buch() - the "
				    "current quadratic order has been deleted");
		return false;
	}

	if (QO->is_subexp_computed())
		return is_principal_subexp(dist);
	else
		return is_principal_buch(dist);
}



//
// qi_class_real::is_principal_subexp()
//
// Task:
//      tests if the equivalence class is the principal class using a variation
//      of the subexponential algorithm due to Abel.
//

bool
qi_class_real::is_principal_subexp () const
{
	debug_handler("qi_class_real", "is_principal_subexp()");

	quadratic_order *QO;
	bool is_prin;
	int num;
	timer t;
	bigfloat dist;
	long oprec = bigfloat::get_precision();

	if (Delta.is_le_zero()) {
		lidia_error_handler("qi_class_real", "is_principal_subexp() - no "
				    "current quadratic order has been defined, or it is not real");
		return false;
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class_real", "is_principal_subexp() - the current "
				    "quadratic order has been deleted");
		return false;
	}

	QO = current_order->get_qo();
	bigfloat::set_precision(QO->prec);

	if (info)
		t.start_timer();

	if (is_one()) {
		is_prin = true;
		dist.assign_zero();
	}
	else {
		num = decimal_length(Delta);
		if ((!QO->is_R_computed() && (num < RSHANKSB)) ||
		    (QO->prin_list.no_of_elements() > 0))
			is_prin = is_principal_buch(dist);
		else {
			QO->regulator(true);
			if (do_verify)
				is_prin = QO->is_in_lattice(*this, dist);
			else
				is_prin = QO->is_in_lattice(*this);
		}
	}

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (is_principal_subexp) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	if (do_verify) {
		if (verify_principal(dist, is_prin)) {
			if (info)
				std::cout << "Principality test is correct." << std::endl;
		}
		else
			std::cout << "Principality test is not correct!" << std::endl;
	}

	bigfloat::set_precision(oprec);
	return is_prin;
}



//
// qi_class_real::is_principal(bigfloat)
//
// Task:
//      tests if the equivalence class is the principal class.  The alogirthm
//      is chosen based on the size of the discriminant.  The distance from
//      the unit ideal to the representative is also computed.
//

bool
qi_class_real::is_principal (bigfloat & dist) const
{
	debug_handler("qi_class_real", "is_principal(bigfloat)");

	bool is_prin;
	int num;

	// FOR TABLE GENERATION
	timer t;
	long t1 = 0;

	bool old = quadratic_order::special;
	quadratic_order::special = false;
	if (old) {
		quadratic_order::special2 = true;
		t.start_timer();
	}

	num = decimal_length(Delta);
	if (num < RSHANKSB)
		is_prin = is_principal_buch(dist);
	else
		is_prin = is_principal_subexp(dist);

	// FOR TABLE GENERATION
	quadratic_order::special2 = false;
	quadratic_order::special = old;
	if (quadratic_order::special) {
		t.stop_timer();
		t1 = t.user_time();
		// ROW:  time for representations    total time
		std::cout << t1 << std::endl;
	}

	return is_prin;
}



//
// qi_class_real::is_principal_buch(bigfloat)
//
// Task:
//      tests if the equivalence class is the principal class using a variation
//      of baby-step giant-step due to Buchmann.  The distance from the unit
//      ideal to the representative is also computed.
//

bool
qi_class_real::is_principal_buch (bigfloat & dist) const
{
	debug_handler("qi_class_real", "is_principal_buch(bigfloat)");

	bigfloat y, usqr, u, temp;
	long upper;
	bigint oa, temp2;
	qi_class_real A, C, BB, D, G, U, *SVAL;
	int flag;
	bool is_prin, done;
	timer t;
	quadratic_order *cur_qo;
	long oprec = bigfloat::get_precision();

	if (Delta.is_le_zero()) {
		lidia_error_handler("qi_class_real", "is_principal_buch(bigfloat) - no "
				    "current quadratic order has been defined, or it is not real");
		return false;
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class_real", "is_principal_buch(bigfloat) - the "
				    "current quadratic order has been deleted");
		return false;
	}

	cur_qo = current_order->get_qo();
	bigfloat::set_precision(cur_qo->prec);

	if ((cur_qo->is_subexp_computed()) &&
	    (cur_qo->prin_list.no_of_elements() == 0)) {
		bigfloat::set_precision(oprec);
		return is_principal_subexp(dist);
	}

	if (info)
		t.start_timer();

	G.assign(*this, 0.0);
	U.assign_one();
	done = false;
	dist.assign_zero();

	if (is_one())
		done = is_prin = true;
	else {
		is_prin = false;

		if (cur_qo->prin_list.no_of_elements() > 0) {
			A = cur_qo->prin_list.last_entry();
			if ((A.is_one()) && (A.d > 0)) {
				// whole cycle has already been computed - table look-up
				SVAL = cur_qo->prin_list.search(G);
				if (SVAL) {
					done = is_prin = true;
					dist = SVAL->get_distance();
				}
				else {
					is_prin = false;
					done = true;
					dist = A.get_distance();
				}
			}
			else {
				SVAL = cur_qo->prin_list.search(G);
				if (SVAL) {
					done = is_prin = true;
					dist = SVAL->get_distance();
				}
				ceil(temp2, A.d);
				if (temp2.is_odd())
					inc(temp2);
				u.assign(temp2);
				y.assign(u);
			}
		}
		else {
			// get hash table size
			temp = ceil(power(bigfloat(Delta), bigfloat(0.25)));
			temp.multiply_by_2();
			multiply(temp, temp, bigfloat(1.33));
			temp.longify(upper);
			cur_qo->prin_list.initialize(static_cast<lidia_size_t>(upper));
			cur_qo->prin_list.set_key_function(&qi_class_real_key);

			add(temp, log(bigfloat(Delta)), 3);
			shift_right(temp, temp, 2);
			ceil(temp, temp);
			temp.bigintify(temp2);
			if (temp2.is_odd())
				inc(temp2);
			u.assign(temp2);
			y.assign(u);

			A.assign_one();
			cur_qo->prin_list.hash(A);
		}
	}

	BB.assign_one();
	D.assign(G);

	if (!done) {
		// oa = (D - b^2) / 4a
		square(oa, A.b);
		subtract(oa, Delta, oa);
		divide(oa, oa, A.a);
		shift_right(oa, oa, 2);
	}

	while (!done) {
		// compute more baby steps
		flag = A.make_list(oa, u, G, U, cur_qo->prin_list);
		if (flag == 1) {
			done = is_prin = true;
			dist.assign(A.get_distance());
		}
		else if (flag == 0) {
			is_prin = false;
			done = true;
			dist.assign(A.get_distance());
		}

		// compute giant steps to u^2
		if (!done) {
			inverse(C, A);
			C.adjust_abs(u);
			square(usqr, u);
		}
		while ((y.compare(usqr) < 0) && (dist.is_zero())) {
			multiply_real(D, D, C);
			D.adjust_abs(y);

			SVAL = cur_qo->prin_list.search(D);
			if (SVAL) {
				// found D in list:  z = y-r
				subtract(dist, D.get_distance(), SVAL->get_distance());
				dist.absolute_value();
				done = is_prin = true;
			}
			else {
				multiply_real(BB, BB, C);
				BB.adjust_abs(y);
				SVAL = cur_qo->prin_list.search(BB);
				if ((SVAL) && (!C.is_one() && !BB.is_one())) {
					// found b in list:  z = y-r
					subtract(dist, BB.get_distance(), SVAL->get_distance());
					dist.absolute_value();
					is_prin = false;
					done = true;
				}
				else
					add(y, y, u);
			}
		}

		u.multiply_by_2();
	}

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (is_principal_buch) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	if (do_verify) {
		if (verify_principal(dist, is_prin)) {
			if (info)
				std::cout << "Principality test is correct." << std::endl;
		}
		else
			std::cout << "Principality test is not correct!" << std::endl;
	}

	if (quadratic_order::special2) {
		std::cout << "0" << std::endl;
		std::cout << "0" << std::endl;
	}

	bigfloat::set_precision(oprec);
	return is_prin;
}



//
// qi_class_real::is_principal_subexp(bigfloat)
//
// Task:
//      tests if the equivalence class is the principal class using a variation
//      of the subexponential algorithm due to Abel.  The distance from the
//      unit ideal to the representative is also computed.
//

bool
qi_class_real::is_principal_subexp (bigfloat & dist) const
{
	debug_handler("qi_class_real", "is_principal_subexp(bigfloat)");

	quadratic_order *QO;
	int num;
	bool is_prin;
	timer t;
	long oprec = bigfloat::get_precision();

	if (Delta.is_le_zero()) {
		lidia_error_handler("qi_class_real", "is_principal_subexp(bigfloat) - no "
				    "current quadratic order has been defined, or it is not real");
		return false;
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class_real", "is_principal_subexp(bigfloat) - the "
				    "current quadratic order has been deleted");
		return false;
	}

	QO = current_order->get_qo();
	bigfloat::set_precision(QO->prec);

	if (info)
		t.start_timer();

	num = decimal_length(Delta);
	if ((!QO->is_R_computed() && (num < RSHANKSB)) ||
	    (QO->prin_list.no_of_elements() > 0)) {
		is_prin = is_principal_buch(dist);
	}
	else {
		QO->regulator(true);
		is_prin = QO->is_in_lattice(*this, dist);
	}

	while (dist.is_lt_zero())
		add(dist, dist, QO->regulator());

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (is_principal_subexp) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	if (do_verify) {
		if (verify_principal(dist, is_prin)) {
			if (info)
				std::cout << "Principality test is correct." << std::endl;
		}
		else
			std::cout << "Principality test is not correct!" << std::endl;
	}

	bigfloat::set_precision(oprec);
	return is_prin;
}



//
// qi_class_real::verify_principal()
//
// Task:
//      verifies whether x is the distance from (1) to the qi_class_real, or
//      the regulator.
//

bool
qi_class_real::verify_principal (const bigfloat & x, const bool is_prin) const
{
	debug_handler("qi_class_real", "verify_principal");

	qi_class_real A, C, U;
	bool OK;
	timer t;
	quadratic_order *cur_qo;
	long oprec = bigfloat::get_precision();

	if (Delta.is_le_zero()) {
		lidia_error_handler("qi_class_real", "verify_principal - no "
				    "current quadratic order has been defined, or it is not real");
		return false;
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class_real", "verify_principal - the "
				    "current quadratic order has been deleted");
		return false;
	}

	cur_qo = current_order->get_qo();
	bigfloat::set_precision(cur_qo->prec);

	A.assign(*this, 0);

	if (info)
		t.start_timer();

	U.assign_one();
	C.assign(nearest(U, x));
	if (info > 1) {
		std::cout << "\nVerifying principality test:" << std::endl;
		std::cout << "  A = " << A << std::endl;
		std::cout << "  x = " << x << std::endl;
		std::cout << "  nearest from (1) = " << C << std::endl;
	}

	if (is_prin)
		OK = C.is_equal(A);
	else
		OK = C.is_one();

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (verify_principal) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	bigfloat::set_precision(oprec);
	return OK;
}



//
// qi_class_real::is_equivalent()
//
// Task:
//      tests if the equivalent classes are equal (reduced representatives are
//      equivalent).  The algorithm is chosen based on the size of the
//      discriminant.
//

bool
qi_class_real::is_equivalent (const qi_class_real & B) const
{
	debug_handler("qi_class_real", "is_equivalent()");

	qi_class_real A;
	bool is_equiv;
	bigfloat dist;

	if (do_verify)
		is_equiv = is_equivalent(B, dist);
	else {
		multiply(A, *this, inverse(B));
		A.assign(A, 0);
		is_equiv = A.is_principal();
	}

	return is_equiv;
}



//
// qi_class_real::is_equivalent_buch()
//
// Task:
//      tests if the equivalent classes are equal (reduced representatives are
//      equivalent) using a variation of baby-step giant-step due to Buchmann.
//

bool
qi_class_real::is_equivalent_buch (const qi_class_real & B) const
{
	debug_handler("qi_class_real", "is_equivalent_buch()");

	qi_class_real A;
	bigfloat dist;
	bool is_equiv;

	if (do_verify)
		is_equiv = is_equivalent_buch(B, dist);
	else {
		multiply(A, *this, inverse(B));
		A.assign(A, 0);
		is_equiv = A.is_principal_buch();
	}

	return is_equiv;
}



//
// qi_class_real::is_equivalent_subexp()
//
// Task:
//      tests if the equivalent classes are equal (reduced representatives are
//      equivalent) using a variation of the subexponential algorithm due to
//      Abel.
//

bool
qi_class_real::is_equivalent_subexp (const qi_class_real & B) const
{
	debug_handler("qi_class_real", "is_equivalent_subexp()");

	qi_class_real A;
	bool is_equiv;
	bigfloat dist;

	if (do_verify)
		is_equiv = is_equivalent_subexp(B, dist);
	else {
		multiply(A, *this, inverse(B));
		A.assign(A, 0);
		is_equiv = A.is_principal_subexp();
	}

	return is_equiv;
}



//
// qi_class_real::is_equivalent(bigfloat)
//
// Task:
//      tests if the equivalent classes are equal (reduced representatives are
//      equivalent).  The algorithm is chosen based on the size of the
//      discriminant.  The distance between the two representatives is also
//      computed.
//

bool
qi_class_real::is_equivalent (const qi_class_real & B, bigfloat & dist) const
{
	debug_handler("qi_class_real", "is_equivalent(bigfloat)");

	bool is_equiv;
	int num;

	num = decimal_length(Delta);
	if (num < RSHANKSB)
		is_equiv = is_equivalent_buch(B, dist);
	else
		is_equiv = is_equivalent_subexp(B, dist);

	return is_equiv;
}



//
// qi_class_real::is_equivalent_buch(bigfloat)
//
// Task:
//      tests if the equivalent classes are equal (reduced representatives are
//      equivalent) using a variation of baby-step giant-step due to Buchmann.
//      The distance between the two representatives is also computed.
//

bool
qi_class_real::is_equivalent_buch (const qi_class_real & B,
				   bigfloat & dist) const
{
	debug_handler("qi_class_real", "is_equivalent_buch(bigfloat)");

	qi_class_real A, BB, C, U;
	bool is_equiv;
	timer t;
	int tv, nrho;
	quadratic_order *cur_qo;
	long oprec = bigfloat::get_precision();

	if (Delta.is_le_zero()) {
		lidia_error_handler("qi_class_real", "is_equivalent_buch(bigfloat) - no "
				    "current quadratic order has been defined, or it is not real");
		return false;
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class_real", "is_equivalent_buch(bigfloat) - the "
				    "current quadratic order has been deleted");
		return false;
	}

	cur_qo = current_order->get_qo();
	bigfloat::set_precision(cur_qo->prec);

	BB.assign(B, 0);

	if (info)
		t.start_timer();

	tv = do_verify;
	do_verify = 0;

	C.assign(*this);
	multiply(A, C, inverse(BB));
	nrho = 0;
	while (A.is_one()) {
		++nrho;
		C.rho();
		multiply(A, C, inverse(BB));
	}
	A.d.assign_zero();

	is_equiv = A.is_principal_buch(dist);

	if (is_equiv) {
		U.assign_one();
		U = nearest(U, dist);
		multiply(U, U, BB);
		dist = U.get_distance();
		while (dist < 0)
			add(dist, dist, cur_qo->regulator());
		while (dist > cur_qo->regulator())
			subtract(dist, dist, cur_qo->regulator());
		A.assign(C, 0);
	}
	else
		A.assign(B, 0);
	dist = A.convert_distance(BB, dist);

	if (is_equiv && nrho) {
		C.assign(nearest(B, dist));
		while (nrho) {
			--nrho;
			C.inverse_rho();
		}
		dist.assign(C.get_distance());
	}

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (is_equivalent_buch) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	do_verify = tv;
	if (do_verify) {
		if (verify_equivalent(B, dist, is_equiv)) {
			if (info)
				std::cout << "Equivalence test is correct." << std::endl;
		}
		else
			std::cout << "Equivalence test is not correct!" << std::endl;
	}

	bigfloat::set_precision(oprec);
	return is_equiv;
}



//
// qi_class_real::is_equivalent_subexp(bigfloat)
//
// Task:
//      tests if the equivalent classes are equal (reduced representatives are
//      equivalent) using a variation of the subexponential algorithm due to
//      Abel.  The distance between the two representatives is also computed.
//

bool
qi_class_real::is_equivalent_subexp (const qi_class_real & B,
				     bigfloat & dist) const
{
	debug_handler("qi_class_real", "is_equivalent_subexp(bigfloat)");

	qi_class_real A, BB, C, U;
	bool is_equiv;
	timer t;
	int tv, nrho;
	quadratic_order *cur_qo;
	long oprec = bigfloat::get_precision();

	if (Delta.is_le_zero()) {
		lidia_error_handler("qi_class_real", "is_equivalent_subexp(bigfloat) - no "
				    "current quadratic order has been defined, or it is not real");
		return false;
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class_real", "is_equivalent_subexp(bigfloat) - the "
				    "current quadratic order has been deleted");
		return false;
	}

	cur_qo = current_order->get_qo();
	bigfloat::set_precision(cur_qo->prec);

	if (info)
		t.start_timer();

	tv = do_verify;
	do_verify = 0;

	BB.assign(B, 0);

	C.assign(*this);
	multiply(A, C, inverse(BB));
	nrho = 0;
	while (A.is_one()) {
		++nrho;
		C.rho();
		multiply(A, C, inverse(BB));
	}
	A.d.assign_zero();

	is_equiv = A.is_principal_subexp(dist);

	if (is_equiv) {
		U.assign_one();
		U = nearest(U, dist);
		multiply(U, U, BB);
		dist = U.get_distance();
		while (dist < 0)
			add(dist, dist, cur_qo->regulator());
		while (dist > cur_qo->regulator())
			subtract(dist, dist, cur_qo->regulator());
		A.assign(C, 0);
	}
	else
		A.assign(B, 0);
	dist = A.convert_distance(BB, dist);

	if (is_equiv && nrho) {
		C.assign(nearest(B, dist));
		while (nrho) {
			--nrho;
			C.inverse_rho();
		}
		dist.assign(C.get_distance());
	}

	while (dist.is_lt_zero())
		add(dist, dist, cur_qo->regulator());

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (is_equivalent_subexp) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	do_verify = tv;
	if (do_verify) {
		if (verify_equivalent(B, dist, is_equiv)) {
			if (info)
				std::cout << "Equivalence test is correct." << std::endl;
		}
		else
			std::cout << "Equivalence test is not correct!" << std::endl;
	}

	bigfloat::set_precision(oprec);
	return is_equiv;
}



//
// qi_class_real::verify_equivalent()
//
// Task:
//      verifies whether x is the distance from B to the qi_class_real, or the
//      regulator.
//

bool
qi_class_real::verify_equivalent (const qi_class_real & B, const bigfloat & x,
				  const bool is_equiv) const
{
	debug_handler("qi_class_real", "verify_equivalent");

	qi_class_real AA, BB, C;
	bool OK;
	timer t;
	quadratic_order *cur_qo;
	long oprec = bigfloat::get_precision();

	if (Delta.is_le_zero()) {
		lidia_error_handler("qi_class_real", "verify_equivalent - no current "
				    "quadratic order has been defined, or it is not real");
		return false;
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class_real", "verify_equivalent - the current "
				    "quadratic order has been deleted");
		return false;
	}

	cur_qo = current_order->get_qo();
	bigfloat::set_precision(cur_qo->prec);

	AA.assign(*this, 0);
	BB.assign(B, 0);

	if (info)
		t.start_timer();

	C.assign(nearest(BB, x));

	if (info > 1) {
		std::cout << "\nVerifying equivalence test:" << std::endl;
		std::cout << "  B = " << BB << std::endl;
		std::cout << "  A = " << AA << std::endl;
		std::cout << "  x = " << x << std::endl;
		std::cout << "  nearest from B = " << C << std::endl;
	}

	if (is_equiv)
		OK = C.is_equal(AA);
	else
		OK = C.is_equal(BB);

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (verify_equivalent) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	bigfloat::set_precision(oprec);
	return OK;
}



//
// qi_class_real::DL()
//
// Task:
//      computes the smallest integer x such that G^x is equivalent to the
//      ideal (if x exists) and the distance from G^x to the ideal.  If x
//      exists, true is returned, otherwise false is returned, x is set to
//      the order of the equivalence class containing G and dist is set to
//      the distance of G^x from (1).  The algorithm uses the DL function of
//      the qi_class class to compute x, and is_principal or is_equivalent
//      to compute dist.
//

bool
qi_class_real::DL (const qi_class_real & G, bigint & x, bigfloat & dist) const
{
	debug_handler("qi_class_real", "DL");

	qi_class A, GG;
	qi_class_real C;
	int tv;
	bool isDL;
	timer t;
	quadratic_order *cur_qo;
	long oprec = bigfloat::get_precision();

	if (Delta.is_le_zero()) {
		lidia_error_handler("qi_class_real", "DL - no "
				    "current quadratic order has been defined, or it is not real");
		return false;
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class_real", "DL - the "
				    "current quadratic order has been deleted");
		return false;
	}

	cur_qo = current_order->get_qo();
	bigfloat::set_precision(cur_qo->prec);


	// FOR TABLE GENERATION
	timer tspecial;
	long t1 = 0;
	bool old = quadratic_order::special;
	quadratic_order::special = false;
	if (old)
		tspecial.start_timer();

	tv = do_verify;
	do_verify = 0;

	if (info)
		t.start_timer();

	// compute x
	A.assign(*this);
	GG.assign(G);
	isDL = A.DL(GG, x);

	// compute distance
	power(C, G, x);
	C.d.assign_zero();

	if (!isDL)
		C.is_principal(dist);
	else
		is_equivalent(C, dist);

	// FOR TABLE GENERATION
	quadratic_order::special = old;
	if (quadratic_order::special) {
		tspecial.stop_timer();
		t1 = tspecial.user_time();
		// ROW:  time for representations    total time
		std::cout << t1 << std::endl;
	}

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (DL) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	do_verify = tv;
	if (do_verify) {
		if (verify_DL(G, x, dist, isDL)) {
			if (info)
				std::cout << "DL is correct." << std::endl;
		}
		else
			std::cout << "DL is not correct!" << std::endl;
	}

	bigfloat::set_precision(oprec);
	return isDL;
}



//
// qi_class_real::verify_DL()
//
// Task:
//      verifies whether x is either the DL of the qi_class to the base G, or
//      the order of G.
//

bool
qi_class_real::verify_DL (const qi_class_real & G, const bigint & x,
			  const bigfloat & dist, const bool is_DL) const
{
	debug_handler("qi_class_real", "verify_DL");

	qi_class_real AA, GG, C;
	bool OK;
	timer t;
	quadratic_order *cur_qo;
	long oprec = bigfloat::get_precision();

	if (Delta.is_le_zero()) {
		lidia_error_handler("qi_class_real", "verify_principal - no "
				    "current quadratic order has been defined, or it is not real");
		return false;
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class_real", "verify_principal - the "
				    "current quadratic order has been deleted");
		return false;
	}

	cur_qo = current_order->get_qo();
	bigfloat::set_precision(cur_qo->prec);

	AA.assign(*this, 0);
	GG.assign(G, 0);

	if (info)
		t.start_timer();


	power(C, G, x);
	C.d.assign_zero();

	if (info > 1) {
		std::cout << "\nVerifying Buchmann-Paulus DL computation:" << std::endl;
		std::cout << "  G = " << GG << std::endl;
		std::cout << "  A = " << AA << std::endl;
		std::cout << "  x = " << x << std::endl;
		std::cout << "  G^x = " << C << std::endl;
		std::cout << "  dist(G^x, A) = " << dist << std::endl;
		if (is_DL)
			std::cout << "  x is log_G of A." << std::endl;
		else
			std::cout << "  x is ord(G)." << std::endl;
	}

	if (is_DL)
		OK = verify_equivalent(C, dist, is_DL);
	else
		OK = C.verify_principal(dist, true);

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (verify_DL) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	bigfloat::set_precision(oprec);
	return OK;
}



//
// subgroup()
//
// Task:
//      computes the structure of the subgroup generated by the ideals in G.
//      The qi_class function of the same name is used to carry out the
//      computation.
//

base_vector< bigint >
subgroup (base_vector< qi_class_real > & G)
{
	debug_handler("qi_class_real", "subgroup");

	lidia_size_t i, l;
	base_vector< qi_class > A;
	base_vector< bigint > S;
	qi_class Amem;

	A.set_mode(EXPAND);
	S.set_mode(EXPAND);
	A.reset();
	S.reset();

	l = G.size();
	A.set_size(l);
	for (i = 0; i < l; ++i)
		A[i] = qi_class(G[i]);

	S = subgroup(A);

	return S;
}



//
// subgroup_BJT()
//
// Task:
//      computes the structure of the subgroup generated by G using a
//      variation of baby-step giant-step due to Buchmann, Jacobson, and Teske
//      which has complexity O(sqrt(H)), where H is the order of the subgroup.
//      The qi_class function of the same name is used to carry out the
//      computation.
//

base_vector< bigint >
subgroup_BJT (base_vector< qi_class_real > & G, base_vector< long > & v)
{
	debug_handler("qi_class_real", "subgroup_BJT");

	lidia_size_t i, l;
	base_vector< qi_class > A;
	base_vector< bigint > S;
	qi_class Amem;

	A.set_mode(EXPAND);
	S.set_mode(EXPAND);
	A.reset();
	S.reset();

	l = G.size();
	A.set_size(l);
	for (i = 0; i < l; ++i)
		A[i] = qi_class(G[i]);

	S = subgroup_BJT(A, v);

	return S;
}



//
// subgroup_h()
//
// Task:
//      computes the structure of the subgroup generated by G assuming that
//      the class number is known.  The class number is computed if it has
//      not been already.  The qi_class function of the same name is used to
//      carry out the computation.
//

base_vector< bigint >
subgroup_h (base_vector< qi_class_real > & G)
{
	debug_handler("qi_class_real", "subgroup_h");

	lidia_size_t i, l;
	base_vector< qi_class > A;
	base_vector< bigint > S;
	qi_class Amem;

	A.set_mode(EXPAND);
	S.set_mode(EXPAND);
	A.reset();
	S.reset();

	l = G.size();
	A.set_size(l);
	for (i = 0; i < l; ++i)
		A[i] = qi_class(G[i]);

	S = subgroup_h(A);

	return S;
}



//
// operator >>
//
// Task:
//      inputs a qi_class_real from the std::istream in.
//

std::istream &
operator >> (std::istream & in, qi_class_real & A)
{
	debug_handler("qi_class_real", "operator >>");

	int n = 0;
	char c;
	bigint ibuf[2];
	bigfloat rbuf;

	long oprec = bigfloat::get_precision();
	bigfloat::set_precision(qi_class_real::current_order->get_qo()->prec);

	in >> c;
	if (c != '(') {
		lidia_error_handler("qi_class_real", "operator >>::(expected");
		return in;
	}

	in >> c;
	while (c != ')' && n != 3) {
		in.putback(c);
		if (n < 2)
			in >> ibuf[n];
		else
			in >> rbuf;
		n++;
		in >> c;
		if (c == ',')
			in >> c;
	}
	A.assign(ibuf[0], ibuf[1], rbuf);

	bigfloat::set_precision(oprec);
	return in;
}



//
// operator <<
//
// Task:
//      outputs a qi_class_real to the std::ostream out.
//

std::ostream &
operator << (std::ostream & out, const qi_class_real & A)
{
	debug_handler("qi_class_real", "operator << ");

	long oprec = bigfloat::get_precision();
	bigfloat::set_precision(qi_class_real::current_order->get_qo()->prec);

	out << "(" << A.a << ", " << A.b << ", " << A.get_distance() << ")";

	bigfloat::set_precision(oprec);
	return out;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
