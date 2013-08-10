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
#include	"LiDIA/qi_class.h"
#include	"LiDIA/number_fields/qo_list.h"
#include	"LiDIA/matrix.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//<TPf>
#if 0
qo_node *qi_class::current_order = NULL;
bigint qi_class::Delta = 0;
bigint qi_class::rootD = 0;
xbigfloat qi_class::rd = 0.0;
bigint qi_class::PEA_L = 0;
bigint qi_class::ordmult = 0;
rational_factorization qi_class::omfact;
int qi_class::info = 0;
int qi_class::do_verify = 0;
long qi_class::xprec = 0;
#endif
//</TPf>


bigint
qi_class_key (const qi_class & G)
{
	debug_handler("qi_class", "qi_class_key(const qi_class&)");
	return G.get_a();
}



//
// qi_class::conductor()
//
// Task:
//      returns the conductor of the ideal.
//

inline bigint
qi_class::conductor () const
{
	debug_handler("qi_class", "conductor");

	bigint temp;

	temp = gcd(a, b);
	return gcd(temp, get_c());
}



//
// qi_class::is_invertible()
//
// Task:
//      returns true if the ideal is invertible.
//

bool
qi_class::is_invertible () const
{
	debug_handler("qi_class", "is_invertible");

	bigint f;

	f = conductor();
	return (f.is_one());
}



//
// qi_class::normalize_imag()
//
// Task:
//      normalizes the imaginary ideal.
//

inline void
qi_class::normalize_imag ()
{
	debug_handler("qi_class", "normalize_imag");

	bigint a2, q, r, temp;

	if (!((b.compare(-a) > 0) && (b.compare(a) <= 0))) {
		// q = floor((a-b) / 2a)
		shift_left(a2, a, 1);
		subtract(temp, a, b);
		div_rem(q, r, temp, a2);
		if ((temp.is_lt_zero()) && (!r.is_zero()))
			dec(q);

		multiply(temp, a2, q);
		add(b, b, temp);
	}
}



//
// qi_class::normalize_real()
//
// Task:
//      normalizes the real ideal.
//

inline void
qi_class::normalize_real ()
{
	debug_handler("qi_class", "normalize_real");

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
		shift_left(a2, a, 1);
		if (a <= rootD) {
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
// qi_class::normalize()
//
// Task:
//      normalizes the ideal.
//

void
qi_class::normalize ()
{
	debug_handler("qi_class", "normalize");

	if (Delta.is_lt_zero())
		normalize_imag();
	else
		normalize_real();
}



//
// qi_class::reduce_imag()
//
// Task:
//      reduces the imaginary ideal.
//

inline void
qi_class::reduce_imag ()
{
	debug_handler("qi_class", "reduce_imag");

	bigint na, oa, a2, t, nt, q, m, r, temp;
	int st, nst;

	normalize_imag();

	t = abs(b);
	if (t.is_zero())
		st = 1;
	else
		st = b.sign();

	shift_left(a2, a, 1);
	div_rem(q, r, t, a2);
	subtract(m, a, r);
	if (m.is_ge_zero()) {
		nt.assign(r);
		nst = -st;
	}
	else {
		add(nt, r, m);
		nst = st;
	}

	// na = (b^2 - D) / 4a
	square(na, b);
	subtract(na, na, Delta);
	divide(na, na, a);
	shift_right(na, na, 2);

	while (a.compare(na) > 0) {
		oa.assign(a);
		a.assign(na);
		t.assign(nt);
		st = nst;

		shift_left(a2, a, 1);
		div_rem(q, r, t, a2);
		subtract(m, a, r);

		// na = oa - q*((r + t) >> 1);
		add(temp, r, t);
		temp.divide_by_2();
		multiply(temp, temp, q);
		subtract(na, oa, temp);

		if (m.is_ge_zero()) {
			nt.assign(r);
			nst = -st;
		}
		else {
			shift_left(nt, m, 1);
			add(nt, nt, r);
			add(na, na, m);
			nst = st;
		}
	}

	b.assign(t);
	if (st < 0)
		b.negate();

	normalize_imag();

	if ((a == na) && (b.is_lt_zero()))
		b.negate();
}



//
// qi_class::reduce_real()
//
// Task:
//      reduces the real ideal.
//

inline void
qi_class::reduce_real ()
{
	debug_handler("qi_class", "reduce_real");

	bigint q, r, a2, oa, na, nb, temp;
	int s;

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

		while (!((abs(temp) < b) && (b.compare(rootD) <= 0))) {
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

			b.assign(nb);
			oa.assign(a);
			a.assign(na);

			shift_left(a2, abs(a), 1);
			subtract(temp, rootD, a2);
			if (temp.is_lt_zero())
				inc(temp);
		}

		a.absolute_value();
	}
}



//
// qi_class::reduce_real(quadratic_number)
//
// Task:
//      reduces the ideal and updates its associated quadratic number.
//

inline void
qi_class::reduce_real (quadratic_number_standard & alpha)
{
	debug_handler("qi_class", "reduce_real");

	bigint q, r, a2, na, oa, nb, OB, BB, NB, temp, beta_x;
	int s;
	quadratic_number_standard beta;

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

		// beta = (2*OB*a + BB*b + BB*rd) / 2
		multiply(a2, OB, a);
		a2.multiply_by_2();
		multiply(beta_x, BB, b);
		add(beta_x, beta_x, a2);
		//<MM>
		beta.assign_order(alpha.get_order());
		//</MM>
		beta.assign(beta_x, BB, 2);

		multiply(alpha, alpha, beta);

		a.absolute_value();
	}
}



//
// qi_class::reduce()
//
// Task:
//      reduces the ideal
//

void
qi_class::reduce ()
{
	debug_handler("qi_class", "reduce");

	if (Delta.is_lt_zero())
		reduce_imag();
	else
		reduce_real();
}



//
// constructor
//

qi_class::qi_class ()
{
	debug_handler("qi_class", "qi_class()");

	if (!current_order)
		current_order = quadratic_order::qo_l().add_last_to_list(current_order);
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

qi_class::qi_class (const bigint & a2, const bigint & b2)
{
	debug_handler("qi_class", "qi_class(bigint, bigint)");

	bigint c, temp;

	if (!current_order)
		current_order = quadratic_order::qo_l().add_last_to_list(current_order);

	a.assign_zero();
	b.assign_zero();
	if ((!a2.is_zero()) && (!Delta.is_zero())) {
		square(c, b2);
		subtract(c, c, Delta);
		shift_left(temp, a2, 2);
		remainder(c, c, temp);
		if (c.is_zero()) {
			a.assign(a2);
			b.assign(b2);
			reduce();

			if (!is_invertible()) {
				a.assign_zero();
				b.assign_zero();
			}
		}
	}
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

qi_class::qi_class (const long a2, const long b2)
{
	debug_handler("qi_class", "qi_class(long, long)");

	bigint c, temp;

	if (!current_order)
		current_order = quadratic_order::qo_l().add_last_to_list(current_order);

	a.assign_zero();
	b.assign_zero();
	if ((a2) && (!Delta.is_zero())) {
		square(c, b2);
		subtract(c, c, Delta);
		shift_left(temp, bigint(a2), 2);
		remainder(c, c, temp);
		if (c.is_zero()) {
			a.assign(a2);
			b.assign(b2);
			reduce();

			if (!is_invertible()) {
				a.assign_zero();
				b.assign_zero();
			}
		}
	}
}



//
// constructor:
//    if qf is regular and primitive, initialize with qf reduced.
//

qi_class::qi_class (const quadratic_form & qf)
{
	debug_handler("qi_class", "qi_class(quadratic_form)");

	if ((!qf.is_regular()) || (!qf.is_primitive())) {
		a.assign_zero();
		b.assign_zero();
	}
	else {
		set_current_order(qf.which_order());
		a.assign(abs(qf.get_a()));
		b.assign(qf.get_b());
		reduce();
	}
}



//
// constructor:
//    if A is invertible, initialize with A reduced.
//

qi_class::qi_class (const quadratic_ideal & A)
{
	debug_handler("qi_class", "qi_class(quadratic_ideal)");

	if (!current_order)
		current_order = quadratic_order::qo_l().add_last_to_list(current_order);

	if (!A.is_invertible()) {
		a.assign_zero();
		b.assign_zero();
	}
	else {
		set_current_order((quadratic_order&)A.which_order());
		a.assign(A.get_a());
		b.assign(A.get_b());
		reduce();
	}
}



//
// constructor:
//    initialize with a copy of A.
//

qi_class::qi_class (const qi_class_real & A)
{
	debug_handler("qi_class", "qi_class(qi_class_real)");

	if (!current_order)
		current_order = quadratic_order::qo_l().add_last_to_list(current_order);

	a.assign(A.get_a());
	b.assign(A.get_b());
}



//
// constructor:
//    initialize with a copy of A.
//

qi_class::qi_class (const qi_class & A)
{
	debug_handler("qi_class", "qi_class(qi_class)");

	if (!current_order)
		current_order = quadratic_order::qo_l().add_last_to_list(current_order);

	a.assign(A.a);
	b.assign(A.b);
}



//
// destructor
//

qi_class::~qi_class()
{
	debug_handler("qi_class", "~qi_class");
}



//
// qi_class::set_current_order()
//
// Task:
//      sets the current quadratic order and the other static variables
//

void
qi_class::set_current_order (quadratic_order & QO)
{
	debug_handler("qi_class", "set_current_order");

	bigfloat temp;

	current_order = quadratic_order::qo_l().add_to_list(current_order, QO);

	Delta.assign(QO.discriminant());
	PEA_L.assign(QO.nu_bound());
	if (!Delta.is_zero()) {
		sqrt(rd, xbigfloat(abs(Delta)), QO.xprec + b_value(Delta) + 11);
		floor(rootD, rd);
	}
	else {
		rd.assign_zero();
		rootD.assign_zero();
	}
	ordmult.assign_zero();
	xprec = QO.xprec;
}



//
// qi_class::verbose()
//
// Task:
//      sets the verbosity of commands.  Currently, the following levels are
//      supported:
//         0 - nothing
//         > 0 - run times
//

void
qi_class::verbose (int state)
{
	debug_handler("qi_class", "verbose");

	if (state <= 0)
		info = 0;
	else
		info = state;
}



//
// qi_class::verification()
//
// Task:
//      sets the level of verifications.  Currently, the following levels are
//      supported:
//         0 - nothing
//         1 - computation is verified (true or false)
//

void
qi_class::verification (int level)
{
	debug_handler("qi_class", "verification");

	if (level <= 0)
		do_verify = 0;
	else
		do_verify = level;
}



//
// qi_class::assign_zero()
//
// Task:
//      set to the zero ideal
//

void
qi_class::assign_zero ()
{
	debug_handler("qi_class", "assign_zero");

	a.assign_zero();
	b.assign_zero();
}



//
// qi_class::assign_one()
//
// Task:
//      set to the unit ideal of the current quadratic_order
//

void
qi_class::assign_one ()
{
	debug_handler("qi_class", "assign_one");

	if (Delta.is_zero()) {
		lidia_error_handler("qi_class", "assign_one - no current quadratic "
				    "order has been defined");
		return;
	}

	a.assign_one();
	b.assign(Delta);
	normalize();
}



//
// qi_class::assign_principal()
//
// Task:
//    set to the reduced principal ideal generated by the element
//
//    (x + (D + sqrt(Delta))/2 y)
//

void
qi_class::assign_principal (const bigint & x, const bigint & y)
{
	debug_handler("qi_class", "assign_principal");

	bigint x2, y2, n, m, k, l, temp, temp2;

	if (Delta.is_zero()) {
		lidia_error_handler("qi_class", "assign_principal - no current quadratic "
				    "order has been defined");
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

	// a = n / m^2
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

	reduce();
}



//
// qi_class::assign(bigint, bigint)
//
// Task:
//    set to the ideal belonging to the current order given by
//
//    [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    If the parameters do not form an ideal, false will be returned and
//    the ideal will not be modified.
//

bool
qi_class::assign (const bigint & a2, const bigint & b2)
{
	debug_handler("qi_class", "assign(bigint, bigint)");

	bigint c, temp;

	a.assign_zero();
	b.assign_zero();

	if (Delta.is_zero()) {
		lidia_error_handler("qi_class", "assign(bigint, bigint) - no current "
				    "quadratic order has been defined");
		return false;
	}

	if (a2.is_zero())
		if (b2.is_zero())
			return true;
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
		return false;
	}
}



//
// qi_class::assign(long, long)
//
// Task:
//    set to the ideal belonging to the current order given by
//
//    [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    If the parameters do not form an ideal, false will be returned and
//    the ideal will not be modified.
//

bool
qi_class::assign (const long a2, const long b2)
{
	debug_handler("qi_class", "assign(long, long)");

	return assign(bigint(a2), bigint(b2));
}



//
// qi_class::assign(quadratic_form)
//
// Task:
//      set to the reduction of qf.  If qf is not regular or not primitive,
//      the error handler will be evoked.
//

void
qi_class::assign (const quadratic_form & qf)
{
	debug_handler("qi_class", "assign(quadratic_form)");

	if (!qf.is_regular()) {
		lidia_error_handler("qi_class", "assign(quadratic_form) - the quadratic "
				    "form must be regular");
		return;
	}

	if (!qf.is_primitive()) {
		lidia_error_handler("qi_class", "assign(quadratic_form) - the quadratic "
				    "form is not primitive");
		return;
	}

	set_current_order(qf.which_order());
	a.assign(abs(qf.get_a()));
	b.assign(qf.get_b());
	reduce();
}



//
// qi_class::assign(quadratic_ideal)
//
// Task:
//      set to the reduction of B.  If B is not invertible, the error handler
//      will be evoked.
//

void
qi_class::assign (const quadratic_ideal & B)
{
	debug_handler("qi_class", "assign(quadratic_ideal)");

	if (!B.is_invertible()) {
		lidia_error_handler("qi_class", "qi_class(quadratic_ideal) - the quadratic "
				    "ideal is not invertible");
		return;
	}

	set_current_order((quadratic_order&)B.which_order());
	a.assign(B.get_a());
	b.assign(B.get_b());
	reduce();
}



//
// qi_class::assign(qi_class_real)
//
// Task:
//      set to a copy of B.
//

void
qi_class::assign (const qi_class_real & B)
{
	debug_handler("qi_class", "assign(qi_class_real)");

	a.assign(B.a);
	b.assign(B.b);
}



//
// qi_class::assign(qi_class)
//
// Task:
//      set to a copy of B.
//

void
qi_class::assign (const qi_class & B)
{
	debug_handler("qi_class", "assign(qi_class)");

	a.assign(B.a);
	b.assign(B.b);
}



//
// operator =
//
// Task:
//      make a copy of an existing qi_class
//

qi_class & qi_class::operator = (const qi_class & A)
{
	debug_handler("qi_class", "operator = ");

	assign(A);
	return *this;
}



//
// qi_class::get_current_order_ptr()
//
// Task:
//      returns a pointer to the current quadratic_order
//

quadratic_order *
qi_class::get_current_order_ptr ()
{
	debug_handler("qi_class", "get_current_order_ptr");

	if (!current_order)
		return NULL;
	else
		return current_order->get_qo();
}



//
// qi_class::get_current_order_ptr()
//
// Task:
//      returns the current quadratic_order
//

quadratic_order &
qi_class::get_current_order ()
{
	debug_handler("qi_class", "get_current_order");

	quadratic_order *QO;

	QO = current_order->get_qo();
	if (!QO) {
		lidia_error_handler("qi_class", "get_current_order - the current "
                                    "quadratic order has been deleted");
		return quadratic_order::zero_QO;
	}

	return *QO;
}



//
// qi_class::discriminant()
//
// Task:
//      returns the discriminant of the current quadratic order
//

bigint
qi_class::discriminant ()
{
	debug_handler("qi_class", "discriminant");

	return Delta;
}



//
// qi_class::get_rootD()
//
// Task:
//      returns floor(sqrt(|Delta|))
//

bigint
qi_class::get_rootD ()
{
	debug_handler("qi_class", "get_rootD");

	return rootD;
}



//
// qi_class::get_rd()
//
// Task:
//      returns sqrt(|Delta|)
//

xbigfloat
qi_class::get_rd ()
{
	debug_handler("qi_class", "get_rd");

	return rd;
}



//
// qi_class::get_a()
//
// Task:
//      returns coefficient a
//

bigint
qi_class::get_a () const
{
	debug_handler("qi_class", "get_a() const");
	return a;
}



//
// qi_class::get_b()
//
// Task:
//      returns coefficient b
//

bigint
qi_class::get_b () const
{
	debug_handler("qi_class", "get_b() const");
	return b;
}



//
// qi_class::get_c()
//
// Task:
//      returns (b^2 - Delta) / 4a
//

bigint
qi_class::get_c () const
{
	debug_handler("qi_class", "get_c() const");

	bigint c, temp;
	if (a.is_zero())
		c.assign_zero();
	else {
		square(c, b);
		subtract(c, c, Delta);
		shift_left(temp, a, 2);
		divide(c, c, temp);
	}
	return c;
}



//
// multiply()
//
// Task:
//      computes the product of A and B (reduced ideal product)
//

void
multiply (qi_class & C, const qi_class & A, const qi_class & B)
{
	debug_handler("qi_class", "multiply");

	if (qi_class::Delta.is_lt_zero())
		multiply_imag(C, A, B);
	else
		multiply_real(C, A, B);
}



//
// multiply_imag()
//
// Task:
//      multiplies two imaginary ideal equivalence classes
//

void
multiply_imag (qi_class & C, const qi_class & A, const qi_class & B)
{
	debug_handler("qi_class", "multiply_imag");

	bigint newa, newb, dpr, v, d, w, ab2, temp;

	// solve dpr = v A.a + w B.a
	dpr.assign(xgcd_left(v, A.a, B.a));

	// C.b = v A.a (B.b - A.b)
	multiply(newb, v, A.a);
	subtract(temp, B.b, A.b);
	multiply(newb, newb, temp);

	// C.a = A.a B.a
	multiply(newa, A.a, B.a);

	if (!dpr.is_one()) {
		add(ab2, A.b, B.b);
		ab2.divide_by_2();
		d.assign(xgcd(v, w, dpr, ab2));

		// C.b = (C.b*v + w(Delta - A.b^2)/2) / d
		multiply(newb, newb, v);

		square(temp, A.b);
		subtract(temp, qi_class::Delta, temp);
		temp.divide_by_2();
		multiply(temp, temp, w);

		add(newb, newb, temp);
		divide(newb, newb, d);

		// C.a = C.a / (d^2)
		square(temp, d);
		divide(newa, newa, temp);

		dpr.assign(d);
	}

	C.a.assign(newa);
	add(newb, newb, A.b);
	shift_left(ab2, C.a, 1);
	remainder(C.b, newb, ab2);

	C.reduce_imag();
}



//
// nucomp()
//
// Task:
//      multiplies two imaginary ideal equivalence classes using the NUCOMP
//      algorithm of Shanks.
//

void
nucomp (qi_class & C, const qi_class & A, const qi_class & B)
{
	debug_handler("qi_class", "nucomp");

	bigint a1, b1, c1, a2, b2, c2, s, n, d1, d;
	bigint u, v, u1, l, aa;
	bigint v1, v3;
	bigint b, e, f, g, q;
	bigint temp;
	bool flag;

	if (A.a < B.a) {
		a1.assign(B.a);
		b1.assign(B.b);
		square(c1, b1);
		subtract(c1, c1, qi_class::Delta);
		shift_left(temp, a1, 2);
		divide(c1, c1, temp);

		a2.assign(A.a);
		b2.assign(A.b);
		square(c2, b2);
		subtract(c2, c2, qi_class::Delta);
		shift_left(temp, a2, 2);
		divide(c2, c2, temp);
	}
	else {
		a1.assign(A.a);
		b1.assign(A.b);
		square(c1, b1);
		subtract(c1, c1, qi_class::Delta);
		shift_left(temp, a1, 2);
		divide(c1, c1, temp);

		a2.assign(B.a);
		b2.assign(B.b);
		square(c2, b2);
		subtract(c2, c2, qi_class::Delta);
		shift_left(temp, a2, 2);
		divide(c2, c2, temp);
	}

	// initialize
	add(s, b1, b2);
	s.divide_by_2();
	subtract(n, b2, s);

	// first Euclidean pass -  solve d = u a1 + v a2
	d.assign(xgcd(u, v, a2, a1));
	if (d.is_one()) {
		multiply(aa, u, n);
		aa.negate();
		flag = false;
	}
	else {
		// second Euclidean pass
		d1.assign(xgcd_left(u1, s, d));

		if (!d1.is_one()) {
			divide(a1, a1, d1);
			divide(a2, a2, d1);
			divide(s, s, d1);
			divide(d, d, d1);
			flag = true;
		}
		else
			flag = false;

		// l = - u1 * ( c1 * u + c2 * v ) mod d */
		remainder(c1, c1, d);
		remainder(temp, c2, d);
		multiply(temp, v, temp);
		multiply(l, u, c1);
		add(l, l, temp);
		multiply(l, l, u1);
		l.negate();
		remainder(l, l, d);
		if (l.is_lt_zero())
			add(l, l, d);

		// aa = - u * ( n / d ) + l * ( a1 / d ) */
		divide(aa, n, d);
		multiply(aa, aa, u);
		divide(temp, a1, d);
		multiply(temp, temp, l);
		subtract(aa, temp, aa);
	}

	// partial reduction
	d.assign(a1);
	v3.assign(aa);
	nugcd(u, d, v1, v3, qi_class::PEA_L);

	if (u.is_zero()) {
		//  u = 0; d = a1; v1 = 1; v3 = aa
		// b = a2 * d + n * u ) / a1  --->    b = a2 (-)
		// e = s * d + c2 * u ) / a1  --->   e = s  (-)

		// f = ( b * v3 + n ) / d
		// b1 = 2 * b * v3 + n

		multiply(temp, a2, v3);
		add(f, temp, n);
		add(b1, temp, f);

		// q = 2 * e * v1 - s       --->   q = s (-)

		if (flag)
			multiply(s, s, d1);

		// a1 = d * b + e * u    --->   a1 = d * b
		multiply(a1, d, a2);

		// b1 = b1 + q
		add(b1, b1, s);
	}
	else {
		// u != 0

		// b = a2 * d + n * u ) / a1
		multiply(temp, a2, d);
		multiply(b, n, u);
		add(b, b, temp);
		divide(b, b, a1);

		// e = s * d + c2 * u ) / a1
		multiply(temp, s, d);
		multiply(e, c2, u);
		add(e, e, temp);
		divide(e, e, a1);

		// f = ( b * v3 + n ) / d
		// b1 = 2 * b * v3 + n
		multiply(q, b, v3);
		add(f, q, n);
		add(b1, q, f);

		// g = ( e * v1 - s ) / u
		// q = 2 * e * v1 - s
		multiply(q, e, v1);
		subtract(g, q, s);
		add(q, q, g);

		if (flag) {
			multiply(v1, v1, d1);
			multiply(u, u, d1);
			multiply(q, q, d1);
		}

		// a1 = d * b + e * u
		multiply(temp, d, b);
		multiply(a1, e, u);
		add(a1, a1, temp);

		// b1 = b1 + q
		add(b1, b1, q);
	}

	C.a.assign(a1);
	C.b.assign(b1);
	C.reduce();
}



//
// multiply_real()
//
// Task:
//      multiplies two real ideal equivalence classes
//

void
multiply_real (qi_class & C, const qi_class & A, const qi_class & B)
{
	debug_handler("qi_class", "multiply_real");

	bigint newa, newb, dpr, v, d, w, ab2, temp;

	// solve dpr = v A.a + w B.a
	dpr.assign(xgcd_left(v, A.a, B.a));

	// C.b = v A.a (B.b - A.b)
	multiply(newb, v, A.a);
	subtract(temp, B.b, A.b);
	multiply(newb, newb, temp);

	// C.a = A.a B.a
	multiply(newa, A.a, B.a);

	if (!dpr.is_one()) {
		add(ab2, A.b, B.b);
		ab2.divide_by_2();
		d.assign(xgcd(v, w, dpr, ab2));

		// C.b = (C.b*v + w(D - A.b^2)/2) / d
		multiply(newb, newb, v);

		square(temp, A.b);
		subtract(temp, qi_class::Delta, temp);
		temp.divide_by_2();
		multiply(temp, temp, w);

		add(newb, newb, temp);
		divide(newb, newb, d);

		// C.a = C.a / (d^2)
		square(temp, d);
		divide(newa, newa, temp);

		dpr.assign(d);
	}

	C.a.assign(newa);
	add(newb, newb, A.b);
	shift_left(ab2, C.a, 1);
	remainder(C.b, newb, ab2);

	C.reduce_real();
}



//
// multiply_real(quadratic_number)
//
// Task:
//      multiplies two real ideal equivalence classes and computes the
//      associated quadratic number
//

void
multiply_real (qi_class & C, const qi_class & A, const qi_class & B,
	       quadratic_number_standard & Cq,
	       const quadratic_number_standard & Aq,
	       const quadratic_number_standard & Bq)
{
	debug_handler("qi_class", "multiply_real(quadratic_number)");

	bigint newa, newb, dpr, v, d, w, ab2, temp;

	// solve dpr = v A.a + w B.a
	dpr.assign(xgcd_left(v, A.a, B.a));

	// C.b = v A.a (B.b - A.b)
	multiply(newb, v, A.a);
	subtract(temp, B.b, A.b);
	multiply(newb, newb, temp);

	// C.a = A.a B.a
	multiply(newa, A.a, B.a);

	if (!dpr.is_one()) {
		add(ab2, A.b, B.b);
		ab2.divide_by_2();
		d.assign(xgcd(v, w, dpr, ab2));

		// C.b = (C.b*v + w(D - A.b^2)/2) / d
		multiply(newb, newb, v);

		square(temp, A.b);
		subtract(temp, qi_class::Delta, temp);
		temp.divide_by_2();
		multiply(temp, temp, w);

		add(newb, newb, temp);
		divide(newb, newb, d);

		// C.a = C.a / (d^2)
		square(temp, d);
		divide(newa, newa, temp);

		dpr.assign(d);
	}

	C.a.assign(newa);
	add(newb, newb, A.b);
	shift_left(ab2, C.a, 1);
	remainder(C.b, newb, ab2);

	multiply(Cq, Aq, Bq);
	divide(Cq, Cq, dpr);

	C.reduce_real(Cq);
}



//
// qi_class::invert()
//
// Task:
//      inverts the qi_class
//

void
qi_class::invert ()
{
	debug_handler("qi_class", "invert");

	if (!((Delta.is_lt_zero()) && (a == get_c()))) {
		b.negate();
		normalize();
	}
}



//
// inverse(qi_class, qi_class)
//
// Task:
//      computes the inverse of B
//

void
inverse (qi_class & A, const qi_class & B)
{
	debug_handler("qi_class", "inverse(qi_class, qi_class)");

	A.assign(B);
	if (!((qi_class::Delta.is_lt_zero()) && (A.a == A.get_c()))) {
		A.b.negate();
		A.normalize();
	}
}



//
// inverse(qi_class)
//
// Task:
//      returns the inverse of B
//

qi_class
inverse (const qi_class & B)
{
	debug_handler("qi_class", "inverse(qi_class)");

	qi_class A;

	A.assign(B);
	if (!((qi_class::Delta.is_lt_zero()) && (A.a == A.get_c()))) {
		A.b.negate();
		A.normalize();
	}

	return A;
}



//
// divide()
//
// Task:
//      computes A*B^-1.
//

void
divide (qi_class & C, const qi_class & A, const qi_class & B)
{
	debug_handler("qi_class", "divide");

	qi_class temp;

	if (B.is_zero()) {
		lidia_error_handler("qi_class", "divide() - zero divisor");
		return;
	}

	temp = inverse(B);
	multiply(C, A, temp);
}



//
// square
//
// Task:
//      multiplies A with itself
//

void
square (qi_class & C, const qi_class & A)
{
	debug_handler("qi_class", "square");

	if (qi_class::Delta.is_lt_zero())
		square_imag(C, A);
	else
		square_real(C, A);
}



//
// square_imag
//
// Task:
//      multiplies the imaginary equivalence class A with itself
//

void
square_imag (qi_class & C, const qi_class & A)
{
	debug_handler("qi_class", "square_imag");

	bigint newb, d, w, temp;

	// solve d = v A.a + w A.b
	d.assign(xgcd_right(w, A.a, A.b));

	// C.b = A.b + w (D - A.b^2) / (2d)
	square(newb, A.b);
	subtract(newb, qi_class::Delta, newb);
	shift_left(temp, d, 1);
	divide(newb, newb, temp);
	multiply(newb, newb, w);
	add(C.b, newb, A.b);

	// C.a = (A.a/d)^2
	divide(temp, A.a, d);
	square(C.a, temp);

	shift_left(temp, C.a, 1);
	remainder(C.b, C.b, temp);

	C.reduce_imag();
}



//
// nudupl
//
// Task:
//      multiplies the imaginary equivalence class A with itself using the
//      NUDUPL algorithm of Shanks.
//

void
nudupl (qi_class & C, const qi_class & A)
{
	debug_handler("qi_class", "nudupl");

	bigint a, b, c;
	bigint temp, u, d1, d, v1, v3, e, g;
	bool flag;

	a.assign(A.a);
	b.assign(A.b);
	square(c, b);
	subtract(c, c, qi_class::Delta);
	shift_left(temp, a, 2);
	divide(c, c, temp);

	// Euclidian pass
	d1.assign(xgcd_left(u, b, a));
	flag = (!d1.is_one());
	if (flag) {
		divide(a, a, d1);
		divide(b, b, d1);
	}

	multiply(v3, u, c);
	v3.negate();

	// partial euclidian algorithm
	d.assign(a);
	nugcd(u, d, v1, v3, qi_class::PEA_L);

	// final squaring
	if (u.is_zero()) {
		if (flag)
			multiply(b, b, d1);

		square(a, d);
		square(c, v3);

		add(temp, d, v3);
		square(temp, temp);
		add(b, b, temp);
		subtract(b, b, a);
		subtract(b, b, c);
	}
	else {
		// u != 0
		multiply(e, b, d);
		multiply(temp, c, u);
		add(e, e, temp);
		divide(e, e, a);

		multiply(temp, e, v1);
		subtract(g, temp, b);
		add(b, temp, g);

		if (flag) {
			multiply(b, b, d1);
			multiply(u, u, d1);
			multiply(v1, v1, d1);
		}

		square(a, d);
		square(c, v3);

		add(temp, d, v3);
		square(temp, temp);
		add(b, b, temp);
		subtract(b, b, a);
		subtract(b, b, c);

		multiply(temp, e, u);
		add(a, a, temp);
	}

	C.a.assign(a);
	C.b.assign(b);
	C.reduce();
}



//
// square_real()
//
// Task:
//      multiplies the real ideal equivalence class with itself
//

void
square_real (qi_class & C, const qi_class & A)
{
	debug_handler("qi_class", "square_real");

	bigint newb, d, w, temp;

	// solve d = v A.a + w A.b
	d.assign(xgcd_right(w, A.a, A.b));

	// C.b = A.b + w (D - A.b^2) / (2d)
	square(newb, A.b);
	subtract(newb, qi_class::Delta, newb);
	shift_left(temp, d, 1);
	divide(newb, newb, temp);
	multiply(newb, newb, w);
	add(C.b, newb, A.b);

	// C.a = (A.a/d)^2
	divide(temp, A.a, d);
	square(C.a, temp);

	shift_left(temp, C.a, 1);
	remainder(C.b, C.b, temp);

	C.reduce_real();
}



//
// power(bigint)
//
// Task:
//      computes A^i
//

void
power (qi_class & C, const qi_class & A, const bigint & i)
{
	debug_handler("qi_class", "power(bigint)");

	if (qi_class::Delta.is_lt_zero())
		power_imag(C, A, i);
	else
		power_real(C, A, i);
}



//
// power(long)
//
// Task:
//      computes A^i

void
power (qi_class & C, const qi_class & A, const long i)
{
	debug_handler("qi_class", "power(long)");

	if (qi_class::Delta.is_lt_zero())
		power_imag(C, A, i);
	else
		power_real(C, A, i);
}



//
// power_imag(bigint)
//
// Task:
//      computes A^i using binary exponentiation
//

void
power_imag (qi_class & C, const qi_class & A, const bigint & i)
{
	debug_handler("qi_class", "power_imag(bigint)");

	qi_class B;
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
			multiply_imag(C, C, B);
		j.divide_by_2();
		if (j.is_gt_zero())
			square_imag(B, B);
	}
}



//
// power_imag(long)
//
// Task:
//      computes A^i using binary exponentiation
//

void
power_imag (qi_class & C, const qi_class & A, const long i)
{
	debug_handler("qi_class", "power_imag(long)");

	qi_class B;
	register long j;

	B.assign(A);
	j = i;
	if (j < 0) {
		B.invert();
		j = -j;
	}
	C.assign_one();
	while (j > 0) {
		if ((j & 1) == 1)
			multiply_imag(C, C, B);
		j >>= 1;
		if (j > 0)
			square_imag(B, B);
	}
}



//
// nupower(bigint)
//
// Task:
//      computes A^i using binary exponentiation and the NUCOMP and NUDUPL
//      algorithms.
//

void
nupower (qi_class & C, const qi_class & A, const bigint & i)
{
	debug_handler("qi_class", "nupower(bigint)");

	qi_class B;
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
			nucomp(C, C, B);
		j.divide_by_2();
		if (j.is_gt_zero())
			nudupl(B, B);
	}
}



//
// nupower(long)
//
// Task:
//      computes A^i using binary exponentiation and the NUCOMP and NUDUPL
//      algorithms.
//

void
nupower (qi_class & C, const qi_class & A, const long i)
{
	debug_handler("qi_class", "nupower(long)");

	qi_class B;
	register long j;

	B.assign(A);
	j = i;
	if (j < 0) {
		B.invert();
		j = -j;
	}
	C.assign_one();
	while (j > 0) {
		if ((j & 1) == 1)
			nucomp(C, C, B);
		j >>= 1;
		if (j > 0)
			nudupl(B, B);
	}
}



//
// power_real(bigint)
//
// Task:
//      computes A^i using binary exponentiation
//

void
power_real (qi_class & C, const qi_class & A, const bigint & i)
{
	debug_handler("qi_class", "power_real(bigint)");

	qi_class B;
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
power_real (qi_class & C, const qi_class & A, const long i)
{
	debug_handler("qi_class", "power_real(long)");

	qi_class B;
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
// operator -
//
// Task:
//      returns A^-1
//

qi_class operator - (const qi_class & A)
{
	debug_handler("qi_class", "operator -");

	qi_class B;

	B = inverse(A);
	return B;
}



//
// operator *
//
// Task:
//      multiplies A and B
//

qi_class operator * (const qi_class & A, const qi_class & B)
{
	debug_handler("qi_class", "operator *");

	qi_class C;

	multiply(C, A, B);
	return C;
}



//
// operator /
//
// Task:
//      multiplies A and B^-1
//

qi_class operator / (const qi_class & A, const qi_class & B)
{
	debug_handler("qi_class", "operator /");

	qi_class C;

	divide(C, A, B);
	return C;
}



//
// operator *=
//
// Task:
//      *this = *this * A
//

qi_class & qi_class::operator *= (const qi_class & A)
{
	debug_handler("qi_class", "operator *= ");

	multiply(*this, *this, A);
	return *this;
}



//
// operator /=
//
// Task:
//      *this = *this / A
//

qi_class & qi_class::operator /= (const qi_class & A)
{
	debug_handler("qi_class", "operator /= ");

	multiply(*this, *this, inverse(A));
	return *this;
}



//
// qi_class::is_zero()
//
// Task:
//      tests if the ideal is the zero ideal
//

bool
qi_class::is_zero () const
{
	debug_handler("qi_class", "is_zero");

	return a.is_zero();
}



//
// qi_class::is_one()
//
// Task:
//      tests if the ideal is the unit ideal
//

bool
qi_class::is_one () const
{
	debug_handler("qi_class", "is_one");

	return a.is_one();
}



//
// qi_class::is_equal()
//
// Task:
//      tests if the ideals are equal
//

bool
qi_class::is_equal (const qi_class & B) const
{
	debug_handler("qi_class", "is_equal");

	return (!(a.compare(B.a)) && !(b.compare(B.b)));
}



//
// operator ==
//
// Task:
//      tests if A and B are exacly equal (same representatives)
//

bool operator == (const qi_class & A, const qi_class & B)
{
	debug_handler("qi_class", "operator == ");

	return (A.is_equal(B));
}



//
// operator !=
//
// Task:
//      tests if A and B are not equal
//

bool operator != (const qi_class & A, const qi_class & B)
{
	debug_handler("qi_class", "operator != ");

	return (!A.is_equal(B));
}



//
// operator !
//
// Task:
//      tests if A is the zero ideal
//

bool operator ! (const qi_class & A)
{
	debug_handler("qi_class", "operator !");

	return (A.is_zero());
}



//
// swap
//
// Task:
//      swaps A and B
//

void
swap (qi_class & A, qi_class & B)
{
	debug_handler("qi_class", "swap");

	qi_class C;

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
generate_prime_ideal (qi_class & A, const bigint & p)
{
	debug_handler("qi_class", "generate_prime_ideal");

	bigint D, b, temp;
	int kro, Dp, ip;

	D.assign(qi_class::discriminant());
	if (D.is_zero()) {
		lidia_error_handler("qi_class", "generate_prime_ideal - no current "
				    "quadratic  order has been defined");
		return false;
	}

	kro = kronecker(D, p);
	if (kro < 0)
		return false;
	else {
		if (p == 2) {
			if (kro == 0) {
				Dp = static_cast<int>(remainder(D, 16L));
				if (Dp < 0)
					Dp += 16;
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
// qi_class::rho()
//
// Task:
//      applies the reduction operator once, resulting in another reduced
//      representative of the equivalence class.  If the current quadratic
//      order is not real, the error handler will be evoked.
//

void
qi_class::rho ()
{
	debug_handler("qi_class", "rho");

	bigint q, r, a2, nb, na, oa, temp;
	int s;

	if (Delta.is_le_zero()) {
		lidia_error_handler("qi_class", "rho() - the current quadratic "
				    "order must be real");
		return;
	}

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

	b.assign(nb);
	a.assign(na);
}



//
// qi_class::rho(quadratic_number)
//
// Task:
//      applies the reduction operator once, resulting in another reduced
//      representative of the equivalence class.  If the current quadratic
//      order is not real, the error handler will be evoked.
//

void
qi_class::rho (bigint & oa, quadratic_number_standard & alpha)
{
	debug_handler("qi_class", "rho(quadratic_number)");

	quadratic_number_standard beta;
	bigint q, r, a2, nb, na, temp;
	int s;

	if (Delta.is_le_zero()) {
		lidia_error_handler("qi_class", "rho() - the current quadratic "
				    "order must be real");
		return;
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

	// beta = (1/a) (nb + rd)/2
	//<MM>
	beta.assign_order(alpha.get_order());
	//</MM>
	beta.assign(nb, 1, a);
	multiply(alpha, alpha, beta);

	b.assign(nb);
	oa.assign(a);
	a.assign(na);
}



//
// apply_rho(qi_class, qi_class)
//
// Task:
//      applies the reduction operator once, resulting in another reduced
//      representative of the equivalence class.  If the current quadratic
//      order is not real, the error handler will be evoked.
//

void
apply_rho (qi_class & A, const qi_class & B)
{
	debug_handler("qi_class", "apply_rho(qi_class, qi_class)");

	A.assign(B);
	A.rho();
}



//
// apply_rho(qi_class)
//
// Task:
//      applies the reduction operator once, resulting in another reduced
//      representative of the equivalence class.  If the current quadratic
//      order is not real, the error handler will be evoked.
//

qi_class
apply_rho (const qi_class & A)
{
	debug_handler("qi_class", "apply_rho(qi_class)");

	qi_class B;

	B.assign(A);
	B.rho();
	return B;
}



//
// qi_class::inverse_rho()
//
// Task:
//      applies the inverse reduction operator once, resulting in another
//      reduced representative of the equivalence class.  If the current
//      quadratic order is not real, the error handler will be evoked.
//

void
qi_class::inverse_rho ()
{
	debug_handler("qi_class", "inverse_rho");

	bigint q, a2, r, temp;

	if (Delta.is_le_zero()) {
		lidia_error_handler("qi_class", "inverse_rho() - the current quadratic "
				    "order must be real");
		return;
	}

	// a = (D - b^2) / 4a
	square(temp, b);
	subtract(temp, Delta, temp);
	divide(temp, temp, a);
	shift_right(a, temp, 2);

	// nq = floor((rootD + b) / 2a)
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
// apply_inverse_rho(qi_class, qi_class)
//
// Task:
//      applies the inverse reduction operator once, resulting in another
//      reduced representative of the equivalence class.  If the current
//      quadratic order is not real, the error handler will be evoked.
//

void
apply_inverse_rho (qi_class & A, const qi_class & B)
{
	debug_handler("qi_class", "apply_inverse_rho(qi_class, qi_class)");

	A.assign(B);
	A.inverse_rho();
}



//
// apply_inverse_rho(qi_class)
//
// Task:
//      applies the inverse reduction operator once, resulting in another
//      reduced representative of the equivalence class.  If the current
//      quadratic order is not real, the error handler will be evoked.
//

qi_class
apply_inverse_rho (const qi_class & A)
{
	debug_handler("qi_class", "apply_inverse_rho(qi_class)");

	qi_class B;

	B.assign(A);
	B.inverse_rho();
	return B;
}



//
// qi_class::make_cycle
//
// Task:
//      computes each reduced representative of the qi_class and adds them to
//      the hash table.  This function is used by the order, DL, subgroup,
//      and class group functions.
//

void
qi_class::make_cycle (long x, hash_table< ideal_node > & HT)
{
	debug_handler("qi_class", "make_cycle(lidia_size_t, hash_table)");

	bigint q, r, a2, nb, na, oa, temp;
	int s;
	qi_class tmp;

	tmp.assign(*this);

	// oa = (D - b^2) / 4a
	square(oa, b);
	subtract(oa, Delta, oa);
	divide(oa, oa, a);
	shift_right(oa, oa, 2);

	while (true) {
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

		b.assign(nb);
		oa.assign(a);
		a.assign(na);

		if (is_equal(tmp))  break;
		HT.hash(ideal_node(*this, x));
	}
}



//
// qi_class::make_cycle
//
// Task:
//      computes each reduced representative of the qi_class and adds them to
//      the hash table.  This function is used by the order, DL, subgroup,
//      and class group functions.
//

void
qi_class::make_cycle (long x, indexed_hash_table< ideal_node > & HT)
{
	debug_handler("qi_class", "make_cycle(long, indexed_hash_table)");

	bigint q, r, a2, nb, na, oa, temp;
	int s;
	qi_class tmp;

	tmp.assign(*this);

	// oa = (D - b^2) / 4a
	square(oa, b);
	subtract(oa, Delta, oa);
	divide(oa, oa, a);
	shift_right(oa, oa, 2);

	while (true) {
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

		b.assign(nb);
		oa.assign(a);
		a.assign(na);

		if (is_equal(tmp))
			break;
		HT.hash(ideal_node(*this, x));
	}
}



//
// qi_class::is_equivalent()
//
// Task:
//      tests if the equivalent classes are equal (reduced representatives are
//      equivalent).  If the current order is real, the qi_class_real function
//      is used.
//

bool
qi_class::is_equivalent (const qi_class & B) const
{
	debug_handler("qi_class", "is_equivalent");

	qi_class_real nA, nB;
	bool equiv;
	quadratic_order *cur_qo;
	long oprec = bigfloat::get_precision();

	cur_qo = current_order->get_qo();
	bigfloat::set_precision(cur_qo->prec);

	if (Delta.is_lt_zero())
		equiv = (this->is_equal(B));
	else {
		nA.assign(*this);
		nB.assign(B);
		equiv = nA.is_equivalent(nB);
	}

	bigfloat::set_precision(oprec);
	return equiv;
}



//
// qi_class::is_principal()
//
// Task:
//      tests if the equivalence class is the principal class.  If the current
//      order is real, the qi_class_real function is used.
//

bool
qi_class::is_principal () const
{
	debug_handler("qi_class", "is_principal");

	qi_class_real nA;
	bool prin;
	quadratic_order *cur_qo;
	long oprec = bigfloat::get_precision();

	cur_qo = current_order->get_qo();
	bigfloat::set_precision(cur_qo->prec);

	if (Delta.is_le_zero())
		prin = is_one();
	else {
		nA.assign(*this);
		prin = nA.is_principal();
	}

	bigfloat::set_precision(oprec);
	return prin;
}



//
// qi_class::order_in_CL()
//
// Task:
//      returns the order of the equivalence class in the class group.  The
//      algorithm is chosen based on the size of the discriminant and
//      whether the class number is known.
//

bigint
qi_class::order_in_CL () const
{
	debug_handler("qi_class", "order_in_CL");

	bigfloat temp;
	long upper, v;
	bigint ord;
	int num;

	if (Delta.is_zero()) {
		lidia_error_handler("qi_class", "order_in_CL() - no current quadratic "
				    "order has been defined");
		return bigint();
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class", "order_in_CL() - the current quadratic "
                                    "order has been deleted");
		return bigint();
	}

	// class number known
	if (current_order->get_qo()->is_h_computed())
		ord = order_h();
	else {
		// class number not known
		num = decimal_length(Delta);

		if (Delta.is_gt_zero()) {
			if (num < OBJT_RB)
				ord = order_BJT(1);
			else if (num < OSHANKS_RB)
				ord = order_shanks();
			else
				ord = order_subexp();
		}
		else {
			if (num < OBJT_IB) {
				temp = ceil(power(bigfloat(-Delta), bigfloat(0.25)));
				temp.longify(upper);
				v = (upper >> 1);
				ord = order_BJT(v);
			}
			else if (num < OSHANKS_IB)
				ord = order_shanks();
			else
				ord = order_subexp();
		}
	}

	return ord;
}



//
// qi_class::order_BJT()
//
// Task:
//      computes the order in the class group using a baby-step giant-step
//      variation due to Buchmann, Jacobson, and Teske which has complexity
//      O(sqrt(x)), where x is the order.
//

bigint
qi_class::order_BJT (long v) const
{
	debug_handler("qi_class", "order_BJT");

	bigint ord;
	timer t;

	if (Delta.is_zero()) {
		lidia_error_handler("qi_class", "order_BJT() - no current quadratic "
                                    "order has been defined");
		return bigint();
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class", "order_BJT() - the current quadratic "
                                    "order has been deleted");
		return bigint();
	}

	if (info)
		t.start_timer();

	if (Delta.is_lt_zero()) {
		if (v < 2)
			v = 2;
		if ((v & 1) == 1)
			++v;
		ord = oBJT_imag(v);
	}
	else
		ord = oBJT_real();

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (order_BJT) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	if (do_verify) {
		if (verify_order(ord)) {
			if (info)
				std::cout << "Order is correct." << std::endl;
		}
		else
			std::cout << "Order is not correct!" << std::endl;
	}


	return ord;
}



//
// qi_class::order_h()
//
// Task:
//      computes the order in the class group assuming that the class number is
//      known.  If not, the class number is computed.
//

bigint
qi_class::order_h () const
{
	debug_handler("qi_class", "order_h");

	bigint h, ord;
	rational_factorization hfact;
	timer t;

	if (Delta.is_zero()) {
		lidia_error_handler("qi_class", "order_h() - no current quadratic "
                                    "order has been defined");
		return bigint();
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class", "order_h() - the current quadratic "
                                    "order has been deleted");
		return bigint();
	}

	if (info)
		t.start_timer();

	h = current_order->get_qo()->class_number();
	hfact = current_order->get_qo()->factor_h();
	if (Delta.is_lt_zero())
		ord = omult_imag(h, hfact);
	else {
		if (current_order->get_qo()->prin_list.no_of_elements() > 0)
			ord = omult_real(h, hfact);
		else
			ord = omult_real_sub(h, hfact);
	}

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (order_h) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	if (do_verify) {
		if (verify_order(ord)) {
			if (info)
				std::cout << "Order is correct." << std::endl;
		}
		else
			std::cout << "Order is not correct!" << std::endl;
	}

	return ord;
}



//
// qi_class::order_mult()
//
// Task:
//      computes the order in the class group given a multiple of the actual
//      order.
//

bigint
qi_class::order_mult (const bigint & h,
		      const rational_factorization & hfact) const
{
	debug_handler("qi_class", "order_mult");

	bigint ord;

	if (Delta.is_zero()) {
		lidia_error_handler("qi_class", "order_mult() - no current quadratic "
                                    "order has been defined");
		return bigint();
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class", "order_mult() - the current quadratic "
                                    "order has been deleted");
		return bigint();
	}

	if (h.is_one())
		ord.assign_one();
	else if (Delta.is_lt_zero())
		ord = omult_imag(h, hfact);
	else
		ord = omult_real(h, hfact);

	if (do_verify) {
		if (verify_order(ord)) {
			if (info)
				std::cout << "Order is correct." << std::endl;
		}
		else
			std::cout << "Order is not correct!" << std::endl;
	}

	return ord;
}



//
// qi_class::order_shanks()
//
// Task:
//      computes the order in the class group using the D^1/5 baby-step
//      giant-step variation due to Shanks.  The regulator will be computed
//      if it has not been already.
//

bigint
qi_class::order_shanks () const
{
	debug_handler("qi_class", "order_shanks");

	long OQ, v;
	bigint hstar, ord;
	bigfloat R, nFI, temp, A, F;
	timer t;

	if (Delta.is_zero()) {
		lidia_error_handler("qi_class", "order_shanks() - no current quadratic "
                                    "order has been defined");
		return bigint();
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class", "order_shanks() - the current quadratic "
                                    "order has been deleted");
		return bigint();
	}

	if (info)
		t.start_timer();

	if (Delta.is_lt_zero()) {
		// compute approximation of h
		OQ = current_order->get_qo()->get_optimal_Q_cnum();
		nFI = current_order->get_qo()->estimate_L1(OQ);
		nFI *= sqrt(bigfloat(-Delta)) / Pi();
		if (Delta == -3)
			nFI *= bigfloat(3.0);
		if (Delta == -4)
			nFI *= bigfloat(2.0);
		temp.assign(round(nFI));
		temp.bigintify(hstar);

		A = current_order->get_qo()->estimate_L1_error(OQ);
		F = exp(A) - bigfloat(1.0);
		temp = bigfloat(1.0) - exp(-A);
		if (temp > F)
			F = temp;
		temp = sqrt(nFI*F + abs(nFI - bigfloat(hstar)));
		temp.longify(v);
		if (v < 2)
			v = 2;
		if ((v & 1) == 1)
			++v;

		ord = os_imag(v, hstar);
	}
	else {
		// compute approximation of h
		R = current_order->get_qo()->regulator();
		OQ = current_order->get_qo()->get_optimal_Q_cnum();
		nFI = current_order->get_qo()->estimate_L1(OQ);
		nFI *= sqrt(bigfloat(Delta)) / (bigfloat(2.0)*R);
		temp.assign(round(nFI));
		temp.bigintify(hstar);

		A = current_order->get_qo()->estimate_L1_error(OQ);
		F = exp(A) - bigfloat(1.0);
		temp = bigfloat(1.0) - exp(-A);
		if (temp > F)
			F = temp;
		temp = sqrt(nFI*F + abs(nFI - bigfloat(hstar)));
		temp.longify(v);
		if (v < 2)
			v = 2;
		if ((v & 1) == 1)
			++v;

		ord = os_real(v, hstar);
	}

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (order_shanks) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	if (do_verify) {
		if (verify_order(ord)) {
			if (info)
				std::cout << "Order is correct." << std::endl;
		}
		else
			std::cout << "Order is not correct!" << std::endl;
	}

	return ord;
}



//
// qi_class::order_subexp()
//
// Task:
//      computes the order in the class group in subexponential time by
//      computing the class number subexponentially, and then calling
//      order_h.
//

bigint
qi_class::order_subexp () const
{
	debug_handler("qi_class", "order_subexp");

	bigint ord;
	timer t;

	if (Delta.is_zero()) {
		lidia_error_handler("qi_class", "order_subexp() - no current quadratic "
                                    "order has been defined");
		return bigint();
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class", "order_subexp() - the current quadratic "
                                    "order has been deleted");
		return bigint();
	}

	if (info)
		t.start_timer();

	ord = order_h();

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (order_subexp) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	return ord;
}



//
// qi_class::verify_order()
//
// Task:
//      verifies whether x is the order of the qi_class in the class group.
//

bool
qi_class::verify_order (const bigint & x) const
{
	debug_handler("qi_class", "verify_order");

	qi_class C;
	bool OK;
	timer t;

	if (info)
		t.start_timer();

	power(C, *this, x);

	if (info > 1) {
		std::cout << "\nVerifying order computation:" << std::endl;
		std::cout << "  A = " << *this << std::endl;
		std::cout << "  x = " << x << std::endl;
		std::cout << "  A^x = " << C << std::endl;
	}

	OK = C.is_principal();

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (verify_order) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	return OK;
}



//
// qi_class::DL()
//
// Task:
//      computes the smallest integer x such that G^x is equivalent to  the
//      ideal, if x exists.  If x exists, true is returned, otherwise false is
//      returned and x is set to the order of the equivalence class containing
//      G.  The algorithm is chosen based on the size of the discriminant
//      and whether the class number is known.
//

bool
qi_class::DL (const qi_class & G, bigint & x) const
{
	debug_handler("qi_class", "DL");

	bool isDL;
	int num;

	// FOR TABLE GENERATION
	timer t;
	long t1 = 0;
	bool old = quadratic_order::special;
	quadratic_order::special = false;
	if (old)
		t.start_timer();

	if (Delta.is_zero()) {
		lidia_error_handler("qi_class", "DL() - no current quadratic "
                                    "order has been defined");
		return false;
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class", "DL() - the current quadratic "
                                    "order has been deleted");
		return false;
	}

	if (current_order->get_qo()->is_h_computed()) {
//    isDL = DL_h(G,x);
		isDL = DL_subexp(G, x);
	}
	else {
		// class number not known
		num = decimal_length(Delta);
		if (Delta.is_gt_zero()) {
			if (num < DLBJT_RB)
				isDL = DL_BJT(G, x);
			else
				isDL = DL_subexp(G, x);
		}
		else {
			if (num < DLBJT_IB)
				isDL = DL_BJT(G, x);
			else
				isDL = DL_subexp(G, x);
		}
	}

	// FOR TABLE GENERATION
	quadratic_order::special = old;
	if (quadratic_order::special) {
		t.stop_timer();
		t1 = t.user_time();
		// ROW:  time for representations    total time
		std::cout << t1 << std::endl;
	}

	return isDL;
}



//
// qi_class::DL_BJT(qi_class,bigint)
//
// Task:
//      computes the discrete log base G using a variation of baby-step
//      giant-step due to Buchmann, Jacobson, and Teske with has complexity
//      O(sqrt(x)), where x is the discrete log.
//

bool
qi_class::DL_BJT (const qi_class & G, bigint & x) const
{
	debug_handler("qi_class", "DL_BJT(qi_class, bigint)");

	bigfloat temp;
	long upper, v;

	if (Delta.is_lt_zero()) {
		temp = ceil(power(bigfloat(-Delta), bigfloat(0.25)));
		temp.longify(upper);
		v = (upper >> 1);
	}
	else
		v = 2;

	return DL_BJT(G, x, v);
}



//
// qi_class::DL_BJT(qi_class,bigint,long)
//
// Task:
//      computes the discrete log base G using a variation of baby-step
//      giant-step due to Buchmann, Jacobson, and Teske with has complexity
//      O(sqrt(x)), where x is the discrete log.
//

bool
qi_class::DL_BJT (const qi_class & G, bigint & x, long v) const
{
	debug_handler("qi_class", "DL_BJT(qi_class, bigint, long)");

	bool isDL;
	timer t;

	if (Delta.is_zero()) {
		lidia_error_handler("qi_class", "DL_BJT() - no current quadratic "
                                    "order has been defined");
		return false;
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class", "DL_BJT() - the current quadratic "
                                    "order has been deleted");
		return false;
	}

	if (info)
		t.start_timer();

	if (v < 2)
		v = 2;
	if ((v & 1) == 1)
		++v;

	if (Delta.is_lt_zero())
		isDL = DLBJT_imag(G, x, v);
	else
		isDL = DLBJT_real(G, x);

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (DL_BJT) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	if (do_verify) {
		if (verify_DL(G, x, isDL)) {
			if (info)
				std::cout << "DL is correct." << std::endl;
		}
		else
			std::cout << "DL is not correct!" << std::endl;
	}

	return isDL;
}



//
// qi_class::DL_h()
//
// Task:
//      computes the discrete log base G assuming the class number is known.
//      The class number is computed if it has not been already.
//

bool
qi_class::DL_h (const qi_class & G, bigint & x) const
{
	debug_handler("qi_class", "DL_h");

	bool isDL;
	timer t;

	if (Delta.is_zero()) {
		lidia_error_handler("qi_class", "DL_h() - no current quadratic "
                                    "order has been defined");
		return false;
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class", "DL_h() - the current quadratic "
                                    "order has been deleted");
		return false;
	}

	if (info)
		t.start_timer();

	if (Delta.is_lt_zero())
		isDL = DLh_imag(G, x);
	else
		isDL = DLh_real(G, x);

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (DL_h) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	if (do_verify) {
		if (verify_DL(G, x, isDL)) {
			if (info)
				std::cout << "DL is correct." << std::endl;
		}
		else
			std::cout << "DL is not correct!" << std::endl;
	}

	return isDL;
}



//
// qi_class::DL_subexp()
//
// Task:
//      computes the discrete log base G using the subexponential algorithm
//      based on that of Buchmann and Duellmann.
//

bool
qi_class::DL_subexp (const qi_class & G, bigint & x) const
{
	debug_handler("qi_class", "DL_subexp");

	bool isDL;
	timer t;

	if (Delta.is_zero()) {
		lidia_error_handler("qi_class", "DL_subexp() - no current quadratic "
                                    "order has been defined");
		return false;
	}

	if (!current_order->get_qo()) {
		lidia_error_handler("qi_class", "DL_subexp() - the current quadratic "
                                    "order has been deleted");
		return false;
	}

	if (info)
		t.start_timer();

	isDL = DLsubexp(G, x);

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (DL_subexp) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	if (do_verify) {
		if (verify_DL(G, x, isDL)) {
			if (info)
				std::cout << "DL is correct." << std::endl;
		}
		else
			std::cout << "DL is not correct!" << std::endl;
	}

	return isDL;
}



//
// qi_class::verify_DL()
//
// Task:
//      verifies whether x is either the DL of the qi_class to the base G, or
//      the order of G.
//

bool
qi_class::verify_DL (const qi_class & G, const bigint & x,
		     const bool is_DL) const
{
	debug_handler("qi_class", "verify_DL");

	qi_class C;
	bool OK;
	timer t;

	if (info)
		t.start_timer();

	power(C, G, x);
	if (info > 1) {
		std::cout << "\nVerifying DL computation:" << std::endl;
		std::cout << "  G = " << G << std::endl;
		std::cout << "  A = " << *this << std::endl;
		std::cout << "  x = " << x << std::endl;
		std::cout << "  G^x = " << C << std::endl;
		if (is_DL)
			std::cout << "  x is log_G of A." << std::endl;
		else
			std::cout << "  x is ord(G)." << std::endl;
	}

	if (is_DL)
		OK = C.is_equivalent(*this);
	else
		OK = C.is_principal();

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (verify_DL) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	return OK;
}



//
// subgroup()
//
// Task:
//      computes the structure of the subgroup generated by the ideals in G.
//      The algorithm is chosen based on the size of the discriminant and
//      whether the class number is known.
//

base_vector< bigint >
subgroup (base_vector< qi_class > & G)
{
	debug_handler("qi_class", "subgroup");

	lidia_size_t i, Gsize;
	base_vector< long > v;
	base_vector< bigint > S;

	v.set_mode(EXPAND);
	v.reset();
	S.set_mode(EXPAND);
	S.reset();

	if (qi_class::get_current_order().is_h_computed()) {
//    S = subgroup_h(G);
		Gsize = G.size();
		v.set_size(Gsize);
		for (i = 0; i < Gsize; ++i)
			v[i] = 2;
		S = subgroup_BJT(G, v);
	}
	else {
		// class number not known
		Gsize = G.size();
		v.set_size(Gsize);
		for (i = 0; i < Gsize; ++i)
			v[i] = 2;
		S = subgroup_BJT(G, v);
	}

	return S;
}



//
// subgroup_BJT()
//
// Task:
//      computes the structure of the subgroup generated by G using a
//      variation of baby-step giant-step due to Buchmann, Jacobson, and Teske
//      which has complexity O(sqrt(H)), where H is the order of the subgroup.
//

base_vector< bigint >
subgroup_BJT (base_vector< qi_class > & G, base_vector< long > & v)
{
	debug_handler("qi_class", "subgroup_BJT");

	base_vector< bigint > S;
	quadratic_order Q;
	timer t;

	Q = qi_class::get_current_order();

	S.set_mode(EXPAND);
	S.reset();

	if (qi_class::info)
		t.start_timer();

	if (qi_class::Delta.is_lt_zero())
		S = subBJT_imag(G, v);
	else
		S = subBJT_real(G);

	if (qi_class::info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (subgroup_BJT) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	return S;
}



//
// subgroup_h()
//
// Task:
//      computes the structure of the subgroup generated by G assuming that
//      the class number is known.  The class number is computed if it has
//      not been already.
//

base_vector< bigint >
subgroup_h (base_vector< qi_class > & G)
{
	debug_handler("qi_class", "subgroup_h");

	base_vector< bigint > S;
	quadratic_order Q;
	timer t;

	Q = qi_class::get_current_order();

	if (qi_class::info)
		t.start_timer();

	S.set_mode(EXPAND);
	S.reset();

	if (qi_class::Delta.is_lt_zero())
		S = subh_imag(G);
	else
		S = subh_real(G);

	if (qi_class::info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (subgroup_h) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	return S;
}



//
// operator >>
//
// Task:
//      inputs a qi_class from the std::istream in.
//

std::istream &
operator >> (std::istream & in, qi_class & A)
{
	debug_handler("qi_class", "operator >>");

	int n = 0;
	char c;
	bigint ibuf[2];

	in >> c;
	if (c != '(') {
		lidia_error_handler("qi_class", "operator >>::(expected");
		return in;
	}

	in >> c;
	while (c != ')' && n != 2) {
		in.putback(c);
		in >> ibuf[n];
		n++;
		in >> c;
		if (c == ',')
			in >> c;
	}
	A.assign(ibuf[0], ibuf[1]);
	return in;
}



//
// operator <<
//
// Task:
//      outputs a qi_class to the std::ostream out.
//

std::ostream &
operator << (std::ostream & out, const qi_class & A)
{
	debug_handler("qi_class", "operator << ");

	out << "(" << A.a << ", " << A.b << ")";
	return out;
}



///////////////////////
// PRIVATE FUNCTIONS //
///////////////////////


//
// qi_class::oBJT_imag()
//
// Task:
//      computes the order in the class group using the BJT algorithm for
//      imaginary orders.  The parameter v is the size of the initial step
//      width.
//

bigint
qi_class::oBJT_imag (long v) const
{
	debug_handler("qi_class", "oBJT_imag");

	bigint y, usqr;
	long upper, u, r, s;
	bigint x;
	bigfloat temp;
	qi_class A, B, C, CINV;
	ideal_node *Inode;
	hash_table< ideal_node > htable;

	// get hash table size
	temp = ceil(power(bigfloat(-Delta), bigfloat(0.25)));
	temp.multiply_by_2();
	multiply(temp, temp, bigfloat(1.33));
	temp.longify(upper);

	htable.initialize(static_cast<lidia_size_t>(upper));
	htable.set_key_function(&ideal_node_key);

	x.assign_zero();
	u = v;
	s = 1;
	power_imag(C, *this, u);
	CINV = inverse(*this);

	A.assign_one();
	B.assign(C);
	y.assign(v);

	while (x.is_zero()) {
		// compute more baby steps
		for (r = s; r <= u; ++r) {
			multiply_imag(A, A, CINV);
			if (s == 1) {
				if (A.is_one()) {
					x.assign(r);
					break;
				}
				else
					htable.hash(ideal_node(A, r));
			}
			else
				htable.hash(ideal_node(A, r));
		}

		// compute giant steps to u^2
		square(usqr, u);
		while ((y.compare(usqr) < 0) && (x.is_zero())) {
			Inode = htable.search(ideal_node(B, 0));
			if (Inode) {
				// found b in list:  x = y+r
				r = Inode->get_index();
				add(x, y, r);
			}
			else {
				// not found, take another giant step
				add(y, y, u);
				multiply_imag(B, B, C);
			}
		}

		if (x.is_zero()) {
			// double u
			s = u+1;
			u <<= 1;
			square_imag(C, C);
		}
	}

	return x;
}



//
// oBJT_real()
//
// Task:
//      computes the order in the class group using an unpublished variation of
//      the BJT algorithm for real orders.  The initial step width is taken to
//      be 1 in order to efficiently distinguish principal ideals.
//      Basically, we store the entire cycle of each baby step if the
//      regulator is small enough.  Otherwise, we store all ideals in the cycle
//      with distance < sqrt(R) and modify the search for each giant step.
//

bigint
qi_class::oBJT_real () const
{
	debug_handler("qi_class", "oBJT_real");

	bigint y, usqr;
	long upper, u, r, s;
	bigint x;
	bigfloat temp, Reg, sqReg, GStepWidth;
	qi_class A, Aprime, B, C, CINV;
	qi_class_real Areal, F, Gstep;
	ideal_node *Inode;
	hash_table< ideal_node > htable;
	quadratic_order *cur_qo;
	long oprec = bigfloat::get_precision();

	cur_qo = current_order->get_qo();
	bigfloat::set_precision(cur_qo->prec);

	// compute regulator
	Reg = cur_qo->regulator();
	sqrt(sqReg, Reg);

	if (cur_qo->is_subexp_computed()) {
		bigfloat::set_precision(oprec);
		return order_h();
	}

	// compute giant step ideal (distance = sqrt(R))
	F.assign_one();
	Gstep.assign(nearest(F, sqReg));
	if (Gstep.get_distance() > sqReg)
		Gstep.inverse_rho();

	// add ideals to principal list, if necessary
	if (cur_qo->prin_list.no_of_elements() > 0)
		F.assign(cur_qo->prin_list.last_entry());
	else {
		F.assign_one();
		F.rho();
		cur_qo->prin_list.hash(F);
	}
	if (Gstep.is_one()) {
		while (!F.is_one()) {
			F.rho();
			cur_qo->prin_list.hash(F);
		}
	}
	else {
		while ((F.get_distance() < Gstep.get_distance()) && (!F.is_one())) {
			F.rho();
			cur_qo->prin_list.hash(F);
		}
		if (!F.is_one()) {
			F.rho();
			cur_qo->prin_list.hash(F);
		}
		if (F.is_one())
			Gstep.assign(F);
	}
	floor(GStepWidth, Gstep.get_distance());

	// get hash table size
	temp = ceil(power(bigfloat(Delta), bigfloat(0.25)));
	temp.multiply_by_2();
	multiply(temp, temp, bigfloat(1.33));
	temp.longify(upper);

	htable.initialize(static_cast<lidia_size_t>(upper));
	htable.set_key_function(&ideal_node_key);

	x.assign_zero();
	s = u = 1;
	C.assign(*this);
	CINV = inverse(*this);

	A.assign_one();
	B.assign(*this);
	y.assign_one();

	// principality test
	if (is_principal())
		x.assign_one();

	while (x.is_zero()) {
		// compute more baby steps
		for (r = s; r <= u; ++r) {
			multiply_real(A, A, CINV);
			htable.hash(ideal_node(A, r));

			if (Gstep.is_one()) {
				// store each cycle of ideals
				A.make_cycle(r, htable);
			}
			else {
				// store ideals with distance < sqReg
				Areal.assign(A);
				Areal.make_list(sqReg, r, htable);
			}
		}

		// compute giant steps to u^2
		square(usqr, u);
		while ((y.compare(usqr) <= 0) && (x.is_zero())) {
			if (Gstep.is_one()) {
				Inode = htable.search(ideal_node(B, 0));
				if (Inode) {
					// found b in list:  x = y+r
					r = Inode->get_index();
					add(x, y, r);
				}
			}
			else {
				F.assign(B, 0.0);
				temp.assign_zero();
				while ((F.get_distance() <= Reg) && (x.is_zero())) {
					Inode = htable.search(ideal_node(qi_class(F), 0));
					if (Inode) {
						// found b in list:  x = y+r
						r = Inode->get_index();
						add(x, y, r);
					}
					else {
						add(temp, temp, GStepWidth);
						multiply_real(F, F, Gstep);
						F.adjust_pos(temp);
					}
				}
			}
			if (x.is_zero()) {
				// not found, take another giant step
				add(y, y, u);
				multiply_real(B, B, C);
				if (y.is_one()) {
					multiply_real(B, B, C);
					inc(y);
				}
			}
		}

		if (x.is_zero()) {
			// double u
			s = u+1;
			u <<= 1;
			square_real(C, C);
		}
	}

	bigfloat::set_precision(oprec);
	return x;
}



//
// qi_class::omult_imag()
//
// Task:
//      computes the order in the class group given a multiple of the actual
//      order.
//

bigint
qi_class::omult_imag (const bigint & h,
		      const rational_factorization & hfact) const
{
	debug_handler("qi_class", "omult_imag");

	register lidia_size_t i;
	register int j, ex, num;
	bigint ord, p, pwr;
	qi_class A;

	ord.assign_one();

	if ((!h.is_one()) && (!this->is_one())) {
		num = hfact.no_of_comp();
		for (i = 0; i < num; ++i) {
			A.assign(*this);
			p.assign(hfact.base(i));
			ex = hfact.exponent(i);
			pwr.assign(h);
			for (j = 0; j < ex; ++j)
				divide(pwr, pwr, p);
			if (pwr > 1)
				power_imag(A, A, pwr);

			while (!A.is_one()) {
				multiply(ord, ord, p);
				power_imag(A, A, p);
				if (ord > h)
					return bigint(0);
			}
		}
	}

	return ord;
}



//
// qi_class::omult_real()
//
// Task:
//      computes the order in the class group given a multiple of the actual
//      order.
//

bigint
qi_class::omult_real (const bigint & h,
		      const rational_factorization & hfact) const
{
	debug_handler("qi_class", "omult_real");

	register lidia_size_t i;
	register int j, ex, num;
	bigint ord, p, pwr;
	qi_class A;

	ord.assign_one();

	if (!h.is_one()) {
		num = hfact.no_of_comp();
		for (i = 0; i < num; ++i) {
			A.assign(*this);
			p.assign(hfact.base(i));
			ex = hfact.exponent(i);
			pwr.assign(h);
			for (j = 0; j < ex; ++j)
				divide(pwr, pwr, p);
			if (pwr > 1)
				power_real(A, A, pwr);

			while (!A.is_principal()) {
				multiply(ord, ord, p);
				power_real(A, A, p);
				if (ord > h)
					return bigint(0);
			}
		}
	}

	return ord;
}



//
// qi_class::omult_real_sub()
//
// Task:
//      computes the order in the class group given a multiple of the actual
//      order, assuming that the class group was computed with a subexponential
//      method
//

bigint
qi_class::omult_real_sub (const bigint & h,
			  const rational_factorization & hfact) const
{
	debug_handler("qi_class", "omult_real_sub");

	register lidia_size_t i, k;
	register int j, ex, num;
	bigint ord, p, pwr;
	qi_class A;
	math_vector< bigint > vec, vecA;
	quadratic_number_standard q;

	ord.assign_one();

	if (!h.is_one()) {
		vecA = current_order->get_qo()->represent_over_FB_sieve(*this, q);

		num = hfact.no_of_comp();
		for (i = 0; i < num; ++i) {
			vec = vecA;
			A.assign(*this);
			p.assign(hfact.base(i));
			ex = hfact.exponent(i);
			pwr.assign(h);
			for (j = 0; j < ex; ++j)
				divide(pwr, pwr, p);
			if (pwr > 1) {
				power_real(A, A, pwr);
				for (k = 0; k < vec.size(); ++k)
					vec[k] = pwr*vec[k];
			}

			while (!current_order->get_qo()->is_in_lattice(vec)) {
				multiply(ord, ord, p);
				if (ord > h)
					return bigint(0);

				power_real(A, A, p);
				for (k = 0; k < vec.size(); ++k)
					vec[k] = p*vec[k];
			}
		}
	}

	return ord;
}



//
// qi_class::os_imag()
//
// Task:
//      computes the order in the class group using Shanks' D^1/5 algorithm
//      for imaginary orders. The parameter v is the size of the initial step
//      width.
//

bigint
qi_class::os_imag (long v, const bigint & hstar) const
{
	debug_handler("qi_class", "os_imag");

	bigint y, usqr;
	long upper, u, r, s;
	bigint x;
	bigfloat temp;
	qi_class A, B, Bcomp, C, CINV, GINV;
	ideal_node *Inode;
	hash_table< ideal_node > htable;

	x.assign_zero();
	u = v;
	y.assign_zero();
	s = 1;

	// check whether this^ordmult = (1)
	if (!ordmult.is_zero()) {
		power_imag(B, *this, ordmult);
		if (B.is_one())
			x = omult_imag(ordmult, omfact);
	}

	if (x.is_zero()) {
		// check whether this^hstar = (1)
		power_imag(B, *this, hstar);
		if (B.is_one()) {
			ordmult.assign(hstar);
			omfact.assign(ordmult);
			omfact.factor();
			x = omult_imag(ordmult, omfact);
		}
		else {
			// get hash table size
			temp = ceil(power(bigfloat(-Delta), bigfloat(0.2)));
			if (bigfloat(v) > temp)
				temp.assign(v);
			temp *= bigfloat(1.33);
			temp.longify(upper);
			htable.initialize(static_cast<lidia_size_t>(upper));
			htable.set_key_function(&ideal_node_key);

			A.assign_one();
			htable.hash(ideal_node(A, 0));
			GINV = inverse(*this);
			power_imag(C, *this, u);
			CINV = inverse(C);
			Bcomp.assign(B);
		}
	}

	while (x.is_zero()) {
		// compute more baby steps
		for (r = s; r <= u; ++r) {
			multiply_imag(A, A, GINV);
			if (A.is_one()) {
				x = r;
				break;
			}
			else
				htable.hash(ideal_node(A, r));
		}

		// compute giant steps to u^2
		square(usqr, u);
		while ((y.compare(usqr) < 0) && (x.is_zero())) {
			Inode = htable.search(ideal_node(Bcomp, 0));
			if (Inode) {
				// found Bcomp in list:  x = hstar - y + r)
				r = Inode->get_index();
				subtract(x, hstar, y);
				add(x, x, r);
				if (x.is_gt_zero()) {
					ordmult.assign(x);
					omfact.assign(x);
					omfact.factor();
					x = omult_imag(ordmult, omfact);
				}
			}

			if (x.is_zero()) {
				Inode = htable.search(ideal_node(B, 0));
				if (Inode) {
					// found B in list:  x = hstar + y +r
					r = Inode->get_index();
					add(ordmult, hstar, y);
					add(ordmult, ordmult, r);
					omfact.assign(ordmult);
					omfact.factor();
					x = omult_imag(ordmult, omfact);
				}
			}

			if (x.is_zero()) {
				// not found, take another giant step
				add(y, y, u);
				multiply_imag(B, B, C);
				if (B.is_one()) {
					add(ordmult, hstar, y);
					omfact.assign(ordmult);
					omfact.factor();
					x = omult_imag(ordmult, omfact);
				}
				else {
					multiply_imag(Bcomp, Bcomp, CINV);
					if (Bcomp.is_one()) {
						subtract(ordmult, hstar, y);
						omfact.assign(ordmult);
						omfact.factor();
						x = omult_imag(ordmult, omfact);
					}
				}
			}
		}

		if (x.is_zero()) {
			// double u
			s = u+1;
			u <<= 1;
			square_imag(C, C);
			CINV = inverse(C);
		}
	}

	return x;
}



//
// qi_class::os_real()
//
// Task:
//      computes the order in the class group using a variation of Shanks'
//      D^1/5 algorithm for real orders.
//      Basically, we store the entire cycle of each baby step if the
//      regulator is small enough.  Otherwise, we store all ideals in the cycle
//      with distance < sqrt(R) and modify the search for each giant step.
//

bigint
qi_class::os_real (long v, const bigint & hstar) const
{
	debug_handler("qi_class", "os_real");

	bigint y, usqr;
	long upper, u, r, s;
	bigint x;
	bigfloat temp, Reg, sqReg, GStepWidth;
	qi_class A, Aprime, B, Bcomp, C, CINV, G, GINV;
	qi_class_real Areal, F, Gstep;
	ideal_node *Inode;
	hash_table< ideal_node > htable;
	quadratic_order *cur_qo;
	long oprec = bigfloat::get_precision();

	cur_qo = current_order->get_qo();
	bigfloat::set_precision(cur_qo->prec);

	x.assign_zero();

	// compute regulator
	Reg = cur_qo->regulator();
	sqrt(sqReg, Reg);

	if (cur_qo->is_subexp_computed()) {
		bigfloat::set_precision(oprec);
		return order_h();
	}

	// compute giant step ideal (distance = sqrt(R))
	F.assign_one();
	Gstep.assign(nearest(F, sqReg));
	if (Gstep.get_distance() > sqReg)
		Gstep.inverse_rho();

	// add ideals to principal list, if necessary
	if (cur_qo->prin_list.no_of_elements() > 0)
		F.assign(cur_qo->prin_list.last_entry());
	else {
		F.assign_one();
		F.rho();
		cur_qo->prin_list.hash(F);
	}
	if (Gstep.is_one()) {
		while (!F.is_one()) {
			F.rho();
			cur_qo->prin_list.hash(F);
		}
	}
	else {
		while ((F.get_distance() < Gstep.get_distance()) && (!F.is_one())) {
			F.rho();
			cur_qo->prin_list.hash(F);
		}
		if (!F.is_one()) {
			F.rho();
			cur_qo->prin_list.hash(F);
		}
		if (F.is_one())
			Gstep.assign(F);
	}
	floor(GStepWidth, Gstep.get_distance());

	// check whether this^ordmult = (1)
	if (!ordmult.is_zero()) {
		power_real(B, *this, ordmult);
		if (B.is_principal())
			x = omult_real(ordmult, omfact);
	}

	s = u = 1;
	y.assign_zero();

	if (x.is_zero()) {
		// check whether this^hstar = (1)
		power_real(B, *this, hstar);
		if (B.is_principal()) {
			ordmult.assign(hstar);
			omfact.assign(ordmult);
			omfact.factor();
			x = omult_real(ordmult, omfact);
		}
		else {
			// get hash table size
			temp = ceil(power(bigfloat(-Delta), bigfloat(0.2)));
			if (bigfloat(v) > temp)
				temp.assign(v);
			temp *= bigfloat(1.33);
			temp.longify(upper);
			htable.initialize(static_cast<lidia_size_t>(upper));
			htable.set_key_function(&ideal_node_key);

			A.assign_one();
			GINV = inverse(*this);
			power_real(C, *this, u);
			CINV = inverse(C);
			Bcomp.assign(B);
		}
	}

	while (x.is_zero()) {
		// compute more baby steps
		for (r = s; r <= u; ++r) {
			multiply_real(A, A, GINV);
			htable.hash(ideal_node(A, r));

			if (Gstep.is_one()) {
				// store each cycle of ideals
				A.make_cycle(r, htable);
			}
			else {
				// store ideals with distance < sqReg
				Areal.assign(A);
				Areal.make_list(sqReg, r, htable);
			}
		}

		// compute giant steps to u^2
		r = -1;
		square(usqr, u);
		while ((y.compare(usqr) < 0) && (x.is_zero())) {
			if (Gstep.is_one()) {
				Inode = htable.search(ideal_node(Bcomp, 0));
				if (Inode) {
					// found Bcomp in list:  x = hstar - y + r)
					r = Inode->get_index();
				}
			}
			else {
				F.assign(Bcomp, 0.0);
				temp.assign_zero();
				while ((F.get_distance() <= Reg) && (x.is_zero())) {
					Inode = htable.search(ideal_node(qi_class(F), 0));
					if (Inode) {
						// found Bcomp in list:  x = hstar - y + r)
						r = Inode->get_index();
					}
					else {
						add(temp, temp, GStepWidth);
						multiply_real(F, F, Gstep);
						F.adjust_pos(temp);
					}
				}
			}
			if (r >= 0) {
				subtract(x, hstar, y);
				add(x, x, r);
				if (x.is_gt_zero()) {
					ordmult.assign(x);
					omfact.assign(x);
					omfact.factor();
					x = omult_real(ordmult, omfact);
				}
			}

			if (x.is_zero()) {
				r = -1;
				if (Gstep.is_one()) {
					Inode = htable.search(ideal_node(B, 0));
					if (Inode) {
						// found B in list:  x = hstar + y + r)
						r = Inode->get_index();
					}
				}
				else {
					F.assign(B, 0.0);
					temp.assign_zero();
					while ((F.get_distance() <= Reg) && (x.is_zero())) {
						Inode = htable.search(ideal_node(qi_class(F), 0));
						if (Inode) {
							// found Bin list:  x = hstar + y + r)
							r = Inode->get_index();
						}
						else {
							add(temp, temp, GStepWidth);
							multiply_real(F, F, Gstep);
							F.adjust_pos(temp);
						}
					}
				}
				if (r >= 0) {
					add(x, hstar, y);
					add(x, x, r);
					ordmult.assign(x);
					omfact.assign(x);
					omfact.factor();
					x = omult_real(ordmult, omfact);
				}
			}

			if (x.is_zero()) {
				// not found, take another giant step
				add(y, y, u);
				multiply_real(B, B, C);
				multiply_real(Bcomp, Bcomp, CINV);
				if (y.is_one()) {
					multiply_real(B, B, C);
					multiply_real(Bcomp, Bcomp, CINV);
					inc(y);
				}
				if (Bcomp.is_principal()) {
					subtract(x, hstar, y);
					ordmult.assign(x);
					omfact.assign(x);
					omfact.factor();
					x = omult_real(ordmult, omfact);
				}
				if (Bcomp.is_principal()) {
					subtract(x, hstar, y);
					ordmult.assign(x);
					omfact.assign(x);
					omfact.factor();
					x = omult_real(ordmult, omfact);
				}
			}
		}

		if (x.is_zero()) {
			// double u
			s = u+1;
			u <<= 1;
			square_real(C, C);
			CINV = inverse(C);
		}
	}

	bigfloat::set_precision(oprec);
	return x;
}



//
// qi_class::DLBJT_imag()
//
// Task:
//      computes the dl base G in the class group using the BJT algorithm for
//      imaginary orders.  The parameter v is the size of the initial step
//      width.
//

bool
qi_class::DLBJT_imag (const qi_class & G, bigint & x, long v) const
{
	debug_handler("qi_class", "DLBJT_imag");

	bigint usqr, y;
	long upper, u, r, s;
	bigfloat temp;
	qi_class A, B, C, CINV, DINV, E;
	ideal_node *Inode;
	hash_table< ideal_node > htable;
	bool DL;

	if (is_one()) {
		x.assign_zero();
		return true;
	}

	if (is_equal(G)) {
		x.assign_one();
		return true;
	}

	// get hash table size
	temp = ceil(power(bigfloat(-Delta), bigfloat(0.25)));
	temp.multiply_by_2();
	multiply(temp, temp, bigfloat(1.33));
	temp.longify(upper);

	htable.initialize(static_cast<lidia_size_t>(upper));
	htable.set_key_function(&ideal_node_key);

	DL = false;
	x.assign_zero();
	u = v;
	s = 1;
	power_imag(C, G, u);
	CINV = inverse(G);
	DINV = inverse(*this);

	A.assign_one();
	B.assign(C);
	y.assign(v);

	while (x.is_zero()) {
		// compute more baby steps
		for (r = s; r <= u; ++r) {
			multiply_imag(A, A, CINV);
			if (s == 1) {
				if (A.is_equal(DINV)) {
					x.assign(r);
					DL = true;
					break;
				}

				if (A.is_one()) {
					x.assign(r);
					break;
				}

				htable.hash(ideal_node(A, r));
			}
			else
				htable.hash(ideal_node(A, r));
		}

		// compute giant steps to u^2
		square(usqr, u);
		while ((y.compare(usqr) < 0) && (x.is_zero())) {
			multiply_imag(E, DINV, B);
			Inode = htable.search(ideal_node(E, 0));
			if (Inode) {
				// found b*d^-1 in list:  x = y+r
				r = Inode->get_index();
				add(x, y, r);
				DL = true;
			}
			else {
				Inode = htable.search(ideal_node(B, 0));
				if (Inode) {
					// found b in list:  x = y+r
					r = Inode->get_index();
					add(x, y, r);
				}
				else {
					// not found, take another giant step
					add(y, y, u);
					multiply_imag(B, B, C);
				}
			}
		}

		// double u
		if (x.is_zero()) {
			s = u+1;
			u <<= 1;
			square_imag(C, C);
		}
	}

	return DL;
}



//
// qi_class::DLBJT_real()
//
// Task:
//      computes the dl in the class group using an unpublished variation of
//      the BJT algorithm for real orders.  The initial step width is taken to
//      be 1 in order to efficiently distinguish principal ideals.
//      Basically, we store the entire cycle of each baby step if the
//      regulator is small enough.  Otherwise, we store all ideals in the cycle
//      with distance < sqrt(R) and modify the search for each giant step.
//

bool
qi_class::DLBJT_real (const qi_class & G, bigint & x) const
{
	debug_handler("qi_class", "DLBJT_real");

	bigint y, usqr;
	long upper, u, r, s;
	bigfloat temp, Reg, sqReg, GStepWidth;
	qi_class A, Aprime, B, C, CINV, DINV, E;
	qi_class_real Areal, F, Gstep;
	ideal_node *Inode;
	hash_table< ideal_node > htable;
	bool is_DL;
	quadratic_order *cur_qo;
	long oprec = bigfloat::get_precision();

	cur_qo = current_order->get_qo();
	bigfloat::set_precision(cur_qo->prec);

	// compute regulator
	Reg = cur_qo->regulator();
	sqrt(sqReg, Reg);

	if (cur_qo->is_subexp_computed()) {
		bigfloat::set_precision(oprec);
		return DLsubexp(G, x);
	}

	// compute giant step ideal (distance = sqrt(R))
	F.assign_one();
	Gstep.assign(nearest(F, sqReg));
	if (Gstep.get_distance() > sqReg)
		Gstep.inverse_rho();

	// add ideals to principal list, if necessary
	if (cur_qo->prin_list.no_of_elements() > 0)
		F.assign(cur_qo->prin_list.last_entry());
	else {
		F.assign_one();
		F.rho();
		cur_qo->prin_list.hash(F);
	}
	if (Gstep.is_one()) {
		while (!F.is_one()) {
			F.rho();
			cur_qo->prin_list.hash(F);
		}
	}
	else {
		while ((F.get_distance() < Gstep.get_distance()) && (!F.is_one())) {
			F.rho();
			cur_qo->prin_list.hash(F);
		}
		if (!F.is_one()) {
			F.rho();
			cur_qo->prin_list.hash(F);
		}
		if (F.is_one())
			Gstep.assign(F);
	}
	floor(GStepWidth, Gstep.get_distance());

	// get hash table size
	temp = ceil(power(bigfloat(Delta), bigfloat(0.25)));
	temp.multiply_by_2();
	multiply(temp, temp, bigfloat(1.33));
	temp.longify(upper);

	htable.initialize(static_cast<lidia_size_t>(upper));
	htable.set_key_function(&ideal_node_key);

	x.assign_zero();
	s = u = 1;
	y.assign_one();
	C.assign(G);
	CINV = inverse(G);
	DINV = inverse(*this);

	A.assign_one();
	B.assign(G);
	is_DL = false;

	// principality test
	if (is_principal()) {
		x.assign_zero();
		is_DL = true;
	}
	else if (G.is_principal())
		x.assign_one();

	while ((x.is_zero()) && (!is_DL)) {
		// compute more baby steps
		for (r = s; r <= u; ++r) {
			multiply_real(A, A, CINV);
			htable.hash(ideal_node(A, r));

			if (Gstep.is_one()) {
				// store each cycle of ideals
				A.make_cycle(r, htable);
			}
			else {
				// store ideals with distance < sqReg
				Areal.assign(A);
				Areal.make_list(sqReg, r, htable);
			}

			if (s == 1) {
				if (Gstep.is_one()) {
					Inode = htable.search(ideal_node(DINV, 0));
					if (Inode) {
						x.assign_one();
						is_DL = true;
					}
				}
				else {
					F.assign(DINV, 0.0);
					temp.assign_zero();
					while ((F.get_distance() <= Reg) && (x.is_zero())) {
						Inode = htable.search(ideal_node(qi_class(F), 0));
						if (Inode) {
							x.assign_one();
							is_DL = true;
						}
						else {
							add(temp, temp, GStepWidth);
							multiply_real(F, F, Gstep);
							F.adjust_pos(temp);
						}
					}
				}
			}
		}

		// compute giant steps to u^2
		square(usqr, u);
		while ((y.compare(usqr) <= 0) && (x.is_zero())) {
			multiply_real(E, B, DINV);
			if (Gstep.is_one()) {
				Inode = htable.search(ideal_node(E, 0));
				if (Inode) {
					// found b*d^-1 in list:  x = y+r
					r = Inode->get_index();
					add(x, y, r);
					is_DL = true;
				}
			}
			else {
				F.assign(E, 0.0);
				temp.assign_zero();
				while ((F.get_distance() <= Reg) && (x.is_zero())) {
					Inode = htable.search(ideal_node(qi_class(F), 0));
					if (Inode) {
						// found b*d^-1 in list:  x = y+r
						r = Inode->get_index();
						add(x, y, r);
						is_DL = true;
					}
					else {
						add(temp, temp, GStepWidth);
						multiply_real(F, F, Gstep);
						F.adjust_pos(temp);
					}
				}
			}

			if (x.is_zero()) {
				if (Gstep.is_one()) {
					Inode = htable.search(ideal_node(B, 0));
					if (Inode) {
						// found b in list:  z = y+r
						r = Inode->get_index();
						add(x, y, r);
					}
				}
				else {
					F.assign(B, 0.0);
					temp.assign_zero();
					while ((F.get_distance() <= Reg) && (x.is_zero())) {
						Inode = htable.search(ideal_node(qi_class(F), 0));
						if (Inode) {
							// found b in list:  z = y+r
							r = Inode->get_index();
							add(x, y, r);
						}
						else {
							add(temp, temp, GStepWidth);
							multiply_real(F, F, Gstep);
							F.adjust_pos(temp);
						}
					}
				}
			}

			if (x.is_zero()) {
				// not found, take another giant step
				add(y, y, u);
				multiply_real(B, B, C);
				if (y.is_one()) {
					multiply_real(B, B, C);
					inc(y);
				}
			}
		}

		if (x.is_zero()) {
			// double u
			s = u+1;
			u <<= 1;
			square_real(C, C);
		}
	}

	bigfloat::set_precision(oprec);
	return is_DL;
}



//
// FIX
//

//
// qi_class::DLh_imag()
//
// Task:
//      computes the dl in the class group assuming that the class number is
//      known.
//

bool
qi_class::DLh_imag (const qi_class & G, bigint & x) const
{
	debug_handler("qi_class", "DLh_imag");

	bigint h;

	h = current_order->get_qo()->class_number();

	return DL_BJT(G, x);
}



//
// FIX
//

//
// qi_class::DLh_real()
//
// Task:
//      computes the dl in the class group assuming that the class number is
//      known.
//

bool
qi_class::DLh_real (const qi_class & G, bigint & x) const
{
	debug_handler("qi_class", "DLh_real");

	bigint h;

	x = 0;

	h = current_order->get_qo()->class_number();

	return DLBJT_real(G, x);
}



//
// qi_class::DLsubexp()
//
// Task:
//      computes the dl in the class group using a subexponential method based
//      on one due to Buchmann and Duellmann.  The class group is computed
//      first.
//

bool
qi_class::DLsubexp (const qi_class & G, bigint & x) const
{
	debug_handler("qi_class", "DLsubexp_imag");

	quadratic_order *curqo;
	qi_class A, B;
	base_vector< bigint > CL;
	math_vector< bigint > vec1, vec2, xvec;
	lidia_size_t numCL, i, j;
	bigint d, m, gdiv, temp, xi, aa, h, c, cj, v;
	bigfloat ftemp;
	long upper;
	bool is_DL;

	is_DL = true;
	curqo = current_order->get_qo();
	CL = curqo->class_group_siqs();
	h = curqo->class_number();
	numCL = CL.size();

	if (!curqo->is_subexp_computed()) {
		if (Delta.is_lt_zero()) {
			ftemp = ceil(power(bigfloat(-Delta), bigfloat(0.25)));
			ftemp.longify(upper);
			return DLBJT_imag(G, x, (upper >> 1));
		}
		else
			return DLBJT_real(G, x);
	}

	if (is_principal()) {
		x.assign_zero();
		return true;
	}

	if (is_equivalent(G)) {
		x.assign_one();
		return true;
	}

	vec1.set_capacity(numCL);
	vec1.set_size(numCL);
	vec2.set_capacity(numCL);
	vec2.set_size(numCL);
	xvec.set_capacity(numCL);
	xvec.set_size(numCL);

	A.assign(G);
	B.assign(*this);
	vec1 = curqo->represent_over_generators_sieve(A);
	vec2 = curqo->represent_over_generators_sieve(B);


#ifdef LIDIA_DEBUG
	std::cout << "G = " << G << std::endl;
	std::cout << " = " << vec1 << std::endl;
	std::cout << "A = " << *this << std::endl;
	std::cout << " = " << vec2 << std::endl;
	std::cout << "CL = " << CL << std::endl;
	std::cout << "gens = " << curqo->generators() << std::endl;
#endif

	// compute xi's and test for solution
	for (i = 0; i < numCL; ++i) {
		gdiv = xgcd_left(aa, vec1[i], CL[i]);
		remainder(temp, vec2[i], gdiv);

#ifdef LIDIA_DEBUG
		std::cout << "i = " << i << ", gdiv = " << gdiv << ", temp = ";
		std::cout << temp << std::endl;
		std::cout << "aa = " << aa << std::endl;
#endif

		if (!temp.is_zero()) {
			is_DL = false;
			break;
		}
		else {
			divide(CL[i], CL[i], gdiv);

			divide(temp, vec2[i], gdiv);
			multiply(xi, aa, temp);
			remainder(xi, xi, CL[i]);
			if (xi.is_lt_zero())
				add(xi, xi, CL[i]);
			xvec[i] = xi;
		}
	}

#ifdef LIDIA_DEBUG
	std::cout << "xvec = " << xvec << std::endl;
	std::cout << "CL = " << CL << std::endl;
#endif

	if (is_DL) {
		x.assign(xvec[0]);
		i = 1;
		m.assign_one();
		while ((i < numCL) && (is_DL)) {
			// test for solution
			for (j = i; j < numCL; ++j) {
				subtract(temp, xvec[i-1], xvec[j]);
				d = gcd(CL[i-1], CL[j]);
				remainder(temp, temp, d);
				if (!temp.is_zero()) {
					is_DL = false;
					break;
				}
			}

			if (is_DL) {
				temp = gcd(m, CL[i-1]);
				multiply(m, m, CL[i-1]);
				divide(m, m, temp);

				// c = x_i - x, d = gcd(m,CL[i])
				subtract(c, xvec[i], x);
				d = gcd(m, CL[i]);

				// solve cj m = d (mod CL[i])
				temp = xgcd_left(cj, m, CL[i]);

				// vj = (cj c) / d (mod CL[i]/d)
				multiply(v, cj, c);
				divide(v, v, d);
				divide(temp, CL[i], d);
				remainder(v, v, temp);
				if (v.is_lt_zero())
					add(v, v, temp);

				// x += vj m
				multiply(temp, v, m);
				add(x, x, temp);

				++i;
			}
		}

		remainder(x, x, h);
		if (x.is_lt_zero())
			add(x, x, h);
	}

#ifdef LIDIA_DEBUG
	std::cout << "x = " << x << std::endl;
	std::cout << "is_DL ? " << is_DL << std::endl;
#endif

	if (!is_DL)
		x = G.order_h();

	return is_DL;
}



//
// subBJT_imag()
//
// Task:
//      computes the structure of the subgroup generated by G using the BJT
//      algorithm for imaginary quadratic orders.  The parameter v is the
//      vector of initial step widths to be used.
//

base_vector< bigint >
subBJT_imag (base_vector< qi_class > & G, base_vector< long > & v)
{
	debug_handler("qi_class", "subBJT_imag");

	bigint usqr, y, det, Bjj;
	long s, u, r, q, upper, Bj, idx;
	lidia_size_t i, j, k, Gsize, rank, numRpr, numQ, curr_index;
	bigfloat temp;
	qi_class Gidl, A, B, C, D, E, H, Gq, GBj;
	matrix< bigint > Bmat, junk, U;
	base_vector< bigint > CL;
	base_vector< long > Rvec, Qvec;
	ideal_node *Inode;
	indexed_hash_table< ideal_node > Q, R;

	Rvec.set_mode(EXPAND);
	Qvec.set_mode(EXPAND);
	CL.set_mode(EXPAND);
	Rvec.reset();
	Qvec.reset();
	CL.reset();
	U.resize(1, 1);

	Gsize = G.size();

	// get maximum set sizes
	temp = ceil(power(bigfloat(-qi_class::discriminant()), bigfloat(0.25)));
	temp.longify(upper);

	Q.initialize(static_cast<lidia_size_t>(upper) << 1);
	Q.set_key_function(&ideal_node_key);
	B.assign_one();
	Q.hash(ideal_node(B, 0));

	R.initialize(static_cast<lidia_size_t>(upper) << 1);
	R.set_key_function(&ideal_node_key);
	R.hash(ideal_node(B, 0));

	rank = 0;
	det.assign_one();
	for (j = 0; j < Gsize; ++j) {
		u = v[j];
		if (u < 2)
			u = 2;
		if ((u & 1) == 1)
			++u;
		y.assign(u);
		s = 1;

		Gidl = G[j];
		A.assign_one();
		power_imag(B, Gidl, u);
		C.assign(B);
		H.assign(inverse(Gidl));

		numRpr = R.no_of_elements();
		curr_index = numRpr;
		numQ = Q.no_of_elements();

		if (Gidl.is_one())
			Bjj.assign_one();
		else
			Bjj.assign_zero();

		// check whether current ideal is in previously generated subgroup
		if ((rank > 0) && (Bjj.is_zero())) {
			for (i = 0; i < numQ; ++i) {
				D.assign(Q[i].get_A());
				multiply_imag(E, D, Gidl);
				Inode = R.search(ideal_node(E, 0));
				if ((E.is_one()) || (Inode)) {
					Bjj.assign_one();
					break;
				}
			}
		}

		while (Bjj.is_zero()) {
			// compute more baby steps
			for (r = s; r <= u; ++r) {
				multiply_imag(A, A, H);
				if ((s == 1) && (r > 1)) {
					for (k = 0; k < numRpr; ++k) {
						D.assign(R[k].get_A());
						multiply_imag(E, D, A);

						// check if relation already found
						Inode = Q.search(ideal_node(E, 0));
						if (Inode) {
							q = Inode->get_index();
							Bjj.assign(r);
							decode_vector(Bmat, Bjj, k, q, Rvec, Qvec, numRpr, numQ);
							break;
						}
						else {
							R.hash(ideal_node(E, curr_index));
							++curr_index;
						}
					}
					if (!Bjj.is_zero())
						break;
				}
				else {
					for (k = 0; k < numRpr; ++k) {
						D.assign(R[k].get_A());
						multiply_imag(E, D, A);
						R.hash(ideal_node(E, curr_index));
						++curr_index;
					}
				}
			}

			// compute giant steps to u^2
			square(usqr, u);
			while ((y.compare(usqr) < 0) && (Bjj.is_zero())) {
				// search
				for (i = 0; i < numQ; ++i) {
					E.assign(Q[i].get_A());
					multiply_imag(D, E, B);
					Inode = R.search(ideal_node(D, 0));
					if (Inode) {
						r = Inode->get_index();
						q = i;
						add(Bjj, y, (r/numRpr));
						if (Bjj > 1) {
							r %= numRpr;
							decode_vector(Bmat, Bjj, r, q, Rvec, Qvec, numRpr, numQ);
						}
						break;
					}
				}

				if (Bjj.is_zero()) {
					// not found, take another giant step
					add(y, y, u);
					multiply_imag(B, B, C);
				}
			}

			// double u
			if (Bjj.is_zero()) {
				s = u+1;
				u <<= 1;
				square_imag(C, C);
			}
		}

		if (!Bjj.is_one())
			++rank;

		if (j < (Gsize-1)) {
			temp = ceil(sqrt(bigfloat(Bjj)));
			temp.longify(Bj);

			// compute new R' (remove entries with too large exponents)
			det *= Bj;
			det.longify(idx);
			numRpr = R.no_of_elements();
			for (i = numRpr-1; i >= idx; --i)
				R.remove_from(i);

			if (!Bjj.is_one()) {
				Rvec[rank-1] = Bj;
				Qvec[rank-1] = 1;

				// compute new Q
				numQ = Q.no_of_elements();
				curr_index = numQ;
				power_imag(GBj, Gidl, Bj);
				Gq.assign(GBj);
				for (i = 1; i < Bj; ++i) {
					if (Bjj.compare(i*Bj) <= 0)
						break;
					for (k = 0; k < numQ; ++k) {
						E.assign(Q[k].get_A());
						multiply_imag(D, E, Gq);
						Q.hash(ideal_node(D, curr_index));
						++curr_index;
					}
					multiply_imag(Gq, Gq, GBj);
					++Qvec[rank-1];
				}
			}
		}
	}

	// compute structure
	if (rank == 0)
		CL[0] = 1;
	else {
		Bmat.snf_havas(U, junk);
		i = 0;
		for (j = 0; j < rank; ++j) {
			Bjj = Bmat.member(j, j);
			if (!Bjj.is_one()) {
				CL[i] = Bjj;
				++i;
			}
		}
	}

	return CL;
}



//
// subBJT_real()
//
// Task:
//      computes the structure of the subgroup generated by G using an
//      unpublished variation of the BJT algorithm for real orders.
//      Basically, we store the entire cycle of each baby step if the
//      regulator is small enough.  Otherwise, we store all ideals in the cycle
//      with distance < sqrt(R) and modify the search for each giant step.
//

base_vector< bigint >
subBJT_real (base_vector< qi_class > & G)
{
	debug_handler("qi_class", "subBJT_real");

	bigint usqr, y, det, Bjj;
	long s, u, r, q, upper, Bj, idx;
	lidia_size_t i, j, k, Gsize, rank, numRpr, numQ, curr_index, numRreps;
	bigfloat temp, Reg, sqReg, GStepWidth;
	qi_class Gidl, A, Aprime, B, C, D, E, H, Gq, GBj;
	qi_class_real Areal, F, Gstep, *FS;
	matrix< bigint > Bmat, junk, U;
	base_vector< bigint > CL;
	base_vector< long > Rvec, Qvec;
	ideal_node *Inode;
	indexed_hash_table< ideal_node > Q, R;
	base_vector< qi_class > Rreps;
	quadratic_order *cur_qo;
	long oprec = bigfloat::get_precision();

	cur_qo = qi_class::current_order->get_qo();
	bigfloat::set_precision(cur_qo->prec);

	U.resize(1, 1);

	// compute regulator
	Reg = cur_qo->regulator_shanks();
	sqrt(sqReg, Reg);

	// compute giant step ideal (distance = sqrt(R))
	F.assign_one();
	Gstep.assign(nearest(F, sqReg));
	if (Gstep.get_distance() > sqReg)
		Gstep.inverse_rho();

	// add ideals to principal list, if necessary
	if (cur_qo->prin_list.no_of_elements() > 0)
		F.assign(cur_qo->prin_list.last_entry());
	else {
		F.assign_one();
		F.rho();
		cur_qo->prin_list.hash(F);
	}

	if (Gstep.is_one()) {
		while (!F.is_one()) {
			F.rho();
			cur_qo->prin_list.hash(F);
		}
	}
	else {
		while ((F.get_distance() < Gstep.get_distance()) && (!F.is_one())) {
			F.rho();
			cur_qo->prin_list.hash(F);
		}
		if (!F.is_one()) {
			F.rho();
			cur_qo->prin_list.hash(F);
		}
		if (F.is_one())
			Gstep.assign(F);
	}
	floor(GStepWidth, Gstep.get_distance());

	Rvec.set_mode(EXPAND);
	Qvec.set_mode(EXPAND);
	CL.set_mode(EXPAND);
	Rreps.set_mode(EXPAND);
	Rvec.reset();
	Qvec.reset();
	CL.reset();
	Rreps.reset();

	Gsize = G.size();

	// get maximum set sizes
	temp = ceil(power(bigfloat(qi_class::discriminant()), bigfloat(0.25)));
	temp.longify(upper);

	Q.initialize(static_cast<lidia_size_t>(upper) << 1);
	Q.set_key_function(&ideal_node_key);
	B.assign_one();
	Q.hash(ideal_node(B, 0));

	R.initialize(static_cast<lidia_size_t>(upper) << 1);
	R.set_key_function(&ideal_node_key);
	R.hash(ideal_node(B, 0));

	Rreps[0] = B;

	rank = 0;
	det.assign_one();
	for (j = 0; j < Gsize; ++j) {
		s = u = 1;
		y.assign_one();

		Gidl = G[j];
		A.assign_one();
		B.assign(Gidl);
		C.assign(B);
		H.assign(inverse(Gidl));

		numQ = Q.no_of_elements();
		curr_index = numRpr = numRreps = Rreps.size();

		Bjj.assign_zero();

		// check whether current ideal is in previously generated subgroup
		if (Bjj.is_zero()) {
			for (i = 0; i < numQ; ++i) {
				D.assign(Q[i].get_A());
				multiply_real(E, D, Gidl);

				if (Gstep.is_one()) {
					FS = cur_qo->prin_list.search(qi_class_real(E));
					if (FS)
						Bjj.assign_one();
					else {
						Inode = R.search(ideal_node(E, 0));
						if (Inode)
							Bjj.assign_one();
					}
				}
				else {
					F.assign(E, 0.0);
					temp.assign_zero();
					while ((F.get_distance() <= Reg) && (Bjj.is_zero())) {
						FS = cur_qo->prin_list.search(F);
						if (FS)
							Bjj.assign_one();
						else {
							Inode = R.search(ideal_node(qi_class(F), 0));
							if (Inode)
								Bjj.assign_one();
							else {
								add(temp, temp, GStepWidth);
								multiply_real(F, F, Gstep);
								F.adjust_pos(temp);
							}
						}
					}
				}

				if (Bjj.is_one())
					break;
			}
		}

		while (Bjj.is_zero()) {
			// compute more baby steps
			for (r = s; r <= u; ++r) {
				multiply_real(A, A, H);
				for (k = 0; k < numRreps; ++k) {
					D.assign(Rreps[k]);
					multiply_real(E, D, A);
					Rreps[curr_index] = E;
					R.hash(ideal_node(E, curr_index));

					if (Gstep.is_one()) {
						// store each cycle of ideals
						E.make_cycle(curr_index, R);
					}
					else {
						// store ideals with distance < sqReg
						Areal.assign(E);
						Areal.make_list(sqReg, curr_index, R);
					}

					++curr_index;
				}
			}

			// compute giant steps to u^2
			square(usqr, u);
			while ((y.compare(usqr) <= 0) && (Bjj.is_zero())) {
				// search
				for (i = 0; i < numQ; ++i) {
					D.assign(Q[i].get_A());
					multiply_real(E, D, B);

					if (Gstep.is_one()) {
						FS = cur_qo->prin_list.search(qi_class_real(E));
						if (FS) {
							r = 0;
							q = i;
							add(Bjj, y, (r/numRpr));
							r %= numRpr;
							decode_vector(Bmat, Bjj, r, q, Rvec, Qvec, numRpr, numQ);
						}
						else {
							Inode = R.search(ideal_node(E, 0));
							if (Inode) {
								r = Inode->get_index();
								q = i;
								add(Bjj, y, (r/numRpr));
								r %= numRpr;
								decode_vector(Bmat, Bjj, r, q, Rvec, Qvec, numRpr, numQ);
							}
						}
					}
					else {
						F.assign(E, 0.0);
						temp.assign_zero();
						while ((F.get_distance() <= Reg) && (Bjj.is_zero())) {
							FS = cur_qo->prin_list.search(F);
							if (FS) {
								r = 0;
								q = i;
								add(Bjj, y, (r/numRpr));
								r %= numRpr;
								decode_vector(Bmat, Bjj, r, q, Rvec, Qvec, numRpr, numQ);
							}
							else {
								Inode = R.search(ideal_node(qi_class(F), 0));
								if (Inode) {
									r = Inode->get_index();
									q = i;
									add(Bjj, y, (r/numRpr));
									r %= numRpr;
									decode_vector(Bmat, Bjj, r, q, Rvec, Qvec, numRpr, numQ);
								}
								else {
									add(temp, temp, GStepWidth);
									multiply_real(F, F, Gstep);
									F.adjust_pos(temp);
								}
							}
						}
					}

					if (Bjj.is_gt_zero())
						break;
				}

				if (Bjj.is_zero()) {
					// not found, take another giant step
					add(y, y, u);
					multiply_real(B, B, C);
					if (y.is_one()) {
						multiply_real(B, B, C);
						inc(y);
					}
				}
			}

			// double u
			if (Bjj.is_zero()) {
				s = u+1;
				u <<= 1;
				square_real(C, C);
			}
		}

		if (!Bjj.is_one())
			++rank;

		if (j < (Gsize-1)) {
			temp = ceil(sqrt(bigfloat(Bjj)));
			temp.longify(Bj);

			// compute new R' (remove entries with too large exponents)
			det *= Bj;
			det.longify(idx);
			numRpr = R.no_of_elements();
			i = numRpr-1;
			while (R[i].get_index() >= idx) {
				R.remove_from(i);
				--i;
			}
			Rreps.set_size(static_cast<lidia_size_t>(idx));

			if (!Bjj.is_one()) {
				Rvec[rank-1] = Bj;
				Qvec[rank-1] = 1;

				// compute new Q
				numQ = Q.no_of_elements();
				curr_index = numQ;
				power_real(GBj, Gidl, Bj);
				Gq.assign(GBj);
				for (i = 1; i < Bj; ++i) {
					if (Bjj.compare(i*Bj) <= 0)
						break;
					for (k = 0; k < numQ; ++k) {
						E.assign(Q[k].get_A());
						multiply_real(D, E, Gq);
						Q.hash(ideal_node(D, curr_index));
						++curr_index;
					}
					multiply_real(Gq, Gq, GBj);
					++Qvec[rank-1];
				}
			}
		}
	}

	// compute structure
	if (rank == 0)
		CL[0] = 1;
	else {
		Bmat.snf_havas(U, junk);
		i = 0;
		for (j = 0; j < rank; ++j) {
			Bjj = Bmat.member(j, j);
			if (!Bjj.is_one()) {
				CL[i] = Bjj;
				++i;
			}
		}
	}

	bigfloat::set_precision(oprec);
	return CL;
}



//
// FIX
//

//
// subh_imag()
//
// Task:
//      computed the structure of the subgroup generated by G assuming the
//      class number is known.
//

base_vector< bigint >
subh_imag (base_vector< qi_class > & G)
{
	debug_handler("qi_class", "subh_imag");

	base_vector< bigint > CL;
	bigint h;

	CL.set_mode(EXPAND);
	CL.reset();

	h = qi_class::current_order->get_qo()->class_number();

	G = G;

	return CL;
}



//
// FIX
//

//
// subh_real()
//
// Task:
//      computed the structure of the subgroup generated by G assuming the
//      class number is known.
//

base_vector< bigint >
subh_real (base_vector< qi_class > & G)
{
	debug_handler("qi_class", "subh_real");

	base_vector< bigint > CL;
	bigint h;

	CL.set_mode(EXPAND);
	CL.reset();

	h = qi_class::current_order->get_qo()->class_number();

	G = G;

	return CL;
}



//
// key functions needed for hash tables
//

bigint
ideal_node_key (const ideal_node & G)
{
	return G.get_A().get_a();
}



void
swap (ideal_node & A, ideal_node & B) {
	ideal_node C;
	C = A;
	A = B;
	B = C;
}



std::istream &
operator >> (std::istream & in, ideal_node & A)
{
	in >> A.A;
	in >> A.index;
	return in;
}



std::ostream &
operator << (std::ostream & out, const ideal_node & A)
{
	out << "[" << A.A << ", index = " << A.index << "]";
	return out;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
