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
#include	"LiDIA/quadratic_form.h"
#include	"LiDIA/nmbrthry_functions.h"
#include	"LiDIA/rational_factorization.h"
#include	"LiDIA/hash_table.h"
#include	"LiDIA/quadratic_order.h"
#include	"LiDIA/qi_class.h"
#include	"LiDIA/qi_class_real.h"
#include	"LiDIA/quadratic_ideal.h"
#include	"LiDIA/number_fields/qo_list.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bigint
quadratic_form_key (const quadratic_form & G)
{
	return G.get_a();
}



//
// qf_floor()
//
// Task:
//      special floor function used by some quadratic form routines
//

inline bigint
qf_floor (bigint x, bigint y)
{
	debug_handler("quadratic_form", "qf_floor");

	bigint d, r;
	div_rem(d, r, x, y);
	if (x.sign()*y.sign() == -1 && !r.is_zero())
		dec(d);
	return d;
}



//
// qf_round()
//
// Task:
//      special round function used by some quadratic form routines
//

inline bigint
qf_round (bigint rn, bigint rd)
{
	debug_handler("quadratic_form", "qf_round");

	rn.multiply_by_2();
	add(rn, rn, rd);
	rd.multiply_by_2();
	return qf_floor(rn, rd);
}



//
// quadratic_form::is_regular_disc(const bigint & Delta)
//
// Task:
//      returns true if Delta is a regular discriminant, false otherwise.
//

inline bool
is_regular_disc (const bigint & Delta)
{
	debug_handler("quadratic_form", "is_regular_disc");

	if (Delta.is_zero() || Delta.is_one())
		return false;

	if (Delta.is_lt_zero()) {
		return true;
	}
	else {
	    // Delta is *irregular* iff Delta is a square, 
	    // i.e. Delta is a perfect power m^k with integral m and 
	    // even k > 0.
	    // Note: is_power returns 0 if Delta is not a perfect power.
	        bigint temp;
		long val = is_power(temp, Delta);
		return (val == 0) || ((val & 1) == 1);
	}
}



//
// quadratic_form::norm_number_pos_def()
//
// Task:
//      returns the normalization number of the positive definite form.
//

inline bigint
quadratic_form::norm_number_pos_def ()
{
	debug_handler("quadratic_form", "norm_number_pos_def");

	bigint a2, temp, q, r;

	// s = floor((a-b) / 2a)
	shift_left(a2, a, 1);
	subtract(temp, a, b);
	div_rem(q, r, temp, a2);
	if ((temp.is_lt_zero()) && (!r.is_zero()))
		dec(q);

	return q;
}



//
// quadratic_form::norm_number_indef()
//
// Task:
//      returns the normalization number of the indefinite form.
//

inline bigint
quadratic_form::norm_number_indef ()
{
	debug_handler("quadratic_form", "norm_number_indef");

#if 0
	bigint a2, absa(a), q, temp;

	absa.absolute_value();
	shift_left(a2, absa, 1);
	if (absa.compare(rootD) <= 0) {
		// q = floor((rootD - b) / 2|a|)
		subtract(temp, rootD, b);
	}
	else {
		// q = floor((|a| - b) / 2|a|)
		subtract(temp, absa, b);
	}
	divide(q, temp, a2);
	if (a.is_lt_zero())
		q.negate();

	return q;
#endif

	bigint abs_a, bn, abs_a_mult_2, tmp, s, _b(b);
	int sign_a;

	_b.negate();

	sign_a = a.sign();
	abs_a.assign(abs(a));

	shift_left(abs_a_mult_2, abs_a, 1);
	if (abs_a.compare(rootD) > 0) {
		s = qf_round(_b, abs_a_mult_2);
		multiply(s, s, sign_a);
	}
	else {
		s = qf_floor(rootD - b , abs_a_mult_2);
		multiply(s, s, sign_a);
	}

	return s;
}



//
// quadratic_form::norm_number()
//
// Task:
//      returns the normalization number of the form.
//

bigint
quadratic_form::norm_number ()
{
	debug_handler("quadratic_form", "norm_number");

	quadratic_form f;

	if (is_pos_definite())
		return norm_number_pos_def();
	else if (is_indefinite())
		return norm_number_indef();
	else if (is_neg_definite()) {
		f.assign(-a, b, -c);
		return f.norm_number_pos_def();
	}
	else
		return bigint(0);
}



//
// quadratic_form::is_reduced_pos_def()
//
// Task:
//      returns true if the positive definite form is reduced.
//

inline bool
quadratic_form::is_reduced_pos_def () const
{
	debug_handler("quadratic_form", "is_reduced_pos_def");

	bool tmp;

	if ((c.compare(a) > 0) || ((c == a) && (b.is_ge_zero())))
		tmp = true;
	else
		tmp = false;

	tmp = (tmp && is_normal());

	return tmp;
}



//
// quadratic_form::is_reduced_indef()
//
// Task:
//      returns true if the indefinite form is reduced.
//

inline bool
quadratic_form:: is_reduced_indef() const
{
	debug_handler("quadratic_form", "is_reduced_indef");

	bigint abs_a, _abs_a, abs_a_2, bn;
	bool tmp;

	abs_a.assign(abs(a));

	negate(_abs_a, abs_a);
	shift_left(abs_a_2, abs_a, 1);
	subtract(bn, rootD, abs_a_2);
	if (bn.is_negative())
		inc(bn);
	bn.absolute_value();
	tmp = ((bn < b) && (b <= rootD));

	return (tmp && is_normal());
}



//
// quadratic_form::is_reduced_irregular()
//
// Task:
//      returns true if the irregular form is reduced.
//

inline bool
quadratic_form::is_reduced_irregular () const
{
	debug_handler("quadratic_form", "is_reduced_irregular");

	bool x;

	if ((a.is_zero()) && (b.is_zero()))
		x = true;
	else if ((a.is_zero()) && (b.is_positive()) && (c.is_ge_zero()) && (c < b))
		x = true;
	else
		x = false;

	return x;
}



//
// quadratic_form::almost_reduce_irregular()
//
// Task:
//      executes the first step in reducing the irregular form.
//

inline void
quadratic_form::almost_reduce_irregular (matrix_GL2Z & U)
{
	debug_handler("quadratic_form", "almost_reduce_irregular");

	matrix_GL2Z V((1), (0), (0), (1));
	bigint d;

	if (!a.is_zero()) {
		if (Delta.is_zero())
			d.assign_zero();
		else {
			d.assign(rootD);
			if (Delta.is_lt_zero())
				d.negate();
		}
		bigint t(a);
		t.multiply_by_2();
		bigint s(b+d);
		s.negate();
		bigint u, v, g;
		g = xgcd(v, u, s, t);
		s = s/g;
		t = t/g;

		u.negate();
		V = matrix_GL2Z (s, u, t, v);

		transform(V);
		multiply(U, U, V);
	}
}



//
// quadratic_form::reduce_irregular(matrix_GL2Z)
//
// Task:
//      reduces the irregular form and computes the corresponding
//      transformation matrix.
//

inline void
quadratic_form::reduce_irregular (matrix_GL2Z & U)
{
	debug_handler("quadratic_form", "reduce_irregular(matrix_GL2Z)");

	matrix_GL2Z W;
	matrix_GL2Z D((0), (-1), (1), (0));

	almost_reduce_irregular(U);

	if (b.is_zero())
		return;

	if (b.is_negative() && c.is_zero()) {
		b.negate();
		W = matrix_GL2Z((0), (-1), (1), (0));
		multiply(U, U, W);
		return;
	}

	if (b.is_negative()) {
		swap(a, c);
		b.negate();
		multiply(U, U, D);
		almost_reduce_irregular(U);
	}
	bigint s(qf_floor(c, b));
	s.negate();
	W = matrix_GL2Z ((1), s, (0), (1));
	transform(W);
	multiply(U, U, W);
}



//
// quadratic_form::reduce_irregular()
//
// Task:
//      reduces the irregular form.
//

inline void
quadratic_form::reduce_irregular ()
{
	debug_handler("quadratic_form", "reduce_irregular");

	matrix_GL2Z U((1), (0), (0), (1));

	reduce_irregular(U);
}



//
// quadratic_form::reduce_pos_def(matrix_GL2Z)
//
// Task:
//      reduces the positive definite form and computes the corresponding
//      transformation matrix.
//

inline void
quadratic_form::reduce_pos_def (matrix_GL2Z & U)
{
	debug_handler("quadratic_form", "reduce_pos_def(matrix_GL2Z)");

	bigint temp;
	matrix_GL2Z V((0), (1), (-1), (0));

	while (!is_reduced_pos_def()) {
		multiply(U, U, V);
		temp.assign(a);
		a.assign(c);
		c.assign(temp);
		b.negate();
		normalize(U);
	}
}



//
// quadratic_form::reduce_pos_def()
//
// Task:
//      reduces the positive definite form.
//

inline void
quadratic_form::reduce_pos_def ()
{
	debug_handler("quadratic_form", "reduce_pos_def");

	bigint oa, a2, t, nt, q, m, r, temp;
	int st, nst;

	normalize_regular();

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

	while (a.compare(c) > 0) {
		oa.assign(a);
		a.assign(c);
		t.assign(nt);
		st = nst;

		shift_left(a2, a, 1);
		div_rem(q, r, t, a2);
		subtract(m, a, r);

		// c = oa - q*((r + t) >> 1);
		add(temp, r, t);
		temp.divide_by_2();
		multiply(temp, temp, q);
		subtract(c, oa, temp);

		if (m.is_ge_zero()) {
			nt.assign(r);
			nst = -st;
		}
		else {
			shift_left(nt, m, 1);
			add(nt, nt, r);
			add(c, c, m);
			nst = st;
		}
	}

	b.assign(t);
	if (st < 0)
		b.negate();

	oa.assign(c);
	normalize_regular();
	c.assign(oa);

	if (!a.compare(c) && (b.is_lt_zero()))
		b.negate();
}



//
// quadratic_form::reduce_indef(matrix_GL2Z)
//
// Task:
//      reduces the indefinite form and computes the corresponding
//      transformation matrix.
//

inline void
quadratic_form::reduce_indef (matrix_GL2Z & U)
{
	debug_handler("quadratic_form", "reduce_indef(matrix_GL2Z)");

	bigint temp;
	matrix_GL2Z V((0), (1), (-1), (0)), W((0), (1), (1), (0));

	while (!is_reduced_indef()) {
		multiply(U, U, W);
		multiply(U, U, V);
		b.negate();
		normalize(U);

		multiply(U, U, W);
		temp.assign(a);
		a.assign(c);
		c.assign(temp);
		normalize(U);
	}
}



//
// quadratic_form::reduce_indef()
//
// Task:
//      reduces the indefinite form.
//

inline void
quadratic_form::reduce_indef ()
{
	debug_handler("quadratic_form", "reduce_indef");

	bigint temp;

	while (!is_reduced_indef()) {
		b.negate();
		normalize_regular();

		temp.assign(a);
		a.assign(c);
		c.assign(temp);
		normalize_regular();
	}
}



//
// quadratic_form::prop_equivalent_pos_def()
//
// Task:
//      returns true if the positive definite forms are properly equivalent.
//

bool
quadratic_form::prop_equivalent_pos_def (const quadratic_form & g) const
{
	debug_handler("quadratic_form", "prop_equivalent_pos_def");

	quadratic_form rf(*this), rg(g);

	rg.reduce();
	rf.reduce();

	return (!rf.compare(rg));
}



//
// quadratic_form::prop_equivalent_neg_def()
//
// Task:
//      returns true if the negative definite forms are properly equivalent.
//

bool
quadratic_form::prop_equivalent_neg_def (const quadratic_form & g) const
{
	debug_handler("quadratic_form", "prop_equivalent_neg_def");

	quadratic_form rf, rg;

	rf.assign(-a, b, -c);
	rg.assign(-g.a, g.b, -g.c);

	rg.reduce();
	rf.reduce();

	return (!rf.compare(rg));
}



//
// quadratic_form::prop_equivalent_indef()
//
// Task:
//      returns true if the indefinite forms are properly equivalent.
//

bool
quadratic_form::prop_equivalent_indef (const quadratic_form & g) const
{
	debug_handler("quadratic_form", "prop_equivalent_indef");

	quadratic_form rf(*this), rg(g), root;

	rf.reduce();
	rg.reduce();

	if (!rf.compare(rg))
		return true;

	root.assign(rg);
	do {
		rg.rho();
		if (rf.is_equal(rg))
			return true;
	} while (root.compare(rg));

	return false;
}



//
// quadratic_form::prop_equivalent_irregular()
//
// Task:
//      returns true if the irregular forms are properly equivalent.
//

bool
quadratic_form::prop_equivalent_irregular (const quadratic_form & g) const
{
	debug_handler("quadratic_form", "prop_equivalent_irregular");

	quadratic_form rf(*this), rg(g);

	rf.reduce();
	rg.reduce();

	return (!rf.compare(rg));
}



//
// quadratic_form::prop_equivalent_pos_def()
//
// Task:
//      returns true if the positive definite forms are properly equivalent
//      and computes the corresponding transformation matrix.
//

bool
quadratic_form::prop_equivalent_pos_def (const quadratic_form & g,
					 matrix_GL2Z & U) const
{
	debug_handler("quadratic_form", "prop_equivalent_pos_def(matrix_GL2Z)");

	quadratic_form rf(*this), rg(g);
	matrix_GL2Z V((1), (0), (0), (1));
	matrix_GL2Z S(V);
	bool equiv;

	rg.reduce(S);
	rf.reduce(V);
	if (rf == rg) {
		V.invert();
		multiply(S, S, V);
		multiply(U, U, S);
		equiv = true;
	}
	else
		equiv = false;

	return equiv;
}



//
// quadratic_form::prop_equivalent_neg_def()
//
// Task:
//      returns true if the negative definite forms are properly equivalent
//      and computes the corresponding transformation matrix.
//

bool
quadratic_form::prop_equivalent_neg_def (const quadratic_form & g,
					 matrix_GL2Z & U) const
{
	debug_handler("quadratic_form", "prop_equivalent_neg_def(matrix_GL2Z)");

	quadratic_form rf(*this), rg(g);
	matrix_GL2Z V((1), (0), (0), (1));
	matrix_GL2Z S(V);
	bool equiv;

	rg.reduce(S);
	rf.reduce(V);
	if (rf == rg) {
		V.invert();
		multiply (S, S, V);
		multiply(U, U, S);
		equiv = true;
	}
	else
		equiv = false;

	return equiv;
}



//
// quadratic_form::prop_equivalent_indef()
//
// Task:
//      returns true if the indefinite forms are properly equivalent
//      and computes the corresponding transformation matrix.
//

bool
quadratic_form::prop_equivalent_indef (const quadratic_form & g,
				       matrix_GL2Z & U) const
{
	debug_handler("quadratic_form", "prop_equivalent_indef(matrxi_GL2Z)");

	quadratic_form rf(*this), rg(g), rg2(g), root;
	matrix_GL2Z V((1), (0), (0), (1));
	matrix_GL2Z S(V), T(V);

	rf.reduce(V);
	rg.reduce(S);

	if (rf == rg) {
		V.invert();
		multiply(S, S, V);
		multiply(U, U, S);
		return true;
	}

	rg2.reduce(T);

	root.assign(rg);
	do {
		rg.rho(S);
		if (!rf.compare(rg)) {
			V.invert();
			multiply(S, S, V);
			multiply(U, U, S);
			return true;
		}

		rg2.inverse_rho(T);
		if (!rf.compare(rg2)) {
			V.invert();
			multiply(T, T, V);
			multiply(U, U, T);
			return true;
		}

	} while (root.compare(rg));

	return false;
}



//
// quadratic_form::prop_equivalent_irregular()
//
// Task:
//      returns true if the irregular forms are properly equivalent
//      and computes the corresponding transformation matrix.
//

bool
quadratic_form::prop_equivalent_irregular (const quadratic_form & g,
					   matrix_GL2Z & U) const
{
	debug_handler("quadratic_form", "prop_equivalent_irregular(matrix_GL2Z)");

	quadratic_form rf(*this), rg(g);
	matrix_GL2Z V ((1), (0), (0), (1));
	matrix_GL2Z W(V);
	bool equiv;

	rf.reduce(V);
	rg.reduce(W);
	if (rf.compare(rg))
		equiv = false;
	else {
		V.invert();
		multiply(W, W, V);
		multiply(U, U, W);
		equiv = true;
	}

	return equiv;
}



//
// quadratic_form::comp_reps_irregular()
//
// Task:
//      computes the representations of N by the irregular form.  If
//      some non-trivial representations exist, true is returned.
//

bool
quadratic_form::comp_reps_irregular (
	sort_vector< pair < bigint, bigint > > & Reps,
	const bigint & N)
{
	debug_handler("quadratic_form", "comp_reps_irregular");

	quadratic_form g;
	bigint NN, temp, temp2, x, y, nx, ny, con;
	lidia_size_t num, j, pos;
	matrix_GL2Z U(1, 0, 0, 1);
	sort_vector< bigint > plist;

	Reps.set_mode(EXPAND);
	Reps.reset();
	plist.set_mode(EXPAND);
	plist.reset();
	g.assign(*this);



	// make form primitive --- content must divide N
	con = content();
	remainder(temp, N, con);
	if (!temp.is_zero()) {
		return false;
	}

	divide(NN, N, con);
	g.assign(a/con, b/con, c/con);


#ifdef LIDIA_DEBUG
	std::cout << "\nIN REPS IRR:  *this = " << *this << ", N = " << N << std::endl;
	std::cout << "g = " << g << std::endl;
	std::cout << "U = \n" << std::endl;
	std::cout << U << std::endl;
#endif

	g.reduce(U);
	num = 0;

#ifdef LIDIA_DEBUG
	std::cout << "g reduced:  " << g << std::endl;
#endif

	//
	// if N = 0, representations are easy
	//

	if (N.is_zero()) {
		Reps[0].assign(0, 0);
		++num;

		y.assign(-1);
		multiply(temp, g.c, y);
		remainder(temp2, temp, g.b);
		while (!temp.is_zero()) {
			dec(y);
			multiply(temp, g.c, y);
			remainder(temp2, temp, g.b);
		}
		divide(x, y, g.b);
		Reps[num].assign(x, y);
		++num;

		return true;
	}


	//
	// if Delta = 0, also easy
	//

	if (Delta.is_zero()) {
		if (!(N % g.c)) {
			divide(y, N, g.c);
			sqrt(temp, y);
			if ((temp*temp) == y) {
				x.assign_zero();
				Reps[num].assign(x, temp);
				++num;

				temp.negate();
				Reps[num].assign(x, temp);
				++num;

				Reps.sort();

				return true;
			}
		}

		return false;
	}


	//
	// general case
	//

	plist = divisors(abs(NN));

#ifdef LIDIA_DEBUG
	std::cout << "plist = " << plist << std::endl;
#endif

	for (j = 0; j < plist.size(); ++j) {
		y.assign(plist[j]);
#ifdef LIDIA_DEBUG
		std::cout << "\nj = " << j << ", trying y = " << y << std::endl;
#endif
		// x = N - Cy^2
		square(temp, y);
		multiply(x, g.c, temp);
		subtract(x, NN, x);

		multiply(temp, g.b, y);
		remainder(temp2, x, temp);
		if (temp2.is_zero()) {
			divide(x, x, temp);
#ifdef LIDIA_DEBUG
			std::cout << "FOUND:  x = " << x << ", y = " << y << std::endl;
#endif
			multiply(temp, U.get_s(), x);
			multiply(temp2, U.get_u(), y);
			add(nx, temp, temp2);
			multiply(temp, U.get_t(), x);
			multiply(temp2, U.get_v(), y);
			add(ny, temp, temp2);
#ifdef LIDIA_DEBUG
			std::cout << "nx = " << nx << ", ny = " << ny << std::endl;
#endif

			if (!Reps.linear_search(pair< bigint, bigint > (nx, ny), pos)) {
				Reps[num].assign(nx, ny);
				++num;
			}
		}
	}

	Reps.sort();

	return (num > 0);
}



//
// constructor:
//    - initialize quadratic order pointer to NULL.
//

quadratic_form::quadratic_form ()
{
	debug_handler("quadratic_form", "quadratic_form()");

	a.assign_zero();
	b.assign_zero();
	c.assign_zero();
	Delta.assign_zero();
	rootD.assign_zero();
	PEA_L.assign_zero();
}



//
// constructor:
//
// Task:
//      set to the form (a, b, c).  If the form is regular, it will point to
//      a quadratic_order of the same discriminant, and such an order will
//      be dynamically created if it doesn't already exist.
//

quadratic_form::quadratic_form (const bigint & a2, const bigint & b2,
				const bigint & c2)
{
	debug_handler("quadratic_form", "quadratic_form(bigint, bigint, bigint)");

	bigint temp;

	a.assign(a2);
	b.assign(b2);
	c.assign(c2);

	square(Delta, b);
	multiply(temp, a, c);
	shift_left(temp, temp, 2);
	subtract(Delta, Delta, temp);

	if (Delta.is_gt_zero())
		sqrt(rootD, Delta);
	else
		rootD.assign_zero();

	newton_root(PEA_L, abs(Delta >> 2), 4);
}



//
// constructor:
//
// Task:
//      set to the form (a, b, c).  If the form is regular, it will point to
//      a quadratic_order of the same discriminant, and such an order will
//      be dynamically created if it doesn't already exist.
//

quadratic_form::quadratic_form (const long a2, const long b2, const long c2)
{
	debug_handler("quadratic_form", "quadratic_form(long, long, long)");

	bigint temp;

	a.assign(a2);
	b.assign(b2);
	c.assign(c2);

	square(Delta, b);
	multiply(temp, a, c);
	shift_left(temp, temp, 2);
	subtract(Delta, Delta, temp);

	if (Delta.is_gt_zero())
		sqrt(rootD, Delta);
	else
		rootD.assign_zero();

	newton_root(PEA_L, abs(Delta >> 2), 4);
}



//
// constructor:
//    initialize with A.
//

quadratic_form::quadratic_form (const qi_class & A)
{
	debug_handler("quadratic_form", "quadratic_form(qi_class)");

	a.assign(A.get_a());
	b.assign(A.get_b());
	c.assign(A.get_c());
	Delta.assign(A.discriminant());
	rootD.assign(A.get_rootD());
	PEA_L.assign(qi_class::PEA_L);
}



//
// constructor:
//    initialize with A.
//

quadratic_form::quadratic_form (const qi_class_real & A)
{
	debug_handler("quadratic_form", "quadratic_form(qi_class_real)");

	a.assign(A.get_a());
	b.assign(A.get_b());
	c.assign(A.get_c());
	Delta.assign(A.discriminant());
	rootD.assign(A.get_rootD());
	PEA_L.assign(qi_class::PEA_L);
}



//
// constructor:
//    initialize with A.
//

quadratic_form::quadratic_form (const quadratic_ideal & A)
{
	debug_handler("quadratic_form", "quadratic_form(quadratic_ideal)");

	a.assign(A.get_a());
	b.assign(A.get_b());
	c.assign(A.get_c());
	Delta.assign(A.discriminant());
	if (Delta.is_gt_zero())
		sqrt(rootD, Delta);
	else
		rootD.assign_zero();

	newton_root(PEA_L, abs(Delta >> 2), 4);
}



//
// constructor:
//    initialize with a copy of f.
//

quadratic_form::quadratic_form (const quadratic_form & f)
{
	debug_handler("quadratic_form", "quadratic_form(quadratic_form)");

	a.assign(f.a);
	b.assign(f.b);
	c.assign(f.c);
	Delta.assign(f.Delta);
	rootD.assign(f.rootD);
	PEA_L.assign(f.PEA_L);
}



//
// destructor:
//    update reference counter of corresponding entry in the list of current
//    quadratic orders.
//

quadratic_form::~quadratic_form()
{
	debug_handler("quadratic_form", "~quadratic_form()");
}



//
// quadratic_form::assign_zero()
//
// Task:
//      set to the zero form (all coefficients 0).
//

void
quadratic_form::assign_zero ()
{
	debug_handler("quadratic_form", "assign_zero");

	a.assign_zero();
	b.assign_zero();
	c.assign_zero();
	Delta.assign_zero();
	rootD.assign_zero();
	PEA_L.assign_zero();
}



//
// quadratic_form::assign_one()
//
// Task:
//      set to the zero form (all coefficients 0).
//

void
quadratic_form::assign_one ()
{
	debug_handler("quadratic_form", "assign_one");

	assign_one(Delta);
}



//
// quadratic_form::assign_one(Delta)
//
// Task:
//    set to the unit form (leading coefficient one) of discriminant Delta.
//    If Delta is not congruent to 0 or 1 modulo 4, the error handler will be
//    called.  If the resulting form is regular, the quadratic_form will
//    point to the quadratic_order of discriminant Delta, and such an order
//    will be dynamically created if it doesn't already exist.
//

void
quadratic_form::assign_one (const bigint & newDelta)
{
	debug_handler("quadratic_form", "assign_one(bigint)");

	bigint tmp(3);
	int temp;

	temp = static_cast<int>(remainder(newDelta, 4));
	if (temp < 0)
		temp += 4;
	if (temp > 1) {
		lidia_error_handler("quadratic_form", "assign_one() - the discriminant must"
				    " be congruent to 0 or 1 modulo 4");
		return;
	}

	if (newDelta.is_zero()) {
		a.assign_one();
		b.assign(2);
		c.assign_one();
		return;
	}

	a.assign_one();
	b = newDelta & tmp;
	square(tmp, b);
	subtract(c, tmp, newDelta);
	shift_right(c, c, 2);

	Delta.assign(newDelta);
	if (Delta.is_gt_zero())
		sqrt(rootD, Delta);
	else
		rootD.assign_zero();

	newton_root(PEA_L, abs(Delta >> 2), 4);

	normalize();
}



//
// quadratic_form::assign(bigint, bigint, bigint)
//
// Task:
//      set to the form (a, b, c).  If the form is regular, it will point to
//      a quadratic_order of the same discriminant, and such an order will
//      be dynamically created if it doesn't already exist.
//

void
quadratic_form::assign (const bigint & a2, const bigint & b2, const bigint & c2)
{
	debug_handler("quadratic_form", "assign(bigint, bigint, bigint)");

	bigint temp;

	a.assign(a2);
	b.assign(b2);
	c.assign(c2);

	square(Delta, b);
	multiply(temp, a, c);
	shift_left(temp, temp, 2);
	subtract(Delta, Delta, temp);

	if (Delta.is_gt_zero())
		sqrt(rootD, Delta);
	else
		rootD.assign_zero();

	newton_root(PEA_L, abs(Delta >> 2), 4);
}



//
// quadratic_form::assign(long, long, long)
//
// Task:
//      set to the form (a, b, c).  If the form is regular, it will point to
//      a quadratic_order of the same discriminant, and such an order will
//      be dynamically created if it doesn't already exist.
//

void
quadratic_form::assign (const long a2, const long b2, const long c2)
{
	debug_handler("quadratic_form", "assign(long, long, long)");

	assign(bigint(a2), bigint(b2), bigint(c2));
}



//
// quadratic_form::assign(qi_class)
//
// Task:
//      set to the form corresponding to A.
//

void
quadratic_form::assign (const qi_class & A)
{
	debug_handler("quadratic_form", "assign(qi_class)");

	a.assign(A.get_a());
	b.assign(A.get_b());
	c.assign(A.get_c());
	Delta.assign(A.discriminant());
	rootD.assign(A.get_rootD());
	PEA_L.assign(qi_class::PEA_L);
}



//
// quadratic_form::assign(qi_class_real)
//
// Task:
//      set to the form corresponding to A.
//

void
quadratic_form::assign (const qi_class_real & A)
{
	debug_handler("quadratic_form", "assign(qi_class_real)");

	a.assign(A.get_a());
	b.assign(A.get_b());
	c.assign(A.get_c());
	Delta.assign(A.discriminant());
	rootD.assign(A.get_rootD());
	PEA_L.assign(qi_class::PEA_L);
}



//
// quadratic_form::assign(quadratic_ideal)
//
// Task:
//      set to the form corresponding to A.
//

void
quadratic_form::assign (const quadratic_ideal & A)
{
	debug_handler("quadratic_form", "assign(quadratic_ideal)");

	a.assign(A.get_a());
	b.assign(A.get_b());
	c.assign(A.get_c());
	Delta.assign(A.discriminant());
	if (Delta.is_gt_zero())
		sqrt(rootD, Delta);
	else
		rootD.assign_zero();

	newton_root(PEA_L, abs(Delta >> 2), 4);
}



//
// quadratic_form::assign(quadratic_form)
//
// Task:
//      set to a copy of g.
//

void
quadratic_form::assign (const quadratic_form & g)
{
	debug_handler("quadratic_form", "assign(quadratic_form)");

	a.assign(g.a);
	b.assign(g.b);
	c.assign(g.c);
	Delta.assign(g.Delta);
	rootD.assign(g.rootD);
	PEA_L.assign(g.PEA_L);
}



//
// operator =
//
// Task:
//    make a copy of an existing quadratic_form
//

quadratic_form & quadratic_form::operator = (const quadratic_form & g)
{
	debug_handler("quadratic_form", "operator = ");

	assign(g);
	return *this;
}



//
// quadratic_form::get_a()
//
// Task:
//      returns coefficient a
//

bigint
quadratic_form::get_a () const
{
	debug_handler("quadratic_form", "get_a() const");
	return a;
}



//
// quadratic_form::get_b()
//
// Task:
//      returns coefficient b
//

bigint
quadratic_form::get_b () const
{
	debug_handler("quadratic_form", "get_b() const");
	return b;
}



//
// quadratic_form::get_c()
//
// Task:
//      returns coefficient c
//

bigint
quadratic_form::get_c () const
{
	debug_handler("quadratic_form", "get_c() const");
	return c;
}



//
// quadratic_form::discriminant()
//
// Task:
//      return the discriminant of the form.
//

bigint
quadratic_form::discriminant () const
{
	debug_handler("quadratic_form", "discriminant() const");
	return Delta;
}



//
// quadratic_form::which_order()
//
// Task:
//      return a pointer to the quadratic_order with the same discriminant
//      (NULL if the form is irregular).
//

quadratic_order &
quadratic_form::which_order () const
{
	debug_handler("quadratic_form", "which_order() const");

	// search for quadratic_order in qo_list
	if (is_regular())
		return *quadratic_order::qo_l().add_to_list(Delta);
	else
		return quadratic_order::zero_QO;
}



//
// which_order(quadratic_form)
//
// Task:
//      return a pointer to the quadratic_order with the same discriminant
//      (NULL if the form is irregular).
//

quadratic_order &
which_order (const quadratic_form & f)
{
	debug_handler("quadratic_form", "which_order(quadratic_form)");

	// search for quadratic_order in qo_list
	if (f.is_regular())
		return *quadratic_order::qo_l().add_to_list(f.Delta);
	else
		return quadratic_order::zero_QO;
}



//
// compose_regular()
//
// Task:
//      compute the composition of two pos def or indef forms.
//

void
compose_regular (quadratic_form & f, const quadratic_form & g1,
		 const quadratic_form & g2)
{
	debug_handler("quadratic_form", "compose_regular");

	bigint Delta, newa, newb, newc, dpr, v, d, w, ab2, temp;

	Delta.assign(g1.Delta);
	dpr.assign(xgcd_left(v, g1.a, g2.a));

	multiply(newb, v, g1.a);
	subtract(temp, g2.b, g1.b);
	multiply(newb, newb, temp);

	multiply(newa, g1.a, g2.a);

	if (!dpr.is_one()) {
		add(ab2, g1.b, g2.b);
		ab2.divide_by_2();
		d.assign(xgcd(v, w, dpr, ab2));

		multiply(newb, newb, v);

		square(temp, g1.b);
		subtract(temp, Delta, temp);
		temp.divide_by_2();
		multiply(temp, temp, w);

		add(newb, newb, temp);
		divide(newb, newb, d);

		square(temp, d);
		divide(newa, newa, temp);

		dpr.assign(d);
	}

	add(newb, newb, g1.b);
	shift_left(ab2, newa, 1);
	remainder(newb, newb, ab2);

	square(newc, newb);
	subtract(newc, newc, Delta);
	shift_left(temp, newa, 2);
	divide(newc, newc, temp);

	f.a.assign(newa);
	f.b.assign(newb);
	f.c.assign(newc);
	f.Delta.assign(Delta);
	f.rootD.assign(g1.rootD);
}



//
// FIX - Check irregular
//

//
// compose()
//
// Task:
//      compute the composition of two forms.  If the forms have different
//      discriminants, the error handler will be evoked.
//

void
compose (quadratic_form & f, const quadratic_form & g1,
	 const quadratic_form & g2)
{
	debug_handler("quadratic_form", "compose");

	bigint Delta, newa, newb, newc, dpr, v, d, w, ab2, temp;

	if (g1.Delta != g2.Delta) {
		lidia_error_handler("quadratic_form", "compose() - forms have different "
				    "discriminants");
		return;
	}

	Delta.assign(g1.Delta);
	if (!g1.is_regular() && ((g1.a == 0) || (g2.a == 0))) {
		quadratic_form h1, h2;

		h1.assign(g1);
		h2.assign(g2);
		h1.reduce();
		h2.reduce();
		multiply(newc, h1.c, h2.c);
		remainder(newc, newc, h1.b);
		newa.assign_zero();
		newb.assign(h1.b);
	}
	else {
		dpr.assign(xgcd_left(v, g1.a, g2.a));

		multiply(newb, v, g1.a);
		subtract(temp, g2.b, g1.b);
		multiply(newb, newb, temp);

		multiply(newa, g1.a, g2.a);

		if (!dpr.is_one()) {
			add(ab2, g1.b, g2.b);
			ab2.divide_by_2();
			d.assign(xgcd(v, w, dpr, ab2));

			multiply(newb, newb, v);

			square(temp, g1.b);
			subtract(temp, Delta, temp);
			temp.divide_by_2();
			multiply(temp, temp, w);

			add(newb, newb, temp);
			divide(newb, newb, d);

			square(temp, d);
			divide(newa, newa, temp);

			dpr.assign(d);
		}

		add(newb, newb, g1.b);
		shift_left(ab2, newa, 1);
		remainder(newb, newb, ab2);

		square(newc, newb);
		subtract(newc, newc, Delta);
		shift_left(temp, newa, 2);
		divide(newc, newc, temp);
	}

	f.a.assign(newa);
	f.b.assign(newb);
	f.c.assign(newc);
	f.Delta.assign(Delta);
	f.rootD.assign(g1.rootD);
}



//
// compose_reduce()
//
// Task:
//      computes a reduced form equivalent to the composition of two forms.
//

void
compose_reduce (quadratic_form & f, const quadratic_form & g1,
		const quadratic_form & g2)
{
	debug_handler("quadratic_form", "compose_reduce");

	compose(f, g1, g2);
	f.reduce();
}



//
// nucomp()
//
// Task:
//      computes a reduced form equivalent to the composition of two forms
//      using the NUCOMP algorithm of Shanks.  If the forms are not positive
//      or negative definite, the error handler will be evoked.
//

void
nucomp (quadratic_form & f, const quadratic_form & g1, const quadratic_form & g2)
{
	debug_handler("quadratic_form", "nucomp");

	bigint a1, b1, c1, a2, b2, c2, s, n, d1, d;
	bigint u, v, u1, l, aa;
	bigint v1, v3;
	bigint b, e, ff, g, q;
	bigint temp;
	bool flag;

	if (g1.Delta != g2.Delta) {
		lidia_error_handler("quadratic_form", "nucomp() - forms have different "
				    "discriminants");
		return;
	}

	if (!g1.is_pos_definite() && !g1.is_neg_definite()) {
		lidia_error_handler("quadratic_form", "nucomp() - forms must be positive "
				    "or negative definite");
		return;
	}

	if (abs(g1.a) < abs(g2.a)) {
		a1.assign(abs(g2.a));
		b1.assign(g2.b);
		c1.assign(abs(g2.c));

		a2.assign(abs(g1.a));
		b2.assign(g1.b);
		c2.assign(abs(g1.c));
	}
	else {
		a1.assign(abs(g1.a));
		b1.assign(g1.b);
		c1.assign(abs(g1.c));

		a2.assign(abs(g2.a));
		b2.assign(g2.b);
		c2.assign(abs(g2.c));
	}

	// initialize
	add(s, b1, b2);
	s.divide_by_2();
	subtract(n, b2, s);

	// first Euclidean pass:  solve d = u a1 + v a2
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
	nugcd(u, d, v1, v3, g1.PEA_L);

	if (u.is_zero()) {
		//  u = 0; d = a1; v1 = 1; v3 = aa
		// b = a2 * d + n * u ) / a1  --->    b = a2 (-)
		// e = s * d + c2 * u ) / a1  --->   e = s  (-)

		// ff = ( b * v3 + n ) / d
		// b1 = 2 * b * v3 + n

		multiply(temp, a2, v3);
		add(ff, temp, n);
		add(b1, temp, ff);
		divide(ff, ff, d);

		// g = ( e * v1 - s ) / u    cannot be computed (u == 0)
		// q = 2 * e * v1 - s       --->   q = s (-)

		// g = ( c2 + s * v3 ) / d

		multiply(g, s, v3);
		add(g, c2, g);
		divide(g, g, d);

		if (flag) {
			multiply(s, s, d1);
			multiply(g, g, d1);
		}

		// a1 = d * b + e * u    --->   a1 = d * b
		multiply(a1, d, a2);

		// c1 = v3 * ff + g * v1  --->   c1 = v3 * ff + g
		multiply(c1, v3, ff);
		add(c1, c1, g);

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

		// ff = ( b * v3 + n ) / d
		// b1 = 2 * b * v3 + n
		multiply(q, b, v3);
		add(ff, q, n);
		add(b1, q, ff);
		divide(ff, ff, d);

		// g = ( e * v1 - s ) / u
		// q = 2 * e * v1 - s
		multiply(q, e, v1);
		subtract(g, q, s);
		add(q, q, g);
		divide(g, g, u);

		if (flag) {
			multiply(v1, v1, d1);
			multiply(u, u, d1);
			multiply(q, q, d1);
		}

		// a1 = d * b + e * u
		multiply(temp, d, b);
		multiply(a1, e, u);
		add(a1, a1, temp);

		// c1 = v3 * ff + g * v1
		multiply(temp, v3, ff);
		multiply(c1, g, v1);
		add(c1, c1, temp);

		// b1 = b1 + q
		add(b1, b1, q);
	}

	if (g1.is_pos_definite()) {
		f.a.assign(a1);
		f.b.assign(b1);
		f.c.assign(c1);
	}
	else {
		negate(f.a, a1);
		f.b.assign(b1);
		negate(f.c, c1);
	}
	f.Delta.assign(g1.Delta);
	f.rootD.assign(g1.rootD);
	f.reduce();
}



//
// quadratic_form::conjugate()
//
// Task:
//      compute the conjugate of the form.
//

void
quadratic_form::conjugate ()
{
	debug_handler("quadratic_form", "conjugate");

	b.negate();
}



//
// get_conjugate(quadratic_form, quadratic_form)
//
// Task:
//      compute the conjugate of g.
//

void
get_conjugate (quadratic_form & f, const quadratic_form & g)
{
	debug_handler("quadratic_form", "get_conjugate(quadratic_form, "
		      "quadratic_form)");

	f.assign(g);
	f.b.negate();
}



//
// get_conjugate(quadratic_form)
//
// Task:
//      return the conjugate of f.
//

quadratic_form
get_conjugate (const quadratic_form & f)
{
	debug_handler("quadratic_form", "get_conjugate(quadratic_form)");

	quadratic_form g;

	g.assign(f);
	g.b.negate();

	return g;
}



//
// divide()
//
// Task:
//      compute the composition of the g1 with the conjugate of g2.
//

void
divide (quadratic_form & f, const quadratic_form & g1, const quadratic_form & g2)
{
	debug_handler("quadratic_form", "divide");

	if (g2.is_zero()) {
		lidia_error_handler("quadratic_form", "divide() - zero divisor");
		return;
	}

	compose(f, g1, get_conjugate(g2));
}



//
// divide_reduce()
//
// Task:
//      compute a reduced form equivalent to the composition of the g1 with
//      the conjugate of g2.
//

void
divide_reduce (quadratic_form & f, const quadratic_form & g1,
	       const quadratic_form & g2)
{
	debug_handler("quadratic_form", "divide_reduce");

	if (g2.is_zero()) {
		lidia_error_handler("quadratic_form", "divide_reduce() - zero divisor");
		return;
	}

	compose_reduce(f, g1, get_conjugate(g2));
}



//
// FIX - check irregular
//

//
// square()
//
// Task:
//      compute the composition of g with itself.
//

void
square (quadratic_form & f, const quadratic_form & g)
{
	debug_handler("quadratic_form", "square");

	bigint Delta, newb, d, w, temp, newa, newc;

	Delta.assign(g.Delta);
	if (!g.is_regular() && (g.a.is_zero())) {
		quadratic_form h;

		h.assign(g);
		h.reduce();
		square(newc, h.c);
		remainder(newc, newc, h.b);
		newa.assign_zero();
		newb.assign(h.b);
	}
	else {
		d.assign(xgcd_right(w, g.a, g.b));

		square(newb, g.b);
		subtract(newb, Delta, newb);
		shift_left(temp, d, 1);
		divide(newb, newb, temp);
		multiply(newb, newb, w);
		add(newb, newb, g.b);

		divide(temp, g.a, d);
		square(newa, temp);

		shift_left(temp, newa, 1);
		remainder(newb, newb, temp);

		square(newc, newb);
		subtract(newc, newc, Delta);
		shift_left(temp, newa, 2);
		divide(newc, newc, temp);
	}

	f.a.assign(newa);
	f.b.assign(newb);
	f.c.assign(newc);
	f.Delta.assign(Delta);
	f.rootD.assign(g.rootD);
}



//
// square_reduce()
//
// Task:
//      compute a reduced form equivalent to the composition of g with itself.
//

void
square_reduce (quadratic_form & f, const quadratic_form & g)
{
	debug_handler("quadratic_form", "square_reduce");

	square(f, g);
	f.reduce();
}



//
// nudupl()
//
// Task:
//      computes a reduced form equivalent to the composition of g with itself
//      using the NUDUPL algorithm of Shanks.  If the form is not positive
//      or negative definite, the error handler will be evoked.
//

void
nudupl (quadratic_form & f, const quadratic_form & g)
{
	debug_handler("quadratic_form", "nudupl");

	bigint a, b, c;
	bigint temp, u, d1, d, v1, v3, e, gg;
	bool flag;

	if (!g.is_pos_definite() && !g.is_neg_definite()) {
		lidia_error_handler("quadratic_form", "nudupl() - form must be positive "
				    "or negative definite");
		return;
	}

	a.assign(abs(g.a));
	b.assign(g.b);
	c.assign(abs(g.c));

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
	nugcd(u, d, v1, v3, g.PEA_L);

	// final squaring
	if (u.is_zero()) {
		multiply(gg, b, v3);
		add(gg, gg, c);

		if (flag) {
			multiply(b, b, d1);
			multiply(gg, gg, d1);
		}

		divide(gg, gg, d);

		square(a, d);
		square(c, v3);

		add(temp, d, v3);
		square(temp, temp);
		add(b, b, temp);
		subtract(b, b, a);
		subtract(b, b, c);

		add(c, c, gg);
	}
	else {
		// u != 0
		multiply(e, b, d);
		multiply(temp, c, u);
		add(e, e, temp);
		divide(e, e, a);

		multiply(temp, e, v1);
		subtract(gg, temp, b);
		add(b, temp, gg);
		divide(gg, gg, u);

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

		multiply(temp, gg, v1);
		add(c, c, temp);
	}

	if (g.is_pos_definite()) {
		f.a.assign(a);
		f.b.assign(b);
		f.c.assign(c);
	}
	else {
		negate(f.a, a);
		f.b.assign(b);
		negate(f.c, c);
	}
	f.Delta.assign(g.Delta);
	f.rootD.assign(g.rootD);
	f.reduce();
}



//
// power()
//
// Task:
//      computes g^i using binary exponentiation.
//

void
power (quadratic_form & f, const quadratic_form & g, const bigint & i)
{
	debug_handler("quadratic_form", "power(bigint)");

	quadratic_form h;
	bigint j;

	h.assign(g);
	j.assign(i);
	if (j.is_lt_zero()) {
		h.conjugate();
		j.absolute_value();
	}
	f.assign_one(g.Delta);
	while (j.is_gt_zero()) {
		if (j.is_odd())
			compose(f, f, h);
		j.divide_by_2();
		if (j.is_gt_zero())
			square(h, h);
	}
}



//
// power()
//
// Task:
//      computes g^i using binary exponentiation
//

void
power (quadratic_form & f, const quadratic_form & g, const long i)
{
	debug_handler("quadratic_form", "power(long)");

	quadratic_form h;
	register long j;

	h.assign(g);
	j = i;
	if (j < 0) {
		h.conjugate();
		j = -j;
	}
	f.assign_one(g.Delta);
	while (j > 0) {
		if ((j & 1) == 1)
			compose(f, f, h);
		j >>= 1;
		if (j > 0)
			square(h, h);
	}
}



//
// power_reduce()
//
// Task:
//      computes a reduced form equivalent to g^i.
//

void
power_reduce (quadratic_form & f, const quadratic_form & g, const bigint & i)
{
	debug_handler("quadratic_form", "power_reduce(bigint)");

	quadratic_form h;
	bigint j;

	h.assign(g);
	j.assign(i);
	if (j.is_lt_zero()) {
		h.conjugate();
		j.absolute_value();
	}
	f.assign_one(g.Delta);
	while (j.is_gt_zero()) {
		if (j.is_odd())
			compose_reduce(f, f, h);
		j.divide_by_2();
		if (j.is_gt_zero())
			square_reduce(h, h);
	}
}



//
// power_reduce()
//
// Task:
//      computes a reduced form equivalent to g^i.
//

void
power_reduce (quadratic_form & f, const quadratic_form & g, const long i)
{
	debug_handler("quadratic_form", "power_reduce(long)");

	quadratic_form h;
	register long j;

	h.assign(g);
	j = i;
	if (j < 0) {
		h.conjugate();
		j = -j;
	}
	f.assign_one(g.Delta);
	while (j > 0) {
		if ((j & 1) == 1)
			compose_reduce(f, f, h);
		j >>= 1;
		if (j > 0)
			square_reduce(h, h);
	}
}



//
// nupower()
//
// Task:
//      computes a reduced form equivalent to g^i using binary exponentiation
//      and the nucomp and nudupl algorithms.  If the form is not positive or
//      negative definite, the error handler will be evoked.
//

void
nupower (quadratic_form & f, const quadratic_form & g, const bigint & i)
{
	debug_handler("quadratic_form", "nupower(qf&, const qf&, const bigint&)");

	quadratic_form h;
	bigint j;

	if (!g.is_pos_definite() && !g.is_neg_definite()) {
		lidia_error_handler("quadratic_form", "nupower() - form must be positive "
				    "or negative definite");
		return;
	}

	h.assign(g);
	j.assign(i);
	if (j.is_lt_zero()) {
		h.conjugate();
		j.absolute_value();
	}
	f.assign_one(g.Delta);
	while (j.is_gt_zero()) {
		if (j.is_odd())
			nucomp(f, f, h);
		j.divide_by_2();
		if (j.is_gt_zero())
			nudupl(h, h);
	}
}



//
// nupower()
//
// Task:
//      computes a reduced form equivalent to g^i using binary exponentiation
//      and the nucomp and nudupl algorithms.  If the form is not positive or
//      negative definite, the error handler will be evoked.
//

void
nupower (quadratic_form & f, const quadratic_form & g, const long i)
{
	debug_handler("quadratic_form", "nupower(long)");

	quadratic_form h;
	register long j;

	if (!g.is_pos_definite() && !g.is_neg_definite()) {
		lidia_error_handler("quadratic_form", "nupower() - form must be positive "
				    "or negative definite");
		return;
	}

	h.assign(g);
	j = i;
	if (j < 0) {
		h.conjugate();
		j = -j;
	}
	f.assign_one(g.Delta);
	while (j > 0) {
		if ((j & 1) == 1)
			nucomp(f, f, h);
		j >>= 1;
		if (j > 0)
			nudupl(h, h);
	}
}



//
// operator -
//
// Task:
//      returns the conjugate of f.
//

quadratic_form operator - (const quadratic_form & f)
{
	debug_handler("quadratic_form", "operator -");

	return get_conjugate(f);
}



//
// operator *
//
// Task:
//      returns the composition of f and g.
//

quadratic_form operator * (const quadratic_form & f, const quadratic_form & g)
{
	debug_handler("quadratic_form", "operator *");

	quadratic_form h;

	compose(h, f, g);
	return h;
}



//
// operator /
//
// Task:
//      returns the composition of f and the conjugate of g.
//

quadratic_form operator / (const quadratic_form & f, const quadratic_form & g)
{
	debug_handler("quadratic_form", "operator /");

	quadratic_form h;

	divide(h, f, g);
	return h;
}



//
// operator *=
//
// Task:
//      *this = *this * f
//

quadratic_form & quadratic_form::operator *= (const quadratic_form & f)
{
	debug_handler("quadratic_form", "operator *= ");

	compose(*this, *this, f);
	return *this;
}



//
// operator /=
//
// Task:
//      *this  = *this / f
//

quadratic_form & quadratic_form::operator /= (const quadratic_form & f)
{
	debug_handler("quadratic_form", "operator /= ");

	compose(*this, *this, get_conjugate(f));
	return *this;
}



//
// quadratic_form::is_zero()
//
// Task:
//      tests if the form is the zero form.
//

bool
quadratic_form::is_zero () const
{
	debug_handler("quadratic_form", "is_zero");

	return ((a.is_zero()) && (b.is_zero()) && (c.is_zero()));
}



//
// quadratic_form::is_one()
//
// Task:
//      tests if the form is the unit form.
//

bool
quadratic_form::is_one () const
{
	debug_handler("quadratic_form", "is_one");

	quadratic_form U;

	U.assign_one(Delta);
	return (is_equal(U));
}



//
// quadratic_form::is_equal()
//
// Task:
//      tests if the forms are equal
//

bool
quadratic_form::is_equal (const quadratic_form & g) const
{
	debug_handler("quadratic_form", "is_equal");

	return (!(a.compare(g.a)) && !(b.compare(g.b)) && !(c.compare(g.c)));
}



//
// quadratic_form::compare()
//
// Task:
//      lexicographically compares the forms
//

int
quadratic_form::compare (const quadratic_form & g) const
{
	debug_handler("quadratic_form", "compare");

	int i = a.compare(g.a);
	if (i == 0)
		i = b.compare(g.b);
	if (i == 0)
		i = c.compare(g.c);

	return i;
}



//
// quadratic_form::abs_compare()
//
// Task:
//      lexicographically compares the forms using the absolute values of the
//      actuall coefficients.
//

int
quadratic_form::abs_compare (const quadratic_form & g) const
{
	debug_handler("quadratic_form", "abs_compare");

	int i = a.abs_compare(g.a);
	if (i == 0)
		i = b.abs_compare(g.b);
	if (i == 0)
		i = c.abs_compare(g.c);

	return i;
}



//
// operator ==
//
// Task:
//      tests if f and g are equal.
//

bool operator == (const quadratic_form & f, const quadratic_form & g)
{
	debug_handler("quadratic_form", "operator == ");

	return f.is_equal(g);
}



//
// operator !=
//
// Task:
//      tests if f and g are not equal.
//

bool operator != (const quadratic_form & f, const quadratic_form & g)
{
	debug_handler("quadratic_form", "operator != ");

	return !f.is_equal(g);
}



//
// operator <=
//
// Task:
//      tests if f is lexicographically less than or equal to g.
//

bool operator <= (const quadratic_form & f, const quadratic_form & g)
{
	debug_handler("quadratic_form", "operator <= ");

	return (f.compare(g) <= 0);
}



//
// operator <
//
// Task:
//      tests if f is lexicographically less than g.
//

bool operator < (const quadratic_form & f, const quadratic_form & g)
{
	debug_handler("quadratic_form", "operator < ");

	return (f.compare(g) < 0);
}



//
// operator >=
//
// Task:
//      tests if f is lexicographically greater than or equal to g.
//

bool operator >= (const quadratic_form & f, const quadratic_form & g)
{
	debug_handler("quadratic_form", "operator >= ");

	return (f.compare(g) >= 0);
}



//
// operator >
//
// Task:
//      tests if f is lexicographically greater than g.
//

bool operator > (const quadratic_form & f, const quadratic_form & g)
{
	debug_handler("quadratic_form", "operator >");


	return (f.compare(g) > 0);
}



//
// operator !
//
// Task:
//      tests if f is the zero form.
//

bool operator ! (const quadratic_form & f)
{
	debug_handler("quadratic_form", "operator !");

	return (f.is_zero());
}



//
// swap()
//
// Task:
//      swap f and g.
//

void
swap (quadratic_form & f, quadratic_form & g)
{
	debug_handler("quadratic_form", "swap");

	quadratic_form h;

	h.assign(f);
	f.assign(g);
	g.assign(h);
}



//
// quadratic_form::definiteness()
//
// Task:
//      returns an integer representing the "definiteness" of the form:
//         2 - irregular
//         1 - positive definite
//         0 - indefinite
//         -1 - negative definite
//

int
quadratic_form::definiteness () const
{
	debug_handler("quadratic_form", "definiteness");

	if (!is_regular())
		return 2;

	if (Delta.is_positive())
		return 0;

	if (a.is_positive())
		return 1;

	return -1;
}



//
// quadratic_form::is_pos_definite()
//
// Task:
//      returns true if the form is positive definite, false otherwise.
//

bool
quadratic_form::is_pos_definite () const
{
	debug_handler("quadratic_form", "is_pos_definite");

	return (Delta.is_negative() && a.is_positive() && is_regular());
}



//
// quadratic_form::is_pos_semidefinite()
//
// Task:
//      returns true if the form is positive semidefinite, false otherwise.
//

bool
quadratic_form::is_pos_semidefinite () const
{
	debug_handler("quadratic_form", "is_pos_semidefinite");

	return (Delta.is_le_zero() && (a.is_positive() || c.is_positive())
		&& is_regular());
}



//
// quadratic_form::is_neg_definite()
//
// Task:
//      returns true if the form is negative definite, false otherwise.
//

bool
quadratic_form::is_neg_definite () const
{
	debug_handler("quadratic_form", "is_neg_definite");

	return (Delta.is_negative() && a.is_negative() && is_regular());
}



//
// quadratic_form::is_neg_semidefinite()
//
// Task:
//      returns true if the form is negative semidefinite, false otherwise.
//

bool
quadratic_form::is_neg_semidefinite () const
{
	debug_handler("quadratic_form", "is_neg_semidefinite");

	return (Delta.is_le_zero() && (a.is_negative() || c.is_negative())
		&& is_regular());
}



//
// quadratic_form::is_indefinite()
//
// Task:
//      returns true if the form is indefinite, false otherwise.
//

bool
quadratic_form::is_indefinite () const
{
	debug_handler("quadratic_form", "is_indefinite");

	return (Delta.is_positive() && is_regular());
}



//
// quadratic_form::is_regular()
//
// Task:
//      returns true if the form is regular, false otherwise.
//

bool
quadratic_form::is_regular () const
{
	debug_handler("quadratic_form", "is_regular");

	return is_regular_disc(Delta);
}



//
// quadratic_form::content()
//
// Task:
//      returns the content of the form (gcd of coefficients)
//

bigint
quadratic_form::content () const
{
	debug_handler("quadratic_form", "content");

	bigint temp;

	temp = gcd(a, b);
	return gcd(temp, c);
}



//
// quadratic_form::is_primitive()
//
// Task:
//      returns true if the form is primitive (content = 1), false otherwise.
//

bool
quadratic_form::is_primitive () const
{
	debug_handler("quadratic_form", "is_primitive");

	bigint f;

	f = content();
	return (f.is_one());
}



//
// quadratic_form::eval()
//
// Task:
//      evaluates the form at (x,y)
//

bigint
quadratic_form::eval (const bigint & x, const bigint & y) const
{
	debug_handler("quadratic_form", "eval");

	bigint temp1, temp2, temp3;

	square(temp2, x);
	multiply(temp3, temp2, a);
	multiply(temp2, x, y);
	multiply(temp2, temp2, b);
	add(temp1, temp2, temp3);
	square(temp2, y);
	multiply(temp3, temp2, c);
	add(temp1, temp1, temp3);

	return temp1;
}



//
// operator ()
//
// Task:
//      evaluates the form at (x,y)
//

bigint quadratic_form::operator () (const bigint &x, const bigint &y)
{
	debug_handler("quadratic_form", "operator ()");

	bigint temp1, temp2, temp3;

	square(temp2, x);
	multiply(temp3, temp2, a);
	multiply(temp2, x, y);
	multiply(temp2, temp2, b);
	add(temp1, temp2, temp3);
	square(temp2, y);
	multiply(temp3, temp2, c);
	add(temp1, temp1, temp3);

	return temp1;
}



//
// quadratic_form::transform(matrix_GL2Z)
//
// Task:
//      transforms the form using the transformation matrix U.
//

void
quadratic_form::transform (const matrix_GL2Z & U)
{
	debug_handler("quadratic_form", "transform(matrix_GL2Z)");

	bigint temp1, temp2, temp3;

	multiply(temp1, U.get_s(), U.get_u());
	multiply(temp1, temp1, a); //asu

	multiply(temp2, U.get_t(), U.get_v());
	multiply(temp2, temp2, c); //ctv
	add(temp1, temp1, temp2);
	temp1.multiply_by_2(); //2(asu+tvc)

	multiply(temp2, U.get_s(), U.get_v());
	multiply(temp3, U.get_t(), U.get_u());
	add(temp2, temp2, temp3);
	multiply(temp2, temp2, b); //b(sv+tu)
	add (temp1, temp1, temp2); //new b


	temp2.assign(eval(U.get_s(), U.get_t())); //new a
	temp3.assign(eval(U.get_u(), U.get_v())); //new c

	matrix_GL2Z V(U);
	if (V.det() == -1) {
		temp2.negate();
		temp1.negate();
		temp3.negate();
	}

	a.assign(temp2);
	b.assign(temp1);
	c.assign(temp3);
}



//
// generate_prime_form()
//
// Task:
//      generates a prime form of the same discriminant
//

bool
generate_prime_form (quadratic_form & f, const bigint & p)
{
	debug_handler("quadratic_form", "generate_prime_form");

	return generate_prime_form(f, p, f.Delta);
}



//
// generate_prime_form()
//
// Task:
//      computes a form of discriminant newDelta with first coefficient p.  If
//      no such form exists, false is returned.  If the form is regular, it
//      will point to a quadratic_order of the same discriminant, and such an
//      order will be dynamically created if it doesn't exist already.
//

bool
generate_prime_form (quadratic_form & f, const bigint & p,
		     const bigint & newDelta)
{
	debug_handler("quadratic_form", "generate_prime_form(bigint)");

	bigint b, c, temp;
	int kro, Dp, ip;

	Dp = static_cast<int>(remainder(newDelta, 4));
	if (Dp < 0)
		Dp += 4;
	if (Dp > 1) {
		lidia_error_handler("quadratic_form", "prime_form() - newDelta must "
				    "be congruent to 0 or 1 modulo 4");
		return false;
	}

	if (newDelta.is_zero()) {
		return false;
	}

	kro = kronecker(newDelta, p);
	if (kro < 0)
		return false;
	else {
		if (p == 2) {
			if (kro == 0) {
				Dp = static_cast<int>(remainder(newDelta, 8));
				if (Dp < 0)
					Dp += 8;
				if (Dp == 0)
					b.assign_zero();
				else
					b.assign(2);
			}
			else
				b.assign(1);
		}
		else {
			if (kro == 0) {
				remainder(temp, newDelta, p*p);
				if (temp.is_zero())
					b.assign(p);
				else
					b.assign_zero();
			}
			else {
				if (is_int(p)) {
					p.intify(ip);
					Dp = static_cast<int>(remainder(newDelta, static_cast<long>(ip)));
					if (Dp < 0)
						Dp += ip;

					b.assign(ressol(Dp, ip));
				}
				else {
					remainder(temp, newDelta, p);
					if (temp.is_lt_zero())
						add(temp, temp, p);
					ressol(b, temp, p);
				}

				if (b.is_lt_zero())
					add(b, b, p);
			}

			if ((newDelta.is_odd()) != (b.is_odd()))
				subtract(b, p, b);
		}

		square(c, b);
		subtract(c, c, newDelta);
		shift_left(temp, p, 2);
		divide(c, c, temp);

		f.a.assign(p);
		f.b.assign(b);
		f.c.assign(c);
		f.Delta.assign(newDelta);
		if (f.Delta.is_gt_zero())
			sqrt(f.rootD, f.Delta);
		else
			f.rootD.assign_zero();

		newton_root(f.PEA_L, abs(f.Delta >> 2), 4);

		return true;
	}
}



//
// quadratic_form::is_normal()
//
// Task:
//      returns true if the form is normal, false otherwise.
//

bool
quadratic_form::is_normal () const
{
	debug_handler("quadratic_form", "is_normal");

	bigint temp;
	int dness;

	dness = definiteness();

	if (dness == 2)
		return is_reduced();
	if (dness == 1)
		return (b.compare(-a) > 0) && (b.compare(a) <= 0);
	else if (dness == -1)
		return (b.compare(-abs(a)) > 0) && (b.compare(abs(a)) <= 0);
	else {
		if (a.compare(rootD) <= 0) {
			shift_left(temp, a, 1);
			subtract(temp, rootD, temp);
			return (temp.compare(b) < 0) && (b.compare(rootD) <= 0);
		}
		else
			return (b.compare(-a) > 0) && (b.compare(a) <= 0);
	}
}



//
// quadratic_form::normalize_regular()
//
// Task:
//      nomalizes the pos def or indef form (no type-test performed).
//

void
quadratic_form::normalize_regular ()
{
	debug_handler("quadratic_form", "normalize_regular()");

	bigint temp, s;

	// test if already normalized
	if (Delta.is_negative()) {
		if ((b.compare(-a) > 0) && (b.compare(a) <= 0))
			return;
		s = norm_number_pos_def();
	}
	else {
		if (a.compare(rootD) <= 0) {
			shift_left(temp, a, 1);
			subtract(temp, rootD, temp);
			if ((temp.compare(b) < 0) && (b.compare(rootD) <= 0))
				return;
		}
		else
			if ((b.compare(-a) > 0) && (b.compare(a) <= 0))
				return;
		s = norm_number_indef();
	}

	multiply(temp, a, s);
	add(temp, temp, b);
	multiply(temp, temp, s);
	add(c, c, temp); //c = as^2 +bs +c

	multiply (temp, s, a);
	temp.multiply_by_2();
	add(b, b, temp); //b = b+2sa
}



//
// quadratic_form::normalize()
//
// Task:
//      nomalizes the form.
//

void
quadratic_form::normalize ()
{
	debug_handler("quadratic_form", "normalize()");

	bigint tmp1, tmp2, s; //tmp1 becomes new b tmp2 new c
	quadratic_form f;

	if (is_normal())
		return;

	if (!is_regular())
		reduce_irregular();
	else if (is_neg_definite()) {
		f.assign(-a, b, -c);
		f.normalize();
		negate(a, f.a);
		b = f.b;
		negate(c, f.c);
	}
	else {
		s = norm_number();
		multiply (tmp1, s, a);
		tmp1.multiply_by_2();
		add(tmp1, b, tmp1); //tmp1 = b+2sa

		multiply(tmp2, a, s);
		add(tmp2, tmp2, b);
		multiply(tmp2, tmp2, s);
		add(tmp2, tmp2, c); //tmp2 = as^2 +bs +c

		b.assign(tmp1);
		c.assign(tmp2);
	}
}



//
// quadratic_form::normalize(matrix_GL2Z)
//
// Task:
//      nomalizes the form and computes the corresponding transformation matrix
//

void
quadratic_form::normalize (matrix_GL2Z & U)
{
	debug_handler("quadratic_form", "normalize(matrix_GL2Z)");

	bigint tmp1, tmp2, s; //tmp1 becomes new b tmp2 new c
	quadratic_form f;

	if (is_normal())
		return;

	if (!is_regular()) {
		matrix_GL2Z V((1), (0), (0), (1));
		reduce_irregular(V);
		multiply(U, U, V);
	}
	else if (is_neg_definite()) {
		matrix_GL2Z V((1), (0), (0), (-1));
		f.assign(-a, b, -c);
		multiply(U, U, V);
		f.normalize(U);
		multiply(U, U, V);
		negate(a, f.a);
		b = f.b;
		negate(c, f.c);
	}
	else {
		s = norm_number();

		matrix_GL2Z S((1), s, (0), (1));

		multiply (tmp1, s, a);
		tmp1.multiply_by_2();
		add(tmp1, b, tmp1); //tmp1 = b+2sa

		multiply(tmp2, a, s);
		add(tmp2, tmp2, b);
		multiply(tmp2, tmp2, s);
		add(tmp2, tmp2, c); //tmp2 = as^2 +bs +c

		b.assign(tmp1);
		c.assign(tmp2);

		multiply(U, U, S);
	}
}



//
// quadratic_form::is_reduced()
//
// Task:
//      returns true if the form is reduced, false otherwise.
//

bool
quadratic_form::is_reduced () const
{
	debug_handler("quadratic_form", "is_reduced");

	int dness;
	quadratic_form f;

	dness = definiteness();

	if (dness == 2)
		return is_reduced_irregular();
	else if (dness == 0)
		return is_reduced_indef();
	else if (dness == 1)
		return is_reduced_pos_def();
	else {
		f.assign(-a, b, -c);
		return f.is_reduced_pos_def();
	}
}



//
// quadratic_form::reduce_regular()
//
// Task:
//      reduces the pos def or indef form (no type test performed)
//

void
quadratic_form::reduce_regular ()
{
	debug_handler("quadratic_form", "reduce_regular()");

	if (Delta.is_lt_zero()) {
		if (!is_reduced_pos_def())
			reduce_pos_def();
	}
	else {
		if (!is_reduced_indef())
			reduce_indef();
	}
}



//
// quadratic_form::reduce()
//
// Task:
//      reduces the form
//

void
quadratic_form::reduce ()
{
	debug_handler("quadratic_form", "reduce()");

	quadratic_form f;

	if (is_reduced())
		return;

	if (!is_regular())
		reduce_irregular();
	else if (is_neg_definite()) {
		a.negate();
		c.negate();
		normalize();
		reduce_pos_def();
		a.negate();
		c.negate();
	}
	else if (is_pos_definite()) {
		normalize();
		reduce_pos_def();
	}
	else {
		normalize();
		reduce_indef();
	}
}



//
// quadratic_form::reduce(matrix_GL2Z)
//
// Task:
//      reduces the form and computes the corresponding transformation matrix.
//

void
quadratic_form::reduce (matrix_GL2Z & U)
{
	debug_handler("quadratic_form", "reduce(matrix_GL2Z)");

	quadratic_form f;

	if (is_reduced())
		return;

	if (!is_regular())
		reduce_irregular(U);
	else if (is_neg_definite()) {
		matrix_GL2Z V((1), (0), (0), (-1));
		a.negate();
		c.negate();
		multiply(U, U, V);
		normalize(U);
		reduce_pos_def(U);
		multiply(U, U, V);
		a.negate();
		c.negate();
	}
	else if (is_pos_definite()) {
		normalize(U);
		reduce_pos_def(U);
	}
	else {
		normalize(U);
		reduce_indef(U);
	}
}



//
// quadratic_form::rho()
//
// Task:
//      applies one reduction step to the form.
//

void
quadratic_form::rho ()
{
	debug_handler("quadratic_form", "rho()");

	bigint temp;

	if (!is_indefinite()) {
		if (is_regular() && !is_reduced()) {
			if (is_neg_definite()) {
				negate(temp, a);
				negate(a, c);
				c.assign(temp);
				b.negate();
				normalize();
				a.negate();
				c.negate();
			}
			else {
				temp.assign(a);
				a.assign(c);
				c.assign(temp);
				b.negate();
				normalize();
			}
		}
	}
	else {
		temp.assign(a);
		a.assign(c);
		c.assign(temp);

		inverse_rho();

		temp.assign(a);
		a.assign(c);
		c.assign(temp);
		normalize();
	}
}



//
// quadratic_form::rho(matrix_GL2Z)
//
// Task:
//      applies one reduction step to the form and computes the corresponding
//      transformation matrix.
//

void
quadratic_form::rho (matrix_GL2Z & U)
{
	debug_handler("quadratic_form", "rho(matrix_GL2Z)");

	bigint temp;
	matrix_GL2Z V((0), (1), (1), (0)), W((1), (0), (0), (-1));
	matrix_GL2Z V1((0), (1), (-1), (0));

	if (!is_indefinite()) {
		if (is_regular() && !is_reduced()) {
			if (is_neg_definite()) {
				multiply(U, U, W);
				a.negate();
				c.negate();

				multiply(U, U, V1);
				temp.assign(a);
				a.assign(c);
				c.assign(temp);
				b.negate();

				normalize(U);

				multiply(U, U, W);
				a.negate();
				c.negate();
			}
			else {
				multiply(U, U, V1);
				temp.assign(a);
				a.assign(c);
				c.assign(temp);
				b.negate();
				normalize(U);
			}
		}
	}
	else {
		multiply(U, U, V);
		temp.assign(a);
		a.assign(c);
		c.assign(temp);

		inverse_rho(U);

		multiply(U, U, V);
		temp.assign(a);
		a.assign(c);
		c.assign(temp);
		normalize(U);
	}
}



//
// quadratic_form::rho_indef()
//
// Task:
//      applies one reduction step to the indefinite form.
//

void
quadratic_form::rho_indef ()
{
	debug_handler("quadratic_form", "rho_indef()");

	bigint temp;

	temp.assign(a);
	a.assign(c);
	c.assign(temp);
	b.negate();
	normalize_regular();
}



//
// quadratic_form::inverse_rho()
//
// Task:
//      applies one inverse reduction step to the form.
//

void
quadratic_form::inverse_rho ()
{
	debug_handler("quadratic_form", "inverse_rho()");

	bigint temp;

	if (is_indefinite() || (is_regular() && !is_reduced())) {
		if (is_neg_definite()) {
			negate(temp, a);
			negate(a, c);
			c.assign(temp);
			b.negate();
			normalize();
			a.negate();
			c.negate();
		}
		else {
			temp.assign(a);
			a.assign(c);
			c.assign(temp);
			b.negate();
			normalize();
		}
	}
}



//
// quadratic_form::inverse_rho(matrix_GL2Z)
//
// Task:
//      applies one inverse reduction step to the form and computes the
//      corresponding transformation matrix.
//

void
quadratic_form::inverse_rho (matrix_GL2Z & U)
{
	debug_handler("quadratic_form", "inverse_rho(matrix_GL2Z)");

	bigint temp;
	matrix_GL2Z V((0), (1), (-1), (0));

	if (is_indefinite() || (is_regular() && !is_reduced())) {
		if (is_neg_definite()) {
			matrix_GL2Z W((1), (0), (0), (-1));

			multiply(U, U, W);
			a.negate();
			c.negate();

			multiply(U, U, V);
			temp.assign(a);
			a.assign(c);
			c.assign(temp);
			b.negate();

			normalize(U);

			multiply(U, U, W);
			a.negate();
			c.negate();
		}
		else {
			multiply(U, U, V);
			temp.assign(a);
			a.assign(c);
			c.assign(temp);
			b.negate();
			normalize(U);
		}
	}
}



//
// quadratic_form::is_equivalent()
//
// Task:
//      returns true if the forms are equivalent and computes the corresponding
//      transformation matrix.
//

bool
quadratic_form::is_equivalent (const quadratic_form & g, matrix_GL2Z & U) const
{
	debug_handler("quadratic_form", "is_equivalent(matrix_GL2Z)");

	int dness;
	quadratic_form hf(*this), hg(g);
	matrix_GL2Z W((-1), (0),
		      (0), (1));
	matrix_GL2Z I((1), (0),
		      (0), (1));

	if (Delta != g.Delta)
		return false;

	dness = definiteness();

	if (is_pos_definite() && g.is_pos_definite())
		return prop_equivalent_pos_def(g, U);

	if (is_neg_definite() && g.is_neg_definite()) {
		matrix_GL2Z H(I);
		hf.transform(W);
		hg.transform(W);
		if (hf.prop_equivalent_pos_def(hg, H)) {
			multiply(H, W, H);
			multiply(H, H, inverse(W));
			multiply(U, U, H);
			return true;
		}
		else
			return false;
	}

	if (is_neg_definite() && g.is_pos_definite()) {
		matrix_GL2Z H(I);
		hf.transform(W);
		if (hf.prop_equivalent_pos_def(hg, H)) {
			multiply(H, W, H);
			multiply(U, U, H);
			return true;
		}
		else
			return false;
	}

	if (is_pos_definite() && g.is_neg_definite()) {
		matrix_GL2Z H(I);
		hg.transform(W);
		if (hf.prop_equivalent_pos_def(hg, H)) {
			multiply(H, H, inverse(W));
			multiply(U, U, H);
			return true;
		}
		else
			return false;
	}


	if (dness == 0) {
		matrix_GL2Z V((1), (0),
			      (0), (1));
		matrix_GL2Z S(V), T(V);

		quadratic_form rf(*this), rg(g), rg2(g), root, rff;

		rf.reduce(S);
		rg.reduce(V);
		rff.assign(-rf.a, rf.b, -rf.c);

		if (rf == rg) {
			S.invert();
			multiply(V, V, S);
			multiply(U, U, V);
			return true;
		}
		if (rff == rg) {
			S.invert();
			multiply(V, V, W);
			multiply(V, V, S);
			multiply(U, U, V);
			return true;
		}

		rg2.reduce(T);

		root.assign(rg);
		do {
			rg.rho(V);
			if (rf == rg) {
				S.invert();
				multiply(V, V, S);
				multiply(U, U, V);
				return true;
			}
			if (rff == rg) {
				S.invert();
				multiply(V, V, W);
				multiply(V, V, S);
				multiply(U, U, V);
				return true;
			}

			rg2.inverse_rho(T);
			if (rf == rg2) {
				S.invert();
				multiply(T, T, S);
				multiply(U, U, T);
				return true;
			}
			if (rff == rg2) {
				S.invert();
				multiply(T, T, W);
				multiply(T, T, S);
				multiply(U, U, T);
				return true;
			}
		} while (root.compare(rg) != 0);

		return false;
	}


	if (dness == 2) {
		matrix_GL2Z Q(I), R(I);
		hf.reduce(Q);
		hg.reduce(R);
		if (!hf.compare(hg)) {
			R.invert();
			multiply(Q, Q, R);
			multiply(U, U, Q);
			return true;
		}
		else {
			hg.transform(W);
			multiply(R, R, W);
			hg.reduce(R);
			if (!hf.compare(hg)) {
				R.invert();
				multiply(Q, Q, R);
				multiply(U, U, Q);
				return true;
			}
			else
				return false;
		}
	}

	return false;
}



//
// quadratic_form::is_prop_equivalent()
//
// Task:
//      returns true if the forms are properly equivalent.
//

bool
quadratic_form::is_prop_equivalent (const quadratic_form & g) const
{
	debug_handler("quadratic_form", "is_prop_equivalent()");

	int dness;

	if (Delta != g.Delta)
		return false;

	dness = definiteness();
	if (dness == 1)
		return prop_equivalent_pos_def(g);
	else if (dness == 2)
		return prop_equivalent_irregular(g);
	else if (dness == 0)
		return prop_equivalent_indef(g);
	else
		return prop_equivalent_neg_def(g);
}



//
// quadratic_form::is_principal()
//
// Task:
//      returns true if the form is equivalent to the unit form.
//

bool
quadratic_form::is_principal () const
{
	debug_handler("quadratic_form", "is_principal()");

	quadratic_form UNIT;

	UNIT.assign_one(Delta);
	if (is_neg_definite()) {
		UNIT.a.negate();
		UNIT.c.negate();
	}
	return is_equivalent(UNIT);
}



//
// quadratic_form::is_equivalent()
//
// Task:
//      returns true if the forms are equivalent.
//

bool
quadratic_form::is_equivalent (const quadratic_form & g) const
{
	debug_handler("quadratic_form", "is_equivalent()");

	quadratic_form hf(*this), hg(g);
	matrix_GL2Z W((-1), (0),
		      (0), (1));
	matrix_GL2Z I((1), (0),
		      (0), (1));

	if (Delta != g.Delta)
		return false;

	if (is_pos_definite() && g.is_pos_definite())
		return prop_equivalent_pos_def(g);

	if (is_neg_definite() && g.is_neg_definite())
		return prop_equivalent_neg_def(g);

	if (is_neg_definite() && g.is_pos_definite()) {
		hf.assign(-a, b, -c);
		return hf.prop_equivalent_pos_def(g);
	}

	if (is_pos_definite() && g.is_neg_definite()) {
		hg.assign(-g.a, g.b, -g.c);
		return prop_equivalent_pos_def(hg);
	}

	if (is_indefinite() && g.is_indefinite()) {
		quadratic_form rf(*this), rg(g), root, rff;

		rf.reduce();
		rg.reduce();
		rff.assign(-rf.a, rf.b, -rf.c);

		if (!rf.compare(rg))
			return true;
		if (!rff.compare(rg))
			return true;

		root.assign(rg);
		do {
			rg.rho();
			if (!rf.compare(rg))
				return true;
			if (!rff.compare(rg))
				return true;
		} while (root.compare(rg));

		return false;
	}

	if (!is_regular() && !g.is_regular()) {
		hf.reduce();
		hg.reduce();
		if (!hf.compare(hg))
			return true;
		else {
			hg.assign(-hg.a, hg.b, -hg.c);
			hg.reduce();
			if (!hf.compare(hg))
				return true;
			return false;
		}
	}

	return false;
}



//
// quadratic_form::is_prop_equivalent()
//
// Task:
//      returns true if the forms are properly equivalent and computes the
//      corresponding transformation matrix.
//

bool
quadratic_form::is_prop_equivalent (const quadratic_form & g,
				    matrix_GL2Z & U) const
{
	debug_handler("quadratic_form", "is_prop_equivalent(matrix_GL2Z)");

	int dness;

	if (Delta != g.Delta)
		return false;

	dness = definiteness();
	if (dness == 1)
		return prop_equivalent_pos_def(g, U);
	else if (dness == 2)
		return prop_equivalent_irregular(g, U);
	else if (dness == 0)
		return prop_equivalent_indef(g, U);
	else
		return prop_equivalent_neg_def(g, U);
}



//
// quadratic_form::is_principal(matrix_GL2Z)
//
// Task:
//      returns true if the form is equivalent to the unit form, and computes
//      the corresponding transformation matrix.
//

bool
quadratic_form::is_principal (matrix_GL2Z & U) const
{
	debug_handler("quadratic_form", "is_principal(matrix_GL2Z)");

	quadratic_form UNIT;

	UNIT.assign_one(Delta);
	if (is_neg_definite()) {
		UNIT.a.negate();
		UNIT.c.negate();
	}
	return is_equivalent(UNIT, U);
}



//
// quadratic_form::order_in_CL()
//
// Task:
//      returns the order in the class group of the equivalence class
//      containing the form.  If the form is regular, the qi_class function of
//      the same name is used for the computation.
//

bigint
quadratic_form::order_in_CL () const
{
	debug_handler("quadratic_form", "order_in_CL");

	qi_class A;
	bigint ord;
	quadratic_form f, g;

	if (is_regular()) {
		ord.assign_zero();
		if (is_primitive()) {
			// search for quadratic_order in qo_list
			qi_class::set_current_order(*quadratic_order::qo_l().add_to_list(Delta));
			A.assign(*this);
			ord.assign(A.order_in_CL());
		}
	}
	else {
		f.assign(*this);
		f.reduce();
		g.assign(f);
		ord.assign_one();
		while (!f.is_one()) {
			compose_reduce(f, f, g);
			inc(ord);
		}
	}

	return ord;
}



//
// quadratic_form::DL()
//
// Task:
//      computes the smallest integer x such that G^x is equivalent to  the
//      form, if x exists.  If x exists, true is returned, otherwise false is
//      returned and x is set to the order of the equivalence class containing
//      G.  If the form is regular, the qi_class function of the same name is
//      used for the computation.
//

bool
quadratic_form::DL (quadratic_form & G, bigint & x) const
{
	debug_handler("quadratic_form", "DL");

	qi_class A, B;
	bool is_DL;

	is_DL = false;
	x.assign_zero();

	if (is_regular()) {
		if (Delta != G.Delta) {
			x.assign_zero();
			return false;
		}

		if ((is_primitive()) && (G.is_primitive())) {
			// search for quadratic_order in qo_list
			qi_class::set_current_order(*quadratic_order::qo_l().add_to_list(Delta));
			A.assign(*this);
			B.assign(G);
			is_DL = A.DL(B, x);
		}
	}

	return is_DL;
}



//
// subgroup()
//
// Task:
//      computes the structure of the subgroup generated by the forms in G.
//      The vector factor_base contains the forms in G which contributed to
//      the subgroup, and U is the transformation matrix used to compute the
//      SNF, i.e., it's inverse represents a map from factor_base to
//      generators of the subgroup.  If the forms are regular, the qi_class
//      function of the same name is used.
//

base_vector< bigint >
subgroup (base_vector< quadratic_form > & G)
{
	debug_handler("quadratic_form", "subgroup");

	lidia_size_t i, j, l;
	base_vector< qi_class > A;
	base_vector< bigint > S;
	qi_class Amem;
	quadratic_form f;

	A.set_mode(EXPAND);
	S.set_mode(EXPAND);
	A.reset();
	S.reset();

	if (G[0].is_regular()) {
		// search for quadratic_order in qo_list
		qi_class::set_current_order(*quadratic_order::qo_l().add_to_list(G[0].Delta));

		l = G.size();
		j = 0;
		for (i = 0; i < l; ++i) {
			if ((qi_class::discriminant() == G[i].Delta) &&
			    (G[i].is_primitive())) {
				A[j] = qi_class(G[i]);
				++j;
			}
		}

		if (A.size() > 0)
			S = subgroup(A);
	}

	return S;
}



//
// quadratic_form::regulator()
//
// Task:
//      returns the regulator of the quadratic order of the same discriminant
//      as the form.  If the form is not irregular, 1 is returned.
//

bigfloat
quadratic_form::regulator ()
{
	debug_handler("quadratic_form", "regulator");

	if (is_indefinite()) {
		return quadratic_order::qo_l().add_to_list(Delta)->regulator();
	}
	else
		return bigfloat(1.0);
}



//
// FIX - class number of irregulars
//

//
// quadratic_form::class_number()
//
// Task:
//      returns the number of equivalence classes of forms of the same
//      discriminant as the form.
//

bigint
quadratic_form::class_number ()
{
	debug_handler("quadratic_form", "class_number");

	if (is_regular()) {
		return quadratic_order::qo_l().add_to_list(Delta)->class_number();
	}
	else
		return bigint(0);
}



//
// FIX - class group of irregulars
//

//
// quadratic_form::class_group()
//
// Task:
//      computes the structure of the group of equivalence classes of forms of
//      the same discriminant as the form.
//

base_vector< bigint >
quadratic_form::class_group ()
{
	debug_handler("quadratic_form", "class_group");

	if (is_regular()) {
		return quadratic_order::qo_l().add_to_list(Delta)->class_group();
	}
	else
		return base_vector< bigint > ();
}



//
// quadratic_form::fundamental_automorphism()
//
// Task:
//      computes the transformation matrix representing the fundamental
//      automorphism.
//

void
quadratic_form::fundamental_automorphism (matrix_GL2Z & U)
{
	debug_handler("quadratic_form", "fundamental_automorphism");

	matrix_GL2Z V((-1), (0),
		      (0), (1));

	if (!is_indefinite())
		return;

	quadratic_form g(*this);
	g.reduce();

	quadratic_form h(g);
	quadratic_form i(g);
	i.transform(V);
	do {
		g.rho(U);
	}while ((g.compare(h)) && (g.compare(i)));
}



//
// quadratic_form::representations()
//
// Task:
//      computes all non-trivial representations of N by the form.  If none
//      exist, false is returned.  If infinitely many exist, the smallest
//      inequivalent ones are computed.
//

bool
quadratic_form::representations (sort_vector< pair < bigint, bigint > > & Reps,
				 const bigint & N)
{
	debug_handler("quadratic_form", "representations");

	quadratic_form g, f, h;
	rational_factorization Nfact;
	bigint con, mult, NN, temp, temp2, ta, tb, tc, tD, x, y, p;
	lidia_size_t num, i, j, k, cands, idx, nfacts, val, v, pos;
	matrix_GL2Z U(1, 0, 0, 1), V;
	base_vector< quadratic_form > plist;
	base_vector< int > elist;
	hash_table< quadratic_form > used_forms;
	quadratic_order *gQO;

	Reps.set_mode(EXPAND);
	Reps.reset();

	// if form is zero, no representations unless N = 0
	if (is_zero()) {
		if (N.is_zero()) {
			Reps[0].assign(0, 0);
			return true;
		}
		else
			return false;
	}


	// if form is irregular, use special function
	if (!is_regular())
		return comp_reps_irregular(Reps, N);



	// if N=0, only one solution
	if (N.is_zero()) {
		Reps[0].assign(0, 0);
		return true;
	}



#ifdef LIDIA_DEBUG
	std::cout << "\nIN REPRESENTATIONS:  *this = " << *this << ", N = ";
	std::cout << N << std::endl;
	std::cout << "Delta:  " << Delta << std::endl;
#endif


	// make form primitive --- content must divide N
	con = content();
	remainder(temp, N, con);
	if (!temp.is_zero()) {
		return false;
	}

	divide(NN, N, con);
	g.assign(a/con, b/con, c/con);

#ifdef LIDIA_DEBUG
	std::cout << "primitive:  g = " << g << ", NN = " << NN << std::endl;
#endif



	// if form is not in a maximal order, compute representations for related
	// form in maximal order first

	gQO = quadratic_order::qo_l().add_to_list(g.Delta);
	con = gQO->conductor();
	if (!gQO->is_maximal() && (gcd(con, NN) > 1)) {
		// compute equivalent form relatively prime to the conductor
		if (gcd(g.a, con) > 1) {
			if (gcd(g.c, con) > 1) {
				ta.assign(g.a + g.b + g.c);
				tb.assign(-g.b - (g.a << 1));
				V = matrix_GL2Z(1, -1, 1, 0);
			}
			else {
				ta.assign(g.c);
				tb.assign(-g.b);
				V = matrix_GL2Z(0, 1, -1, 0);
			}
		}
		else {
			ta.assign(g.a);
			tb.assign(g.b);
			V = U;
		}

		square(tc, tb);
		subtract(tc, tc, g.Delta);
		divide(tc, tc, ta);
		shift_right(tc, tc, 2);
		f.assign(ta, tb, tc);

#ifdef LIDIA_DEBUG
		std::cout << "NON-MAX:  con = " << con << std::endl;
		std::cout << "f = " << f << std::endl;
		std::cout << "V = " << V << std::endl;
#endif


		// compute form with same norm in the maximal order
		xgcd(temp, temp2, con, f.a);
		tD = (f.Delta/(con*con)) % 2;
		tb.assign(f.b*temp + f.a*tD*temp2);

		square(tc, tb);
		subtract(tc, tc, f.Delta/(con*con));
		divide(tc, tc, f.a);
		shift_right(tc, tc, 2);
		h.assign(f.a, tb, tc);

#ifdef LIDIA_DEBUG
		std::cout << "u = " << temp << ", v = " << temp2 << std::endl;
		std::cout << "u con + v a = " << temp*con + temp2*f.a << std::endl;
		std::cout << "h = " << h << std::endl;
#endif

		if (!h.representations(Reps, NN))
			return false;
		else {
			// compute representations for original form

#ifdef LIDIA_DEBUG
			std::cout << "\nComputing representations for orignial form:" << std::endl;
			std::cout << "con = " << con << std::endl;
#endif

			j = 0;
			for (i = 0; i < Reps.size(); ++i) {
				if ((Reps[i].right() % con) == 0) {
					divide(temp, Reps[i].right(), con);

					ta.assign(f.a);
					multiply(tb, f.b, temp);
					multiply(tc, f.c, temp*temp);
					subtract(tc, tc, NN);

					square(tD, tb);
					subtract(tD, tD, (ta*tc) << 2);
					sqrt(tD, abs(tD));

#ifdef LIDIA_DEBUG
					std::cout << "\nRep[" << i << "] = " << Reps[i] << std::endl;
					std::cout << "ty = " << temp << std::endl;
					std::cout << "f:  " << ta << " X^2 + " << tb << " X + " << tc << std::endl;
					std::cout << "D = " << tb*tb - 4*ta*tc << std::endl;
					std::cout << "sqrt(D) = " << tD << std::endl;
#endif

					ta.multiply_by_2();
					subtract(temp2, -tb, tD);
					if ((temp2 % ta) == 0) {
						divide(temp2, temp2, ta);
						Reps[j].assign(temp2, temp);

						temp = V.get_s()*Reps[j].left() + V.get_u()*Reps[j].right();
						temp2 = V.get_t()*Reps[j].left() + V.get_v()*Reps[j].right();
						Reps[j].assign(temp, temp2);
						++j;
					}
					else {
						add(temp2, -tb, tD);
						if ((temp2 % ta) == 0) {
							divide(temp2, temp2, ta);
							Reps[j].assign(temp2, temp);

							temp = V.get_s()*Reps[j].left() + V.get_u()*Reps[j].right();
							temp2 = V.get_t()*Reps[j].left() + V.get_v()*Reps[j].right();
							Reps[j].assign(temp, temp2);
							++j;
						}
#ifdef LIDIA_DEBUG
						else
							std::cout << "ERROR:  can't find X!!!" << std::endl;
#endif
					}
				}
			}

			Reps.set_size(j);

			return true;
		}
	}



	g.reduce(U);

	Nfact.assign(abs(NN));
	Nfact.factor();

	plist.set_mode(EXPAND);
	plist.reset();
	elist.set_mode(EXPAND);
	elist.reset();
	num = 0;

#ifdef LIDIA_DEBUG
	std::cout << "reduced:  g = " << g << std::endl;
	std::cout << "U = \n" << U << std::endl;
	std::cout << "Delta = " << g.Delta << std::endl;
	std::cout << "Nfact = " << Nfact << std::endl;
#endif



	// find all prime forms with norm dividing NN
	cands = 0;
	nfacts = Nfact.no_of_comp();
	f.assign_one(g.Delta);
	mult.assign_one();
	for (i = 0; i < nfacts; ++i) {
		k = Nfact.exponent(i);

		if (!generate_prime_form(h, Nfact.base(i), g.Delta)) {
			while (k > 1) {
				k -= 2;
				multiply(mult, mult, Nfact.base(i));
			}
			if (k)
				return false;
		}
		else {
			if (g.is_neg_definite()) {
				h.a.negate();
				h.c.negate();
			}

			if ((g.Delta % Nfact.base(i)) == 0) {
				while (k > 1) {
					k -= 2;
					multiply(mult, mult, Nfact.base(i));
				}
			}

			while (k > 0) {
				compose(f, f, h);
				plist[cands] = h;
				elist[cands] = 1;
				++cands;
				--k;
			}
		}
	}
	divide(NN, NN, mult*mult);



#ifdef LIDIA_DEBUG
	std::cout << "cands = " << cands << ", plist.size = " << plist.size() << "\n";
	std::cout << plist << std::endl;
	std::cout << elist << std::endl;
	std::cout << "NN = " << NN << std::endl;
	std::cout << "mult = " << mult << std::endl;
#endif


	used_forms.initialize(1 << (cands));
	used_forms.set_key_function(&quadratic_form_key);


	// find all forms with X^2 coefficient equal to NN
	idx = (1 << (cands));
	for (j = 1; j <= idx; ++j) {
		if (!used_forms.search(f)) {
			used_forms.hash(f);

			if (NN.sign() != f.a.sign()) {
				f.a.negate();
				f.c.negate();
			}

#ifdef LIDIA_DEBUG
			std::cout << "\nj = " << j << ", f = " << f << std::endl;
			std::cout << elist << std::endl;
#endif

			V = U;
			if (f.is_prop_equivalent(g, V)) {
				sqrt(temp, abs(NN/f.a));
				multiply(x, V.get_s(), mult*temp);
				multiply(y, V.get_t(), mult*temp);
				if (x.is_lt_zero() && y.is_lt_zero()) {
					x.negate();
					y.negate();
				}

#ifdef LIDIA_DEBUG
				std::cout << "MAYBE FOUND REP:  (" << x << ", " << y << ")" << std::endl;
				std::cout << "MAYBE:  x_1 = " << V.get_s() << ", y_1 = " << V.get_t() << std::endl;
				std::cout << "MAYBE:  temp = " << temp << std::endl;
				std::cout << "MAYBE:  mult = " << mult << std::endl;
				std::cout << "MAYBE:  NN = " << NN << ", f.a = " << f.a << std::endl;
				h = *this;
				h.transform(V);
				std::cout << "this(V) = " << h << std::endl;
				std::cout << "eval(x, y) = " << eval(x, y) << std::endl;
#endif

				if (eval(x, y) == N) {
#ifdef LIDIA_DEBUG
					std::cout << "FOUND REP:  (" << x << ", " << y << ")" << std::endl;
#endif
					if (!Reps.linear_search(pair< bigint, bigint > (x, y), pos)) {
						Reps[num].assign(x, y);
						++num;
					}
				}
			}

#ifdef LIDIA_DEBUG
			f.reduce();
			std::cout << "f reduced:  " << f << std::endl;
#endif
		}


		// compute next ideal (exp vectors arranged in Grey code)
		if (j < idx) {
			val = j;
			v = 0;
			while (!(val & 1)) {
				val >>= 1;
				++v;
			}

			elist[v] = -elist[v];

			f.assign_one(g.Delta);
			for (i = 0; i < cands; ++i) {
				if (elist[i] == 1)
					compose(f, f, plist[i]);
				else
					compose(f, f, get_conjugate(plist[i]));
			}
		}
	}

	// sort the representations
	Reps.sort();

	return (num > 0);
}



//
// operator >>
//
// Task:
//      inputs a quadratic_form from the std::istream in.
//

std::istream &
operator >> (std::istream & in, quadratic_form & f)
{
	debug_handler("quadratic_form", "operator >>");

	int n = 0;
	char c;
	bigint ibuf[3];

	in >> c;
	if (c != '(') {
		lidia_error_handler("quadratic_form", "operator >>::(expected");
		return in;
	}

	in >> c;
	while (c != ')' && n != 3) {
		in.putback(c);
		in >> ibuf[n];
		n++;
		in >> c;
		if (c == ',')
			in >> c;
	}
	f.assign(ibuf[0], ibuf[1], ibuf[2]);
	return in;
}



//
// operator <<
//
// Task:
//      outputs a quadratic_form to the std::ostream out.
//

std::ostream &
operator << (std::ostream & out, const quadratic_form & f)
{
	debug_handler("quadratic_form", "operator << ");

	out << "(" << f.a << ", " << f.b << ", " << f.c << ")";
	return out;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
