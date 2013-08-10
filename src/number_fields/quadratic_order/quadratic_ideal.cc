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
#include	"LiDIA/number_fields/qo_list.h"
#include	"LiDIA/quadratic_order.h"
#include	"LiDIA/quadratic_number_standard.h"
#include	"LiDIA/quadratic_number_power_product.h"
#include	"LiDIA/quadratic_ideal.h"
#include	"LiDIA/qi_class.h"
#include	"LiDIA/qi_class_real.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#define LIDIA_CAST_TO_QO(a) (*(const_cast<quadratic_order*>(&(a))))


bigint
quadratic_ideal_key(const quadratic_ideal & G)
{
	return G.get_a();
}


//
// Note: If reduced_flag == true, then q = 1/a. To save storage
//       q is set to one in this case. Therefore, each function
//       that depends on the value of q has to check whether the
//       reduced_flag is true.
//
// MM
//


//
//   Normalization
//

//
// quadratic_ideal::normalize()
//
// Task:
//      normalizes the quadratic_ideal, i.e.
//      normalizes the form (a,b,c).
//

void quadratic_ideal
::normalize()
{
	debug_handler("quadratic_ideal", "normalize");

	if (discriminant().is_lt_zero())
		normalize_imag();
	else
		normalize_real();
}


//
// quadratic_ideal::normalize_imag()
//
// Task:
//      normalizes the imaginary quadratic_ideal by
//      normalizing the form (a,b,c) where c = (b^2 - D)/4a,
//      q is unchanged.
//      (a,b,c) normal iff -a < b <= a.
//

void quadratic_ideal
::normalize_imag()
{
	debug_handler("quadratic_ideal", "normalize_imag");

	bigint a2, nq, r, temp;

	if (!((b.compare(-a) > 0) && (b.compare(a) <= 0))) {
		// q = floor((a-b) / 2a)
		shift_left(a2, a, 1);
		subtract(temp, a, b);
		div_rem(nq, r, temp, a2);
		if ((temp.is_lt_zero()) && (!r.is_zero()))
			dec(nq);

		multiply(temp, a2, nq);
		add(b, b, temp);
	}
}



//
// quadratic_ideal::normalize_real()
//
// Task:
//      normalizes the real quadratic_ideal by
//      normalizing the form (a,b,c) where c = (b^2 - D)/4a,
//      q is unchanged.
//      (a,b,c) normal iff  |a| >  sqrt(D): -|a| < b <= |a|
//                          |a| <= sqrt(D): sqrt(D) - 2|a| < b < sqrt(D)
//
//      Note that for the ideals, a is always positive, i.e., |a| = a,
//      because Z a = Z -a.
//

void quadratic_ideal
::normalize_real()
{
	debug_handler("quadratic_ideal", "normalize_real");

	bigint a2, nq, r, temp, rootD;
	bool is_normal;

	a.absolute_value(); // MM

	sqrt(rootD, discriminant());
	shift_left(a2, a, 1);
	if (a.compare(rootD) <= 0) {
		subtract(temp, rootD, a2);
		is_normal = ((temp.compare(b) < 0) && (b.compare(rootD) <= 0));
	}
	else
		is_normal = ((b.compare(-a) > 0) && (b.compare(a) <= 0));

	if (!is_normal) {
		// shift_left(a2,a,1); MM
		if (a <= rootD) {
			// q = floor((rootD - b) / 2a)
			subtract(temp, rootD, b);
		}
		else {
			// q = floor((a - b) / 2a)
			subtract(temp, a, b);
		}
		div_rem(nq, r, temp, a2);
		if ((temp.is_lt_zero()) && (!r.is_zero()))
			dec(nq);

		multiply(temp, a2, nq);
		add(b, b, temp);
	}
}


//
// quadratic_ideal::is_normalized()
//
// Task:
//      returns true, if the ideal is normalized;
//      false otherwise
//
// Note:
//      This function is for debbuging purposes only,
//      because the ideals should always be in
//      normalized representation.
//

bool quadratic_ideal
::is_normalized() const
{
	debug_handler("quadratic_ideal", "is_normalized()");

	if (discriminant().is_lt_zero())
		return this->is_normalized_imag();
	else
		return this->is_normalized_real();
}

bool quadratic_ideal
::is_normalized_imag() const
{
	debug_handler("quadratic_ideal", "is_normalized_imag()");

	return ((b.compare(-a) > 0) && (b.compare(a) <= 0));
}

bool quadratic_ideal
::is_normalized_real() const
{
	debug_handler("quadratic_ideal", "is_normalized_real()");

	bigint a1, a2, temp, rootD;

	a1.assign(abs(a));
	sqrt(rootD, discriminant());
	shift_left(a2, a1, 1);

	if (a1.compare(rootD) <= 0) {
		subtract(temp, rootD, a2);
		return ((temp.compare(b) < 0) && (b.compare(rootD) <= 0));
	}
	else
		return ((b.compare(-a1) > 0) && (b.compare(a1) <= 0));
}



//
//  Reduction
//

//
// quadratic_ideal::reduce()
//
// Task:
//      reduces the ideal. Let I = q (Z a + Z (b+sqrt(D))/2). Then
//      I is reduced by reducing the form (a,b,c) and setting q to 1/a'.
//      That means, the reduced ideals I are exactly those in which 1 is
//      a minimum.
//      Note that for storage reasons, instead of setting
//      q = 1/a, q is set to 1 as well as the reduced_flag is set to 1.
//

void quadratic_ideal
::reduce()
{
	debug_handler("quadratic_ideal", "reduce");
	quadratic_number_standard g;

	if (!this->is_zero()) {
	    if (discriminant().is_lt_zero())
		reduce_imag(g, false);
	    else
		reduce_real(g, false);
	}
}


//
// quadratic_ideal::reduce(quadratic_number_standard & g)
//
// Task:
//      reduces the ideal.
//
//      The reducing number g such that I/g is the reduced ideal is
//      returned.
//

void quadratic_ideal
::reduce(quadratic_number_standard & g)
{
	debug_handler("quadratic_ideal", "reduce(quadratic_number_standard &g)");

	if (discriminant().is_lt_zero())
		reduce_imag(g, true);
	else
		reduce_real(g, true);
}


//
// quadratic_ideal::reduce_imag()
//
// Task:
//      The function reduces the imaginary
//      quadratic ideal using the following algorithm:
//
//      Write the ideal in the form q * (Z + Z * (b+sqrt(Delta))/(2a)).
//      It is reduced, iff q = 1 and the form (a,b,c) is reduced.
//
//      To reduce the form (a,b,c) use the following algorithm:
//
//      -) normalize the form (a,b,c), c = (b^2 - Delta) / 4a
//      -) while (a > c)
//          { swap (a,c); b.negate(); normalize(a,b,c); }
//      -) if (a == c && b.is_negative()) b.negate();
//
//      I.e, the form (a,b,c) is reduced and q is unchanged.
//      Note that for computing the new c value during normalization
//      the formula
//         b = q 2a + r, -a < r <= a
//         c = c - q * ((b+r)/2)
//      is used. (for example, see Cohen, A course in computational
//      number theory, Algorithm 5.4.2).
//
//      If the flag "compute_reducing_number" is true, the reducing
//      number g, with this/g is reduced, will be computed.
//      Note that normalization does not change the ideal and swapping
//	a and c followed by a negation of b is the same as dividing the
//	ideal by (b+sqrt(D)) / 2c.
//
// Author: Mike Jacobson (MJJ), initial version
//         Markus Maurer (MM),  modified structure and added computation
//                              of reducing number
//

void quadratic_ideal
::reduce_imag(quadratic_number_standard & g, bool compute_reducing_number)
{
	debug_handler("quadratic_ideal", "reduce_imag(quadratic_number&, bool)");

	if (reduced_flag) {
		if (compute_reducing_number) {
			g.assign_order(this->which_order());
			g.assign_one();
		}
		return;
	}

	//
	// reduce the form (a,b,c)
	//

	bigint Delta, c, a2, qq, r, temp;
	quadratic_number_standard h;
	bigint one;

	Delta.assign(discriminant());

	if (compute_reducing_number) {
		g.assign_order(this->which_order());
		g.assign_one();
		h.assign_order(this->which_order());
		one.assign_one();
	}

	normalize_imag();

	// c = (b^2 - D) / 4a
	square(c, b);
	subtract(c, c, Delta);
	divide(c, c, a);
	shift_right(c, c, 2);

	// while (a > c)
	while (a.compare(c) > 0) {

		// a2 = 2*c (= 2*a after swap)
		shift_left(a2, c, 1);

		// multiply reducing number by
		// (b+sqrt(Delta)) / 2c
		if (compute_reducing_number) {
			h.assign(b, one, a2);
			multiply(g, g, h);
		}

		// exchange a and c, negate b
		swap (a, c);
		b.negate();

		// normalize, i.e.,
		// b = qq 2a + r, -a < r <= a,
		div_rem(qq, r, b, a2);

		if (r.compare(a) > 0)  // r > a
		{
			subtract(r, r, a2);
			inc(qq);
		}
		else {
			negate(temp, a);
			if (r.compare(temp) <= 0) // r <= -a
			{
				add(r, r, a2);
				dec(qq);
			}
		}

		// determine new c
		// c = c - qq*((b + r) / 2);
		add(temp, b, r);
		temp.divide_by_2();
		multiply(temp, temp, qq);
		subtract(c, c, temp);

		// assign b
		b.assign(r);
	}

	// for a == c, b > 0 required
	if ((a == c) && b.is_lt_zero()) {
		if (compute_reducing_number) {
			shift_left(a2, c, 1);
			h.assign(b, one, a2);
			multiply(g, g, h);
		}
		b.negate();
	}

	// Let I be the ideal before the above reduction has been applied.
	// We have this = I/g = q(I) * (Za + Z (b+sqrt(D))/2 with (a,b,c) reduced.
	// To complete the reduction we have to compute
	// this = this/(q*a), g *= (q*a)
	// to obtain a reduced this = 1/a (Za + Z (b+sqrt(D))/2.

	if (compute_reducing_number) {
		multiply (q, q, a);
		multiply (g, g, q);
	}

	q.assign_one();
	reduced_flag = true;
}




//
// quadratic_ideal::reduce_real()
//
// Task:
//      Reduces the real quadratic_ideal by computing
//      the continued fraction expansion of (b+sqrt(D))/2a
//      using Tenners algorithm (see Williams, Mollin,
//      Computation of the class number of a real quadratic field,
//      p. 13). It also computes a reducing number g such that
//      I / g is the reduced ideal, if the flag compute_reducing_number
//      is true.
//
// Note:
//      The ideal I = q * (Z + Z * (b+sqrt(Delta))/(2a)) is reduced, iff
//      q = 1 and the form (a,b,c) is reduced.
//
//      The indefinite form (a,b,c) is reduced, iff
//      |sqrt(Delta)-2|a|| < b < sqrt(Delta).
//
// Authors: Michael Jacobson, initial version
//          Markus Maurer, added computation of reducing number using
//                         Mikes implementation in qi_class_real.
//

void quadratic_ideal
::reduce_real(quadratic_number_standard & g, bool compute_reducing_number)
{
	debug_handler("quadratic_ideal", "reduce_real(quadratic_number, bool)");

	if (reduced_flag) {
		if (compute_reducing_number) {
			g.assign_order(this->which_order());
			g.assign_one();
		}
		return;
	}

	//
	// Reduce the form (a,b,c)
	//

	bigint Delta, rootD;
	bigint nq, r, a2, oa, na, nb, temp;
	int s;

	//<MM>
	bigint a0, OB, BB, NB;
	quadratic_number_standard h;
	//</MM>

	Delta.assign(discriminant());
	sqrt(rootD, Delta);

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

		//<MM>
		if (compute_reducing_number) {
			OB.assign_one();
			BB.assign_zero();
			shift_left(a0, a, 1);
		}
		//</MM>

		while (!((abs(temp) < b) && (b.compare(rootD) <= 0))) {
			s = a.is_gt_zero() ? 0 : 1;

			// (rootD+b+s) = (2a)*q + r
			shift_left(a2, a, 1);
			add(temp, rootD, b);
			if (s)
				inc(temp);
			div_rem(nq, r, temp, a2);

			//<MM>
			if (compute_reducing_number) {
				multiply(temp, nq, BB);
				add(NB, temp, OB);
				OB.assign(BB);
				BB.assign(NB);
			}
			//</MM>

			// nb = rootD + s - r;
			subtract(nb, rootD, r);
			if (s)
				inc(nb);

			// na = oa -q*((nb - b) >> 1);
			subtract(temp, nb, b);
			temp.divide_by_2();
			multiply(temp, temp, nq);
			subtract(na, oa, temp);

			b.assign(nb);
			oa.assign(a);
			a.assign(na);

			shift_left(a2, abs(a), 1);
			subtract(temp, rootD, a2);
			if (temp.is_lt_zero())
				inc(temp);
		}

		//<MM>
		if (compute_reducing_number) {
			// temp = 2*OB*a + BB*b
			multiply(a2, OB, a);
			a2.multiply_by_2();
			multiply(temp, BB, b);
			add(temp, temp, a2);

			// reducing number g with 1/g = (2*OB*a + BB*b + BB*sqrt(Delta)) /a0
			g.assign_order(this->which_order());
			g.set(temp, BB, a0);
			g.invert();
			g.absolute_value();
		}
		//</MM>

		a.absolute_value();
	}
	else // form (a, b, c) already reduced
	{
		g.assign_order(this->which_order());
		g.assign_one();
	}

	// Let I be the ideal before the above reduction has been applied.
	// We have this = I/g = q(I) * (Za + Z (b+sqrt(D)))/2 with (a,b,c)
	// reduced. To complete the reduction we have to compute
	// this = this/(q*a), g *= (q*a) to obtain a reduced
	// this = 1/a (Za + Z (b+sqrt(D)))/2.

	if (compute_reducing_number) {
		multiply (q, q, a);
		multiply (g, g, q);
	}
	q.assign_one();
	reduced_flag = true;
}




//
//  rho operator
//



//
// quadratic_ideal::rho()
//
// Task:
//      applies the reduction operator to the ideal.
//

void
quadratic_ideal::rho()
{
	debug_handler("quadratic_ideal", "rho()");
	quadratic_number_standard g;

	if (discriminant().is_lt_zero()) {
		if (!is_reduced())
			rho_imag(g);
	}
	else
		rho_real(g);
}


//
// quadratic_ideal::rho()
//
// Task:
//      applies the reduction operator rho to the ideal I.
//      rho(I) = I / g.
//

void
quadratic_ideal::rho(quadratic_number_standard & g)
{
	debug_handler("quadratic_ideal", "rho(quadratic_number&)");

	if (discriminant().is_lt_zero()) {
		if (!is_reduced())
			rho_imag(g);
		else {
			g.assign_order(this->which_order());
			g.assign_one();
		}
	}
	else
		rho_real(g);
}



//
// apply_rho(quadratic_ideal, quadratic_ideal)
//
// Task:
//      applies the reduction operator to the ideal.
//

void
apply_rho(quadratic_ideal & A, const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal",
		      "apply_rho(quadratic_ideal, quadratic_ideal)");

	A.assign(B);
	A.rho();
}



//
// apply_rho(quadratic_ideal)
//
// Task:
//      returns the result of applying the reduction operator to the ideal.
//

quadratic_ideal
apply_rho(const quadratic_ideal & A)
{
	debug_handler("quadratic_ideal",
		      "apply_rho(quadratic_ideal)");

	quadratic_ideal B;

	B.assign(A);
	B.rho();
	return B;
}


//
// quadratic_ideal::rho_imag()
//
// Task:
//     Apply the reduction operator once to the imaginary quadratic_ideal I.
//     The reduction operator is defined as
//
//      rho(I) = I/g, g = sigma(normalization((-b+sqrt(Delta))/(2a))).
//
//     and we also have
//
//      rho(I) = sigma(rho_inverse(normalization(sigma(I)))).
//
//     If I is reduced, than the rho operator is the inverse of the
//     inverse_rho operator.
//
// Condition:
//     I must be normalized.
//

void quadratic_ideal
::rho_imag(quadratic_number_standard & g)
{
	debug_handler("quadratic_ideal", "rho_imag(quadratic_number&)");

	bigint Delta, a2, temp;

	Delta.assign(discriminant());

	if (!reduced_flag)
		multiply(q, q, a);

	//
	// 1) normalization(sigma(I)), i.e., normalization(a,-b,..)
	//

	// -a < -b + 2sa <= a
	// Note: We assume that I is normalized, i.e., -a < b <= a.
	if (a != b)
		b.negate();


	// g = sigma(-b+sqrt(Delta))/(2a)
	// Note: Instead of computing the conjugate, we compute the negation
	// of it which does not change I/g.
	shift_left(a2, a, 1);
	g.assign_order(this->which_order());
	g.set(-b, bigint(1), a2);

	//
	// 2) sigma(rho_inverse()) = normalize(c,b,..)
	//

	// a = (b^2 - D) / 4a
	square(temp, b);
	subtract(temp, temp, Delta);
	divide(temp, temp, a);
	shift_right(a, temp, 2);

	normalize_imag();

	if (!reduced_flag)
		divide (q, q, a);
}


//
// quadratic_ideal::rho_real()
//
// Task:
//     Apply the reduction operator once to the imaginary quadratic_ideal I.
//     The reduction operator is defined as
//
//      rho(I) = I/g, g = sigma(normalization((-b+sqrt(Delta))/(2a))).
//
//     and we also have
//
//      rho(I) = sigma(rho_inverse(normalization(sigma(I)))).
//
//     If I is reduced, than the rho operator is the inverse of the
//     inverse_rho operator.
//
// Condition:
//     I must be normalized.
//


void quadratic_ideal
::rho_real(quadratic_number_standard & g)
{
	debug_handler("quadratic_ideal", "rho_real(quadratic_number&)");

	bigint Delta, rootD;
	bigint s, r, c, a2, temp;

	Delta.assign(discriminant());
	sqrt(rootD, Delta);

	if (!reduced_flag)
		multiply(q, q, a);

	//
	// 1) normalization(sigma(I)), i.e., normalization(a,-b,..)
	//

	// s = sign(a) floor((|a|+b)/(2|a|),  if |a| >= sqrt(Delta)
	// s = sign(a) floor(rootD+b)/(2|a|), if |a| <  sqrt(Delta)
	//
	// Then normalization(a,-b,..) = (a,-b+2sa,..)
	//
	// Furthermore, I is normal iff
	//     -|a| < b <= |a|,                        if |a| >= sqrt(Delta)
	//     |sqrt(Delta) - 2|a|| < b < sqrt(Delta), if |a| <  sqrt(Delta)
	//
	// Note: a > 0 and I is normal.

	shift_left(a2, a, 1);

	if (a > rootD) {
		if (a != b)
			b.negate();
	}
	else {
		add(temp, rootD, b); // temp > 0
		div_rem(s, r, temp, a2);
		multiply(temp, s, a2);
		subtract(b, temp, b);
	}

	// g = sigma(b+sqrt(Delta))/(2a)
	// Note: Instead of computing the conjugate, we compute the negation
	// of it which does not change I/g, but yields g > 0, if the ideal was
	// reduced.
	g.assign_order(this->which_order());
	g.set(-b, bigint(1), a2);

	//
	// 2) sigma(rho_inverse()) = normalize(c,b,..)
	//

	// c = (D - b^2) / 4a
	square(temp, b);
	subtract(c, Delta, temp);
	divide(c, c, a2);
	shift_right(a, c, 1);

	a.absolute_value();
	normalize_real();

	if (!reduced_flag)
		divide (q, q, a);
}




//
//  inverse rho operator
//

//
// quadratic_ideal::inverse_rho()
//
// Task:
//      applies the inverse reduction operator to the ideal.
//

void
quadratic_ideal::inverse_rho()
{
	debug_handler("quadratic_ideal", "inverse_rho()");
	quadratic_number_standard g;

	if (discriminant().is_lt_zero()) {
		if (!is_reduced())
			inverse_rho_imag(g);
	}
	else
		inverse_rho_real(g);
}

//
// quadratic_ideal::rho()
//
// Task:
//      applies the inverse reduction operator inverse_rho
//      to the ideal I.
//      inverse_rho(I) = I / g.
//

void
quadratic_ideal::inverse_rho(quadratic_number_standard & g)
{
	debug_handler("quadratic_ideal",
		      "inverse_rho(quadratic_number&)");

	if (discriminant().is_lt_zero()) {
		if (!is_reduced())
			inverse_rho_imag(g);
		else {
			g.assign_order(this->which_order());
			g.assign_one();
		}
	}
	else
		inverse_rho_real(g);
}


//
// apply_inverse_rho(quadratic_ideal, quadratic_ideal)
//
// Task:
//      applies the inverse reduction operator to the ideal.
//

void
apply_inverse_rho(quadratic_ideal & A, const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "apply_inverse_rho(quadratic_ideal, "
		      "quadratic_ideal)");

	A.assign(B);
	A.inverse_rho();
}



//
// apply_inverse_rho(quadratic_ideal)
//
// Task:
//      returns the result of applying the inverse reduction operator to the
//      ideal.
//

quadratic_ideal
apply_inverse_rho(const quadratic_ideal & A)
{
	debug_handler("quadratic_ideal", "apply_inverse_rho(quadratic_ideal)");

	quadratic_ideal B;

	B.assign(A);
	B.inverse_rho();
	return B;
}

//
// quadratic_ideal::inverse_rho_imag()
//
// Task:
//      Apply the inverse reduction operator once to the imaginary
//      quadratic_ideal I.
//
//      The inverse rho operator is defined as
//
//        inverse_rho(I) = I/g, g = (b+sqrt(D))/2a;
//
//      If I is reduced than |g| > 1 and g is right neighbor of 1 on I.
//      We also have
//
//       inverse_rho(I) = normalization(q * a/c * (Z c + Z (-b+sqrt(D))/2)).
//
// Condition:
//      I must be normalized.
//

void quadratic_ideal
::inverse_rho_imag(quadratic_number_standard & g)
{
	debug_handler("quadratic_ideal",
		      "inverse_rho_imag(quadratic_number&)");

	bigint Delta, nq, a2, r, temp;

	Delta.assign(discriminant());

	if (!reduced_flag)
		multiply (q, q, a);

	// a2 = 2a
	shift_left(a2, a, 1);

	// g = (b+sqrt(D)) / 2a
	g.assign_order(this->which_order());
	g.set(b, bigint(1), a2);

	//
	// Normalize (c,-b,..)
	//

	// a = (b^2 - D) / 4a
	square(temp, b);
	subtract(temp, temp, Delta);
	divide(temp, temp, a);
	shift_right(a, temp, 2);

	if (!reduced_flag)
		divide (q, q, a);

	// a2 = 2*a
	shift_left(a2, a, 1);

	// q = floor((a+b) / 2a)
	add(temp, a, b);
	div_rem(nq, r, temp, a2);
	if ((temp.is_lt_zero()) && (!r.is_zero()))
		dec(nq);

	// b = 2a(nq) - b
	multiply(temp, a2, nq);
	subtract(b, temp, b);
}



//
// quadratic_ideal::inverse_rho_real()
//
// Task:
//      Apply the inverse reduction operator once to the imaginary
//      quadratic_ideal I.
//
//      The inverse rho operator is defined as
//
//        inverse_rho(I) = I/g, g = (b+sqrt(D))/2a;
//
//      If I is reduced than |g| > 1 and g is right neighbor of 1 on I.
//      We also have
//
//       inverse_rho(I) = normalization(q * a/c * (Z c + Z (-b+sqrt(D))/2)).
//
// Condition:
//      I must be normalized.
//

void quadratic_ideal
::inverse_rho_real(quadratic_number_standard & g)
{
	debug_handler("quadratic_ideal",
		      "inverse_rho_real(quadratic_number&)");

	bigint Delta, rootD;
	bigint nq, a2, r, temp;

	Delta.assign(discriminant());
	sqrt(rootD, Delta);

	if (!reduced_flag)
		multiply (q, q, a);

	// a2 = 2*a
	shift_left(a2, a, 1);

	// g = (b+sqrt(D)) / 2a
	g.assign_order(this->which_order());
	g.set(b, bigint(1), a2);

	// a = (D - b^2) / 4a
	square(temp, b);
	subtract(temp, Delta, temp);
	divide(temp, temp, a);
	shift_right(a, temp, 2);
	a.absolute_value();

	if (!reduced_flag)
		divide (q, q, a);

	// a2 = 2*a
	shift_left(a2, a, 1);

	if (a <= rootD) {
		// q = floor((rootD + b) / 2a)
		add(temp, rootD, b);
	}
	else {
		// q = floor((a + b) / 2a)
		add(temp, a, b);
	}

	div_rem(nq, r, temp, a2);
	if ((temp.is_lt_zero()) && (!r.is_zero()))
		dec(nq);

	// b = 2a(nq) - b
	multiply(temp, a2, nq);
	subtract(b, temp, b);
}



//
//  Multiplication
//

//
// quadratic_ideal::multiply_imag()
//
// Task:
//      compute the product of two imaginary quadratic ideals.
//

void multiply_imag(quadratic_ideal & C,
		   const quadratic_ideal & A,
		   const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "multiply_imag(qi, qi, qi)");

	bigint Delta, newa, newb, dpr, v, d, w, ab2, temp;
	bigrational newq;

	Delta.assign(A.discriminant());

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
		subtract(temp, Delta, temp);
		temp.divide_by_2();
		multiply(temp, temp, w);

		add(newb, newb, temp);
		divide(newb, newb, d);

		// C.a = C.a / (d^2)
		square(temp, d);
		divide(newa, newa, temp);

		dpr.assign(d);
	}

	add(newb, newb, A.b);
	shift_left(ab2, newa, 1);
	remainder(newb, newb, ab2);

	multiply(newq, A.get_q(), B.get_q());
	multiply(newq, newq, bigrational(dpr));

	C.a = newa;
	C.b = newb;
	C.q = newq;
	C.QO = quadratic_order::qo_l().add_to_list(C.QO, A.QO);
	C.normalize_imag();
	C.reduced_flag = false;
}



//
// quadratic_ideal::multiply_real()
//
// Task:
//      compute the product of two real quadratic ideals.
//

void multiply_real(quadratic_ideal & C,
		   const quadratic_ideal & A,
		   const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "multiply_real(qi, qi, qi)");

	bigint Delta, newa, newb, dpr, v, d, w, ab2, temp;
	bigrational newq;

	Delta.assign(A.discriminant());

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
		subtract(temp, Delta, temp);
		temp.divide_by_2();
		multiply(temp, temp, w);

		add(newb, newb, temp);
		divide(newb, newb, d);

		// C.a = C.a / (d^2)
		square(temp, d);
		divide(newa, newa, temp);

		dpr.assign(d);
	}

	add(newb, newb, A.b);
	shift_left(ab2, newa, 1);
	remainder(newb, newb, ab2);

	multiply(newq, A.get_q(), B.get_q());
	multiply(newq, newq, bigrational(dpr));

	C.a = newa;
	C.b = newb;
	C.q = newq;
	C.QO = quadratic_order::qo_l().add_to_list(C.QO, A.QO);
	C.normalize_real();
	C.reduced_flag = false;
}



//
// quadratic_ideal::square_imag()
//
// Task:
//      compute the square of an imaginary quadratic ideal.
//

void square_imag(quadratic_ideal & C, const quadratic_ideal & A)
{
	debug_handler("quadratic_ideal", "square_imag(qi, qi)");

	bigint Delta, newa, newb, d, w, temp;
	bigrational newq;

	Delta.assign(A.discriminant());

	// solve d = v A.a + w A.b
	d.assign(xgcd_right(w, A.a, A.b));

	// C.b = A.b + w (D - A.b^2) / (2d)
	square(newb, A.b);
	subtract(newb, Delta, newb);
	shift_left(temp, d, 1);
	divide(newb, newb, temp);
	multiply(newb, newb, w);
	add(newb, newb, A.b);

	// C.a = (A.a/d)^2
	divide(temp, A.a, d);
	square(newa, temp);

	shift_left(temp, newa, 1);
	remainder(newb, newb, temp);

	square(newq, A.get_q());
	multiply(newq, newq, d);

	C.a = newa;
	C.b = newb;
	C.q = newq;
	C.QO = quadratic_order::qo_l().add_to_list(C.QO, A.QO);
	C.normalize_imag();
	C.reduced_flag = false;
}



//
// quadratic_ideal::square_real()
//
// Task:
//      compute the square of a real quadratic ideal.
//

void square_real(quadratic_ideal & C, const quadratic_ideal & A)
{
	debug_handler("quadratic_ideal", "square_real(qi, qi)");

	bigint Delta, newa, newb, d, w, temp;
	bigrational newq;

	Delta.assign(A.discriminant());

	// solve d = v A.a + w A.b
	d.assign(xgcd_right(w, A.a, A.b));

	// C.b = A.b + w (D - A.b^2) / (2d)
	square(newb, A.b);
	subtract(newb, Delta, newb);
	shift_left(temp, d, 1);
	divide(newb, newb, temp);
	multiply(newb, newb, w);
	add(newb, newb, A.b);

	// C.a = (A.a/d)^2
	divide(temp, A.a, d);
	square(newa, temp);

	shift_left(temp, newa, 1);
	remainder(newb, newb, temp);

	square(newq, A.get_q());
	multiply(newq, newq, d);

	C.a = newa;
	C.b = newb;
	C.q = newq;
	C.QO = quadratic_order::qo_l().add_to_list(C.QO, A.QO);
	C.normalize_real();
	C.reduced_flag = false;
}



//
//  constructors / destructor
//

//
// constructor:
//    - initialize quadratic order pointer to NULL.
//

quadratic_ideal::quadratic_ideal()
{
	debug_handler("quadratic_ideal", "quadratic_ideal()");

	QO = NULL;
	a.assign_zero();
	b.assign_zero();
	q.assign_zero();
	QO = quadratic_order::qo_l().add_to_list(QO, quadratic_order::zero_QO);
	reduced_flag = false;
}



//
// constructor
//
// Task:
//    set to the ideal belonging to the current order given by
//
//    q2 [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    If the parameters do not form an ideal of the current quadratic order,
//    or the current quadratic order is zero, the ideal will be set to the
//    zero ideal.
//

quadratic_ideal::quadratic_ideal(const bigint & a2, const bigint & b2,
                                 const bigrational & q2)
{
	debug_handler("quadratic_ideal", "quadratic_ideal(bigint, bigint, bigrational)");

	bigint c, temp;

	QO = NULL;
	a.assign_zero();
	b.assign_zero();
	q.assign(q2);
	QO = quadratic_order::qo_l().add_last_to_list(QO);

	if ((!a2.is_zero()) && (!discriminant().is_zero())) {
		square(c, b2);
		subtract(c, c, discriminant());
		shift_left(temp, a2, 2);
		remainder(c, c, temp);

		// (b^2 - Delta) must be divisible by 4a
		if (c.is_zero()) {
			a.assign(a2);
			b.assign(b2);
			normalize();
		}
	}
	reduced_flag = false;
}



//
// constructor
//
// Task:
//    set to the ideal belonging to the current order given by
//
//    q2 [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    If the parameters do not form an ideal of the current quadratic order,
//    or the current quadratic order is zero, the ideal will be set to the
//    zero ideal.
//

quadratic_ideal::quadratic_ideal(const long a2, const long b2,
                                 const bigrational & q2)
{
	debug_handler("quadratic_ideal", "quadratic_ideal(long, long, bigrational)");

	bigint c, temp;

	QO = NULL;
	a.assign_zero();
	b.assign_zero();
	q.assign(q2);
	QO = quadratic_order::qo_l().add_last_to_list(QO);

	if ((a2) && (!discriminant().is_zero())) {
		square(c, b2);
		subtract(c, c, discriminant());
		shift_left(temp, bigint(a2), 2);
		remainder(c, c, temp);

		// (b^2 - Delta) must be divisible by 4a
		if (c.is_zero()) {
			a.assign(a2);
			b.assign(b2);
			normalize();
		}
	}
	reduced_flag = false;
}



//
// constructor
//
// Task:
//    set to the ideal belonging to the current order given by
//
//    q2 [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    If the parameters do not form an ideal of the given quadratic order,
//    or the given quadratic order is zero, the ideal will be set to the
//    zero ideal.
//

quadratic_ideal::quadratic_ideal(const bigint & a2, const bigint & b2,
                                 const bigrational & q2, const quadratic_order & QO2)
{
	debug_handler("quadratic_ideal", "quadratic_ideal(bigint, bigint, bigrational, "
		      "quadratic_order)");

	bigint c, temp;

	QO = NULL;
	a.assign_zero();
	b.assign_zero();
	q.assign(q2);
	QO = quadratic_order::qo_l().add_to_list(QO, LIDIA_CAST_TO_QO(QO2));

	if ((!a2.is_zero()) && (!QO2.is_zero())) {
		square(c, b2);
		subtract(c, c, QO2.discriminant());
		shift_left(temp, a2, 2);
		remainder(c, c, temp);

		// (b^2 - Delta) must be divisible by 4a
		if (c.is_zero()) {
			a.assign(a2);
			b.assign(b2);
			normalize();
		}
	}
	reduced_flag = false;
}



//
// constructor
//
// Task:
//    set to the ideal belonging to the current order given by
//
//    q2 [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    If the parameters do not form an ideal of the given quadratic order,
//    or the given quadratic order is zero, the ideal will be set to the
//    zero ideal.
//

quadratic_ideal::quadratic_ideal(const long a2, const long b2,
                                 const bigrational & q2, const quadratic_order & QO2)
{
	debug_handler("quadratic_ideal", "quadratic_ideal(long, long, bigrational, "
		      "quadratic_order)");

	bigint c, temp;

	QO = NULL;
	a.assign_zero();
	b.assign_zero();
	q.assign(q2);
	QO = quadratic_order::qo_l().add_to_list(QO, LIDIA_CAST_TO_QO(QO2));

	if ((a2) && (!QO2.is_zero())) {
		square(c, b2);
		subtract(c, c, QO2.discriminant());
		shift_left(temp, bigint(a2), 2);
		remainder(c, c, temp);

		// (b^2 - Delta) must be divisible by 4a
		if (c.is_zero()) {
			a.assign(a2);
			b.assign(b2);
			normalize();
		}
	}
	reduced_flag = false;
}




//
// constructor:
//    if qf is regular, initialize with qf.
//

quadratic_ideal::quadratic_ideal(const quadratic_form & qf)
{
	debug_handler("quadratic_ideal", "quadratic_ideal(quadratic_form)");

	QO = NULL;
	if (!qf.is_regular()) {
		a.assign_zero();
		b.assign_zero();
		q.assign_zero();
		QO = quadratic_order::qo_l().add_to_list(QO, quadratic_order::zero_QO);
	}
	else {
		a.assign(abs(qf.get_a()));
		b.assign(qf.get_b());
		q.assign_one();
		QO = quadratic_order::qo_l().add_to_list(QO, qf.which_order());
		normalize();
	}
	reduced_flag = false;
}



//
// constructor:
//    initialize with the qi_class B.
//

quadratic_ideal::quadratic_ideal(const qi_class & B)
{
	debug_handler("quadratic_ideal", "quadratic_ideal(qi_class)");

	QO = NULL;
	a.assign(B.get_a());
	b.assign(B.get_b());
	q.assign_one();
	QO = quadratic_order::qo_l().add_to_list(QO, qi_class::get_current_order());
	reduced_flag = true;
}



//
// constructor:
//    initialize with the qi_class_real B.
//

quadratic_ideal::quadratic_ideal(const qi_class_real & B)
{
	debug_handler("quadratic_ideal", "quadratic_ideal(qi_class_real)");

	QO = NULL;
	a.assign(B.get_a());
	b.assign(B.get_b());
	q.assign_one();
	QO = quadratic_order::qo_l().add_to_list(QO, qi_class::get_current_order());
	reduced_flag = true;
}



//
// constructor:
//    initialize with a copy of the quadratic_ideal B.
//

quadratic_ideal::quadratic_ideal(const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "quadratic_ideal(quadratic_ideal)");

	QO = NULL;
	a.assign(B.a);
	b.assign(B.b);
	q.assign(B.q);
	QO = quadratic_order::qo_l().add_to_list(QO, B.QO);
	reduced_flag = B.reduced_flag;
}


//
// constructor:
//    initialize with the quadratic_order QO2.
//

quadratic_ideal
::quadratic_ideal(const quadratic_order & QO2)
{
	debug_handler("quadratic_ideal", "quadratic_ideal(quadratic_order)");
	QO = NULL;
	this->assign_one(QO2);
}


//
// destructor:
//    update reference counter of corresponding entry in the list of current
//    quadratic orders.
//

quadratic_ideal::~quadratic_ideal()
{
	debug_handler("quadratic_ideal", "~quadratic_ideal()");
	quadratic_order::qo_l().clear(QO);
}



//
//  Assignments
//

//
// quadratic_ideal::assign_zero()
//
// Task:
//    set to the zero ideal (all coefficients zero).  The ideal will belong
//    to the quadratic order which was most recently referenced.  If no
//    quadratic orders have been defined, the error_handler will be called.
//

void
quadratic_ideal::assign_zero()
{
	debug_handler("quadratic_ideal", "assign_zero()");

	assign_zero(*quadratic_order::qo_l().last());
	reduced_flag = false;
}



//
// quadratic_ideal::assign_one()
//
// Task:
//    set to the unit ideal (leading coefficient one).  The ideal will belong
//    to the quadratic order which was most recently referenced.  If no
//    quadratic orders have been defined, the error_handler will be called.
//

void
quadratic_ideal::assign_one()
{
	debug_handler("quadratic_ideal", "assign_one()");

	assign_one(*quadratic_order::qo_l().last());
	reduced_flag = true;
}



//
// quadratic_ideal::assign_principal()
//
// Task:
//    set to the principal ideal generated by the element
//
//    (x + (D + sqrt(Delta))/2 y)
//
//    The ideal will belong to the quadratic order which was most recently
//    referenced.  If no quadratic orders have been defined, the error_handler
//    will be called.
//

void
quadratic_ideal::assign_principal(const bigint & x, const bigint & y)
{
	debug_handler("quadratic_ideal", "assign_principal(bigint, bigint)");

	assign_principal(x, y, *quadratic_order::qo_l().last());
	reduced_flag = false;
}



//
// quadratic_ideal::assign(bigint, bigint, bigrational)
//
// Task:
//    set to the ideal given by
//
//    q2 [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    The ideal will belong to the quadratic order which was most recently
//    referenced.  If no quadratic orders have been defined, the error_handler
//    will be called.
//
//    If the parameters do not form an ideal, false will be returned and
//    the ideal will not be modified.
//

bool
quadratic_ideal::assign(const bigint & a2, const bigint & b2,
                        const bigrational & q2)
{
	debug_handler("quadratic_ideal", "assign(bigint, bigint, bigrational)");

	reduced_flag = false;
	return assign(a2, b2, q2, *quadratic_order::qo_l().last());
}



//
// quadratic_ideal::assign(long, long, bigrational)
//
// Task:
//    set to the ideal given by
//
//    q2 [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    The ideal will belong to the quadratic order which was most recently
//    referenced.  If no quadratic orders have been defined, the error_handler
//    will be called.
//
//    If the parameters do not form an ideal, false will be returned and
//    the ideal will not be modified.
//

bool
quadratic_ideal::assign(const long a2, const long b2, const bigrational & q2)
{
	debug_handler("quadratic_ideal", "assign(long, long, bigrational)");

	reduced_flag = false;
	return assign(a2, b2, q2, *quadratic_order::qo_l().last());
}



//
// quadratic_ideal::assign_zero(quadratic_order)
//
// Task:
//    set to the zero ideal (all coefficients zero) belonging to QO2.
//

void
quadratic_ideal::assign_zero(const quadratic_order & QO2)
{
	debug_handler("quadratic_ideal", "assign_zero(quadratic_order)");

	a.assign_zero();
	b.assign_zero();
	q.assign_zero();
	QO = quadratic_order::qo_l().add_to_list(QO, LIDIA_CAST_TO_QO(QO2));
	reduced_flag = false;
}



//
// quadratic_ideal::assign_one(quadratic_order)
//
// Task:
//    set to the unit ideal (leading coefficient one) belonging to QO2.
//

void
quadratic_ideal::assign_one(const quadratic_order & QO2)
{
	debug_handler("quadratic_ideal", "assign_one(quadratic_order)");

	QO = quadratic_order::qo_l().add_to_list(QO, LIDIA_CAST_TO_QO(QO2));
	if (QO2.is_zero()) {
		a.assign_zero();
		b.assign_zero();
		q.assign_zero();
	}
	else {
		a.assign_one();
		b.assign(QO2.discriminant());
		q.assign_one();
		normalize();
	}
	reduced_flag = true;
}



//
// quadratic_ideal::assign(quadratic_order)
//
// Task:
//    set to the unit ideal (leading coefficient one) belonging to QO2.
//

void
quadratic_ideal::assign(const quadratic_order & QO2)
{
	assign_one(QO2);
}



//
// quadratic_ideal::assign_principal(bigint, bigint, quadratic_order)
//
// Task:
//    set to the principal ideal belonging to QO2 generated by the element
//
//    (x + (Delta + sqrt(Delta))/2 y)
//

void
quadratic_ideal::assign_principal(const bigint & x, const bigint & y,
                                  const quadratic_order & QO2)
{
	debug_handler("quadratic_ideal", "assign_principal(bigint, bigint, "
		      "quadratic_order)");

	bigint x2, y2, Delta, n, m, k, l, temp, temp2;

	QO = quadratic_order::qo_l().add_to_list(QO, LIDIA_CAST_TO_QO(QO2));
	Delta = QO2.discriminant();
	if (Delta.is_zero()) {
		a.assign_zero();
		b.assign_zero();
		q.assign_zero();
	}
	else {
		// use (x2 + y2 \sqrt(Delta))/2 form
		y2.assign(y);
		shift_left(x2, x, 1);
		multiply(temp2, y2, Delta);
		add(x2, x2, temp2);

		// compute norm of (x2 + y2 sqrt(Delta))/2
		square(n, x2);
		square(temp, y2);
		multiply(temp, temp, Delta);
		subtract(n, n, temp);
		shift_right(n, n, 2);
		n.absolute_value();

		// solve m = k y2 + l (x2 + y2 Delta) / 2
		add(temp2, temp2, x2);
		temp2.divide_by_2();
		m.assign(xgcd(k, l, y2, temp2));

		// a = n / m^2
		square(temp, m);
		divide(a, n, temp);

		// b = (kx2 + l(x2+y2)D/2 ) / m
		add(temp2, x2, y2);
		multiply(temp2, temp2, l);
		multiply(b, temp2, Delta);
		b.divide_by_2();
		multiply(temp, k, x2);
		add(b, b, temp);
		divide(b, b, m);
		shift_left(temp, a, 1);
		remainder(b, b, temp);

		// q = m
		q.assign(m);

		normalize();
	}
	reduced_flag = false;
}




//
// quadratic_ideal::assign(bigint, bigint, bigrational, quadratic_order)
//
// Task:
//    set to the ideal belonging to QO2 given by
//
//    q2 [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//    If the parameters do not form an ideal in QO2, false will be returned and
//    the ideal will not be modified.
//

bool
quadratic_ideal::assign(const bigint & a2, const bigint & b2,
                        const bigrational & q2, const quadratic_order & QO2)
{
	debug_handler("quadratic_ideal", "assign(bigint, bigint, bigrational, "
		      "quadratic_order)");

	bigint c, temp;

	reduced_flag = false;
	a.assign_zero();
	b.assign_zero();
	q.assign(q2);
	QO = quadratic_order::qo_l().add_to_list(QO, LIDIA_CAST_TO_QO(QO2));

	if (a2.is_zero())
		if (b2.is_zero())
			return true;
		else
			return false;
	else {
		square(c, b2);
		subtract(c, c, QO2.discriminant());
		shift_left(temp, a2, 2);
		remainder(c, c, temp);

		// (b^2 - Delta) must be divisible by 4a
		if (c.is_zero()) {
			a.assign(a2);
			b.assign(b2);
			normalize();
			return true;
		}
		else
			return false;
	}
}



//
// quadratic_ideal::assign(long, long, bigrational, quadratic_order)
//
// Task:
//  set to the ideal belonging to QO2 given by
//
//  q2 [ a2 Z + (b2 + sqrt(Delta))/2 Z ]
//
//  If the parameters do not form an ideal in QO2, false will be returned
//  and the ideal will not be modified.
//

bool
quadratic_ideal::assign(const long a2, const long b2,
			const bigrational & q2, const quadratic_order & QO2)
{
	debug_handler("quadratic_ideal", "assign(long, long, bigrational, "
		      "quadratic_order)");
	reduced_flag = false;
	return assign(bigint(a2), bigint(b2), q2, QO2);
}



//
// quadratic_ideal::assign(quadratic_form)
//
// Task:
//    set to the normalization of qf.  If qf is not regular, the error_handler
//    will be called.
//

void
quadratic_ideal::assign(const quadratic_form & qf)
{
	debug_handler("quadratic_ideal", "assign(quadratic_form)");

	if (!qf.is_regular()) {
		lidia_error_handler("quadratic_ideal", "assign(quadratic_form) - "
				    "the quadratic form is not regular");
		return;
	}

	reduced_flag = false;
	a.assign(abs(qf.get_a()));
	b.assign(qf.get_b());
	q.assign_one();
	QO = quadratic_order::qo_l().add_to_list(QO, qf.which_order());
	normalize();
}



//
// quadratic_ideal::assign(qi_class)
//
// Task:
//    sets to the reduced ideal B.
//

void
quadratic_ideal::assign(const qi_class & B)
{
	debug_handler("quadratic_ideal", "assign(qi_class)");

	a.assign(B.get_a());
	b.assign(B.get_b());
	q.assign_one();
	QO = quadratic_order::qo_l().add_to_list(QO, qi_class::get_current_order());
	reduced_flag = true;
}



//
// quadratic_ideal::assign(qi_class_real)
//
// Task:
//    sets to the reduced ideal B.
//

void
quadratic_ideal::assign(const qi_class_real & B)
{
	debug_handler("quadratic_ideal", "assign(qi_class_real)");

	a.assign(B.get_a());
	b.assign(B.get_b());
	q.assign_one();
	QO = quadratic_order::qo_l().add_to_list(QO, qi_class::get_current_order());
	reduced_flag = true;
}



//
// quadratic_ideal::assign(quadratic_ideal)
//
// Task:
//    sets to a copy of B.
//

void
quadratic_ideal::assign(const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "assign(quadratic_ideal)");

	a.assign(B.a);
	b.assign(B.b);
	q.assign(B.q);
	reduced_flag = B.reduced_flag;
	QO = quadratic_order::qo_l().add_to_list(QO, B.QO);
}



//
// quadratic_ideal::assign_order()
//
// Task:
//  change the quadratic order to which the ideal belongs. If the ideal is
//  not an ideal in the new order, nothing is changed and false is returned.
//  Otherwise, true is returned.
//


bool
quadratic_ideal::assign_order(const quadratic_order & QO2)
{
	debug_handler("quadratic_ideal", "assign_order");

	bigint c, temp;

	if (is_zero()) {
		QO = quadratic_order::qo_l().add_to_list(QO, LIDIA_CAST_TO_QO(QO2));
		return true;
	}
	else {
		square(c, b);
		subtract(c, c, QO2.discriminant());
		shift_left(temp, a, 2);
		remainder(c, c, temp);

		// (b^2 - Delta) must be divisible by 4a
		if (c.is_zero()) {
			QO = quadratic_order::qo_l().add_to_list(QO, LIDIA_CAST_TO_QO(QO2));

			if (reduced_flag)
				divide (q, q, a);
			reduced_flag = false;

			return true;
		}
		else
			return false;
	}
}

//
// quadratic_ideal::set_order()
//
// Task:
//      same as assign_order
//

bool
quadratic_ideal::set_order(const quadratic_order & QO2)
{
	debug_handler("quadratic_ideal", "set_order");
	return assign_order(QO2);
}



//
// operator =
//
// Task:
//    make a copy of an existing quadratic_ideal
//


quadratic_ideal & quadratic_ideal::operator = (const quadratic_ideal & A)
{
	debug_handler("quadratic_ideal", "operator = ");

	assign(A);
	return *this;
}



//
// quadratic_ideal::get_a()
//
// Task:
//      returns coefficient a
//

bigint
quadratic_ideal::get_a() const
{
	return a;
}



//
// quadratic_ideal::get_b()
//
// Task:
//      returns coefficient b
//

bigint
quadratic_ideal::get_b() const
{
	return b;
}



//
// quadratic_ideal::get_q()
//
// Task:
//      returns coefficient q
//

bigrational
quadratic_ideal::get_q() const
{
	if (reduced_flag) {
		bigrational qq;
		qq.assign(bigint(1), a);
		return qq;
	}
	else
		return q;
}



//
// quadratic_ideal::get_c()
//
// Task:
//    return (b^2 - Delta)/(4a)
//


bigint
quadratic_ideal::get_c() const
{
	debug_handler("quadratic_ideal", "get_c");

	bigint Delta, c, temp;

	Delta.assign(discriminant());

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
// quadratic_ideal::which_order()
//
// Task:
//      return a pointer to the quadratic_order to which the ideal belongs.
//

const quadratic_order &
quadratic_ideal::which_order() const
{
	debug_handler("quadratic_ideal", "which_order()");

	if (!QO->get_qo()) {
		lidia_error_handler("quadratic_ideal", "which_order() - the quadratic "
				    "order of this ideal has been deleted");
		return quadratic_order::zero_QO;
	}

	return *QO->get_qo();
}


//
// quadratic_ideal::get_order()
//
// Task:
//      same as which_order
//

const quadratic_order &
quadratic_ideal::get_order() const
{
	debug_handler("quadratic_ideal", "get_order()");
	return which_order();
}



//
// quadratic_ideal::_order(quadratic_ideal)
//
// Task:
//      return a pointer to the quadratic_order to which the ideal belongs.
//

const quadratic_order &
which_order(const quadratic_ideal & A)
{
	debug_handler("quadratic_ideal", "which_order(quadratic_ideal)");

	if (!A.QO->get_qo()) {
		lidia_error_handler("quadratic_ideal", "which_order() - the quadratic "
				    "order of this ideal has been deleted");
		return quadratic_order::zero_QO;
	}

	return *A.QO->get_qo();
}



//
// quadratic_ideal::discriminant()
//
// Task:
//      return the discriminant of the quadratic order to which the ideal
//      belongs.  If no quadratic order has been defined for the ideal, or if
//      its quadratic order has been deleted, the error_handler is called.
//

bigint
quadratic_ideal::discriminant() const
{
	debug_handler("quadratic_ideal", "discriminant");

	if (!QO->get_qo()) {
		lidia_error_handler("quadratic_ideal", "discriminant() - the quadratic "
				    "order of this ideal has been deleted");
		return bigint();
	}

	return QO->get_qo()->discriminant();
}


//
// quadratic_ideal::add()
//
// Task:
//      computes the sum of two quadratic_ideals.  If the ideals belong to
//      different quadratic orders, the error_handler is called.
//
// Author: Markus Maurer
//
// Note: Not efficient, should be improved. (FIX ME)
//

void add(quadratic_ideal & C, const quadratic_ideal & A,
	 const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "add");

	if (A.discriminant() != B.discriminant()) {
		lidia_error_handler("quadratic_ideal::add",
				    "Ideals are associated with "
				    "different quadratic orders.");
		return;
	}

	bigint d1, d2, m1, m2;

	// q1 = m1 / d1
	//
	if (A.reduced_flag) {
		m1.assign_one();
		d1.assign(A.a);
	}
	else {
		m1.assign(A.q.numerator());
		d1.assign(A.q.denominator());
	}

	// q2 = m2 / d2
	//
	if (B.reduced_flag) {
		m2.assign_one();
		d2.assign(B.a);
	}
	else {
		m2.assign(B.q.numerator());
		d2.assign(B.q.denominator());
	}

	// A + B = 1/(d1d2) *
	//       = < m1*d2*a1, m2*d1*a2, m1*d2 (b1 + sqD)/2, m2*d1 (b2+sqD)/2 >
	//

	// <g> = < m1*d2*a1, m2*d1*a2 >
	//
	bigint m1d2, m2d1, g;

	multiply(m1d2, m1, d2);
	multiply(m2d1, m2, d1);

	g = gcd(m1d2*A.a, m2d1*B.a);

	// I = q * <a, (b+sqD)/2> = < m1*d2 (b1 + sqD)/2, m2*d1 (b2+sqD)/2 >
	//
	quadratic_number_standard alpha, beta;
	quadratic_ideal I;

	alpha.assign_order(A.get_order());
	beta.assign_order (A.get_order());

	alpha.assign(m1d2*A.b, m1d2, bigint(2));
	beta.assign (m2d1*B.b, m2d1, bigint(2));

#ifdef QI_ADD_DEBUG
	std::cout << "First assignment:" << std::endl;
	std::cout << "alpha = " << alpha << std::endl;
	std::cout << "beta = " << beta << std::endl;
#endif

	I.assign_module_of(bigrational(1), alpha, beta);
#ifdef QI_ADD_DEBUG
	std::cout << "I = " << I << std::endl;
#endif

	// Now: A + B = 1/(d1d2) * ( Z g + I )
	//            = 1/(d1d2) * < Z g + Z m a + Z m(b+sqD)/2 >
	//
	bigint m = I.get_q().numerator();

	g = gcd(g, m * I.a);

#ifdef QI_ADD_DEBUG
	std::cout << "Second assignment:" << std::endl;
	std::cout << "alpha = " << alpha << std::endl;
	std::cout << "beta = " << beta << std::endl;
#endif

	alpha.assign(g, bigint(0), bigint(1));
	beta.assign (m*I.b, m, bigint(2));

	C.assign_module_of(bigrational(1), alpha, beta);
	divide (C.q, C.q, d1*d2);
	C.QO = quadratic_order::qo_l().add_to_list(C.QO, A.QO);
	C.reduced_flag = false;
	C.normalize();
}



//
// quadratic_ideal::intersect()
//
// Task:
//      computes the intersection of two quadratic_ideals.  If the ideals
//      belong to different quadratic orders, the error_handler is called.
//

void
intersect(quadratic_ideal & C, const quadratic_ideal & A,
          const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "intersect");

	if (A.discriminant() != B.discriminant()) {
		lidia_error_handler("quadratic_ideal", "intersect() - ideals are "
				    "associated with different quadratic orders");
		return;
	}

	multiply(C, A, B);
	C /= (A+B);
}




//
// quadratic_ideal::multiply()
//
// Task:
//    computes the product of two quadratic_ideals.  If the ideals belong to
//    different quadratic orders, the error_handler is called.
//

void
multiply(quadratic_ideal & C, const quadratic_ideal & A,
         const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "multiply");

	if (A.discriminant() != B.discriminant()) {
		lidia_error_handler("quadratic_ideal", "multiply() - ideals are "
				    "associated with different quadratic orders");
		return;
	}

	if (A.discriminant() < 0)
		multiply_imag(C, A, B);
	else
		multiply_real(C, A, B);
}




//
// quadratic_ideal::divide()
//
// Task:
//   computes the quotient of two quadratic_ideals.  If the ideals belong to
//   different quadratic orders, the error_handler is called.  The quotient
//   is defined as A*B^-1.  If B is not invertible, its conjugate is used.
//

void
divide(quadratic_ideal & C, const quadratic_ideal & A,
       const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "divide");

	quadratic_ideal temp;

	if (B.is_zero()) {
		lidia_error_handler("quadratic_ideal", "divide() - zero divisor");
		return;
	}

	temp = inverse(B);
	multiply(C, A, temp);
}



//
// quadratic_ideal::conjugate()
//
// Task:
//      computes the conjugate of the ideal.
//

void
quadratic_ideal::conjugate()
{
	debug_handler("quadratic_ideal", "conjugate");

	if (!((discriminant().is_lt_zero()) && (a == get_c()))) {
		b.negate();
		normalize();
	}
}



//
// get_conjugate(quadratic_ideal, quadratic_ideal)
//
// Task:
//      computes the conjugate of B.
//

void
get_conjugate(quadratic_ideal & A, const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "get_conjugate(quadratic_ideal, "
		      "quadratic_ideal)");

	A.assign(B);
	if (!((A.discriminant().is_lt_zero()) && (A.a == A.get_c()))) {
		A.b.negate();
		A.normalize();
	}
}



//
// get_conjugate(quadratic_ideal)
//
// Task:
//      returns the conjugate of B.
//

quadratic_ideal
get_conjugate(const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "get_conjugate(quadratic_ideal)");

	quadratic_ideal A;

	A.assign(B);
	if (!((A.discriminant().is_lt_zero()) && (A.a == A.get_c()))) {
		A.b.negate();
		A.normalize();
	}

	return A;
}



//
// quadratic_ideal::invert()
//
// Task:
//      computes the inverse of B if B is invertible, otherwise computes its
//      conjugate.
//

void
quadratic_ideal::invert()
{
	debug_handler("quadratic_ideal", "invert");

	bigrational n;

	if (!is_invertible()) {
		if (!((discriminant().is_lt_zero()) && (a == get_c()))) {
			b.negate();
			normalize();
		}
	}
	else {

		if (reduced_flag)
			divide (q, q, a);
		reduced_flag = false;

		n.assign(norm());
		b.negate();
		normalize();
		divide(q, q, n);
	}
}



//
// inverse(quadratic_ideal, quadratic_ideal)
//
// Task:
//      computes the inverse of B if B is invertible, otherwise computes its
//      conjugate.
//

void
inverse(quadratic_ideal & A, const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "inverse(quadratic_ideal, quadratic_ideal)");

	bigrational n;

	if (B.is_zero())
		A.assign_zero(B.which_order());
	else {
		if (!B.is_invertible()) {
			A.assign(B);
			if (!((A.discriminant().is_lt_zero()) && (A.a == A.get_c()))) {
				A.b.negate();
				A.normalize();
			}
		}
		else {
			n.assign(B.norm());
			A.assign(B);

			if (A.reduced_flag)
				divide (A.q, A.q, A.a);
			A.reduced_flag = false;

			A.b.negate();
			A.normalize();
			divide(A.q, A.q, n);
		}
	}
}



//
// inverse(quadratic_ideal)
//
// Task:
//      returnss the inverse of B if B is invertible, otherwise returns its
//      conjugate.
//

quadratic_ideal
inverse(const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "inverse(quadratic_ideal)");

	quadratic_ideal A;
	bigrational n;

	if (B.is_zero())
		A.assign_zero(B.which_order());
	else {
		if (!B.is_invertible()) {
			A.assign(B);
			if (!((A.discriminant().is_lt_zero()) && (A.a == A.get_c()))) {
				A.b.negate();
				A.normalize();
			}
		}
		else {

			n.assign(B.norm());
			A.assign(B);

			if (A.reduced_flag)
				divide (A.q, A.q, A.a);
			A.reduced_flag = false;

			A.b.negate();
			A.normalize();
			divide(A.q, A.q, n);
		}
	}

	return A;
}



//
// quadratic_ideal::square()
//
// Task:
//      computes the square of A.
//

void
square(quadratic_ideal & C, const quadratic_ideal & A)
{
	debug_handler("quadratic_ideal", "square");

	if (A.discriminant().is_lt_zero())
		square_imag(C, A);
	else
		square_real(C, A);
}



//
// quadratic_ideal::power()
//
// Task:
//      computes A^i using binary exponentiation.
//

void power(quadratic_ideal & C, const quadratic_ideal & A, const bigint & i)
{
	debug_handler("quadratic_ideal", "power(quadratic_ideal, quadratic_ideal, "
		      "bigint)");
	quadratic_ideal B;
	bigint j;

	B.assign(A);
	j.assign(i);
	if (j.is_lt_zero()) {
		B.invert();
		j.absolute_value();
	}
	C.assign_one(A.which_order());
	while (j.is_gt_zero()) {
		if (j.is_odd())
			multiply(C, C, B);
		j.divide_by_2();
		if (j.is_gt_zero())
			square(B, B);
	}
}



//
// quadratic_ideal::power()
//
// Task:
//      computes A^i using binary exponentiation.
//

void
power(quadratic_ideal & C, const quadratic_ideal & A, const long i)
{
	debug_handler("quadratic_ideal", "power(quadratic_ideal, "
		      "quadratic_ideal, long)");

	quadratic_ideal B;
	register long j;

	B.assign(A);
	j = i;
	if (j < 0) {
		B.invert();
		j = -j;
	}
	C.assign_one(A.which_order());
	while (j > 0) {
		if ((j & 1) == 1)
			multiply(C, C, B);
		j >>= 1;
		if (j > 0)
			square(B, B);
	}
}



//
// operator -
//
// Task:
//      returns the inverse of A.
//

quadratic_ideal operator - (const quadratic_ideal & A)
{
	debug_handler("quadratic_ideal", "operator -");

	return inverse(A);
}



//
// operator +
//
// Task:
//      adds A and B.
//

quadratic_ideal operator + (const quadratic_ideal & A,
			    const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "operator +");

	quadratic_ideal C;
	add(C, A, B);
	return C;
}



//
// operator &
//
// Task:
//      intersection of A and B.
//

quadratic_ideal operator & (const quadratic_ideal & A, const quadratic_ideal &B)
{
	debug_handler("quadratic_ideal", "operator &");

	quadratic_ideal C;

	intersect(C, A, B);
	return C;
}



//
// operator *
//
// Task:
//      multiplies A and B.
//

quadratic_ideal operator * (const quadratic_ideal & A, const quadratic_ideal &B)
{
	debug_handler("quadratic_ideal", "operator *");

	quadratic_ideal C;

	multiply(C, A, B);
	return C;
}



//
// operator /
//
// Task:
//      divides A by B.
//

quadratic_ideal operator / (const quadratic_ideal & A, const quadratic_ideal &B)
{
	debug_handler("quadratic_ideal", "operator /");

	quadratic_ideal C;

	divide(C, A, B);
	return C;
}



//
// operator +=
//
// Task:
//      *this = *this + A.
//

quadratic_ideal & quadratic_ideal::operator += (const quadratic_ideal & A)
{
	debug_handler("quadratic_ideal", "operator += ");

	add(*this, *this, A);
	return *this;
}



//
// operator &=
//
// Task:
//      *this = *this & A
//

quadratic_ideal & quadratic_ideal::operator &= (const quadratic_ideal & A)
{
	debug_handler("quadratic_ideal", "operator &= ");

	intersect(*this, *this, A);
	return *this;
}



//
// operator *=
//
// Task:
//      *this = *this * A
//

quadratic_ideal & quadratic_ideal::operator *= (const quadratic_ideal & A)
{
	debug_handler("quadratic_ideal", "operator *= ");

	multiply(*this, *this, A);
	return *this;
}



//
// operator /=
//
// Task:
//      *this = *this / A.
//

quadratic_ideal & quadratic_ideal::operator /= (const quadratic_ideal & A)
{
	debug_handler("quadratic_ideal", "operator /= ");

	multiply(*this, *this, inverse(A));
	return *this;
}




//
// quadratic_ideal::is_zero()
//
// Task:
//      tests if the ideal is the zero ideal.
//

bool
quadratic_ideal::is_zero() const
{
	debug_handler("quadratic_ideal", "is_zero");

	return ((a.is_zero()) && (q.is_zero()));
}



//
// quadratic_ideal::is_one()
//
// Task:
//      tests if the ideal is the unit ideal.
//

bool
quadratic_ideal::is_one() const
{
	debug_handler("quadratic_ideal", "is_one");

	return ((a.is_one()) && (q.is_one()));
}



//
// quadratic_ideal::is_equal()
//
// Task:
//      tests if the ideal is equal to B.
//

bool quadratic_ideal
::is_equal(const quadratic_ideal & B) const
{
	debug_handler("quadratic_ideal", "is_equal(quadratic_ideal)");

	if (discriminant() != B.discriminant())
		return false;
	else if (a.compare(B.a) || b.compare(B.b))
		return false;
	else if (reduced_flag == B.reduced_flag) {
		if (reduced_flag)
			return true;
		else
			return (!q.compare(B.q));
	}
	else if (reduced_flag) {
		bigrational h;
		multiply (h, a, B.q);
		return h.is_one();
	}
	else {
		bigrational h;
		multiply (h, B.a, q);
		return h.is_one();
	}
}


//
// quadratic_ideal::is_subset()
//
// Task:
//      tests if the ideal is a subset of B.
//

bool
quadratic_ideal::is_subset(const quadratic_ideal & B) const
{
	debug_handler("quadratic_ideal", "is_subset");

	bool is_sub;
	bigrational temp;

	if (discriminant() != B.discriminant())
		return false;

	divide(temp, norm(), B.norm());
	is_sub = (is_bigint(temp));

	return is_sub;
}



//
// quadratic_ideal::is_proper_subset()
//
// Task:
//      tests if the ideal is a proper subset of B.
//

bool
quadratic_ideal::is_proper_subset(const quadratic_ideal & B) const
{
	debug_handler("quadratic_ideal", "is_proper_subset");

	return ((this->is_subset(B)) && (!this->is_equal(B)));
}



//
// operator ==
//
// Task:
//      tests if A and B are equal.
//

bool operator == (const quadratic_ideal & A, const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "operator == ");

	return A.is_equal(B);
}



//
// operator !=
//
// Task:
//      tests if A and B are not equal.
//

bool operator != (const quadratic_ideal & A, const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "operator != ");

	return !A.is_equal(B);
}



//
// operator <=
//
// Task:
//      tests if A is a subset of B.
//

bool operator <= (const quadratic_ideal & A, const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "operator <= ");

	return (A.is_subset(B));
}



//
// operator <
//
// Task:
//      tests if A is a proper subset of B.
//

bool operator < (const quadratic_ideal & A, const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "operator < ");

	return (A.is_proper_subset(B));
}



//
// operator >=
//
// Task:
//      tests if B is a subset of A.
//

bool operator >= (const quadratic_ideal & A, const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "operator >= ");

	return (B.is_subset(A));
}



//
// operator >
//
// Task:
//      tests if B is a proper subset of A.
//

bool operator > (const quadratic_ideal & A, const quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "operator >");

	return (B.is_proper_subset(A));
}



//
// operator !
//
// Task:
//      tests if A is the zero ideal.
//

bool operator ! (const quadratic_ideal & A)
{
	debug_handler("quadratic_ideal", "operator !");

	return (A.is_zero());
}



//
// swap()
//
// Task:
//      swap A and B.
//

void
swap(quadratic_ideal & A, quadratic_ideal & B)
{
	debug_handler("quadratic_ideal", "swap");

	quadratic_ideal C;

	C.assign(A);
	A.assign(B);
	B.assign(C);
}



//
// quadratic_ideal::smallest_rational()
//
// Task:
//      return the smallest rational number in the ideal.
//

bigrational
quadratic_ideal::smallest_rational() const
{
	debug_handler("quadratic_ideal", "smallest_rational");

	bigrational temp;

	if (reduced_flag)
		temp.assign_one();
	else
		multiply(temp, q, a);
	return temp;
}



//
// quadratic_ideal::is_integral()
//
// Task:
//      returns true if the ideal is integral, false otherwise.
//

bool
quadratic_ideal::is_integral() const
{
	debug_handler("quadratic_ideal", "is_integral");

	if (reduced_flag)
		return a.is_one();
	else
		return is_bigint(q);
}



//
// quadratic_ideal::conductor()
//
// Task:
//      returns the conductor of the ideal.
//

bigint
quadratic_ideal::conductor() const
{
	debug_handler("quadratic_ideal", "conductor");

	bigint temp;

	temp = gcd(a, b);
	return gcd(temp, get_c());
}



//
// quadratic_ideal::is_invertible()
//
// Task:
//      returns true if the ideal is invertible, false otherwise.
//

bool
quadratic_ideal::is_invertible() const
{
	debug_handler("quadratic_ideal", "is_invertible");

	bigint f;

	f = conductor();
	return (f.is_one());
}



//
// quadratic_ideal::norm()
//
// Task:
//      returns the norm of the ideal.
//

bigrational
quadratic_ideal::norm() const
{
	debug_handler("quadratic_ideal", "norm");

	bigrational nm;
	bigint temp;

	if (reduced_flag)
		nm.assign(conductor(), a);
	else {
		multiply(temp, a, conductor());
		square(nm, q);
		multiply(nm, nm, temp);
	}

	return nm;
}



//
// quadratic_ideal::ring_of_multipliers()
//
// Task:
//      returns the ring of multipliers of the ideal.
//

void
quadratic_ideal::ring_of_multipliers(quadratic_order & rm)
{
	debug_handler("quadratic_ideal", "ring_of_multipliers");

	bigint f, newDelta;

	f = conductor();
	newDelta = (QO->get_qo()->discriminant()) / (f*f);

	rm.assign(newDelta);
}


//
// quadratic_ideal::denominator()
//
// Task:
//      Returns the denominator d of the ideal I, i.e.,
//      d is the smallest positive integer with d*I integral
//      with respect to the order over which I is defined.
//
// Author:
//      Markus Maurer
//

bigint quadratic_ideal
::denominator () const
{
	debug_handler("quadratic_ideal", "denominator()");

	if (reduced_flag)
		return a;
	else
		return q.denominator();
}



//
// quadratic_ideal::multiply_by_denominator()
//
// Task:
//      Multiplies the ideal I by its denominator d, i.e.,
//      d is the smallest positive integer with d*I integral
//      with respect to the order over which I is defined.
//
// Author:
//      Markus Maurer
//

void quadratic_ideal
::multiply_by_denominator()
{
	debug_handler("quadratic_ideal", "multiply_by_denominator()");

	if (reduced_flag) {
		reduced_flag = false;
		q.assign_one();
	}
	else
		q.multiply_by_denominator();
}



//
// generate_prime_ideal(quadratic_ideal, bigint);
//
// Task:
//   computes a prime ideal lying over p.  The ideal will belong
//   to the quadratic order which was most recently referenced.  If no
//   quadratic orders have been defined, the error_handler will be called.
//   If no such ideal exists, false is returned.
//

bool
generate_prime_ideal(quadratic_ideal & A, const bigint & p)
{
	debug_handler("quadratic_ideal", "generate_prime_ideal(quadratic_ideal, "
		      "bigint)");

	return generate_prime_ideal(A, p, *quadratic_order::qo_l().last());
}



//
// generate_prime_ideal(quadratic_ideal, bigint, quadratic_order);
//
// Task:
//   computes a prime ideal in QO2 lying over p.  If such an ideal does not
//   exist, false is returned.
//

bool
generate_prime_ideal(quadratic_ideal & A, const bigint & p,
                     const quadratic_order & QO2)
{
	debug_handler("quadratic_ideal", "generate_prime_ideal(quadratic_ideal, "
		      "bigint, quadratic_order)");

	bigint Delta, b, temp;
	int kro, Dp, ip;

	Delta.assign(QO2.discriminant());
	if (Delta.is_zero()) {
		A.assign_zero(QO2);
		return false;
	}

	kro = kronecker(Delta, p);
	if (kro < 0)
		return false;
	else {
		if (p == 2) {
			if (kro == 0) {
				Dp = static_cast<int>(remainder(Delta, 16));
				if (Dp == 8)
					b.assign_zero();
				else
					b.assign(p);
			}
			else
				b.assign(1);
		}
		else {
			if (kro == 0) {
				remainder(temp, Delta, p*p);
				if (temp.is_zero())
					b.assign(p);
				else
					b.assign_zero();
			}
			else {
				if (is_int(p)) {
					p.intify(ip);
					Dp = static_cast<int>(remainder(Delta, static_cast<long>(ip)));
					if (Dp < 0)
						Dp += ip;

					b.assign(ressol(Dp, ip));
				}
				else {
					remainder(temp, Delta, p);
					if (temp.is_lt_zero())
						add(temp, temp, p);
					ressol(b, temp, p);
				}

				if (b.is_lt_zero())
					add(b, b, p);
			}

			if ((Delta.is_odd()) != (b.is_odd()))
				subtract(b, p, b);
		}

		A.assign(p, b, 1, QO2);
		return true;
	}
}



//
// quadratic_ideal::is_reduced()
//
// Task:
//      returns true if the ideal is reduced, false otherwise.
//
//      The ideal is reduced, if it is of the form
//      I = 1/a (Z a + Z (b+sqrt(D))/2) with (a,b,c) reduced, which is
//      equivalent to the fact that 1 is a minimum in I.
//

bool quadratic_ideal
::is_reduced() const
{
	debug_handler("quadratic_ideal", "is_reduced");

	if (reduced_flag)
		return true;

	bigint Delta, rootD, c, temp;

	Delta.assign(discriminant());

	if (a.is_zero())
		return true;

	if (Delta.is_lt_zero()) {
		c.assign(get_c());
		return (a.compare(c) < 0) || ((!a.compare(c)) && (b.is_ge_zero()));
	}
	else {
		sqrt(rootD, Delta);
		shift_left(temp, a, 1);
		subtract(temp, rootD, temp);
		if (temp.is_negative())
			inc(temp);
		temp.absolute_value();
		return (temp.compare(b) < 0) && (b.compare(rootD) <= 0);
	}
}



//
// quadratic_ideal::is_equivalent()
//
// Task:
//      returns true if the ideal is equivalent to B, false otherwise.
//

bool
quadratic_ideal::is_equivalent(const quadratic_ideal & B) const
{
	debug_handler("quadratic_ideal", "is_equivalent");

	bool equiv;

	if (discriminant() != B.discriminant())
		return false;

	equiv = false;
	if ((is_invertible()) && (B.is_invertible())) {
		qi_class A, Bc;
		qi_class::set_current_order(LIDIA_CAST_TO_QO(which_order()));
		A.assign(*this);
		Bc.assign(B);
		equiv = A.is_equivalent(Bc);
	}
	else {
		quadratic_form A, Bc;
		A.assign(*this);
		Bc.assign(B);
		equiv = A.is_equivalent(Bc);
	}

	return equiv;
}



//
// quadratic_ideal::is_principal()
//
// Task:
//      returns true if the ideal is principal, false otherwise.
//

bool
quadratic_ideal::is_principal() const
{
	debug_handler("quadratic_ideal", "is_principal");

	qi_class A;
	bool prin;

	prin = false;
	if (is_invertible()) {
		qi_class::set_current_order(LIDIA_CAST_TO_QO(which_order()));
		A.assign(*this);
		prin = A.is_principal();
	}

	return prin;
}



//
// quadratic_ideal::order_in_CL()
//
// Task:
//   returns the order to the equivalence class containing the ideal.
//   If the ideal is not invertible, 0 is returned. The qi_class function
//   of the same name is used for the computation.
//

bigint
quadratic_ideal::order_in_CL() const
{
	debug_handler("quadratic_ideal", "order_in_CL");

	qi_class A;
	bigint ord;

	ord.assign_zero();

	if (is_invertible()) {
		qi_class::set_current_order(LIDIA_CAST_TO_QO(which_order()));
		A.assign(*this);
		ord.assign(A.order_in_CL());
	}

	return ord;
}



//
// quadratic_ideal::DL()
//
// Task:
//  computes the smallest integer x such that G^x is equivalent to the
//  ideal, if x exists.  If x exists, true is returned, otherwise false is
//  returned and x is set to the order of the equivalence class containing
//  G. The qi_class function of the same name is used for the computation.
//

bool
quadratic_ideal::DL(quadratic_ideal & G, bigint & x) const
{
	debug_handler("quadratic_ideal", "DL");

	qi_class A, B;
	bool is_DL;

	if (discriminant() != G.discriminant()) {
		x.assign_zero();
		return false;
	}

	is_DL = false;
	x.assign_zero();

	if ((is_invertible()) && (G.is_invertible())) {
		qi_class::set_current_order(LIDIA_CAST_TO_QO(which_order()));
		A.assign(*this);
		B.assign(G);
		is_DL = A.DL(B, x);
	}

	return is_DL;
}



//
// subgroup()
//
// Task:
//    computes the structure of the subgroup generated by the ideals in G.
//    The qi_class function of the same name is used for the computation.
//

base_vector< bigint >
subgroup(base_vector< quadratic_ideal > & G)
{
	debug_handler("quadratic_ideal", "subgroup");

	lidia_size_t i, j, l;
	base_vector< qi_class > A;
	base_vector< bigint > S;
	qi_class Amem;

	A.set_mode(EXPAND);
	S.set_mode(EXPAND);
	A.reset();
	S.reset();

	qi_class::set_current_order(LIDIA_CAST_TO_QO(G[0].which_order()));

	l = G.size();
	j = 0;
	for (i = 0; i < l; ++i) {
		if ((qi_class::discriminant() == G[i].discriminant()) &&
		    (G[i].is_invertible())) {
			A[j] = qi_class(G[i]);
			++j;
		}
	}

	if (A.size() > 0)
		S = subgroup(A);

	return S;
}



//
// quadratic_ideal::regulator()
//
// Task:
//      returns the regulator of the quadratic order to which the ideal
//      belongs.  The error_handler is called if no quadratic order has been
//      assigned to this ideal or if the quadratic order no longer exists.
//

bigfloat
quadratic_ideal::regulator()
{
	debug_handler("quadratic_ideal", "regulator");

	if (!QO->get_qo()) {
		lidia_error_handler("quadratic_ideal", "regulator() - the quadratic "
				    "order of this ideal has been deleted");
		return bigfloat();
	}

	return QO->get_qo()->regulator();
}



//
// quadratic_ideal::class_number()
//
// Task:
//      returns the class number of the quadratic order to which the ideal
//      belongs.  The error_handler is called if no quadratic order has been
//      assigned to this ideal or if the quadratic order no longer exists.
//

bigint
quadratic_ideal::class_number()
{
	debug_handler("quadratic_ideal", "class_number");

	if (!QO->get_qo()) {
		lidia_error_handler("quadratic_ideal", "class_number() - the quadratic "
				    "order of this ideal has been deleted");
		return bigint();
	}

	return QO->get_qo()->class_number();
}



//
// quadratic_ideal::class_group()
//
// Task:
//      returns the class group of the quadratic order to which the ideal
//      belongs.  The error_handler is called if no quadratic order has been
//      assigned to this ideal or if the quadratic order no longer exists.
//

base_vector< bigint >
quadratic_ideal::class_group()
{
	debug_handler("quadratic_ideal", "class_group");

	if (!QO->get_qo()) {
		lidia_error_handler("quadratic_ideal", "class_group() - the quadratic "
				    "order of this ideal has been deleted");
		return base_vector< bigint > ();
	}

	return QO->get_qo()->class_group();
}




//
// operator >>
//
// Task:
//      inputs a quadratic_ideal from the std::istream in.
//

std::istream & operator >> (std::istream & in, quadratic_ideal & A)
{
	debug_handler("quadratic_ideal", "operator >>");

	int n = 0;
	char c;
	bigint ibuf[2];
	bigrational qbuf;

	in >> c;
	if (c != '(') {
		lidia_error_handler("quadratic_ideal", "operator >>::(expected");
		return in;
	}

	in >> c;
	while (c != ')' && n != 3) {
		in.putback(c);
		if (n == 0)
			in >> qbuf;
		else if (n < 3)
			in >> ibuf[n-1];

		n++;
		in >> c;
	}
	A.assign(ibuf[0], ibuf[1], qbuf);

	return in;
}



//
// operator <<
//
// Task:
//      outputs a quadratic_ideal to the std::ostream out.
//

std::ostream & operator << (std::ostream & out, const quadratic_ideal & A)
{
	debug_handler("quadratic_ideal", "operator << ");

	if (A.reduced_flag) {
		bigrational qq;
		qq.assign(bigint(1), A.a);
		out << "(" << qq << " " << A.a << " " << A.b << ")";
	}
	else
		out << "(" << A.q << " " << A.a << " " << A.b << ")";

	return out;
}



//
// assign_module_of(nq,q1,q2)
//
// Task:
//      Computes q * (Za + Z (b+sqD)/2) = nq * (Z q1 + Z q2),
//
//      with a,b in Z, q in Q.
//
// Note: The function does NOT assign an order to the object.
//       It only initializes the coefficients (q,a,b).
//
// Author: Markus Maurer
//

void quadratic_ideal
::assign_module_of (const bigrational & nq,
		    quadratic_number_standard q1,
		    quadratic_number_standard q2)
{
	debug_handler ("quadratic_ideal",
		       "assign_module_of(const bigrational&, "
		       "const quadratic_number&, const quadratic_number_standard&, "
		       "const quadratic_order&)");

#ifdef DEBUG
	std::cout << "assign_module_of::(bigrational q, qn q1, qn q2, qo O)" << std::endl;
	std::cout << "q = " << nq << std::endl;
	std::cout << "q1 = " << q1 << std::endl;
	std::cout << "q2 = " << q2 << std::endl;
	std::cout << "O = " << O << std::endl;
#endif

	bigint x, y, d, d1, d2;
	quadratic_number_standard h1, h2, h3;

	// q1 = (a1 + b1 sqrt(D)) / d1
	// q2 = (a2 + b2 sqrt(D)) / d2
	//
	// I = q ( Z q1 + Z q2 )

	//
	// Compute a transformation T such that
	// the basis of I/(q/(d1*d2)) will be of the form
	//
	// I/(q/(d1*d2)) = Z a3 + Z (a4 + b4 sqrt(Delta))
	//
	// with a3, a4, b4 in Z.
	//

	// Write I = q ( Z q1 + Z q2 ) with denominator
	// of q1 and q2 equal to one.
	d1 = q1.get_d();
	d2 = q2.get_d();

	divide (q, nq, d1);
	divide (q, q, d2);

	q1.multiply_by_denominator();
	q2.multiply_by_denominator();
	multiply (q1, q1, d2);
	multiply (q2, q2, d1);

	// determine the xgcd of the sqrt(Delta) coefficients
	d1 = q1.get_b();
	d2 = q2.get_b();

#ifdef DEBUG
	std::cout << "d1 = " << d1 << std::endl;
	std::cout << "d2 = " << d2 << std::endl;
#endif

	d = xgcd(x, y, d1, d2);

	// apply unimodular transformation matrix  ( d2/d  x)
	//                                         (-d1/d  y)
	divide (d2, d2, d);
	divide (d1, d1, d);
	d2.negate();

#ifdef DEBUG
	std::cout << "trafo matrix: (" << d2 << "  " << x << ")" << std::endl;
	std::cout << "              (" << d1 << "  " << y << ")" << std::endl;
#endif

	multiply(h2, d2, q1);
	multiply(h3, d1, q2);
	add(h1, h2, h3);

	multiply(h2, x, q1);
	multiply(h3, y, q2);
	add(h2, h2, h3);

#ifdef DEBUG
	std::cout << "transformation matrix applied:" << std::endl;
	std::cout << "h1 = " << h1 << std::endl;
	std::cout << "h2 = " << h2 << std::endl;
#endif

	// build the standard representation; we should have
	// I = q * (Z h1 + Z h2)
	// with
	// h1 = a3
	// h2 = a4 + b4 sqrt(D)
	//
	// and b4 | a4 and (2 b4) | a3.
	//
	// Hence: I = Q ( Z a + Z (b+sqrt(Delta))/2 )
	//        with
	//        Q = q * (2 b4), a = a3 / (2 b4), b = a4 / b4.

	if (!h1.get_b().is_zero())
		lidia_error_handler("quadratic_ideal::assign_module_of("
				    "const bigrational&, "
				    "const quadratic_number&, const quadratic_number&,)"
				    "const quadratic_order&",
				    "h1 not in Z.");

	// determine b = a4 / b4
	div_rem(b, y, h2.get_a(), h2.get_b());
	if (!y.is_zero())
		lidia_error_handler("quadratic_ideal::assign_module_of("
				    "const bigrational&, "
				    "const quadratic_number&, const quadratic_number&,)"
				    "const quadratic_order&",
				    "b4 does not divide a4.");

	// determine a = a3 / (2 b4)
	shift_left (d1, h2.get_b(), 1);
	div_rem(a, y, h1.get_a(), d1);
	a.absolute_value();

	if (!y.is_zero())
		lidia_error_handler("quadratic_ideal::assign_module_of("
				    "const bigrational&, "
				    "const quadratic_number&, const quadratic_number&,)"
				    "const quadratic_order&",
				    "2 b4 does not divide a3.");

	// determine q = q / (2 b4)
	multiply(q, q, d1);
	q.absolute_value();

#ifdef DEBUG
	std::cout << "O = " << O << std::endl;
	std::cout << "q = " << q << std::endl;
	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << b << std::endl;
#endif
}




//
// assign(bigrational,quadratic_number_standard,quadratic_number_standard)
//
// Task:
//      If q * (Z q1 + Z q2) is a quadratic ideal,
//      i.e., a Z module of rank 2, the function returns true, and
//      assigns the standard representation. Otherwise the function
//      returns false, and initializes the ideal with the order.
//
// Author: Markus Maurer
//

bool quadratic_ideal
::assign (const bigrational & nq,
	  const quadratic_number_standard & q1,
	  const quadratic_number_standard & q2)
{
	debug_handler ("quadratic_ideal",
		       "assign(const bigrational&, "
		       "const quadratic_number&, "
		       "const quadratic_number_standard &, "
		       "const quadratic_order&)");

	this->assign_module_of(nq, q1, q2);

	if ((this->a).is_zero()) {
		this->assign_one(q1.get_order());
		return false;
	}
	else {
		this->reduced_flag = false;
		QO = quadratic_order::qo_l().add_to_list(
			QO, LIDIA_CAST_TO_QO(q1.get_order()));
		this->normalize();
		return true;
	}
}


//
// multiply(J,I,g)
//
// Task:
//       J = I * g
//
// Author: Markus Maurer
//

void multiply (quadratic_ideal & J,
	       const quadratic_ideal & I,
	       const quadratic_number_standard & g)
{
	debug_handler ("quadratic_ideal",
		       "multiply(quadratic_ideal&, const quadratic_ideal&, "
		       "const quadratic_number_standard &)");

	quadratic_number_standard h;
	bigint a2;

	//
	// I = q * (Z a + Z (b+sqrt(D))/2).
	//   = q*a*(Z + Z (b+sqrt(D))/2a)
	//
	// Set J.q = q * a
	//     h   = g * (b+sqrt(D))/2a
	//
	// and J = J.q * (Z g + Z h).
	//

	if (I.reduced_flag)
		J.q.assign_one();
	else
		multiply (J.q, I.q, I.a);

	shift_left(a2, I.a, 1);
	h.assign_order(I.which_order());
	h.set(I.b, bigint(1), a2);
	multiply(h, h, g);

	J.assign(J.q, g, h);
	J.reduced_flag = false;
}


//
// divide(J,I,g)
//
// Task:
//       J = I / g
//
// Author: Markus Maurer
//

void divide (quadratic_ideal & J,
	     const quadratic_ideal & I,
	     const quadratic_number_standard & g)
{
	debug_handler ("quadratic_ideal",
		       "divide(quadratic_ideal&, const quadratic_ideal&, "
		       "const quadratic_number_standard &)");

	quadratic_number_standard h;
	inverse(h, g);
	multiply(J, I, g);
}



//
// local_close
//
// Let I = *this. The function transforms I into the reduced ideal
// J = I / alpha, where alpha is a minimum in I that is
// close to t with regard to k-1. It also computes an absolute
// k-approximation "a" to Ln alpha. The function uses reduction, rho, and
// inverse_rho operations to determine J.
//
// Author: Markus Maurer
//

void quadratic_ideal
::local_close(quadratic_number_standard & alpha,
	      xbigfloat & a,
	      xbigfloat t,
	      long k)
{
	if (this->discriminant() < 0) {
		lidia_error_handler("quadratic_ideal::local_close",
				    "No real quadratic ideal.");
		return;
	}

	quadratic_ideal B, C;
	quadratic_number_standard beta, gamma, mu;
	xbigfloat b, c;

	C = *this;

	// Reduce I by division by gamma
	// and approximate Ln gamma by c.
	C.reduce(gamma);
	gamma.get_absolute_Ln_approximation(c, k);

	// Choose direction and find b <= t < c
	if (c <= t) {
		while (c <= t) {
			// Copy C to B
			B = C; beta = gamma; b = c;

			// Move forward with C
			C.inverse_rho(mu);
			gamma *= mu;
			gamma.get_absolute_Ln_approximation(c, k);
		}
	}
	else {
		B = C; beta = gamma; b = c;
		while (b > t) {
			// Copy B to C
			C = B; gamma = beta; c = b;

			// Move backward with B
			B.rho(mu);
			beta *= mu;
			beta.get_absolute_Ln_approximation(b, k);
		}
	}

	// Choose closest
	if (t - b < c - t) {
		*this = B;
		alpha = beta;
		a = b;
	}
	else {
		*this = C;
		alpha = gamma;
		a = c;
	}
}


//
// order_close
//
// Let I = *this be a quadratic order. The function transforms I
// into the reduced ideal J = I / alpha, where alpha is a minimum
// in I that is  close to t with regard to k-1. It also computes an
// absolute k-approximation "a" to Ln alpha. The function uses
// repeated squaring to determine J.
//
// Author: Markus Maurer
//

void quadratic_ideal
::order_close (quadratic_number_power_product & alpha,
	       xbigfloat & a,
	       xbigfloat t,
	       long k)
{
	if (this->discriminant() < 0) {
		lidia_error_handler("quadratic_ideal::order_close",
				    "No real quadratic ideal.");
		return;
	}

	long n, j;
	quadratic_number_standard beta;
	xbigfloat b;

	if (!this->is_one())
		lidia_error_handler("quadratic_ideal::order_close",
				    "Ideal is not an order.");

	alpha.assign_one(this->get_order());
	a.assign_zero();

	n = t.b_value() + 1;

	if (n <= 0) {
		if (k >= 2)
		{ this->local_close (beta, b, t, k); alpha = beta; a = b; }
	}
	else {
		shift_right(t, t, n);

		for (j = 1; j <= n; j++) {
			square(*this, *this);
			alpha.square(alpha);
			a.multiply_by_2();
			t.multiply_by_2();

			this->local_close (beta, b, t-a, k+1+2*(n-j));
			alpha.multiply(alpha, beta);
			a += b;
		}
	}
}


//
// close
//
// Let I = *this. The function transforms I into the reduced ideal
// J = I / alpha, where alpha is a minimum in I that is
// close to t with regard to k-1. It also computes an absolute
// k-approximation "a" to Ln alpha. It uses order_close and local_close.
//
// Author: Markus Maurer
//

void quadratic_ideal
::close (quadratic_number_power_product & alpha,
	 xbigfloat & a,
	 xbigfloat t,
	 long k)
{
	if (this->discriminant() < 0) {
		lidia_error_handler("quadratic_ideal::close",
				    "No real quadratic ideal.");
		return;
	}

	quadratic_ideal J;
	quadratic_number_power_product beta;
	quadratic_number_standard gamma;
	xbigfloat b, c;

	J.assign_one(this->get_order());
	J.order_close(beta, b, t, k+2);

	*this *= J;
	this->local_close (gamma, c, t-b, k+1);

	alpha.multiply(beta, gamma);
	add(a, b, c);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
