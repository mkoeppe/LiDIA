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
//	Author	: Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/gf2n_poly_modulus.h"
#include	"LiDIA/finite_fields/Fp_polynomial_util.h"
#include	"LiDIA/udigit.h"
#include        "LiDIA/arith.inl"
#include	"crossover.tbl"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



const int MAX_TABLE = 500;


gf2n_polynomial gf2n_poly_modulus::help_pol;

//******** constructor *******************************************

gf2n_poly_modulus::gf2n_poly_modulus(const gf2n_poly_modulus & p)
{
	int i;
	used_table_entries = p.used_table_entries;
	size_tab = p.size_tab;
	mod.assign(p.mod);

	tab = new gf2n_polynomial [size_tab];
	for (i = 0; i < size_tab; i++)
		tab[i].assign(p.tab[i]);
}



void gf2n_poly_modulus::assign(const gf2n_poly_modulus & a)
{
	if (this != &a) {
		int i;

		if (tab != NULL)
			delete [] tab;

		used_table_entries = a.used_table_entries;
		size_tab = a.size_tab;
		tab = new gf2n_polynomial[size_tab];
		for (i = 0; i < size_tab; i++)
			tab[i].assign(a.tab[i]);
		mod.assign(a.mod);
	}
}



//************ build a gf2n_poly_modulus, initialization ************

void gf2n_poly_modulus::build(const gf2n_polynomial & modulus)
{
	if (modulus.deg < 2)
		lidia_error_handler("gf2n_poly_modulus", "set_modulus::degree of "
				    "modulus < 2");

	if (!modulus.is_monic())
		lidia_error_handler("gf2n_poly_modulus", "set_modulus::modulus not monic");

	int d = modulus.deg;

	help_pol.set_size(2*d);
	mod.assign(modulus);

	//  tab[i] = recip_{i+d} mod X^{d+i} for i = 0,...,d-2
	//  the table is computed on the fly, i.e. if a value is needed, then
	//  it is computed and stored in the table.
	//  There is a bound MAX_TABLE for the maximal number of entries allowed
	//  in the table. If the size is bigger than this bound, then only "the
	//  hig end" values are stored (i.e. tab[size_tab-MAX_TAB], ... ,
	//  tab[size_tab - 1]. We use the heuristic that computation of the
	//  reciprocal for a polynomial of small degree is faster than large degree
	//  reciprocals.

	if (tab != NULL)
		delete[] tab;

	used_table_entries = 0;
	size_tab = d-1;
	tab = new gf2n_polynomial[size_tab];

	gf2n_polynomial::resize_stack_for_degree(2*d);

	int i;

	for (i = 0; i < size_tab; i++)
		tab[i].assign_zero();
}



//---------------------------------------------------------
// plain_rem means that we do not use the table !!

void plain_rem(gf2n_polynomial & x, const gf2n_polynomial & a,
	       gf2n_poly_modulus & F)
{
	if (a.deg < F.mod.deg) {
		x.assign(a);
		return;
	}

	remainder(x, a, F.mod);
}



//------------------------------------------------------------------

void kara_rem(gf2n_polynomial & x, const gf2n_polynomial & a,
	      gf2n_poly_modulus & F)
{
	if ((a.deg <= F.mod.deg) ||
	    (a.deg > 2*F.mod.deg-2)) {
		remainder(x, a, F.mod);
		return;
	}

	lidia_size_t deg_a = a.degree(), deg_b = F.mod.degree();
	lidia_size_t deg_q = deg_a - deg_b;

	gf2n_polynomial P1(a.degree()), P2(a.degree()), P3(a.degree());

	if (F.tab[deg_q].is_zero())  // not yet computed !!
	{
		copy_reverse(P3, F.mod, 0, deg_b);
		invert(P2, P3, deg_q + 1);

		if ((F.size_tab > MAX_TABLE) && (deg_q < F.size_tab - MAX_TABLE))
			copy_reverse(P1, P2, 0, deg_q);
		else {
			copy_reverse(F.tab[deg_q], P2, 0, deg_q);
			P1.assign(F.tab[deg_q]);
			F.used_table_entries ++;
		}
	}
	else
		P1.assign(F.tab[deg_q]);

	floor(P3, a, deg_b);
	multiply(P2, P3, P1);
	trunc(P2, P2, 2*deg_q+1);
	floor(P3, P2, deg_q);
	multiply(P2, P3, F.mod);
	add(x, a, P2);
}



void remainder(gf2n_polynomial & x, const gf2n_polynomial& a,
	       gf2n_poly_modulus & F)
{
	if (F.mod.degree() > 4)
		kara_rem(x, a, F);
	else
		plain_rem(x, a, F);
}



//*************  functions with reduction: multiply, square ************


void multiply(gf2n_polynomial & c, const gf2n_polynomial & a,
	      const gf2n_polynomial & b, gf2n_poly_modulus &F)
{
	if (a.deg + b.deg < F.mod.deg)
		multiply(c, a, b);
	else {
		multiply(gf2n_poly_modulus::help_pol, a, b);
		remainder(c, gf2n_poly_modulus::help_pol, F);
	}
}



//-----------------------------------------------------------------

void square(gf2n_polynomial & c, const gf2n_polynomial & b,
	    gf2n_poly_modulus & F)
{
	if ((b.deg << 1) < F.mod.deg)
		square(c, b);
	else {
		square(gf2n_poly_modulus::help_pol, b);
		if (F.mod.degree() < CROSSOVER_REDUCE_SQUARE)
			remainder(c, gf2n_poly_modulus::help_pol, F.mod);
		else
			kara_rem(c, gf2n_poly_modulus::help_pol, F);
	}
}



//****************************************************************

void shift_left(gf2n_polynomial & erg, const gf2n_polynomial & a,
		gf2n_poly_modulus & F, unsigned int d)
{
	shift_left(erg, a, d);
	remainder(erg, erg, F);
}



//*********************************************************************
//  uses squaring trick to compute xq^i mod F.mod for i = 0, ..., d -1
//  no deallocation of memory !!

gf2n_polynomial * power_table(const gf2n_polynomial & xq,
			      gf2n_poly_modulus & F, int d)
{
	gf2n_polynomial * table;
	int i, j;

	if (d < 2)
		lidia_error_handler("gf2n_poly_modulus", "power_table::index d < 2");

	table = new gf2n_polynomial [d];
	table[0].assign_one();
	table[1].assign(xq);

	if (d == 3) {
		square(table[2], table[1], F);
		return table;
	}

	i = 2; j = 3;

	while (j < d) {
		while (i < d) {
			square(table[i], table[i/2], F);
			i <<= 1;
		}
		i = j; j += 2;
		multiply(table[i], table[i-1], xq, F);
		i <<= 1;
	}
	return table;
}



//****** high level functions used in factoring *************************
//
// Function determines erg = X^q mod F.mod. In addition to the usual
// power function for general polynomials, we also implemented a blocked
// version. First we determine the optimal blocking factor with the following
// cost function (experimental result for interesting parameters):
//
// If blocking factor d is chosen, then the expected cost are
//   (q % d) + max(d - bit_length(deg(p)), 0) + deg(pol)/2  square_mods
//   (deg(pol)-1)/2 mult_mod
//   floor (p/d) reductions
//
// Since blocking strategy needs computation of a new reduction table, it
// can be worse than power. If the optimal blocking factor is therefore 1,
// then we use the power function; otherwise we use the blocking version.

const float cost_red = 3.3;
const float cost_mult_red = 1.6; // precomputed values
const float cost_square_red = 1.0;
const int MAX_BLOCK_FACTOR = 100;


void Xq (gf2n_polynomial & erg, gf2n_poly_modulus & F, unsigned int d)
{
	int dp = F.mod.deg;
	int block_factor = 1, i, j, s;
	float minimum = d*cost_square_red, cost = 0.0;

	//-------- find optimal blocking factor first ---------------

	j = 2;
	i = static_cast<int>(LiDIA::log2(static_cast<double>(dp)));

	while (j <= MAX_BLOCK_FACTOR) {
		cost = (dp / 2 + (d % j))
			* cost_square_red + (dp-1)/2 * cost_mult_red
			+ (d / j) * cost_red;
		if (j > i)
			cost += (j-i) * cost_square_red;

		if (cost < minimum) {
			block_factor = j; minimum = cost;
		}
		j++;
	}

	if (block_factor == 1)  // don't use blocking strategy,
	{
		gf2n_polynomial x(1);
		power(erg, x, bigint(1) << d, F);
		return;
	}

	// now we use blocking: block_factor many mults are done before reduction
	// we need a power_table for this !!

	gf2n *h2;
	h2 = new gf2n[dp];

	for (i = 1; i < dp; i++)
		h2[i].assign_zero();

	erg.set_size(dp);

	if (static_cast<unsigned int>(1 << block_factor) < static_cast<unsigned int>(dp) &&
	    block_factor < static_cast<int>(CHAR_BIT * SIZEOF_INT))
		erg.assign(gf2n_polynomial(1 << block_factor));
	else {
		j = static_cast<int>(LiDIA::log2(static_cast<double>(dp)));
		erg.assign(gf2n_polynomial(1 << j));

		j = block_factor - j;
		while (j > 0) {
			square(erg, erg, F);
			j--;
		}
	}

	// compute a power_table for faster reduction

	gf2n_polynomial * table = power_table (erg, F, dp);
	gf2n tt, ht;
	gf2n_polynomial *ptr;

	for (j = 2; j <= static_cast<int>(d / block_factor); j++) {
		// the outer loop
		tt.assign(erg.coeff[0]); // the constant term
		if (!tt.is_zero())
			for (i = 0; i < block_factor; i++)
				square(tt, tt);
		h2[0].assign(tt);

		for (s = 1; s < dp; s++) {
			// all other terms
			tt.assign(erg.coeff[s]);
			if (tt.is_zero())
				continue;

			for (i = 0; i < block_factor; i++)
				square(tt, tt);

			ptr = table+s;
			int end = dp-1;

			if (table[s].deg < dp-1)
				end = table[s].deg;

			for (i = 0; i <= end; i++) {
				multiply(ht, tt, ptr->coeff[i]);
				add(h2[i], h2[i], ht);
			}
		}

		for (i = 0; i < dp; i++) {
			// copying and reinitialization
			erg.coeff[i].assign(h2[i]);
			h2[i].assign_zero();
		}
	} // of outer loop

	erg.deg = dp-1; // new degree computation
	while (erg.coeff[erg.deg].is_zero() && erg.deg > 0)
		erg.deg--;

	if (erg.deg == 0 && erg.coeff[0].is_zero())
		erg.deg = -1;

	if (d % block_factor)
		for (j = d % block_factor; j; j--)
			square(erg, erg, F);

	delete[] h2;
	delete[] table;
}



//******************************************************************

void power(gf2n_polynomial & c, const gf2n_polynomial & aa, const bigint & bb,
	   gf2n_poly_modulus & F)
{
	gf2n_polynomial a(aa);
	bigint b(bb);

	if (b.is_negative()) {
		negate(b, b);
		if (!xgcd_left(a, a, F.mod).is_one()) {
			lidia_error_handler("gf2n_poly_modulus", "power::negative exponent"
					    " for non-invertible polynomial");
			return;
		}
	}

	if (b.is_zero() || a.is_one())
		c.assign_one();
	else
		if (b.is_one())
			c.assign(a);
		else {
			gf2n_polynomial p(a);
			int i;

			for (i = b.bit_length()-2; i >= 0; i--) {
				square(p, p, F);
				if (b.bit(i))
					multiply(p, a, p, F);
			}
			c.assign(p);
		}
}



//*********** Functions for Atkin: compute degree **********************
// determines res = f(X) for polynomial X, uses table with
// powers X^i mod if possible (i.e. with pointer != NULL)
// Alias res, f is allowed
//
//  1st version: uses table
//  2nd version: uses gf2n_poly_modulus for reduction and Horner

void compose(gf2n_polynomial & res, const gf2n_polynomial & f,
             const gf2n_polynomial & x, const gf2n_polynomial *table)
{
	register int i, fd = f.degree();
	gf2n_polynomial r(0), help;
	gf2n h2;

	if (f.is_zero()) {
		res.assign_zero();
		return;
	}

	if (fd == 0) {
		res.assign(f);
		return;
	}

	if (x != table[1])
		lidia_error_handler("gf2n_polymod", "compose::Input table differs"
				    " from input polynomial");

	r.set_coefficient(f.const_term(), 0);
	for (i = 1; i <= fd; i++) {
		f.get_coefficient(h2, i);
		multiply_by_scalar(help, h2, table[i]);
		add(r, r, help);
	}
	res.assign(r);
}



void compose(gf2n_polynomial & res, const gf2n_polynomial & f,
             const gf2n_polynomial & x, gf2n_poly_modulus & F)
{
	register int i, fd = f.degree();
	gf2n_polynomial r(0);
	gf2n h2;

	if (f.is_zero()) {
		res.assign_zero();
		return;
	}

	if (fd == 0) {
		res.assign(f);
		return;
	}

	multiply(r, f.lead_coeff(), x, F);
	i = fd - 1;
	while (i > 0) {
		add(r, r, f.coeff[i]);
		i--;
		multiply(r, r, x, F);
	}
	add(res, r, f.const_term());
}



//
// Task:	The common degree of the irreducible factors of F.mod
//              are computed.
//
// Conditions:	F.mod is assumed to be an "equal degree" polynomial,
//		Xq = x^q mod F.mod
//              d is multiple of common degree, if -1, no multiple is known
//
//              Usually we use an early abort strategy. Note that d MUST
//              be a multiply of common degree or -1.


lidia_size_t compute_degree(const gf2n_polynomial & Xq,
			    gf2n_poly_modulus & F,
			    lidia_size_t d)
{
	gf2n_polynomial help;

	if (d == -1)
		d = F.mod.degree();

	lidia_size_t r, stop_comp;
	gf2n_polynomial *table;

	table = power_table(Xq, F, F.mod.deg);

	help.assign(Xq);

	r = 2; // compute the maximal non-trivial divisor of d
	while (d % r != 0)
		r++;
	stop_comp = d / r;

#ifdef DEBUG
	stop_comp = d;
#endif

	for (r = 2; r <= stop_comp; r++)
	{                   // now compute h = X^(q^r) mod F.mod
		compose(help, help, Xq, table);

		if (help.is_x()) {
			delete[] table;
			return r;
		}

#ifdef DEBUG
		gf2n_polynomial check;

		gcd(check, help + gf2n_polynomial (1), F.mod);
		if (check.degree() > 0) {
			gf2n_polynomial hx;
			remainder(check, F.mod, check);
			assert(check.is_zero());
			lidia_error_handler("compute_degree", "Factor found -->no (d...d)"
					    " splitting type");
		}
#endif

	}

#ifdef DEBUG
	lidia_error_handler("gf2n_polynomial", "compute_degree::no degree found");
#endif

	delete[] table;
	return d;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
