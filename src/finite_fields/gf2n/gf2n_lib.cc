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
#include	"LiDIA/gf2n.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//******* static variables for fast solve_quadratic *************

gf2n * gf2n::table_solve_quadratic;
bigint * gf2n::table_solve_quadratic_masks = NULL;
int gf2n::table_solve_quadratic_number_masks = 0;
bool gf2n::table_solve_quadratic_init = false;


//**** relative degree of element over prime field **************

unsigned int relative_degree(const gf2n & a)
{
	return a.relative_degree();
}



//***** absolute degree over prime field *************************

unsigned int get_absolute_degree(const gf2n & x)
{
	gf2n a(x);

	return gf2n::get_absolute_degree();
}



//****** solve_quadratic ****************************************
// Input: polynomial X^2 + a1*X + a0
// Output: true and set root, if polynomial splits,
//         false otherwise
//***************************************************************

bool solve_quadratic(gf2n & root, const gf2n & a1, const gf2n & a0)
{
	gf2n c, h;
	gf2n c0, c1; // coefficients of c1*z + c0 mod (z^2+z+c)
	gf2n res0, res1; // coefficients of res1*z + res0 mod (z^2+z+c) 
	int i, degree;

	// check trivial cases

	if (a1.is_zero()) {
		sqrt(root, a0);
		return true;
	}

	if (a0.is_zero()) {
		root.assign_zero();
		return true;
	}

	// variable transformation: solve z^2 + z + c = 0 for c= f0/f1^2
	// -->root = z*f1 is solution for original equation

	square(c, a1);
	c.invert();
	multiply(c, c, a0);

	if (c.trace() != 0)  // polynomials splits iff trace = 0
		return false;

	degree = gf2n::get_absolute_degree();

	int d = 1, counter = 0;
	c1.randomize(degree);
	res1.assign(c1);
	c0.randomize(degree);
	res0.assign(c0);

	do {
		// start trace computation for random polynomial

		// c1*z + 0 modulo (z^2+z+c) 
		counter ++;

		for (i = 1; i < degree; i++) {
			// (c1*z+c0) < -- (c1*z+c0)^2
			square(c1, c1);
			square(c0, c0);
			multiply(h, c1, c);
			add(c0, c0, h);
			add(res1, res1, c1);
			add(res0, res0, c0);
		}
		square(c1, c1);
		square(c0, c0);
		multiply(h, c1, c);
		add(c0, c0, h);

		if (!res1.is_zero()) {
			// degree one->check whether root found

			if (!res1.is_one())      // make monic
				multiply(res0, res0, inverse(res1));

			// z + c0 should be factor of z^2+z+c, check this

			square(h, res0);
			add(h, h, res0);
			add(h, h, c);
			if (h.is_zero())  // root 'c0' found
				break;
		}

		if (counter >= 5) {
			counter = 0;
			if (d == degree)
				d = 1;
			else
				d++;

			while (degree % d != 0 && d < degree)
				d++;
		}

		c1.randomize(degree / d);
		res1.assign(c1);
		c0.randomize(degree);
		res0.assign(c0);
	}
	while (true);

	multiply(root, res0, a1);
	return (true);
}



//****** characteristic *******************************************

bigint characteristic(const gf2n & a)
{
	gf2n b(a);

	return bigint(2);
}



//***** number_of_elements ****************************************

bigint number_of_elements(const gf2n & a)
{
	return (bigint(1) << a.extension_degree());
}



//--------------------------------------------------------------------
// The following functions are worthwhile  for solving many quadratic
// polynomial equations only (as in point counting). You have to first
// initialize a table with the next function.
//
// Basic idea: we determine roots for polynomials
//
//   y^2 + y + x^j = 0  where x is the main variable of the polynomial
//                      basis chosen for GF(2^n)
// These are then combined.

void gf2n::initialize_table_for_solve_quadratic()
{
	if (gf2n::table_solve_quadratic_init == true)
		return;

	bigint h;
	gf2n gf2n_one, z;
	int *trace_1;
	unsigned int number_trace_1 = 0;
	unsigned int i, j;

	gf2n_one.assign_one();

	trace_1 = new int[gf2n::degree + 1];

	gf2n::table_solve_quadratic = new gf2n[gf2n::degree + 40];
	gf2n::table_solve_quadratic_init = true;

	for (i = 0; i < gf2n::degree; i++) {
		shift_left(h, bigint(1), i);
		z.assign(h);
		if (!solve_quadratic(gf2n::table_solve_quadratic[i], gf2n_one, z)) {
			// y^2 + y + x^i has no solution. This is stored as solution 0.
			trace_1[number_trace_1] = i;
			number_trace_1 ++;
			gf2n::table_solve_quadratic[i].assign_zero();
		}
	}

	// if there are more than 2 j such that y^2 + y + x^j has no root,
	// we try combinations of these (y^2 + y + (x^{j1} + x^{j2}). These
	// must have a root, since it's trace is 0!!

	if (number_trace_1 >= 2) {
		// there is more than 1 irreducible polynmomial
		int d = number_trace_1 * (number_trace_1 - 1) / 2;

		if (d >= 40) {
			// enlarge table and copy already computed entries
			gf2n * t;

			t = new gf2n [gf2n::degree + d + 1];

			for (i = 0; i < gf2n::degree; i++)
				t[i].assign(gf2n::table_solve_quadratic[i]);

			delete[] gf2n::table_solve_quadratic;
			gf2n::table_solve_quadratic = t;
		}

		gf2n::table_solve_quadratic_masks = new bigint[d];
		gf2n::table_solve_quadratic_number_masks = 0;

		for (i = 0; i < number_trace_1; i++)  // now try combinations
			for (j = i+1; j < number_trace_1; j++) {
				shift_left(h, bigint(1), trace_1[i]);
				add(h, h, bigint(1) << trace_1[j]);
				z.assign(h);
				if (solve_quadratic(z, gf2n_one, z)) {
					// has a solution
					gf2n::table_solve_quadratic[gf2n::degree +
								   gf2n::table_solve_quadratic_number_masks].assign(z);
					gf2n::table_solve_quadratic_masks[gf2n::table_solve_quadratic_number_masks].assign(h);
					gf2n::table_solve_quadratic_number_masks++;
				}
			}
		delete[] trace_1;
		return;
	}
	delete[] trace_1;
}



//--------------------------------------------------------------------
// This function deletes the static tables to release the memory. They
// must be called per hand.

void gf2n::delete_table_for_solve_quadratic()
{
	if (gf2n::table_solve_quadratic_init == true) {
		delete[] gf2n::table_solve_quadratic;
		if (gf2n::table_solve_quadratic_masks != NULL)
			delete[] gf2n::table_solve_quadratic_masks;
		gf2n::table_solve_quadratic_masks = NULL;
		gf2n::table_solve_quadratic_init = false;
		gf2n::table_solve_quadratic_number_masks = 0;
	}
}



//-------------------------------------------------------------------------
// The main function foir solving quadratic polynomials with tables.
//
//  Tries to find a root of Y^2 + a1*Y + a0 = 0. If this polynomial is
//  irreducible, then 'false' is returned; otherwise the functions returns
//  tru and sets 'root' to one found root. If the bolean variable
//  'sure_has_root' is true, then the existence of a root is not checked, but
//   assumed.
//
//  The algorithm uses the information computed in the tables.

bool solve_quadratic_with_table(gf2n & root, const gf2n & a1, const gf2n & a0,
				bool sure_has_root)
{
	gf2n c, res0; // coefficients of res1*z + res0 mod (z^2+z+c) 


	//   check trivial cases first

	if (! gf2n::table_solve_quadratic_init)  // no tables initialized -->d
		return solve_quadratic(root, a1, a0); // standar version

	if (a1.is_zero()) {
		sqrt(root, a0);
		return true;
	}

	if (a0.is_zero()) {
		root.assign_zero();
		return true;
	}

	// variable transformation: solve Z^2 + Z + c = 0 for c = f0/f1^2
	// -->root = z*f1 is solution for original equation

	square(c, a1);
	c.invert();
	multiply(c, c, a0);

	if (!sure_has_root)
		if (c.trace() != 0)  // polynomials splits if and only if trace = 0
			return false;

	res0.assign_zero();
	if (!c.is_reduced())
		partial_reduce2[gf2n::invsel](c.element);

	gf2n_word mask; // now we check the bits in the polynomial rep. of c
	bigint h(0);
	int bad_bits = 0, j;
	unsigned int i;

	for (i = 0; i < gf2n::degree; i++) {
		mask = 1 << (i % GF2N_WORDSIZE);
		if (c.element[i / GF2N_WORDSIZE] & mask)
			if (!gf2n::table_solve_quadratic[i].is_zero())   // this x^i factorizes
				add(res0, res0, gf2n::table_solve_quadratic[i]);
			else {
				if ((bad_bits & 1) == 0)       // there is no pair yet
					shift_left(h, bigint(1), i);
				else {
					add(h, h, bigint(1) << i);

					// otherwise we have two indices stored in h
					// find a precomputed value that has these two bits set

					for (j = 0; j < gf2n::table_solve_quadratic_number_masks; j++)
						if (h == gf2n::table_solve_quadratic_masks[j]) {
							add(res0, res0, gf2n::table_solve_quadratic [gf2n::degree
												    + j ]);
							break;
						}
					h.assign_zero();
				}
				bad_bits ++;
			}
	}

	multiply(root, res0, a1); // compute the solution for original equation
	return (true);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
