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
#include	"LiDIA/eco_gf2n.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//---------------------------------------------------------------------
// compute res0 + res * Y = (res0 + res1*Y)^2 mod f
//

void eco_gf2n::square_y_term(ff_pol & res0, ff_pol & res1,
			     const ff_pol & curve, ff_polmod & f)
{
	ff_pol h;

	square(h, res1, f);
	shift_left(res1, h, f, 1);
	multiply(h, h, curve, f);
	square(res0, res0, f);
	add(res0, res0, h);
}



void eco_gf2n::Ytop_f_plain (ff_pol & res0, ff_pol & res1, ff_polmod & f)
{
	int d = A6.relative_degree(); // if A6 is in a smaller field, the
	// computation can be done in subfield !!
	int i;
	ff_pol curve(3);

	curve.set_coefficient(A6, 0);
	remainder(curve, curve, f.modulus());

	res1.assign_x();
	res0.assign(curve);

	for (i = 2; i <= d; i++)
		square_y_term(res0, res1, curve, f);
}



//======================================================================
// special functions for block version
//
// compute g^(2^r) with using precomputed table.
// table[i] = x^(2^r)^i mod f.

void eco_gf2n::power_block(ff_pol & g, ff_pol * table, int r)
{
	gf2n tt;
	ff_pol h(0);
	int dp(g.degree()), i, s;

	tt.assign(g.const_term());
	if (!tt.is_zero())
		for (i = 0; i < r; i++)
			square(tt, tt);

	h.set_coefficient(tt, 0);

	for (s = 1; s <= dp; s++) {
		// all other terms
		g.get_coefficient(tt, s);
		if (tt.is_zero())
			continue;

		for (i = 0; i < r; i++)
			square(tt, tt);

		add(h, h, tt*table[s]);
	}
	g.assign(h);
}



//====================================================================

const float cost_red = 3.3;
const float cost_mult_red = 1.6; // precomputed values
const float cost_square_red = 1.0;
const int MAX_BLOCK_FACTOR = 100;



void eco_gf2n::Ytop_f (ff_pol & res0, ff_pol & res1, ff_polmod & f)
{
	int d = A6.relative_degree(); // if A6 is in a smaller field, the
	// computation can be done in subfield !!
	int dp(f.modulus().degree());
	int block_factor = 1, i, j, s;
	float minimum = 2*(d-1)*cost_square_red, cost = 0.0;

	// compute optimal block factor

	j = 2;
	s = 1 << j;

	while (j <= MAX_BLOCK_FACTOR) {
		cost = cost_square_red * (j + dp/2 + 2 * (d % j))
			+ cost_mult_red * ((dp-1)/2 + 2 * d/j)
			+ cost_red * d / j;

		if (cost < minimum) {
			block_factor = j; minimum = cost;
		}
		j++; s = 1 << j;
	}

	// now do the fast exponentiation

	if (block_factor == 1) {
		Ytop_f_plain (res0, res1, f);
		return;
	}

	ff_pol curve(3);
	ff_pol x_block, y0_block, y1_block;

	curve.set_coefficient(A6, 0);
	remainder(curve, curve, f.modulus());

	x_block.assign_x(); // x^(2^bf)
	y0_block.assign_zero(); // Y^(2^bf) = y0_block + Y * y1_block
	y1_block.assign_one();

	for (i = 1; i <= block_factor; i++) {
		square(x_block, x_block, f);
		square_y_term(y0_block, y1_block, curve, f);
	}

	ff_pol *table = power_table (x_block, f, dp);
	ff_pol h;

	res1.assign_one();
	res0.assign_zero();

	for(j = 1; j <=  (d / block_factor); j++) {
		// the outer loop
		power_block(res0, table, block_factor);
		power_block(res1, table, block_factor);

		multiply(h, res1, y0_block, f);
		add(res0, res0, h);
		multiply(res1, res1, y1_block, f);
	}

	if (d % block_factor)
		for (j = d % block_factor; j; j--)
			square_y_term(res0, res1, curve, f);


	delete[] table;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
