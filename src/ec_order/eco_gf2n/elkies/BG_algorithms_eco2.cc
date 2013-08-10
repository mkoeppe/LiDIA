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
#include	"LiDIA/eco_gf2n.h"
#include	"LiDIA/weco2_rat_function.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//-----------------------------------------------------------------
//  Let M be the matrix whose columns are given as a(x) * x^i mod F
//  for i = 0, ..., F.degree()-1. The function determines the first
//  row of this matrix and returns it as vector v.

void eco_gf2n::
mult_matrix(base_vector< ff_element > & v, const ff_pol & aa,
	    ff_polmod & F)
{
	int i, d;
	ff_pol a(aa);

	d = F.modulus().degree();
	v[0].assign(a.const_term());

	for (i = 1; i < d; i++) {
		multiply_by_x_mod(a, a, F.modulus());
		v[i].assign(a.const_term());
	}
}



void eco_gf2n::inner_product(ff_element & x, const base_vector< ff_element > & a,
			     const ff_pol & b)
{
	lidia_size_t n = comparator< lidia_size_t >::min(a.size(), b.degree() + 1);
	lidia_size_t i;
	gf2n t;

	x.assign_zero();

	for (i = 0; i < n; i++) {
		b.get_coefficient(t, i);
		multiply(t, a[i], t);
		add(x, x, t);
	}
}



//******************************************************************
// This functions implements the fast search from the paper
// Vmueller and Mmaurer: "Fast ev search". The function inner_product
// is already available in Fp_polynomial.

int eco_gf2n::search_in_table(const ff_rat & f,
			      ff_rat * table, ff_polmod & F,
			      int size)
{
	int d = F.modulus().degree();
	base_vector< ff_element > v1(d, d), v2(d, d);
	ff_element r1, r2;
	ff_element res1, res2;

	if (f.is_zero() && table[0].is_zero())
		return 0;

	mult_matrix(v1, f.numerator(), F);
	mult_matrix(v2, f.denominator(), F);

	for (int i = 0; i < size; i++) {
		inner_product(r1, v1, table[i].denominator());
		inner_product(r2, v2, table[i].numerator());

		res1.assign(r1);
		res2.assign(r2);

#ifdef DEBUG
		if (equal(f, table[i], F)) {
			std::cout << "\nDEBUG: found match for i = " << i << std::flush;

			if (res1 == res2)
				std::cout << " --> OK " << std::flush;
			else
				std::cout << "\n --> ERROR " << std::flush;
		}
#endif

		if (res1 == res2) {
			// pseudo_match found
			if (equal(f, table[i], F))
				return i;
		}
	}
	return -1;
}



//---------------------------------------------------------------------
//  the Babystep Giantstep version of schoofpart

bool eco_gf2n::schoofpart_BG (lidia_size_t & ev, const ff_pol & fC,
			      lidia_size_t* klist)
{

	if (l < lower_bound_for_BG) {
		// don't use BG for such small primes !!
		char strat = ev_strategy;
		bool rc;

		ev_strategy = eco_gf2n::EV_DIVISION_POLYNOMIAL;
		rc = schoofpart(ev, fC, klist, false);
		ev_strategy = strat;
		return rc;
	}

	lidia_size_t  i, j, d;
	lidia_size_t  number_BG;
	ff_pol  xp, yp0, yp1;
	ff_polmod ffC;

	ffC.build (fC);

	//--------------------- compute number of babysteps/giantsteps first ----

	d = klist[klist[0]]; // the maximal element
	if (l - klist[1] > static_cast<unsigned int>(d))
		// take care of the sign !!
		d = l - klist[1];

	number_BG = static_cast<lidia_size_t>(std::sqrt(static_cast<double>(d) / 2.0)) + 1;

	while(2*number_BG*(number_BG + 1) < d)
		number_BG ++;

	if (info) {
		std::cout<<"\nTesting "<< d;
		if (d == 1) {
			std::cout <<" possibility using Babystep Giantstep Algorithm\n(";
			std::cout <<number_BG<<" Baby- / Giantsteps) :  "<<std::flush;
		}
		else {
			std::cout <<" possibilities using Babystep Giantstep Algorithm\n(";
			std::cout <<number_BG<<" Baby- / Giantsteps) :  "<<std::flush;
		}
	}

#ifdef TIMING
	timer t;
	t.set_print_mode();
	t.start_timer();
#endif

	Xq (xp, ffC, A6.relative_degree()); // compute frob(P) = (X, Y)^q
	Ytop_f (yp0, yp1, ffC);

#ifdef TIMING
	t.stop_timer();
	if (info)
		std::cout << "\nComputation of (X, Y)^q needs time : " << t << std::flush;
#endif

	if (xp.is_x()) {
		// match for k = +/- 1.
		if (yp0.is_zero()) {
			ev = 1;
			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is 1" << std::flush;
			return true;
		}

		if (yp0.is_x()) {
			ev = -1;
			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is -1" << std::flush;
			return true;
		}
	}

	// initialize the formal point computations

	weco2_rat_function::initialize(eco_gf2n::A6, ffC);

	weco2_rat_function P, frob_P, kP, P2, bgP;
	ff_rat *table;
	ff_rat frat;

	frob_P.assign(xp, yp0, yp1);

#ifdef DEBUG
	weco2_rat_function frob(frob_P), H;
#endif

	P.assign_xy();
	kP.assign_xy();

	// start preparing the babystep table

	table = new ff_rat [number_BG]; // table[i] = (i+1) * P
	table[0].assign(P.get_y0());
	if (info)
		std::cout << "B " << std::flush;

	multiply_by_2(kP, P);
	P2.assign(kP);
	table[1].assign(kP.get_y0());
	if (info)
		std::cout << "B " << std::flush;

	if (equal(frob_P.get_x(), kP.get_x(), ffC)) {
		// + / - 2
		if (equal(frob_P.get_y0(), kP.get_y0(), ffC)) {
			ev = 2;
			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is " << ev << std::flush;
			delete[] table;
			return true;
		}
		else
			if (equal(frob_P.get_y1(), kP.get_y1(), ffC)) {
				ev = -2;
				if (info)
					std::cout << "\n\nEigenvalue of Frobenius is " << ev << std::flush;
				delete[] table;
				return true;
			}
	}


	// Babystep Table is computed with doubling trick

	i = 4; j = 3;
	while(j <= number_BG) {
		while(i <= number_BG) {
			multiply_by_2(kP, kP);
			if (i == number_BG - 1)
				bgP.assign(kP);

			if (equal(frob_P.get_x(), kP.get_x(), ffC)) {
				if (equal(frob_P.get_y0(), kP.get_y0(), ffC)) {
					ev = i;
					if (info)
						std::cout << "\n\nEigenvalue of Frobenius is " << ev << std::flush;
					delete[] table;
					return true;
				}
				else
					if (equal(frob_P.get_y1(), kP.get_y1(), ffC)) {
						ev = -i;
						if (info)
							std::cout << "\n\nEigenvalue of Frobenius is " << ev << std::flush;
						delete[] table;
						return true;
					}
			}

#ifdef DEBUG
			H.assign_xy();
			multiply(H, i, H);
			assert(kP == H);
#endif

			table[i-1].assign(kP.get_y0());
			i <<= 1;
			if (info)
				std::cout << "B " << std::flush;
		}
		i = j; j += 2;
		add(kP, P, P2);

		if (i == number_BG - 1)
			bgP.assign(kP);

		if (equal(frob_P.get_x(), kP.get_x(), ffC)) {
			if (equal(frob_P.get_y0(), kP.get_y0(), ffC)) {
				ev = i;
				if (info)
					std::cout << "\n\nEigenvalue of Frobenius is " << ev << std::flush;
				delete[] table;
				return true;
			}
			else
				if (equal(frob_P.get_y1(), kP.get_y1(), ffC)) {
					ev = -i;
					if (info)
						std::cout << "\n\nEigenvalue of Frobenius is " << ev << std::flush;
					delete[] table;
					return true;
				}
		}
#ifdef DEBUG
		H.assign_xy();
		multiply(H, i, H);
		assert(kP == H);
#endif
		P.assign(kP);
		table[i-1].assign(kP.get_y0());

		if (info)
			std::cout << "B " << std::flush;
		i <<= 1;
	}               // end of table computation

	if (info)
		std::cout << "  ";

	kP.assign_xy();
	add(bgP, bgP, kP); // = number_BG * (X, Y)
	subtract(frob_P, frob_P, bgP);
	multiply_by_2(kP, bgP);
	negate(kP, kP);

	for (j = 0; j <= number_BG; j++) {
		// iteration j:  frob_P = (X, Y)^p - (2*j+1)
		// * number_BG * (X, Y)
		if (info)
			std::cout << "G " << std::flush;

		if (frob_P.is_zero()) {
			ev = (2*j+1) * number_BG;

#ifdef DEBUG
			P.assign_xy();
			multiply(H, ev, P);
			assert(frob == H);
#endif

			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is " << ev << std::flush;
			delete[] table;
			return true;
		}

		i = search_in_table(frob_P.get_y0(), table, ffC, number_BG);

		if (i != -1) {
			// match found: frob_P = (i+1) * (X, Y)
			ev = (2*j+1) * number_BG + i + 1;

#ifdef DEBUG
			P.assign_xy();
			multiply(H, ev, P);
			assert(frob == H);
#endif

			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is " << ev << std::flush;
			delete[] table;
			return true;
		}

		add(frat, frob_P.get_y0(), frob_P.get_x(), ffC);

		i = search_in_table(frat, table, ffC, number_BG);

		if (i != -1) {
			// match found: frob_P = -(i+1) * (X, Y)
			ev = (2*j+1) * number_BG - i - 1;

#ifdef DEBUG
			P.assign_xy();
			multiply(H, ev, P);
			assert(frob == H);
#endif

			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is " << ev << std::flush;
			delete[] table;
			return true;
		}

		add(frob_P, frob_P, kP);
	}

	delete[] table;
	lidia_error_handler("eco_gf2n", "schoofpart_BG::found no solution");
	return false;
}



//********************************************************************
//
// Now all the functions for the funny Babystep Giantstep idea are coming.
//


//--------------------------------------------------------------------
// The following function doubles the point with x-coordinate x1 and
// returns only the x-coordinate of the result (using formulas for
// projective coordinates)

void eco_gf2n::double_point_x_only(ff_rat & x3, const ff_rat & x1,
				   ff_polmod & F)
{
	ff_pol x(x1.numerator()), z(x1.denominator()), h;
	ff_element sqrt4_a6;

	sqrt(sqrt4_a6, A6); sqrt(sqrt4_a6, sqrt4_a6);
	multiply(h, x, z, F);
	square(h, h, F); // denominator(x3)

	multiply_by_scalar(z, z, sqrt4_a6);
	add(x, x, z);
	square(x, x, F);
	square(x, x, F);
	x3.assign(x, h);
}



//------------------------------------------------------------------
// Since the cost for Babysteps are lower than the cost for Giantsteps
// (double_point_x_only needs 3 squarings and 1 mult --> weight 1.9,
// Giantsteps with divpol:  5.5 mults + 1 square --> weight 5.8)
// we use a weighted formula for computing costs.
// Heuristik: square = 30 % of multiplication

void eco_gf2n::
find_number_of_FBG_steps(lidia_size_t & baby, lidia_size_t & giant,
			 double & average,
			 lidia_size_t *klist)
{
	int i, j;
	float min;
	int b_tmp, g_tmp;
	lidia_size_t alpha;
	const float weight_B = 1.9;
	const float weight_G = 5.8;

	baby = giant = 0;

	for (i = 1; i <= klist[0]; i++) {
		// compute maximal exponent 2^j*frob
		// if alpha = klist[i], set b_tmp
		min = weight_G * klist[i];
		alpha = klist[i];
		b_tmp = 0;
		g_tmp = alpha;

		for (j = 1; j < static_cast<int>(l); j ++) {
			alpha  = (alpha << 1) % l;
			alpha = comparator< lidia_size_t >::min(alpha, l-alpha);
			if (weight_B * j + weight_G * alpha  < min)    // 2^j * ev = alpha mod l
			{
				min = weight_B * j + weight_G * alpha;
				b_tmp = j;
				g_tmp = alpha;
			}
		}

		if (b_tmp > baby)
			baby = b_tmp;

		if (g_tmp > giant)
			giant = g_tmp;
	}

	average = 0.0;
	giant = 0;

	for (i = 1; i <= klist[0]; i++) {
		alpha = klist[i];
		g_tmp = alpha;

		for (j = 1; j <= baby; j++) {
			alpha = (alpha << 1) % l;
			alpha = comparator< lidia_size_t >::min(alpha, l-alpha);
			if (alpha < g_tmp)
				g_tmp = alpha;
		}
		if (g_tmp > giant)
			giant = g_tmp;
		average += static_cast<double>(weight_B * baby + weight_G * g_tmp);
	}
	average /= static_cast<double>(klist[0]);
	average /= static_cast<double>(weight_G);
}



//----------------------------------------------------------------------

bool eco_gf2n::schoofpart_FBG (lidia_size_t & ev, const ff_pol & fC,
			       lidia_size_t* klist, bool c_sign)
{
	lidia_size_t   i, j;
	lidia_size_t   number_baby, number_giant;
	double average;
	ff_pol  xp_pol; // X^q mod fC
	ff_rat  xp; // and as rational function
	ff_polmod ffC;

	if (l < lower_bound_for_FBG) {
		char strat = ev_strategy;
		bool rc;

		ev_strategy = eco_gf2n::EV_DIVISION_POLYNOMIAL;
		rc = schoofpart(ev, fC, klist, false);
		ev_strategy = strat;
		return rc;
	}

	ffC.build (fC);
	find_number_of_FBG_steps(number_baby, number_giant, average, klist);

	if (info) {
		std::cout << "\nTesting " << klist[0];
		if (klist[0] == 1)
			std::cout << " possibility using Funny Babystep Giantstep Algorithm\n(";
		else
			std::cout << " possibilities using Funny Babystep Giantstep Algorithm\n(";
		std::cout << number_baby << " Babysteps, " << number_giant << " Giantsteps, ";
		std::cout << average << " expected point additions) :   " << std::flush;
	}

#ifdef TIMING
	timer t;
	t.set_print_mode();
	t.start_timer();
#endif

	Xq (xp_pol, ffC, A6.relative_degree());

#ifdef TIMING
	t.stop_timer();
	if (info)
		std::cout << "\nComputing X^q mod f_C(X) needs time : " << t << std::flush;
#endif

	xp.assign(xp_pol);

	if (xp_pol.is_x()) {
		// match for k = +/- 1.
		ev = 1;
		if (info)
			std::cout << "\n\nEigenvalue of Frobenius is +/- 1" << std::flush;
		if (c_sign)
			compute_sign(ev, fC);
		return true;
	}

	PsiPowers* psi_pow;
	ff_rat *table;
	ff_pol num, den;
	int top;

	init_psi_comp (psi_pow, top, number_giant + 1, ffC);

	// now I'm testing first whether (X^q, .) = alpha *(X, Y) with
	// alpha even, 2 <= alpha <= number_giant

	for (i = 2; i <= number_giant; i += 2) {
		while ((i+1) > top)
			next_psi(psi_pow, top, ffC);

		shift_left (num, *(psi_pow[i].pow2), ffC);
		multiply (den, *(psi_pow[i+1].pow1), *(psi_pow[i-1].pow1),
			  ffC);
		add (num, num, den);
		multiply(den, xp_pol, *(psi_pow[i].pow2), ffC);

		if (num == den) {
			ev = i;
			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is +/- " << ev << std::flush;
			free_psi(psi_pow, 0, number_giant + 4);

			if (c_sign)
				compute_sign(i, fC);
			return true;
		}
	}

	// next we determine a table of babysteps as
	// table[i] = x-coord(2^j * (X^q, .))

	table = new ff_rat [number_baby+1];
	table[0].assign(xp);

	for (i = 1; i <= number_baby; i++) {
		double_point_x_only(xp, xp, ffC);

		table[i].assign(xp);
		if (info)
			std::cout << "B " << std::flush;
	}

	// now do Giantsteps, but necessary only for odd multipliers !!

	udigit uinv;

	if (info)
		std::cout << "  ";

	for (j = 1; j <= number_giant; j += 2) {
		while ((j+1) > top)
			next_psi(psi_pow, top, ffC);

		if (info)
			std::cout << "G " << std::flush;

		shift_left (num, *(psi_pow[j].pow2), ffC);
		multiply (den, *(psi_pow[j+1].pow1), *(psi_pow[j-1].pow1),
			  ffC);
		add (num, num, den);

		i = search_in_table(ff_rat(num, *(psi_pow[j].pow2)), table,
				    ffC, number_baby + 1);

		if (i != -1) {
			// match found
			uinv = power_mod(static_cast<udigit>(2),
						static_cast<udigit>(i),
						static_cast<udigit>(l));
			uinv = invert_mod(uinv, static_cast<udigit>(l));
			ev = static_cast<int>(uinv * j) % l;
			delete[] table;
			free_psi(psi_pow, 0, number_giant + 4);
			if (info)
				std::cout << "\n\nEigenvalue of Frobenius is +/- " << ev << std::flush;

			if (c_sign)
				compute_sign(ev, fC);
			return true;
		}
	}

	delete[] table;
	free_psi(psi_pow, 0, number_giant + 4);
	lidia_error_handler("eco_gf2n", "FBG::nothing found");
	return false;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
